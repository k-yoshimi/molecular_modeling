"""Calculators for the relaxation stage.

Machine-learned backends need PyTorch; it is imported only inside the functions that build those
calculators, so that the parent process of a run never loads it (see ``relax.WorkerPool``).
"""

from __future__ import annotations

import importlib.util
import multiprocessing
import os
import sys
import time
from importlib.metadata import PackageNotFoundError, version

from ase.calculators.calculator import Calculator

from .config import CONFIG_RELAX_METHODS

#: Methods that can be chosen in the config file or on the command line.
CONFIG_METHODS: tuple[str, ...] = CONFIG_RELAX_METHODS
#: Methods for tests only: ``run_csp(..., method=...)`` accepts them, the config file does not.
TEST_METHODS: tuple[str, ...] = ("EMT", "EMT-CRASH-ON-AU", "EMT-BROKEN", "EMT-SLOW")
#: Every name accepted by :func:`make_calculator`.
METHODS: tuple[str, ...] = CONFIG_METHODS + TEST_METHODS
#: Machine-learned methods: they need the ``csp-ml`` extra and neutral singlet ``atoms.info``.
ML_METHODS: frozenset[str] = frozenset({"OrbMol-v2"})

#: Methods checked by :func:`preflight` before a run touches its output directory.
PREFLIGHT_METHODS: frozenset[str] = ML_METHODS

ML_INSTALL_HINT = 'install the machine-learned relaxation dependencies with: pip install -e ".[csp-ml]"'


class RelaxBackendError(RuntimeError):
    """The calculator of the chosen relaxation method cannot be built or used."""


def make_calculator(method: str = "GFN2-xTB") -> Calculator:
    """Return a new ASE calculator for ``method`` (one of :data:`METHODS`).

    Raises
    ------
    ValueError
        If ``method`` is not supported.
    """
    if method in ("GFN2-xTB", "GFN1-xTB"):
        from tblite.ase import TBLite

        return TBLite(method=method, verbosity=0)
    if method == "OrbMol-v2":
        quiet_warp()
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.inference.calculator import ORBCalculator

        model, adapter = pretrained.orbmol_v2(device="cpu")
        return ORBCalculator(model, adapter, device="cpu")
    if method in TEST_METHODS:
        return _test_calculator(method)
    raise ValueError(f"unsupported method {method!r}; use one of {', '.join(CONFIG_METHODS)}")


#: Calculators of the machine-learned methods built in this (worker) process, by method name.
_CACHE: dict[str, Calculator] = {}


def get_calculator(method: str) -> Calculator:
    """Calculator for one relaxation in this process.

    Machine-learned methods (:data:`ML_METHODS`) load their model once per process and are reused
    (after ``reset()``); every other method gets a new calculator each time, as before, so that no
    state carries over from one candidate to the next.
    """
    if method not in ML_METHODS:
        return make_calculator(method)
    calc = _CACHE.get(method)
    if calc is None:
        calc = _CACHE[method] = make_calculator(method)
    calc.reset()
    return calc


def _test_calculator(method: str) -> Calculator:
    from ase.calculators.emt import EMT

    if method == "EMT-CRASH-ON-AU":  # kills the worker process when gold is present

        class _CrashOnAu(EMT):
            def calculate(self, atoms=None, properties=("energy",), system_changes=("positions",)):
                if atoms is not None and 79 in atoms.numbers:
                    os._exit(1)
                super().calculate(atoms, properties, system_changes)

        return _CrashOnAu()
    if method == "EMT-BROKEN":  # stands for a backend whose model cannot be loaded
        raise RuntimeError("test backend is broken")
    if method == "EMT-SLOW":  # stands for a backend that hangs while loading
        time.sleep(30)
    return EMT()


def quiet_warp() -> None:
    """Show only warnings and errors of NVIDIA Warp (used by orb-models), not its start-up messages."""
    if "warp" not in sys.modules and importlib.util.find_spec("warp") is None:
        return
    import warp

    warp.config.log_level = max(warp.config.log_level, warp.LOG_WARNING)


def weights_cache_dir() -> str:
    """Directory in which orb-models keeps downloaded weights.

    As in ``cached_path``: ``CACHED_PATH_CACHE_ROOT`` whenever it is set (an empty value means the
    current directory), else ``~/.cache/cached_path``.
    """
    root = os.environ.get("CACHED_PATH_CACHE_ROOT")
    if root is None:
        return os.path.expanduser(os.path.join("~", ".cache", "cached_path"))
    return os.path.abspath(root)


def _weights_info(method: str) -> dict[str, str]:
    """Weights file and its local directory, for the methods that download weights."""
    if method != "OrbMol-v2":
        return {}
    import inspect

    from cached_path import cached_path
    from orb_models.forcefield import pretrained

    url = str(inspect.signature(pretrained.orbmol_v2).parameters["weights_path"].default)
    local = str(cached_path(url))
    return {"weights": url.rsplit("/", 1)[-1], "cache": os.path.dirname(local)}


def exercise(calc: Calculator) -> None:
    """Compute energy, forces and stress of one water molecule in a periodic cell with ``calc``.

    The first call of a machine-learned calculator is slow (the model compiles itself); this makes it
    happen here instead of in the first relaxation.
    """
    from ase.build import molecule

    atoms = molecule("H2O")
    atoms.set_cell([10.0, 10.0, 10.0])
    atoms.center()
    atoms.pbc = True
    atoms.info.update(charge=0, spin=1)
    atoms.calc = calc
    atoms.get_potential_energy()
    atoms.get_forces()
    atoms.get_stress()
    atoms.calc = None


def _preflight_child(method: str, conn) -> None:
    """Build the calculator and evaluate one water molecule in a periodic cell; send the outcome."""
    try:
        exercise(make_calculator(method))
        conn.send(("ok", _weights_info(method)))
    except Exception as exc:  # noqa: BLE001 - reported to the parent as text
        conn.send(("error", f"{type(exc).__name__}: {exc}"))
    finally:
        conn.close()


def preflight(method: str, timeout_s: float = 600.0) -> dict[str, str]:
    """Check in a separate (spawned) process that ``method`` can build its calculator and compute.

    Building the calculator downloads the weights on first use. The child process ends after the
    check, so no model stays in this process. Returns ``{"weights": ..., "cache": ...}`` for methods
    that download weights, else an empty dict.

    Raises
    ------
    RelaxBackendError
        One line naming the method, the cause, the weights cache and what to do without network
        access, if the check fails, the child dies, or no answer comes within ``timeout_s`` seconds.
    """
    ctx = multiprocessing.get_context("spawn")
    recv, send = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=_preflight_child, args=(method, send), daemon=True)
    proc.start()
    send.close()
    try:
        if recv.poll(timeout_s):
            try:
                status, payload = recv.recv()
            except EOFError:
                proc.join(5)
                status, payload = "error", f"the check process died (exit code {proc.exitcode})"
        else:
            status, payload = "error", f"TimeoutError: no answer within {timeout_s:g} s"
    finally:
        if proc.is_alive():
            proc.terminate()
        proc.join(5)
        recv.close()
    if status == "ok":
        return dict(payload)
    cause = " ".join(str(payload).split())
    cache = weights_cache_dir()
    raise RelaxBackendError(
        f"relaxation method {method} is not usable: {cause}; weights are cached in {cache}; if the weights "
        'could not be downloaded because this computer has no network access, run python -c "from '
        "orb_models.forcefield import pretrained; "
        "pretrained.orbmol_v2(device='cpu')\" on one with access and copy the "
        f"{cache} directory"
    )


def backend_info(method: str) -> str:
    """Package and version behind ``method``, e.g. ``"tblite 0.5.0"`` (read from metadata, not imported)."""
    if method in TEST_METHODS:
        return "test calculator"
    package = "orb-models" if method in ML_METHODS else "tblite"
    try:
        return f"{package} {version(package)}"
    except PackageNotFoundError:
        return f"{package} (version unknown)"
