"""Calculators for the relaxation stage.

Machine-learned backends need PyTorch; it is imported only inside the functions that build those
calculators, so that the parent process of a run never loads it (see ``relax.WorkerPool``).
"""

from __future__ import annotations

import os
import time

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
