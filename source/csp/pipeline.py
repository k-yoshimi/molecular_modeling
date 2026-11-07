"""End-to-end crystal structure prediction: generate, screen, relax, deduplicate, write."""

from __future__ import annotations

import csv
import importlib.util
import logging
import math
import os
import shutil
import sys
import time
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write

from .backends import METHODS, ML_INSTALL_HINT, ML_METHODS, PREFLIGHT_METHODS, backend_info, preflight
from .candidate import MAX_CIF, Candidate
from .config import SEED_GENERATOR_WEIGHT, CSPConfig, SeedSpec
from .dedupe import deduplicate
from .deps import require
from .gen_random import default_z, generate_random
from .gen_seed import generate_from_seed
from .gen_voronoi import orient_candidate
from .molecule import load_molecule, molecule_volume
from .outdir import INCOMPLETE_NAME, check_single_process, hold_lock, mark_incomplete
from .relax import WorkerPool, relax_many
from .screen import check_supported_elements, prescreen

logger = logging.getLogger(__name__)

CSV_HEADER = (
    "rank,id,energy_per_mol_eV,rel_energy_kJmol,density_gcm3,space_group_in,space_group_out,origin,max_mol_rmsd,status,"
    "relax_method"
)
EV_TO_KJMOL = 96.485
SAFETY_MARGIN = 0.9


#: Packages of the ``csp`` extra needed by every run; ``tblite`` is needed in addition for GFN methods.
REQUIRED_PACKAGES = ("pyxtal", "pymatgen", "spglib")
#: Packages of the ``csp-ml`` extra needed by the machine-learned methods; only their presence is checked,
#: so that the parent process never imports PyTorch.
ML_PACKAGES = ("orb_models", "torch")


def check_dependencies(method: str = "GFN2-xTB") -> None:
    """Raise ``ImportError`` naming every missing package needed for ``method``.

    An unknown ``method`` raises ``ValueError``.

    For the GFN methods ``tblite.ase`` (which loads the compiled tblite library) must import and provide
    ``TBLite``; importing the root ``tblite`` package alone does not prove the calculator works. For the
    machine-learned methods the packages of the ``csp-ml`` extra are only looked up, not imported;
    whether they work is checked by :func:`~molecular_modeling.csp.backends.preflight` in a child process.
    """
    if method not in METHODS:
        raise ValueError(f"unsupported method {method!r}; use one of {', '.join(METHODS)}")
    if method in ML_METHODS:
        require(list(REQUIRED_PACKAGES), purpose="crystal structure prediction")
        missing = [name for name in ML_PACKAGES if importlib.util.find_spec(name) is None]
        if missing:
            py = "" if sys.version_info >= (3, 12) else " (it needs Python 3.12 or newer)"
            raise ImportError(
                f"relaxation method {method}{py} needs {', '.join(missing)}, which cannot be imported; "
                f"{ML_INSTALL_HINT}"
            )
        return
    names = list(REQUIRED_PACKAGES) + (["tblite", "tblite.ase"] if method.startswith("GFN") else [])
    require(names, {"tblite.ase": "TBLite"}, purpose="crystal structure prediction")


@contextmanager
def _run_log(path: Path) -> Iterator[logging.Handler]:
    """Send ``molecular_modeling.csp`` log records to ``path`` and stderr for the duration of the run.

    Yields the file handler, so that records can be written to the file alone.
    """
    root = logging.getLogger("molecular_modeling.csp")
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handlers: list[logging.Handler] = [file_handler, logging.StreamHandler()]
    old_level = root.level
    for h in handlers:
        h.setFormatter(fmt)
        h.setLevel(logging.INFO)
        root.addHandler(h)
    root.setLevel(logging.INFO)
    try:
        yield file_handler
    finally:
        for h in handlers:
            root.removeHandler(h)
            h.close()
        root.setLevel(old_level)


def _tag(origin: str) -> str:
    return origin.replace(":", "_").replace("/", "_")


def _spacegroup_out(atoms: Atoms, symprec: float = 0.1) -> int:
    """Space group number of ``atoms`` found by spglib (0 if detection fails)."""
    import spglib

    ds = spglib.get_symmetry_dataset(
        (atoms.get_cell()[:], atoms.get_scaled_positions(), atoms.get_atomic_numbers()), symprec=symprec
    )
    if ds is None:
        return 0
    return int(ds["number"] if isinstance(ds, dict) else ds.number)


def _split_counts(cfg: CSPConfig) -> tuple[int, int]:
    """Per-space-group counts for the random and Voronoi generators."""
    n = cfg.n_candidates_per_sg
    wr = float(cfg.generators.get("random", 0.0))
    wv = float(cfg.generators.get("voronoi", 0.0))
    if wr + wv <= 0:
        return 0, 0
    n_rand = round(n * wr / (wr + wv))
    return n_rand, n - n_rand


def _seed_candidates(cfg: CSPConfig, mol: Atoms) -> list[tuple[SeedSpec, list[Candidate]]]:
    """Candidates from every seed structure (none if the seed generator weight is 0).

    Raises ``ValueError`` naming the seed file if it cannot be read or does not fit the molecule.
    """
    if float(cfg.generators.get("seed", SEED_GENERATOR_WEIGHT)) <= 0:
        return []
    out = []
    for seed in cfg.seed_structures:
        try:
            out.append((seed, generate_from_seed(mol, seed)))
        except (OSError, ValueError, AssertionError, StopIteration, IndexError, KeyError) as exc:
            raise ValueError(f"seed structure {seed.cif}: {type(exc).__name__}: {exc}") from exc
    return out


def _generate(
    cfg: CSPConfig,
    mol: Atoms,
    rng: np.random.Generator,
    seeded: list[tuple[SeedSpec, list[Candidate]]] | None = None,
) -> list[Candidate]:
    """Generate all candidates (random, Voronoi-oriented and seed-based).

    ``random`` and ``voronoi`` weights split ``n_candidates_per_sg`` between those generators. The
    ``seed`` weight only switches the seed generator on (> 0, the default when absent) or off (0):
    a seed always yields one candidate per substitution choice. ``seeded`` are seed candidates already
    made by :func:`_seed_candidates` (made here if ``None``).
    """
    n_rand, n_vor = _split_counts(cfg)
    out: list[Candidate] = []
    for sg in cfg.space_groups:
        z = cfg.z[sg] if cfg.z is not None and sg in cfg.z else default_z(sg, cfg.z_prime)
        rand: list[Candidate] = []
        if n_rand:
            rand = generate_random(mol, sg, z, n_rand, rng)
            if len(rand) < n_rand:
                logger.info("space group %d: random generator returned %d of %d requested", sg, len(rand), n_rand)
        vor: list[Candidate] = []
        if n_vor:
            # each skeleton gives two orientations, or one when the turned copy would be identical (a
            # symmetric molecule); draw more skeletons until the share is filled, at most n_vor of them
            n_skel = 0
            while len(vor) < n_vor and n_skel < n_vor:
                skeletons = generate_random(mol, sg, z, min(math.ceil((n_vor - len(vor)) / 2), n_vor - n_skel), rng)
                if not skeletons:
                    break
                n_skel += len(skeletons)
                for sk in skeletons:
                    vor.extend(orient_candidate(sk))
            vor = vor[:n_vor]
            if len(vor) < n_vor:
                logger.info("space group %d: voronoi generator returned %d of %d requested", sg, len(vor), n_vor)
        logger.info("space group %d (Z=%d): random %d, voronoi %d", sg, z, len(rand), len(vor))
        out.extend(rand + vor)
    seed_weight = float(cfg.generators.get("seed", SEED_GENERATOR_WEIGHT))
    if cfg.seed_structures and seed_weight <= 0:
        logger.info("seed generator weight is %g: %d seed structure(s) not used", seed_weight, len(cfg.seed_structures))
    if seeded is None:
        seeded = _seed_candidates(cfg, mol)
    for seed, seed_cands in seeded:
        logger.info("seed %s: %d candidates", seed.cif, len(seed_cands))
        out.extend(seed_cands)
    counters: dict[tuple[str, int], int] = defaultdict(int)
    for c in out:
        key = (_tag(c.origin), c.space_group)
        c.meta["id"] = f"{key[0]}-{key[1]}-{counters[key]}"
        counters[key] += 1
    return out


def _relax_with_budget(
    cands: list[Candidate], mol: Atoms, cfg: CSPConfig, method: str, t_start: float
) -> list[Candidate]:
    """Relax ``cands`` in rounds of ``cfg.n_procs``, stopping when the next round would exceed the budget.

    All rounds run on one :class:`~molecular_modeling.csp.relax.WorkerPool`, made ready (workers started,
    calculators loaded) before the first round and again after a worker died. That time counts against
    the budget but not as part of a round. The first round always runs. Before each further round the
    remaining time is compared with the longest round so far; if it is shorter, the remaining candidates
    are dropped.
    """
    budget = cfg.time_budget_min * 60 * SAFETY_MARGIN
    n_round = max(cfg.n_procs, 1)
    done: list[Candidate] = []
    rest = list(cands)
    longest = 0.0  # wall time of the slowest round so far
    n_done_rounds = 0
    with WorkerPool(n_round, method) as pool:
        while rest:
            pool.ensure()
            if n_done_rounds:
                remaining = budget - (time.perf_counter() - t_start)
                if remaining < longest:
                    logger.info(
                        "time budget: %.1f s per candidate, %.0f s left; relaxing 0 of %d remaining candidates",
                        longest / n_round,
                        max(remaining, 0.0),
                        len(rest),
                    )
                    break
            batch, rest = rest[:n_round], rest[n_round:]
            t0 = time.perf_counter()
            done = done + relax_many(batch, mol, cfg, method, pool=pool)
            longest = max(longest, time.perf_counter() - t0)
            n_done_rounds += 1
    return done


#: Directory (inside ``out_dir``) in which the results are written before they are moved into place.
STAGING_NAME = ".csp_staging"


def _remove_results(out_dir: Path) -> None:
    """Delete the results of an earlier run: ``candidates.csv`` first, then ranked CIFs and ``validation.md``.

    Callers mark the directory incomplete before calling this, so a partly cleaned directory is never
    taken for a finished run; deleting the CSV first also keeps it from listing CIFs that are gone.
    """
    (out_dir / "candidates.csv").unlink(missing_ok=True)
    for old in [*out_dir.glob("[0-9][0-9][0-9]_*.cif"), out_dir / "validation.md"]:
        old.unlink(missing_ok=True)
    shutil.rmtree(out_dir / STAGING_NAME, ignore_errors=True)


def _mark_complete(out_dir: Path) -> None:
    """Remove the incomplete marker: the very last step of publishing a run's results."""
    (out_dir / INCOMPLETE_NAME).unlink(missing_ok=True)


def _write_outputs(cfg: CSPConfig, ranked: list[Candidate], method: str) -> None:
    """Write ``candidates.csv`` and CIFs of the top :data:`MAX_CIF` candidates.

    Every row of ``candidates.csv`` names the relaxation ``method`` in its last column.

    The directory is marked incomplete (:data:`~molecular_modeling.csp.outdir.INCOMPLETE_NAME`)
    first; then results of an earlier run (``candidates.csv``, ranked CIFs ``NNN_<id>.cif``,
    ``validation.md``) are removed, so the directory only holds files of this run. Everything is
    written into a staging directory; the CIFs are then moved into place, ``candidates.csv`` comes
    next, and the marker is removed last. If any step fails, the marker stays and validation refuses
    the directory.
    """
    out_dir = Path(cfg.out_dir)
    mark_incomplete(out_dir)
    _remove_results(out_dir)
    stage = out_dir / STAGING_NAME
    stage.mkdir()
    emin = ranked[0].energy_per_mol if ranked else 0.0
    rows: list[list[str]] = []
    for rank, c in enumerate(ranked, start=1):
        assert c.energy_per_mol is not None and emin is not None
        rel = (c.energy_per_mol - emin) * EV_TO_KJMOL
        dens = "" if c.density is None else f"{c.density:.4f}"
        rmsd = c.meta.get("max_mol_rmsd")
        rows.append(
            [
                str(rank),
                c.meta["id"],
                f"{c.energy_per_mol:.6f}",
                f"{rel:.3f}",
                dens,
                str(c.space_group),
                str(_spacegroup_out(c.atoms)),
                c.origin,
                "" if rmsd is None else f"{rmsd:.4f}",
                c.status,
                method,
            ]
        )
        if rank <= MAX_CIF:
            write(stage / f"{rank:03d}_{c.meta['id']}.cif", c.atoms)
    # csv.writer quotes fields that contain commas or quotes (e.g. from seed file names)
    with open(stage / "candidates.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(CSV_HEADER.split(","))
        writer.writerows(rows)
    for cif in sorted(stage.glob("[0-9][0-9][0-9]_*.cif")):
        os.replace(cif, out_dir / cif.name)
    os.replace(stage / "candidates.csv", out_dir / "candidates.csv")
    stage.rmdir()
    _mark_complete(out_dir)


def _log_reasons(relaxed: list[Candidate]) -> None:
    """Log, per status and reason, how many candidates failed or were rejected, with the first message."""
    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for c in relaxed:
        if c.status == "failed":
            msg = str(c.meta.get("error", "unknown error")).splitlines()[0] if c.meta.get("error") else "unknown error"
            groups[("failed", msg)].append(msg)
        elif c.status == "rejected_deformed":
            msg = str(c.meta.get("reject_reason", "unknown reason"))
            groups[("rejected_deformed", msg.split(" (")[0])].append(msg)
    for (status, _key), msgs in groups.items():
        logger.info("  %d %s: %s", len(msgs), status, msgs[0])


def _run_stages(
    cfg: CSPConfig,
    mol: Atoms,
    rng: np.random.Generator,
    seeded: list[tuple[SeedSpec, list[Candidate]]],
    method: str,
    out_dir: Path,
) -> list[Candidate]:
    """Generate, prescreen, relax, deduplicate and write; return the ranked candidates."""
    t_start = time.perf_counter()
    t = time.perf_counter()
    cands = _generate(cfg, mol, rng, seeded)
    by_origin: dict[str, int] = defaultdict(int)
    for c in cands:
        by_origin[c.origin.split(":")[0]] += 1
    logger.info("generated %d candidates %s in %.1f s", len(cands), dict(by_origin), time.perf_counter() - t)

    t = time.perf_counter()
    kept, counts = prescreen(cands, molecule_volume(mol), cfg)
    logger.info("prescreen %s in %.1f s", counts, time.perf_counter() - t)

    t = time.perf_counter()
    relaxed = _relax_with_budget(kept, mol, cfg, method, t_start)
    n_ok = sum(c.status == "relaxed" for c in relaxed)
    n_rej = sum(c.status == "rejected_deformed" for c in relaxed)
    n_fail = sum(c.status == "failed" for c in relaxed)
    logger.info(
        "relax: %d relaxed, %d rejected_deformed, %d failed in %.1f s",
        n_ok,
        n_rej,
        n_fail,
        time.perf_counter() - t,
    )
    _log_reasons(relaxed)

    t = time.perf_counter()
    unique = deduplicate([c for c in relaxed if c.status == "relaxed" and c.energy_per_mol is not None])
    ranked = sorted(unique, key=lambda c: c.energy_per_mol if c.energy_per_mol is not None else math.inf)
    ranked = [replace(c) for c in ranked]
    logger.info("dedupe: %d -> %d in %.1f s", n_ok, len(ranked), time.perf_counter() - t)

    _write_outputs(cfg, ranked, method)
    logger.info("wrote %d rows to %s; total %.1f s", len(ranked), out_dir, time.perf_counter() - t_start)
    return ranked


def resolve_method(cfg: CSPConfig, method: str | None) -> tuple[str, str]:
    """The relaxation method of a run and where it came from.

    ``method`` (an argument of :func:`run_csp`) overrides ``cfg.relax_method``; the source is then
    ``"argument"``, otherwise ``cfg.relax_method_source``.
    """
    if method is None:
        return cfg.relax_method, cfg.relax_method_source
    return method, "argument"


def run_csp(cfg: CSPConfig, method: str | None = None) -> list[Candidate]:
    """Run the full CSP pipeline and write ``candidates.csv``, CIF files and ``run.log``.

    Parameters
    ----------
    cfg : CSPConfig
        Run settings. Output goes to ``cfg.out_dir``.
    method : str, optional
        Relaxation method; overrides ``cfg.relax_method`` (``"EMT"`` is available for fast tests).

    Returns
    -------
    list of Candidate
        Deduplicated relaxed candidates, lowest energy first.

    Raises
    ------
    ImportError
        Before anything is written, if a package of the ``csp`` extra needed for ``method`` is missing.
    ValueError, FileNotFoundError
        Before anything is written: under MPI with more than one process, for an unknown ``method``, an unreadable or unsupported molecule, or
        a seed structure that cannot be read or used.
    FileExistsError
        If another run holds ``out_dir`` (lock file ``.csp_run.lock``).
    RelaxBackendError
        Before anything is written, if a machine-learned ``method`` cannot build its calculator (e.g. the
        weights cannot be downloaded); also later if the relaxation workers cannot start.

    Notes
    -----
    The results of an earlier run in ``out_dir`` (``candidates.csv``, ranked CIFs, ``validation.md``)
    are deleted when the run starts. If the run stops with an error or is interrupted, the error is
    written to ``run.log`` and no results are left, so validation cannot mistake old results for new ones.
    """
    check_single_process()
    method, source = resolve_method(cfg, method)
    check_dependencies(method)
    # everything that can fail on bad input is checked before out_dir is touched, so a failed rerun
    # leaves the previous (complete) results and their run.log as they were
    mol = load_molecule(cfg.molecule)
    check_supported_elements(mol)
    rng = np.random.default_rng(cfg.rng)
    seeded = _seed_candidates(cfg, mol)
    # last, since it takes tens of seconds: mistakes in the input are reported without that wait
    backend = preflight(method) if method in PREFLIGHT_METHODS else {}

    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with hold_lock(out_dir, "molecular-modeling-csp"), _run_log(out_dir / "run.log") as file_log:
        # warnings from loading the config were logged before run.log existed: repeat them in the file only
        for msg in cfg.config_warnings:
            file_log.handle(logger.makeRecord(logger.name, logging.WARNING, __file__, 0, "config: %s", (msg,), None))
        weights = f"; weights {backend['weights']} in {backend['cache']}" if backend else ""
        logger.info("relax method %s (from %s); %s%s", method, source, backend_info(method), weights)
        try:
            mark_incomplete(out_dir)  # before anything of the earlier results is changed
            _remove_results(out_dir)
            ranked = _run_stages(cfg, mol, rng, seeded, method, out_dir)
        except BaseException as exc:
            logger.error(
                "run stopped (%s: %s); the results are incomplete (marked by %s): rerun",
                type(exc).__name__,
                exc,
                INCOMPLETE_NAME,
            )
            shutil.rmtree(out_dir / STAGING_NAME, ignore_errors=True)
            raise
    return ranked
