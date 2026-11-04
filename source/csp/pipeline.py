"""End-to-end crystal structure prediction: generate, screen, relax, deduplicate, write."""

from __future__ import annotations

import csv
import logging
import math
import time
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write

from .candidate import MAX_CIF, Candidate
from .config import SEED_GENERATOR_WEIGHT, CSPConfig
from .dedupe import deduplicate
from .deps import require
from .gen_random import default_z, generate_random
from .gen_seed import generate_from_seed
from .gen_voronoi import orient_candidate
from .molecule import load_molecule, molecule_volume
from .relax import relax_many
from .screen import check_supported_elements, prescreen

logger = logging.getLogger(__name__)

CSV_HEADER = (
    "rank,id,energy_per_mol_eV,rel_energy_kJmol,density_gcm3,space_group_in,space_group_out,origin,max_mol_rmsd,status"
)
EV_TO_KJMOL = 96.485
SAFETY_MARGIN = 0.9


#: Packages of the ``csp`` extra needed by every run; ``tblite`` is needed in addition for GFN methods.
REQUIRED_PACKAGES = ("pyxtal", "pymatgen", "spglib")


def check_dependencies(method: str = "GFN2-xTB") -> None:
    """Raise ``ImportError`` naming every missing package of the ``csp`` extra needed for ``method``.

    For the GFN methods ``tblite.ase`` (which loads the compiled tblite library) must import and provide
    ``TBLite``; importing the root ``tblite`` package alone does not prove the calculator works.
    """
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


def _generate(cfg: CSPConfig, mol: Atoms, rng: np.random.Generator) -> list[Candidate]:
    """Generate all candidates (random, Voronoi-oriented and seed-based).

    ``random`` and ``voronoi`` weights split ``n_candidates_per_sg`` between those generators. The
    ``seed`` weight only switches the seed generator on (> 0, the default when absent) or off (0):
    a seed always yields one candidate per substitution choice.
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
            skeletons = generate_random(mol, sg, z, math.ceil(n_vor / 2), rng)
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
    for seed in cfg.seed_structures if seed_weight > 0 else []:
        seeded = generate_from_seed(mol, seed)
        logger.info("seed %s: %d candidates", seed.cif, len(seeded))
        out.extend(seeded)
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

    The first round always runs. Before each further round the remaining time is compared with the
    longest round so far; if it is shorter, the remaining candidates are dropped.
    """
    budget = cfg.time_budget_min * 60 * SAFETY_MARGIN
    n_round = max(cfg.n_procs, 1)
    done: list[Candidate] = []
    rest = list(cands)
    longest = 0.0  # wall time of the slowest round so far
    n_done_rounds = 0
    while rest:
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
        done = done + relax_many(batch, mol, cfg, method)
        longest = max(longest, time.perf_counter() - t0)
        n_done_rounds += 1
    return done


def _write_outputs(cfg: CSPConfig, ranked: list[Candidate]) -> None:
    """Write ``candidates.csv`` and CIFs of the top :data:`MAX_CIF` candidates.

    Ranked CIFs (``NNN_<id>.cif``) left in ``cfg.out_dir`` by an earlier run are removed first, so
    the directory only holds files of this run.
    """
    out_dir = Path(cfg.out_dir)
    for old in out_dir.glob("[0-9][0-9][0-9]_*.cif"):
        old.unlink()
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
            ]
        )
        if rank <= MAX_CIF:
            write(out_dir / f"{rank:03d}_{c.meta['id']}.cif", c.atoms)
    # csv.writer quotes fields that contain commas or quotes (e.g. from seed file names)
    with open(out_dir / "candidates.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(CSV_HEADER.split(","))
        writer.writerows(rows)


def run_csp(cfg: CSPConfig, method: str = "GFN2-xTB") -> list[Candidate]:
    """Run the full CSP pipeline and write ``candidates.csv``, CIF files and ``run.log``.

    Parameters
    ----------
    cfg : CSPConfig
        Run settings. Output goes to ``cfg.out_dir``.
    method : str
        Calculator for relaxation (``"EMT"`` is available for fast tests).

    Returns
    -------
    list of Candidate
        Deduplicated relaxed candidates, lowest energy first.

    Raises
    ------
    ImportError
        Before anything is written, if a package of the ``csp`` extra needed for ``method`` is missing.
    """
    check_dependencies(method)
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with _run_log(out_dir / "run.log") as file_log:
        # warnings from loading the config were logged before run.log existed: repeat them in the file only
        for msg in cfg.config_warnings:
            file_log.handle(logger.makeRecord(logger.name, logging.WARNING, __file__, 0, "config: %s", (msg,), None))
        t_start = time.perf_counter()
        rng = np.random.default_rng(cfg.rng)
        mol = load_molecule(cfg.molecule)
        check_supported_elements(mol)

        t = time.perf_counter()
        cands = _generate(cfg, mol, rng)
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

        t = time.perf_counter()
        unique = deduplicate([c for c in relaxed if c.status == "relaxed" and c.energy_per_mol is not None])
        ranked = sorted(unique, key=lambda c: c.energy_per_mol if c.energy_per_mol is not None else math.inf)
        ranked = [replace(c) for c in ranked]
        logger.info("dedupe: %d -> %d in %.1f s", n_ok, len(ranked), time.perf_counter() - t)

        _write_outputs(cfg, ranked)
        logger.info("wrote %d rows to %s; total %.1f s", len(ranked), out_dir, time.perf_counter() - t_start)
    return ranked
