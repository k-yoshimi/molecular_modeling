"""End-to-end crystal structure prediction: generate, screen, relax, deduplicate, write."""
from __future__ import annotations

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

from .candidate import Candidate
from .config import CSPConfig
from .dedupe import deduplicate
from .gen_random import default_z, generate_random
from .gen_seed import generate_from_seed
from .gen_voronoi import orient_candidate
from .molecule import load_molecule, molecule_volume
from .relax import relax_many
from .screen import prescreen

logger = logging.getLogger(__name__)

CSV_HEADER = (
    "rank,id,energy_per_mol_eV,rel_energy_kJmol,density_gcm3,space_group_in,space_group_out,origin,max_mol_rmsd,status"
)
EV_TO_KJMOL = 96.485
MAX_CIF = 50
SAFETY_MARGIN = 0.9


@contextmanager
def _run_log(path: Path) -> Iterator[None]:
    """Send ``molecular_modeling.csp`` log records to ``path`` and stderr for the duration of the run."""
    root = logging.getLogger("molecular_modeling.csp")
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    handlers: list[logging.Handler] = [logging.FileHandler(path, mode="w", encoding="utf-8"), logging.StreamHandler()]
    old_level = root.level
    for h in handlers:
        h.setFormatter(fmt)
        h.setLevel(logging.INFO)
        root.addHandler(h)
    root.setLevel(logging.INFO)
    try:
        yield
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
    """Generate all candidates (random, Voronoi-oriented and seed-based)."""
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
    for seed in cfg.seed_structures:
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
    """Relax ``cands``, truncating the list if the time budget would be exceeded."""
    if not cands:
        return []
    n_first = min(cfg.n_procs, len(cands))
    t0 = time.perf_counter()
    done = relax_many(cands[:n_first], mol, cfg, method)
    round_time = time.perf_counter() - t0  # wall time of one round of ``n_first`` parallel relaxations
    per_cand = round_time / n_first
    rest = cands[n_first:]
    if rest:
        remaining = cfg.time_budget_min * 60 * SAFETY_MARGIN - (time.perf_counter() - t_start)
        # Only whole rounds count: a partly filled last round still takes as long as a full one.
        fit = max(0, int(remaining // round_time)) * n_first if round_time > 0 else len(rest)
        if fit < len(rest):
            logger.info(
                "time budget: %.1f s per candidate, %.0f s left; relaxing %d of %d remaining candidates",
                per_cand, max(remaining, 0.0), fit, len(rest),
            )
            rest = rest[:fit]
    if rest:
        done = done + relax_many(rest, mol, cfg, method)
    return done


def _write_outputs(cfg: CSPConfig, ranked: list[Candidate]) -> None:
    out_dir = Path(cfg.out_dir)
    emin = ranked[0].energy_per_mol if ranked else 0.0
    lines = [CSV_HEADER]
    for rank, c in enumerate(ranked, start=1):
        assert c.energy_per_mol is not None and emin is not None
        rel = (c.energy_per_mol - emin) * EV_TO_KJMOL
        dens = "" if c.density is None else f"{c.density:.4f}"
        rmsd = c.meta.get("max_mol_rmsd")
        lines.append(
            f"{rank},{c.meta['id']},{c.energy_per_mol:.6f},{rel:.3f},{dens},{c.space_group},"
            f"{_spacegroup_out(c.atoms)},{c.origin},{'' if rmsd is None else f'{rmsd:.4f}'},{c.status}"
        )
        if rank <= MAX_CIF:
            write(out_dir / f"{rank:03d}_{c.meta['id']}.cif", c.atoms)
    (out_dir / "candidates.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    """
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with _run_log(out_dir / "run.log"):
        t_start = time.perf_counter()
        rng = np.random.default_rng(cfg.rng)
        mol = load_molecule(cfg.molecule)

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
            n_ok, n_rej, n_fail, time.perf_counter() - t,
        )

        t = time.perf_counter()
        unique = deduplicate([c for c in relaxed if c.status == "relaxed" and c.energy_per_mol is not None])
        ranked = sorted(unique, key=lambda c: c.energy_per_mol if c.energy_per_mol is not None else math.inf)
        ranked = [replace(c) for c in ranked]
        logger.info("dedupe: %d -> %d in %.1f s", n_ok, len(ranked), time.perf_counter() - t)

        _write_outputs(cfg, ranked)
        logger.info("wrote %d rows to %s; total %.1f s", len(ranked), out_dir, time.perf_counter() - t_start)
    return ranked
