"""Cell relaxation of crystal candidates with a tight-binding (or test) calculator, in parallel."""

from __future__ import annotations

import os
from concurrent.futures import Future, ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import replace

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.filters import FrechetCellFilter
from ase.neighborlist import natural_cutoffs, neighbor_list
from ase.optimize import FIRE
from molecular_modeling.csp.candidate import Candidate
from molecular_modeling.csp.config import CSPConfig
from molecular_modeling.csp.molecule import molecule_rmsd
from molecular_modeling.csp.screen import crystal_density, whole_molecules
from molecular_modeling.molecule_unwrap import bond_pairs, is_connected, unwrap_molecule


def make_calculator(method: str = "GFN2-xTB") -> Calculator:
    """Return an ASE calculator for ``method`` (``"GFN2-xTB"``, ``"GFN1-xTB"``, or ``"EMT"`` for tests).

    Raises
    ------
    ValueError
        If ``method`` is not supported.
    """
    if method in ("GFN2-xTB", "GFN1-xTB"):
        from tblite.ase import TBLite

        return TBLite(method=method, verbosity=0)
    if method == "EMT":
        from ase.calculators.emt import EMT

        return EMT()
    if method == "EMT-CRASH-ON-AU":  # test-only: kills the worker process when gold is present
        from ase.calculators.emt import EMT

        class _CrashOnAu(EMT):
            def calculate(self, atoms=None, properties=("energy",), system_changes=("positions",)):
                if atoms is not None and 79 in atoms.numbers:
                    os._exit(1)
                super().calculate(atoms, properties, system_changes)

        return _CrashOnAu()
    raise ValueError(f"unsupported method {method!r}; use 'GFN2-xTB', 'GFN1-xTB' or 'EMT'")


def _mirror(mol: Atoms) -> Atoms:
    out = mol.copy()
    out.positions = -out.positions
    return out


def max_deformation(atoms: Atoms, ref_mol: Atoms) -> float:
    """Largest molecular RMSD (Å, heavy atoms) over the molecules of a crystal.

    Each molecule is unwrapped across the periodic boundary and compared with both the
    reference and its mirror image (symmetry operations may invert the molecule); the smaller
    value counts. Returns ``inf`` if a molecule is not a single connected piece.
    """
    mirror = _mirror(ref_mol)
    mol_id = atoms.arrays["mol_id"]
    worst = 0.0
    for m in np.unique(mol_id):
        idx = np.nonzero(mol_id == m)[0]
        syms = [atoms.symbols[i] for i in idx]
        pos = unwrap_molecule(atoms.positions[idx], syms, np.asarray(atoms.cell), atoms.pbc)
        if not is_connected(pos, syms):
            return float("inf")
        mol = Atoms(syms, positions=pos)
        worst = max(worst, min(molecule_rmsd(ref_mol, mol), molecule_rmsd(mirror, mol)))
    return worst


def molecules_bonded_together(atoms: Atoms) -> bool:
    """True if any molecule is covalently bonded to another molecule or to a periodic image of itself.

    Every ``mol_id`` molecule is first made whole (unwrapped across the cell boundary). Bonds are then
    looked for with the same criterion as the molecule split (1.2 x natural covalent radii, periodic
    images included); a bond counts as intermolecular unless it joins two atoms of the same molecule
    at the same lattice shift. This catches molecules bonded to each other and molecules that have
    become infinite chains through their own images, which a per-molecule shape check cannot see.

    Raises
    ------
    ValueError
        If a molecule is not one connected piece.
    """
    whole = whole_molecules(atoms)
    i, j, shift = neighbor_list("ijS", whole, natural_cutoffs(whole, mult=1.2))
    mol = whole.arrays["mol_id"]
    return bool(np.any((mol[i] != mol[j]) | np.any(shift != 0, axis=1)))


def intramolecular_bonds(atoms: Atoms) -> dict[int, set[tuple[int, int]]]:
    """Bonded atom pairs of every ``mol_id`` molecule, after unwrapping it across the cell boundary.

    Pairs are indices into the molecule's own atoms (in the order they appear in ``atoms``) and use the
    same criterion as the molecule split (1.2 x natural covalent radii).
    """
    mol_id = atoms.arrays["mol_id"]
    cell = np.asarray(atoms.cell)
    out: dict[int, set[tuple[int, int]]] = {}
    for m in np.unique(mol_id):
        idx = np.nonzero(mol_id == m)[0]
        syms = [atoms.symbols[i] for i in idx]
        pos = unwrap_molecule(atoms.positions[idx], syms, cell, atoms.pbc)
        out[int(m)] = bond_pairs(pos, syms)
    return out


def _bond_changes(before: dict[int, set[tuple[int, int]]], after: dict[int, set[tuple[int, int]]]) -> tuple[int, int]:
    """Numbers of intramolecular bonds added and removed between two :func:`intramolecular_bonds` results."""
    added = sum(len(after.get(m, set()) - b) for m, b in before.items())
    removed = sum(len(b - after.get(m, set())) for m, b in before.items())
    return added, removed


def relax_candidate(
    cand: Candidate, ref_mol: Atoms, method: str, fmax: float, steps: int, max_mol_rmsd: float
) -> Candidate:
    """Relax positions and cell, then fill energy, density and deformation metadata.

    Sets ``status`` to ``"relaxed"``, ``"rejected_deformed"`` (a molecule deviates more than
    ``max_mol_rmsd`` from the reference/its mirror image, or is broken, or has gained or lost a bond
    within itself (e.g. a ring closed or a hydrogen moved to another atom), or molecules have bonded to
    each other or to their own periodic images; the reason is in ``meta["reject_reason"]``) or ``"failed"``
    (exception text in ``meta["error"]``). Exceptions are never raised.
    The input candidate is mutated in place and returned.
    """
    try:
        atoms = cand.atoms.copy()
        bonds_before = intramolecular_bonds(atoms)
        atoms.calc = make_calculator(method)
        opt = FIRE(FrechetCellFilter(atoms), logfile=None)
        converged = bool(opt.run(fmax=fmax, steps=steps))
        n_mol = len(np.unique(atoms.arrays["mol_id"]))
        energy = float(atoms.get_potential_energy()) / n_mol
        atoms.calc = None
        dev = max_deformation(atoms, ref_mol)
        cand.atoms = atoms
        cand.energy_per_mol = energy
        cand.density = crystal_density(atoms)
        cand.meta["max_mol_rmsd"] = dev
        cand.meta["converged"] = converged
        cand.meta.pop("reject_reason", None)
        try:
            bonded = molecules_bonded_together(atoms)
        except ValueError:  # a broken molecule; max_deformation reports it as infinite deformation
            bonded = False
        added, removed = _bond_changes(bonds_before, intramolecular_bonds(atoms))
        if bonded:
            cand.meta["reject_reason"] = "molecules bonded to each other or to their own periodic images"
        elif added or removed:
            cand.meta["reject_reason"] = f"intramolecular bonds changed ({added} added, {removed} removed)"
        elif not dev <= max_mol_rmsd:
            cand.meta["reject_reason"] = f"molecule deformed (RMSD {dev:.3f} > {max_mol_rmsd} Å)"
        cand.status = "rejected_deformed" if "reject_reason" in cand.meta else "relaxed"
    except Exception as exc:  # noqa: BLE001 - a single failure must not stop the run
        cand.status = "failed"
        cand.meta["error"] = f"{type(exc).__name__}: {exc}"
    return cand


def _init_worker() -> None:
    os.environ["OMP_NUM_THREADS"] = "1"


def _worker(args: tuple) -> Candidate:
    cand, ref_mol, method, fmax, steps, max_mol_rmsd = args
    os.environ["OMP_NUM_THREADS"] = "1"
    return relax_candidate(cand, ref_mol, method, fmax, steps, max_mol_rmsd)


def _failed(cand: Candidate, exc: BaseException) -> Candidate:
    cand.status = "failed"
    cand.meta["error"] = f"{type(exc).__name__}: {exc}"
    return cand


def _run_parallel(jobs: list[tuple], keys: list[int], workers: int, results: list[Candidate | None]) -> list[int]:
    """Run ``jobs[k]`` for ``k`` in ``keys`` on one pool; fill ``results`` and return unresolved keys.

    Unresolved are jobs whose result was lost because the pool broke (a worker died), whether they
    had been submitted or not; which of the submitted ones crashed is unknown at this point.
    """
    unresolved: list[int] = []
    futures: dict[int, Future] = {}
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as ex:
        for n, k in enumerate(keys):
            try:
                futures[k] = ex.submit(_worker, jobs[k])
            except BrokenProcessPool:  # a worker already died; the rest is never submitted
                unresolved.extend(keys[n:])
                break
        for k, fut in futures.items():
            try:
                results[k] = fut.result()
            except BrokenProcessPool:
                unresolved.append(k)
            except Exception as exc:  # noqa: BLE001 - e.g. result unpickling errors
                results[k] = _failed(jobs[k][0], exc)
    return sorted(unresolved)


def _run_isolated(jobs: list[tuple], keys: list[int], results: list[Candidate | None]) -> None:
    """Run jobs one at a time in a single worker process, so a crash is pinned on the job that caused it.

    A job is blamed only if it was running when its worker died. If a fresh pool cannot accept a job,
    that job fails with a message saying so, and the others are still tried.
    """
    queue = list(keys)
    ex: ProcessPoolExecutor | None = None
    fresh = False
    try:
        while queue:
            k = queue[0]
            if ex is None:
                ex = ProcessPoolExecutor(max_workers=1, initializer=_init_worker)
                fresh = True
            try:
                fut = ex.submit(_worker, jobs[k])
            except BrokenProcessPool as exc:
                ex.shutdown(wait=True)
                ex = None
                if fresh:  # a new pool that cannot take a job: give up on this job only
                    results[k] = _failed(jobs[k][0], RuntimeError(f"worker process could not be started: {exc}"))
                    queue.pop(0)
                continue
            fresh = False
            queue.pop(0)
            try:
                results[k] = fut.result()
            except BrokenProcessPool:
                results[k] = _failed(jobs[k][0], RuntimeError("worker process died"))
                ex.shutdown(wait=True)
                ex = None
            except Exception as exc:  # noqa: BLE001 - e.g. result unpickling errors
                results[k] = _failed(jobs[k][0], exc)
    finally:
        if ex is not None:
            ex.shutdown(wait=True)


def relax_many(cands: list[Candidate], ref_mol: Atoms, cfg: CSPConfig, method: str = "GFN2-xTB") -> list[Candidate]:
    """Relax candidates in parallel with ``cfg.n_procs`` processes; results keep the input order.

    ``meta["pyxtal"]`` is dropped before sending to workers (it may not be picklable). If a worker
    process dies, the jobs whose results were lost (submitted or not) are re-run one at a time in a
    single worker, so only the job that crashes its worker is marked ``"failed"``.
    """
    if not cands:
        return []
    jobs = []
    for c in cands:
        meta = {k: v for k, v in c.meta.items() if k != "pyxtal"}
        jobs.append((replace(c, meta=meta), ref_mol, method, cfg.relax_fmax, cfg.relax_steps, cfg.max_mol_rmsd))
    results: list[Candidate | None] = [None] * len(jobs)
    keys = list(range(len(jobs)))
    if cfg.n_procs > 1:
        keys = _run_parallel(jobs, keys, cfg.n_procs, results)
    _run_isolated(jobs, keys, results)
    return [r for r in results if r is not None]
