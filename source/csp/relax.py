"""Cell relaxation of crystal candidates in parallel, with the calculator of the chosen method."""

from __future__ import annotations

import logging
import multiprocessing
import os
import sys
import time
from concurrent.futures import Future, ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import replace

import numpy as np
from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.neighborlist import natural_cutoffs, neighbor_list
from ase.optimize import FIRE
from molecular_modeling.csp.backends import (  # noqa: F401 - METHODS and make_calculator are re-exported
    METHODS,
    ML_METHODS,
    RelaxBackendError,
    exercise,
    get_calculator,
    make_calculator,
)
from molecular_modeling.csp.candidate import Candidate
from molecular_modeling.csp.config import CSPConfig
from molecular_modeling.csp.molecule import molecule_rmsd
from molecular_modeling.csp.screen import crystal_density, whole_molecules
from molecular_modeling.molecule_unwrap import bond_pairs, is_connected, unwrap_molecule

logger = logging.getLogger(__name__)


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
        if method in ML_METHODS:  # the machine-learned models need charge and spin: neutral singlet
            atoms.info.update(charge=0, spin=1)
        atoms.calc = get_calculator(method)
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


#: Seconds a warm-up job waits for the other workers of its pool to start.
WARMUP_TIMEOUT_S = 600.0

#: Barrier shared by the workers of a pool that is warming up (set by :func:`_init_worker`).
_BARRIER = None


def _init_worker(method: str, barrier=None) -> None:
    """Prepare a worker process: one thread, and the calculator of ``method`` built (its model loaded) once.

    For machine-learned methods the calculator also evaluates a small structure, so that the slow first
    call (the model compiles itself) happens here and not in the first relaxation round.
    """
    global _BARRIER
    _BARRIER = barrier
    os.environ["OMP_NUM_THREADS"] = "1"
    calc = get_calculator(method)
    torch = sys.modules.get("torch")
    if torch is not None:
        torch.set_num_threads(1)
    if method in ML_METHODS:
        exercise(calc)


def _warmup() -> int:
    """Warm-up job: wait until every worker of the pool is running (and so has loaded its model).

    Without the wait, a worker that started early could take every warm-up job while the others are
    still loading, and their loading time would fall into the first round.
    """
    if _BARRIER is not None:
        _BARRIER.wait(WARMUP_TIMEOUT_S)
    return os.getpid()


def _new_executor(workers: int, method: str, barrier=None) -> ProcessPoolExecutor:
    """Process pool for ``method``: workers are started with spawn (never fork) and run :func:`_init_worker`.

    ``barrier`` (a spawn-context ``Barrier`` for ``workers`` parties) is handed to every worker for the
    warm-up jobs of :class:`WorkerPool`.
    """
    return ProcessPoolExecutor(
        max_workers=workers,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_init_worker,
        initargs=(method,) if barrier is None else (method, barrier),
    )


class WorkerPool:
    """A process pool kept open across several :func:`relax_many` calls (one relaxation stage).

    :meth:`ensure` starts the pool when there is none (at first, and after a worker died) and waits
    until every worker has loaded its calculator; a broken pool is dropped by :meth:`mark_broken`.
    """

    def __init__(self, workers: int, method: str) -> None:
        self.workers = workers
        self.method = method
        self._ex: ProcessPoolExecutor | None = None
        #: Process ids that answered the last warm-up (one per worker).
        self.warm_pids: list[int] = []

    def ensure(self) -> float:
        """Start the pool if needed and return the seconds that took (0 if it was ready).

        Raises
        ------
        RelaxBackendError
            If the workers cannot start, e.g. because the calculator cannot be built.
        """
        if self._ex is not None:
            return 0.0
        t0 = time.perf_counter()
        ex = _new_executor(self.workers, self.method, multiprocessing.get_context("spawn").Barrier(self.workers))
        try:
            self.warm_pids = [fut.result() for fut in [ex.submit(_warmup) for _ in range(self.workers)]]
        except Exception as exc:  # noqa: BLE001 - BrokenProcessPool or an error from the initializer
            ex.shutdown(wait=True, cancel_futures=True)
            raise RelaxBackendError(
                f"relaxation workers could not start ({self.method}): {type(exc).__name__}: {exc}"
            ) from exc
        except BaseException:  # e.g. Ctrl-C while the models load: do not leave the workers running
            ex.shutdown(wait=False, cancel_futures=True)
            raise
        self._ex = ex
        dt = time.perf_counter() - t0
        logger.info("worker pool: %d workers ready in %.1f s (%s)", self.workers, dt, self.method)
        return dt

    @property
    def executor(self) -> ProcessPoolExecutor:
        """The running pool; :meth:`ensure` must have been called."""
        if self._ex is None:
            raise RuntimeError("WorkerPool.ensure() must be called first")
        return self._ex

    def mark_broken(self) -> None:
        """Shut a pool down whose worker died; the next :meth:`ensure` starts a new one."""
        self.close()

    def close(self) -> None:
        if self._ex is not None:
            self._ex.shutdown(wait=True, cancel_futures=True)
            self._ex = None

    def __enter__(self) -> WorkerPool:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


def _worker(args: tuple) -> Candidate:
    cand, ref_mol, method, fmax, steps, max_mol_rmsd = args
    os.environ["OMP_NUM_THREADS"] = "1"
    return relax_candidate(cand, ref_mol, method, fmax, steps, max_mol_rmsd)


def _failed(cand: Candidate, exc: BaseException) -> Candidate:
    cand.status = "failed"
    cand.meta["error"] = f"{type(exc).__name__}: {exc}"
    return cand


def _run_parallel(
    jobs: list[tuple],
    keys: list[int],
    workers: int,
    results: list[Candidate | None],
    ex: ProcessPoolExecutor | None = None,
) -> list[int]:
    """Run ``jobs[k]`` for ``k`` in ``keys`` on one pool; fill ``results`` and return unresolved keys.

    The pool is ``ex`` (left open) or, if ``None``, a new one that is closed at the end. Unresolved are
    jobs whose result was lost because the pool broke (a worker died), whether they had been submitted
    or not; which of the submitted ones crashed is unknown at this point.
    """
    unresolved: list[int] = []
    futures: dict[int, Future] = {}
    own = ex is None
    pool = _new_executor(workers, jobs[keys[0]][2]) if ex is None else ex
    try:
        for n, k in enumerate(keys):
            try:
                futures[k] = pool.submit(_worker, jobs[k])
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
    finally:
        if own:
            pool.shutdown(wait=True)
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
                ex = _new_executor(1, jobs[k][2])
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


def relax_many(
    cands: list[Candidate],
    ref_mol: Atoms,
    cfg: CSPConfig,
    method: str = "GFN2-xTB",
    pool: WorkerPool | None = None,
) -> list[Candidate]:
    """Relax candidates in parallel with ``cfg.n_procs`` processes; results keep the input order.

    With ``pool`` (already made ready by :meth:`WorkerPool.ensure`) the jobs run on its workers, whatever
    ``cfg.n_procs`` is, and the pool is left open; if one of its workers dies the pool is marked broken.
    Without it a pool is made for this call (only if ``cfg.n_procs > 1``) and closed again.

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
    if pool is not None:
        keys = _run_parallel(jobs, keys, pool.workers, results, pool.executor)
        if keys:
            pool.mark_broken()
    elif cfg.n_procs > 1:
        keys = _run_parallel(jobs, keys, cfg.n_procs, results)
    _run_isolated(jobs, keys, results)
    return [r for r in results if r is not None]
