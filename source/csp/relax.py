"""Cell relaxation of crystal candidates with a tight-binding (or test) calculator, in parallel."""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import replace

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE
from molecular_modeling.csp.candidate import Candidate
from molecular_modeling.csp.config import CSPConfig
from molecular_modeling.csp.molecule import molecule_rmsd
from molecular_modeling.csp.screen import crystal_density
from molecular_modeling.molecule_unwrap import is_connected, unwrap_molecule


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


def relax_candidate(
    cand: Candidate, ref_mol: Atoms, method: str, fmax: float, steps: int, max_mol_rmsd: float
) -> Candidate:
    """Relax positions and cell, then fill energy, density and deformation metadata.

    Sets ``status`` to ``"relaxed"``, ``"rejected_deformed"`` (a molecule deviates more than
    ``max_mol_rmsd`` from the reference/its mirror image, or is broken) or ``"failed"``
    (exception text in ``meta["error"]``). Exceptions are never raised.
    The input candidate is mutated in place and returned.
    """
    try:
        atoms = cand.atoms.copy()
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
        cand.status = "relaxed" if dev <= max_mol_rmsd else "rejected_deformed"
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


def relax_many(cands: list[Candidate], ref_mol: Atoms, cfg: CSPConfig, method: str = "GFN2-xTB") -> list[Candidate]:
    """Relax candidates in parallel with ``cfg.n_procs`` processes; results keep the input order.

    ``meta["pyxtal"]`` is dropped before sending to workers (it may not be picklable).
    """
    if not cands:
        return []
    jobs = []
    for c in cands:
        meta = {k: v for k, v in c.meta.items() if k != "pyxtal"}
        jobs.append((replace(c, meta=meta), ref_mol, method, cfg.relax_fmax, cfg.relax_steps, cfg.max_mol_rmsd))
    results: list[Candidate | None] = [None] * len(jobs)
    pending = list(range(len(jobs)))
    workers = cfg.n_procs
    while pending:
        broken: list[int] = []
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as ex:
            futures = {k: ex.submit(_worker, jobs[k]) for k in pending}
            for k, fut in futures.items():
                try:
                    results[k] = fut.result()
                except BrokenProcessPool:
                    broken.append(k)
                except Exception as exc:  # noqa: BLE001 - e.g. result unpickling errors
                    results[k] = _failed(jobs[k][0], exc)
        if workers == 1 and broken:
            # with a single worker the first unresolved job in submission order is the one that crashed
            results[broken[0]] = _failed(jobs[broken[0]][0], RuntimeError("worker process died"))
            broken = broken[1:]
        pending = broken
        workers = 1  # re-run unresolved jobs one process at a time to isolate the crashing one
    return [r for r in results if r is not None]
