"""Molecule loading and rigid-body alignment with atom matching."""
from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read
from molecular_modeling.molecule_unwrap import is_connected
from scipy.optimize import linear_sum_assignment


@dataclass
class AlignResult:
    """Result of :func:`match_and_align`.

    ``rotation @ mobile[perm[i]] + translation`` approximates ``target[i]``.
    ``rmsd`` is over heavy atoms only when ``heavy_only`` was set.
    """

    rotation: np.ndarray
    translation: np.ndarray
    perm: np.ndarray
    rmsd: float


def load_molecule(path: str | Path) -> Atoms:
    """Load a molecule from ``.toml`` (``data/*.toml`` format) or ``.xyz``.

    The simple (unweighted) centroid is moved to the origin.

    Raises
    ------
    ValueError
        If the molecule is fragmented or the file type is unsupported.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".toml":
        with open(path, "rb") as f:
            d = tomllib.load(f)
        pos = np.asarray(d["positions"], dtype=float)
        if d.get("use_fractional", False):
            pos = pos @ np.asarray(d["lattice_vectors"], dtype=float)
        atoms = Atoms(symbols=list(d["symbols"]), positions=pos)
    elif suffix == ".xyz":
        atoms = read(str(path), format="xyz")  # type: ignore[assignment]
    else:
        raise ValueError(f"unsupported molecule file type: {path}")
    atoms.set_pbc(False)
    atoms.positions -= atoms.positions.mean(axis=0)
    if not is_connected(atoms.get_positions(), atoms.get_chemical_symbols()):
        raise ValueError(f"molecule {path} is fragmented (bond graph is not connected)")
    return atoms


def _kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t minimising |R m + t - x| for paired rows."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, tc - rot @ mc


def _assign(mob_pos: np.ndarray, tgt_pos: np.ndarray, mob_sym: np.ndarray, tgt_sym: np.ndarray) -> np.ndarray:
    """Per-element optimal assignment; returns perm with perm[i] = mobile index for target i."""
    perm = np.empty(len(tgt_pos), dtype=int)
    for el in np.unique(tgt_sym):
        ti = np.nonzero(tgt_sym == el)[0]
        mi = np.nonzero(mob_sym == el)[0]
        cost = np.linalg.norm(tgt_pos[ti][:, None, :] - mob_pos[mi][None, :, :], axis=2) ** 2
        r, c = linear_sum_assignment(cost)
        perm[ti[r]] = mi[c]
    return perm


def _initial_rotations(mobile: np.ndarray, target: np.ndarray) -> list[np.ndarray]:
    """Four proper rotations aligning principal axes (sign combinations with det=+1)."""
    _, vm = np.linalg.eigh(np.cov((mobile - mobile.mean(0)).T, bias=True))
    _, vt = np.linalg.eigh(np.cov((target - target.mean(0)).T, bias=True))
    for v in (vm, vt):  # make the axis frames right-handed
        if np.linalg.det(v) < 0:
            v[:, 2] *= -1
    return [vt @ np.diag(s) @ vm.T for s in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]


def match_and_align(
    mobile: Atoms,
    target: Atoms,
    heavy_only: bool = True,
    element_map: dict[str, str] | None = None,
) -> AlignResult:
    """Find the rigid transform and atom permutation that best map ``mobile`` onto ``target``.

    Starts from four principal-axis alignments; for each, alternates per-element
    linear assignment and Kabsch superposition (up to 20 iterations) and keeps
    the lowest RMSD.

    Parameters
    ----------
    mobile, target : ase.Atoms
        Molecules with the same composition (after ``element_map``).
    heavy_only : bool
        Fit and report RMSD on non-hydrogen atoms; hydrogens are then matched
        under the final transform and included in ``perm``.
    element_map : dict, optional
        Elements treated as identical, e.g. ``{"S": "O"}`` renames S to O in
        both molecules before matching.

    Returns
    -------
    AlignResult

    Raises
    ------
    ValueError
        If the compositions differ.
    """
    emap = element_map or {}
    msym = np.array([emap.get(s, s) for s in mobile.get_chemical_symbols()])
    tsym = np.array([emap.get(s, s) for s in target.get_chemical_symbols()])
    if sorted(msym) != sorted(tsym):
        raise ValueError("mobile and target have different compositions")
    mpos = mobile.get_positions()
    tpos = target.get_positions()

    fit = (tsym != "H") if heavy_only and np.any(tsym != "H") else np.ones(len(tsym), dtype=bool)
    m_fit = np.nonzero(msym != "H")[0] if fit.sum() < len(fit) else np.arange(len(msym))
    t_fit = np.nonzero(fit)[0]

    best: tuple[float, np.ndarray, np.ndarray, np.ndarray] | None = None
    for r0 in _initial_rotations(mpos[m_fit], tpos[t_fit]):
        rot = r0
        trans = tpos[t_fit].mean(0) - rot @ mpos[m_fit].mean(0)
        prev = np.inf
        for _ in range(20):
            moved = mpos[m_fit] @ rot.T + trans
            sub = _assign(moved, tpos[t_fit], msym[m_fit], tsym[t_fit])
            rot, trans = _kabsch(mpos[m_fit][sub], tpos[t_fit])
            cur = float(np.sqrt(np.mean(np.sum((mpos[m_fit][sub] @ rot.T + trans - tpos[t_fit]) ** 2, axis=1))))
            if abs(prev - cur) < 1e-10:
                break
            prev = cur
        if best is None or cur < best[0]:
            best = (cur, rot, trans, sub)
    assert best is not None
    rmsd, rot, trans, sub = best

    perm = np.empty(len(tsym), dtype=int)
    perm[t_fit] = m_fit[sub]
    rest = np.setdiff1d(np.arange(len(tsym)), t_fit)
    if len(rest):
        rest_m = np.setdiff1d(np.arange(len(msym)), m_fit)
        moved = mpos[rest_m] @ rot.T + trans
        p = _assign(moved, tpos[rest], msym[rest_m], tsym[rest])
        perm[rest] = rest_m[p]
    return AlignResult(rotation=rot, translation=trans, perm=perm, rmsd=rmsd)


def molecule_rmsd(ref: Atoms, mol: Atoms) -> float:
    """RMSD (Å, heavy atoms) between two molecules after optimal matching and superposition."""
    return match_and_align(mol, ref).rmsd
