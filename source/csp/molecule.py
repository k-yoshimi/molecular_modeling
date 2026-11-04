"""Molecule loading and rigid-body alignment with atom matching."""

from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from ase.io import read
from molecular_modeling.build_system_from_toml import expand_symbols
from molecular_modeling.molecule_unwrap import is_connected
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation

N_RANDOM_STARTS = 16


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


#: Atoms of the input molecule closer than this (Å) are duplicates or overlapping.
MIN_ATOM_DISTANCE = 0.5

_TOKEN = re.compile(r"(?:[A-Z][a-z]?\d*)+")


def _toml_atoms(path: Path) -> Atoms:
    try:
        with open(path, "rb") as f:
            d = tomllib.load(f)
    except (tomllib.TOMLDecodeError, UnicodeDecodeError) as exc:
        raise ValueError(f"{path}: invalid TOML: {exc}") from None
    for key in ("symbols", "positions"):
        if key not in d:
            raise ValueError(f"{path}: '{key}' is required")
    tokens = d["symbols"]
    if not isinstance(tokens, list) or not all(isinstance(t, str) for t in tokens):
        raise ValueError(f"{path}: 'symbols' must be a list of element symbols, got {tokens!r}")
    for t in tokens:
        if not _TOKEN.fullmatch(t):
            raise ValueError(f"{path}: symbol {t!r} is not an element symbol (write e.g. 'C', 'Cl' or 'C6H6')")
    symbols = expand_symbols(tokens)
    unknown = sorted({s for s in symbols if s not in atomic_numbers})
    if unknown:
        raise ValueError(f"{path}: unknown element(s) {unknown}")
    try:
        pos = np.asarray(d["positions"], dtype=float)
    except (TypeError, ValueError):
        raise ValueError(f"{path}: 'positions' must be a list of [x, y, z] coordinates") from None
    if pos.size == 0:
        pos = pos.reshape(0, 3)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"{path}: 'positions' must be a list of [x, y, z] coordinates")
    if len(symbols) != len(pos):
        raise ValueError(f"{path}: {len(symbols)} symbols but {len(pos)} positions")
    if d.get("use_fractional", False):
        if "lattice_vectors" not in d:
            raise ValueError(f"{path}: use_fractional = true needs 'lattice_vectors'")
        lat = np.asarray(d["lattice_vectors"], dtype=float)
        if lat.shape != (3, 3):
            raise ValueError(f"{path}: 'lattice_vectors' must be three [x, y, z] vectors")
        pos = pos @ lat
    return Atoms(symbols=symbols, positions=pos)


def _xyz_atoms(path: Path) -> Atoms:
    try:
        # extxyz reads plain XYZ too and follows a Properties= header (pos may not be the first columns)
        frames = read(str(path), index=":", format="extxyz")
    except Exception as exc:  # ASE raises many kinds of errors for malformed files
        raise ValueError(f"{path}: cannot be read as XYZ: {type(exc).__name__}: {exc}") from None
    frames = frames if isinstance(frames, list) else [frames]
    if len(frames) != 1:
        raise ValueError(f"{path}: holds {len(frames)} frames (structures); give a file with exactly one molecule")
    return frames[0]


def load_molecule(path: str | Path) -> Atoms:
    """Load a molecule from ``.toml`` (``data/*.toml`` format) or ``.xyz`` (exactly one frame).

    The simple (unweighted) centroid is moved to the origin.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file cannot be read or is unsupported, a required key is missing, symbols are not element
        symbols, symbols and positions differ in number, the molecule is empty, two atoms are closer than
        :data:`MIN_ATOM_DISTANCE`, or the molecule is fragmented. The message names the file.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"molecule file not found: {path}")
    suffix = path.suffix.lower()
    if suffix == ".toml":
        atoms = _toml_atoms(path)
    elif suffix == ".xyz":
        atoms = _xyz_atoms(path)
    else:
        raise ValueError(f"unsupported molecule file type: {path} (use .toml or .xyz)")
    if len(atoms) == 0:
        raise ValueError(f"{path}: the molecule has no atoms")
    pos = atoms.get_positions()
    if not np.all(np.isfinite(pos)):
        raise ValueError(f"{path}: positions must be finite numbers")
    if len(atoms) > 1:
        dist = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2)
        np.fill_diagonal(dist, np.inf)
        i, j = np.unravel_index(np.argmin(dist), dist.shape)
        if dist[i, j] < MIN_ATOM_DISTANCE:
            raise ValueError(
                f"{path}: atoms {min(i, j)} and {max(i, j)} are only {dist[i, j]:.3f} Å apart "
                "(duplicate or overlapping atoms)"
            )
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
    """Four principal-axis alignments plus ``N_RANDOM_STARTS`` deterministic random rotations."""
    _, vm = np.linalg.eigh(np.cov((mobile - mobile.mean(0)).T, bias=True))
    _, vt = np.linalg.eigh(np.cov((target - target.mean(0)).T, bias=True))
    for v in (vm, vt):  # make the axis frames right-handed
        if np.linalg.det(v) < 0:
            v[:, 2] *= -1
    rots = [vt @ np.diag(s) @ vm.T for s in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]
    # Principal axes are arbitrary for degenerate inertia (cubes, spheres, ...): add fixed-seed random starts.
    rots.extend(Rotation.random(N_RANDOM_STARTS, random_state=0).as_matrix())
    return rots


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


#: vdW radius (Å) used for elements that ASE has no value for.
FALLBACK_VDW_RADIUS = 2.0
#: Grid points per Å^3 for :func:`molecule_volume`.
VOLUME_POINTS_PER_A3 = 50.0


def molecule_volume(mol: Atoms) -> float:
    """Van der Waals volume (Å^3): union of vdW spheres, estimated on a regular grid.

    The grid has about 50 points per Å^3 and is placed deterministically relative to the
    molecule, so the result does not depend on any random state.
    """
    from ase.data import vdw_radii

    pos = mol.get_positions()
    radii = np.array([vdw_radii[z] for z in mol.get_atomic_numbers()], dtype=float)
    radii[~np.isfinite(radii)] = FALLBACK_VDW_RADIUS
    h = VOLUME_POINTS_PER_A3 ** (-1.0 / 3.0)
    lo = (pos - radii[:, None]).min(axis=0) - h
    hi = (pos + radii[:, None]).max(axis=0) + h
    axes = [lo[d] + h * (np.arange(int(np.ceil((hi[d] - lo[d]) / h))) + 0.5) for d in range(3)]
    inside = np.zeros([len(a) for a in axes], dtype=bool)
    for p, r in zip(pos, radii):
        sl = []
        for d in range(3):
            i0 = max(int(np.floor((p[d] - r - lo[d]) / h)), 0)
            i1 = min(int(np.ceil((p[d] + r - lo[d]) / h)) + 1, len(axes[d]))
            sl.append(slice(i0, i1))
        gx, gy, gz = np.meshgrid(axes[0][sl[0]], axes[1][sl[1]], axes[2][sl[2]], indexing="ij")
        d2 = (gx - p[0]) ** 2 + (gy - p[1]) ** 2 + (gz - p[2]) ** 2
        inside[sl[0], sl[1], sl[2]] |= d2 <= r * r
    return float(inside.sum() * h**3)
