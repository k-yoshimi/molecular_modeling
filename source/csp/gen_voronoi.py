"""Voronoi-guided molecular orientation.

Each molecule's own shape (the convex hull of points on its atoms' van der Waals spheres) is
rotated to best fill its Voronoi cell, the cell built from the molecular centres of the crystal.
Only rotations are applied to atoms; the shape is scaled for the overlap comparison only.
"""

from __future__ import annotations

import copy
import logging
from typing import Any

import numpy as np
from ase import Atoms
from ase.geometry import minkowski_reduce
from molecular_modeling.voronoi_geom import _build_system_seeds_cart, _voronoi_cell_vertices_halfspaces_cart
from scipy.optimize import linprog, minimize
from scipy.spatial import ConvexHull, HalfspaceIntersection
from scipy.spatial.transform import Rotation

from .candidate import Candidate
from .molecule import FALLBACK_VDW_RADIUS

logger = logging.getLogger(__name__)

#: Nelder–Mead settings for the refinement after the Euler-angle grid search.
NM_OPTIONS = {"xatol": 1e-2, "fatol": 1e-4, "maxiter": 400, "maxfev": 800, "disp": False}


def _sphere_directions(n: int) -> np.ndarray:
    """``n`` deterministic unit vectors spread over the sphere, in antipodal pairs.

    A Fibonacci lattice of ``ceil(n / 2)`` points on the upper hemisphere plus their negations, cut to ``n``. For even
    ``n`` the set is centrosymmetric, so its mean is exactly zero.
    """
    m = (n + 1) // 2
    z = (np.arange(m) + 0.5) / m  # upper hemisphere, so negations never coincide
    phi = np.pi * (1.0 + 5.0**0.5) * np.arange(m)
    r = np.sqrt(1.0 - z**2)
    half = np.column_stack([r * np.cos(phi), r * np.sin(phi), z])
    return np.vstack([half, -half])[:n]


def molecule_shape_points(mol: Atoms, n_sphere: int = 26) -> np.ndarray:
    """Point cloud describing the molecule's shape.

    ``n_sphere`` points are placed on the van der Waals sphere of every atom; the molecule's
    shape is the convex hull of these points. Coordinates are relative to the molecule's
    simple (unweighted) centroid of atoms.

    Parameters
    ----------
    mol : ase.Atoms
        The molecule.
    n_sphere : int
        Points per atom.

    Returns
    -------
    ndarray
        ``(len(mol) * n_sphere, 3)`` points in Å.
    """
    from ase.data import vdw_radii

    pos = mol.get_positions()
    pos = pos - pos.mean(axis=0)
    radii = np.array([vdw_radii[z] for z in mol.get_atomic_numbers()], dtype=float)
    radii[~np.isfinite(radii)] = FALLBACK_VDW_RADIUS
    dirs = _sphere_directions(int(n_sphere))
    return (pos[:, None, :] + radii[:, None, None] * dirs[None, :, :]).reshape(-1, 3)


def site_voronoi_vertices(cell: np.ndarray, coms: np.ndarray, index: int, rng: int = 1) -> np.ndarray:
    """Vertices of the Voronoi cell of molecular centre ``index`` (not clipped to the unit cell).

    The lattice is first Minkowski-reduced and every centre is moved by a lattice vector to its image
    nearest to ``coms[index]`` (in the reduced basis); the seeds are these centres replicated over
    ``±rng`` translations of the reduced basis. Without the reduction a skewed but equivalent basis
    would leave near neighbours out of the ``±rng`` block and give a far too large cell.

    Returns
    -------
    ndarray
        ``(M, 3)`` vertices relative to ``coms[index]``.
    """
    cell = np.asarray(cell, float)
    coms = np.asarray(coms, float)
    rcell, _ = minkowski_reduce(cell)
    rcell = np.asarray(rcell, float)
    centre = coms[index]
    frac = (coms - centre) @ np.linalg.inv(rcell)
    near = centre + (frac - np.round(frac)) @ rcell
    near[index] = centre
    seeds = _build_system_seeds_cart(rcell, near, rng=rng)
    verts = _voronoi_cell_vertices_halfspaces_cart(centre, seeds)
    return np.asarray(verts, float) - centre


def _interior_point(halfspaces: np.ndarray) -> np.ndarray | None:
    """Chebyshev centre of ``{x : A x + b <= 0}``; None if the region has no interior."""
    a, b = halfspaces[:, :3], halfspaces[:, 3]
    norm = np.linalg.norm(a, axis=1)
    c = np.array([0.0, 0.0, 0.0, -1.0])
    res = linprog(c, A_ub=np.column_stack([a, norm]), b_ub=-b, bounds=[(None, None)] * 3 + [(0, None)])
    if not res.success or res.x[3] < 1e-9:
        return None
    return np.asarray(res.x[:3], float)


def _intersection_volume(halfspaces: np.ndarray, interior: np.ndarray | None = None) -> float:
    """Volume of ``{x : A x + b <= 0}``; 0.0 if it is empty, flat, or Qhull fails."""
    if interior is None or np.any(halfspaces[:, :3] @ interior + halfspaces[:, 3] >= -1e-12):
        interior = _interior_point(halfspaces)
        if interior is None:
            return 0.0
    try:
        hs = HalfspaceIntersection(halfspaces, interior)
        return float(ConvexHull(hs.intersections).volume)
    except Exception:  # Qhull errors for degenerate regions
        return 0.0


def convex_iou(a_pts: np.ndarray, b_pts: np.ndarray) -> float:
    """Intersection over union of the convex hulls of two point clouds.

    The intersection is the halfspace intersection of both hulls' facet equations, using the
    origin as interior point when it lies strictly inside both (otherwise the Chebyshev centre).
    Returns 0.0 if the hulls have no common interior.
    """
    ha, hb = ConvexHull(np.asarray(a_pts, float)), ConvexHull(np.asarray(b_pts, float))
    inter = _intersection_volume(np.vstack([ha.equations, hb.equations]), np.zeros(3))
    union = ha.volume + hb.volume - inter
    return float(inter / union) if union > 0 else 0.0


class _ShapeInCell:
    """Precomputed hulls for repeated IoU evaluation of a rotated, scaled shape inside a fixed cell.

    Rotating a hull rotates its facet normals and keeps its offsets, so neither hull is rebuilt.
    """

    def __init__(self, shape_pts: np.ndarray, cell_vertices: np.ndarray) -> None:
        shape = ConvexHull(np.asarray(shape_pts, float))
        cell = ConvexHull(np.asarray(cell_vertices, float))
        self.v_cell = float(cell.volume)
        self.scale = (self.v_cell / shape.volume) ** (1.0 / 3.0)
        self.v_shape = float(shape.volume) * self.scale**3
        self.shape_normals = shape.equations[:, :3]
        self.shape_offsets = shape.equations[:, 3] * self.scale
        self.cell_eq = cell.equations

    def iou(self, rot: Rotation) -> float:
        normals = rot.apply(self.shape_normals)
        hs = np.vstack([np.column_stack([normals, self.shape_offsets]), self.cell_eq])
        inter = _intersection_volume(hs, np.zeros(3))
        union = self.v_shape + self.v_cell - inter
        return float(inter / union) if union > 0 else 0.0


def best_orientation(
    shape_pts: np.ndarray, cell_vertices: np.ndarray, grid_num: tuple[int, int, int] = (9, 9, 9)
) -> tuple[Rotation, float]:
    """Rotation of a shape that maximises the IoU of its convex hull with a Voronoi cell.

    For the comparison only, the shape is scaled by ``(V_cell / V_shape)^(1/3)`` so both have
    the same volume. A grid over XYZ Euler angles in ``[-180, 180)`` is searched first, then
    Nelder–Mead refines from the best grid point (as in ``optimize_voronoi_overlap``). The
    search is deterministic.

    Parameters
    ----------
    shape_pts : ndarray
        Shape point cloud, centred on the rotation centre (e.g. :func:`molecule_shape_points`).
    cell_vertices : ndarray
        Cell vertices relative to the same centre (e.g. :func:`site_voronoi_vertices`).
    grid_num : tuple of int
        Grid points per Euler angle.

    Returns
    -------
    (Rotation, float)
        The rotation to apply to the shape (about the origin) and the IoU it reaches.
    """
    match = _ShapeInCell(shape_pts, cell_vertices)

    def objective(x: np.ndarray) -> float:
        return 1.0 - match.iou(Rotation.from_euler("xyz", x, degrees=True))

    axes = [np.linspace(-180.0, 180.0, int(n), endpoint=False) for n in grid_num]
    best_x, best_f = np.zeros(3), objective(np.zeros(3))
    for x in _euler_grid(*axes):
        f = objective(x)
        if f < best_f:
            best_x, best_f = x, f
    res = minimize(objective, x0=best_x, method="Nelder-Mead", options=NM_OPTIONS)
    if float(res.fun) < best_f:
        best_x, best_f = np.asarray(res.x, float), float(res.fun)
    return Rotation.from_euler("xyz", best_x, degrees=True), 1.0 - best_f


def _euler_grid(*axes: np.ndarray) -> np.ndarray:
    """All combinations of the given 1-D arrays as rows."""
    grids = np.meshgrid(*axes, indexing="ij")
    return np.column_stack([g.ravel() for g in grids])


def _plane_normal(pos: np.ndarray, symbols: list[str]) -> np.ndarray:
    """Unit vector along which the molecule is thinnest (the normal of a planar molecule).

    The principal axis with the smallest eigenvalue of the gyration tensor of the heavy atoms (all
    atoms if there are fewer than two heavy atoms). A 180° turn about it swaps the two ends of an
    elongated molecule and keeps its plane.
    """
    heavy = np.array([s != "H" for s in symbols])
    p = pos[heavy] if heavy.sum() >= 2 else pos
    p = p - p.mean(axis=0)
    _, vecs = np.linalg.eigh(p.T @ p)
    return vecs[:, 0]


def _site_representatives(cand: Candidate, n_mol: int) -> list[int]:
    """Molecule index whose orientation is optimised for each independent site.

    With a PyXtal object, sites are emitted in order, each contributing ``wp.multiplicity``
    molecules; the first molecule of each site is its representative. Otherwise every
    molecule is independent.
    """
    px = cand.meta.get("pyxtal")
    if px is None:
        return list(range(n_mol))
    reps, start = [], 0
    for site in px.mol_sites:
        reps.append(start)
        start += int(site.wp.multiplicity)
    if start != n_mol:
        raise ValueError(f"PyXtal sites give {start} molecules but the atoms have {n_mol}")
    return reps


def orient_candidate(cand: Candidate, grid_num: tuple[int, int, int] = (9, 9, 9)) -> list[Candidate]:
    """Orient the molecules of a skeleton so each one's shape best fills its Voronoi cell.

    For every independent molecule the rotation Δ about its centroid maximising the shape/cell
    IoU is found with :func:`best_orientation`. With ``meta["pyxtal"]`` (random skeletons) one
    representative molecule per Wyckoff site is optimised and the site's orientation becomes
    ``Δ @ old``, so all symmetry copies follow; otherwise (seed skeletons, P1) every molecule is
    rotated independently about its own centroid.

    Two candidates are returned: Δ, and Δ followed by a 180° rotation of each optimised
    molecule about its centroid and the axis normal to its plane (the principal axis of the
    heavy-atom gyration tensor with the smallest extent). That turn exchanges the two ends of an
    elongated molecule (head and tail, e.g. the O and S ends) and keeps the molecular plane, which
    the shape/cell IoU cannot tell apart. Both have ``origin = "voronoi"`` and ``meta["voronoi_iou"]`` (mean over optimised
    molecules), ``meta["voronoi_iou_per_site"]``, ``meta["voronoi_variant"]`` (``"best"`` or
    ``"flip"``) and ``meta["parent_origin"]``. The input candidate is not modified.
    """
    atoms = cand.atoms
    ids = np.asarray(atoms.arrays["mol_id"], int)
    n_mol = int(ids.max()) + 1
    pos = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    cell = atoms.get_cell().array
    coms = np.array([pos[ids == k].mean(axis=0) for k in range(n_mol)])

    reps = _site_representatives(cand, n_mol)
    rots: dict[str, list[Rotation]] = {"best": [], "flip": []}
    ious: dict[str, list[float]] = {"best": [], "flip": []}
    for k in reps:
        sel = ids == k
        rel = pos[sel] - coms[k]
        mol_syms = [s for s, m in zip(symbols, sel) if m]
        shape = molecule_shape_points(Atoms(mol_syms, positions=rel))
        verts = site_voronoi_vertices(cell, coms, k)
        rot, iou = best_orientation(shape, verts, grid_num)
        axis = _plane_normal(rot.apply(rel), mol_syms)
        flip = Rotation.from_rotvec(np.pi * axis) * rot
        rots["best"].append(rot)
        rots["flip"].append(flip)
        ious["best"].append(iou)
        ious["flip"].append(_ShapeInCell(shape, verts).iou(flip))

    out: list[Candidate] = []
    for variant in ("best", "flip"):
        meta: dict[str, Any] = {k: v for k, v in cand.meta.items() if k != "pyxtal"}
        if "pyxtal" in cand.meta:
            new_atoms = _apply_to_pyxtal(cand, rots[variant], meta, len(atoms) // n_mol, n_mol)
            if new_atoms is None:
                logger.warning("Voronoi orientation (%s) broke a molecule apart; variant skipped", variant)
                continue
        else:
            new_atoms = atoms.copy()
            new_pos = pos.copy()
            for k, rot in zip(reps, rots[variant]):
                sel = ids == k
                new_pos[sel] = coms[k] + rot.apply(pos[sel] - coms[k])
            new_atoms.set_positions(new_pos)
        meta.update(
            voronoi_iou=float(np.mean(ious[variant])),
            voronoi_iou_per_site=[float(x) for x in ious[variant]],
            voronoi_variant=variant,
            parent_origin=cand.origin,
        )
        out.append(Candidate(atoms=new_atoms, origin="voronoi", space_group=cand.space_group, meta=meta))
    return out


def _apply_to_pyxtal(
    cand: Candidate, rots: list[Rotation], meta: dict[str, Any], n_per_mol: int, n_mol: int
) -> Atoms | None:
    """Rotate each PyXtal site by its Δ on a copy of the crystal; store the copy in ``meta``."""
    from .gen_random import _to_atoms

    crystal = copy.deepcopy(cand.meta["pyxtal"])
    for site, rot in zip(crystal.mol_sites, rots):
        site.orientation.matrix = rot.as_matrix() @ site.orientation.matrix
    new_atoms = _to_atoms(crystal, n_per_mol, n_mol)
    meta["pyxtal"] = crystal
    return new_atoms


def generate_voronoi(skeletons: list[Candidate], grid_num: tuple[int, int, int] = (9, 9, 9)) -> list[Candidate]:
    """Apply :func:`orient_candidate` to every skeleton and concatenate the results."""
    out: list[Candidate] = []
    for skel in skeletons:
        out.extend(orient_candidate(skel, grid_num))
    return out
