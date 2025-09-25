# voronoi_geom.py
# Geometry helpers (Cartesian-first, union-seed Voronoi).
from __future__ import annotations
from typing import List, Sequence, Optional
import numpy as np
from scipy.spatial import Voronoi, HalfspaceIntersection, ConvexHull
import trimesh
from ase.atoms import Atoms
from scipy.spatial.transform import Rotation as R


# =========================
# Basic shapes & converters
# =========================

class Polyhedron:
    """Minimal placeholder for a polyhedron shape."""
    def __init__(self, verts: np.ndarray, faces: np.ndarray):
        self.verts = np.asarray(verts, float)  # (V, 3)
        self.faces = np.asarray(faces, int)    # (F, 3 or 4, ...)


def _build_unit_cell_box_mesh(cell: np.ndarray) -> trimesh.Trimesh:
    """Create a convex mesh of the *unit* cell parallelepiped defined by `cell` (3x3, rows=a,b,c)."""
    corners_frac = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=float)
    corners_cart = corners_frac @ cell  # (8,3)
    hull = trimesh.convex.convex_hull(corners_cart)
    return hull


def _trimesh_to_polyhedron(mesh: trimesh.Trimesh) -> "Polyhedron":
    """Convert a trimesh.Trimesh (triangulated) to Polyhedron."""
    verts = mesh.vertices.view(np.ndarray)
    faces = mesh.faces.view(np.ndarray)   # triangles (F,3)
    return Polyhedron(verts=verts, faces=faces)


def _triangulate_faces(faces: np.ndarray) -> np.ndarray:
    """Ensure faces are triangles; fan-triangulate n-gons."""
    faces = np.asarray(faces, int)
    if faces.ndim != 2 or faces.shape[0] == 0:
        return np.zeros((0, 3), dtype=int)
    if faces.shape[1] == 3:
        return faces
    tris = []
    for face in faces:
        f = np.asarray(face, int)
        for k in range(1, len(f) - 1):
            tris.append([f[0], f[k], f[k + 1]])
    return np.asarray(tris, int)


def polyhedron_to_trimesh(poly: Polyhedron) -> "trimesh.Trimesh":
    """Convert Polyhedron to a processed trimesh.Trimesh."""
    verts = np.asarray(poly.verts, float)
    faces = _triangulate_faces(np.asarray(poly.faces, int))
    return trimesh.Trimesh(vertices=verts, faces=faces, process=True)


# =======================================
# System Voronoi via periodic tiling (±1)
# =======================================

def compute_voronoi_cells_from_coms(cell: np.ndarray, coms_cart: np.ndarray) -> List["Polyhedron"]:
    """
    Build Voronoi cells for COM (Center of Mass) points under periodicity by tiling {-1,0,1}^3 and
    clipping to the unit cell defined by `cell` (3x3, rows=a,b,c).
    Returns a list of Polyhedron in the same order as `coms_cart`.
    """
    if coms_cart.size == 0:
        return []
    cell = np.asarray(cell, float).reshape(3, 3)
    inv_cell = np.linalg.inv(cell)

    coms_frac = coms_cart @ inv_cell
    offsets = np.array([[i, j, k] for i in (-1, 0, 1)
                                 for j in (-1, 0, 1)
                                 for k in (-1, 0, 1)], dtype=float)  # (27,3)
    points_frac = (coms_frac[None, :, :] + offsets[:, None, :]).reshape(-1, 3)
    points_cart = points_frac @ cell

    center_off_idx = int(np.where((offsets == 0).all(axis=1))[0][0])
    M = coms_frac.shape[0]
    center_indices = (center_off_idx * M) + np.arange(M)

    vor = Voronoi(points_cart)
    cell_box = _build_unit_cell_box_mesh(cell)

    cells: List[Polyhedron] = []
    for m, pt_idx in enumerate(center_indices):
        region_idx = vor.point_region[pt_idx]
        vert_ids = vor.regions[region_idx]
        if (vert_ids is None) or (len(vert_ids) == 0) or (-1 in vert_ids):
            # fallback: small cube around COM (robust against unbounded regions)
            c = coms_cart[m]
            span = 0.5
            verts = c + np.array([
                [+span, +span, +span], [+span, +span, -span], [+span, -span, +span], [+span, -span, -span],
                [-span, +span, +span], [-span, +span, -span], [-span, -span, +span], [-span, -span, -span],
            ])
            faces = np.array([
                [0, 1, 3], [0, 3, 2], [4, 5, 7], [4, 7, 6],
                [0, 1, 5], [0, 5, 4], [2, 3, 7], [2, 7, 6],
                [0, 2, 6], [0, 6, 4], [1, 3, 7], [1, 7, 5],
            ], dtype=int)
            cells.append(Polyhedron(verts, faces))
            continue

        region_verts = vor.vertices[np.asarray(vert_ids, dtype=int)]
        try:
            region_hull = trimesh.convex.convex_hull(region_verts)
        except Exception:
            eps = 1e-8
            region_hull = trimesh.convex.convex_hull(region_verts + eps * np.random.randn(*region_verts.shape))
        try:
            clipped = region_hull.intersection(cell_box, engine="scad")
            if clipped.is_empty:
                clipped = region_hull.intersection(cell_box)
        except Exception:
            clipped = region_hull
        cells.append(_trimesh_to_polyhedron(clipped))
    return cells


# ============================
# Rigid motions & small helpers
# ============================

def rotate_and_place_polyhedron(
    poly: Polyhedron,
    euler_deg: Sequence[float],
    src_center: np.ndarray,
    dst_center: np.ndarray,
) -> Polyhedron:
    """Rotate poly about `src_center` (dictionary COM), then translate to `dst_center` (site COM)."""
    rot = R.from_euler('xyz', np.asarray(euler_deg, float), degrees=True)
    rotated = rot.apply(poly.verts - src_center) + dst_center
    return Polyhedron(rotated, poly.faces.copy())


def compute_system_voronoi_cells(system: Atoms, groups: List[List[int]]) -> List["Polyhedron"]:
    if len(groups) == 0:
        return []
    cell = system.get_cell().array
    coms_cart = np.array([system.positions[idxs].mean(axis=0) for idxs in groups], dtype=float)
    return compute_voronoi_cells_from_coms(cell, coms_cart)


def rotate_polyhedron(poly: Polyhedron, euler_deg: Sequence[float], center: np.ndarray) -> Polyhedron:
    """Rotate polyhedron (XYZ Euler degrees) about a given center."""
    rot = R.from_euler('xyz', np.asarray(euler_deg, float), degrees=True)
    rotated = rot.apply(poly.verts - center) + center
    return Polyhedron(rotated, poly.faces.copy())


def wrap_angles_deg(x: np.ndarray) -> np.ndarray:
    """Wrap Euler angles to [-180, 180) for stability."""
    y = (np.asarray(x, float) + 180.0) % 360.0 - 180.0
    y[np.isclose(y, 180.0)] = -180.0
    return y


# ======================
# Generic volume helper
# ======================

def polyhedron_volume(poly: Polyhedron) -> float:
    """Compute the (absolute) volume of a closed triangular/polygonal mesh via tetrahedra sum."""
    v = np.asarray(poly.verts, float)
    f = np.asarray(poly.faces, int)
    if f.ndim != 2 or f.shape[0] == 0:
        return 0.0
    if f.shape[1] != 3:
        # fan triangulation
        tri_list = []
        for face in f:
            face = np.asarray(face, int)
            for k in range(1, len(face) - 1):
                tri_list.append([face[0], face[k], face[k + 1]])
        if not tri_list:
            return 0.0
        tri = np.asarray(tri_list, int)
    else:
        tri = f
    cross = np.cross(v[tri[:, 1]], v[tri[:, 2]])
    vol6 = np.einsum("ij,ij->i", v[tri[:, 0]], cross).sum()
    return float(abs(vol6) / 6.0)


# ==========================================
# Cartesian-only Voronoi union machinery
# ==========================================

def _lattice_translations_cart(A: np.ndarray, rng: int = 1) -> np.ndarray:
    """
    Build translation vectors for a +/- rng supercell using Cartesian lattice A.
    A: (3,3) with rows a,b,c in Å (Cartesian)
    Returns: ((2*rng+1)^3, 3)
    """
    A = np.asarray(A, float)
    a, b, c = A
    shifts = [i*a + j*b + k*c for i in range(-rng, rng+1)
                               for j in range(-rng, rng+1)
                               for k in range(-rng, rng+1)]
    return np.asarray(shifts, float)


def _unique_rows(x: np.ndarray, atol: float = 1e-8) -> np.ndarray:
    """Remove near-duplicate rows by rounding to a tolerance grid."""
    x = np.asarray(x, float)
    grid = np.round(x / atol).astype(np.int64)
    _, idx = np.unique(grid, axis=0, return_index=True)
    return x[np.sort(idx)]


def _voronoi_cell_vertices_halfspaces_cart(p_ref: np.ndarray, others: np.ndarray) -> np.ndarray:
    """
    Compute vertices of the Voronoi cell (Cartesian Å) for p_ref given 'others'
    via intersection of perpendicular-bisector halfspaces.
    """
    p = np.asarray(p_ref, float)
    Q = np.asarray(others, float)

    halfspaces = []
    for q in Q:
        n = q - p
        if float(n @ n) < 1e-20:
            continue  # skip duplicates
        c = 0.5 * float(n @ (p + q))  # plane offset
        a = n
        b = -c
        halfspaces.append(np.r_[a, b])

    if len(halfspaces) == 0:
        return np.empty((0, 3), float)

    halfspaces = np.asarray(halfspaces, float)
    interior = p + 1e-10 * np.array([0.1234, -0.5678, 0.91011], float)  # strict interior point
    hs = HalfspaceIntersection(halfspaces, interior)
    return hs.intersections  # (M,3)


def _build_system_seeds_cart(A_sys: np.ndarray, coms_sys_cart: np.ndarray, rng: int = 1) -> np.ndarray:
    """System seeds = all system COMs replicated by translations of A_sys over +/- rng."""
    T_sys = _lattice_translations_cart(A_sys, rng=rng)
    seeds_sys = (coms_sys_cart[None, :, :] + T_sys[:, None, :]).reshape(-1, 3)
    return seeds_sys


def _build_dict_seeds_cart(
    A_dict: np.ndarray,
    coms_dict_cart: np.ndarray,
    angles_deg: Sequence[float],
    p_ref_cart: np.ndarray,
    rng: int = 1,
    scale: float = 1.0,
) -> np.ndarray:
    """
    Dictionary seeds (Cartesian) with uniform scaling:
      1) recenter dict COMs so that dict[0] is origin,
      2) rotate by R(angles),
      3) uniformly scale by `scale` (applied to both basis COM offsets and lattice),
      4) translate so dict[0] coincides with p_ref_cart,
      5) replicate by translations of the rotated+scaled lattice.
    """
    coms_dict_cart = np.asarray(coms_dict_cart, float)
    p0 = coms_dict_cart[0].copy()
    base = coms_dict_cart - p0  # (K,3)

    Rmat = R.from_euler('xyz', np.asarray(angles_deg, float), degrees=True).as_matrix()
    base_rot = (Rmat @ base.T).T
    base_rot_scaled = scale * base_rot
    placed = base_rot_scaled + np.asarray(p_ref_cart, float)

    A_dict = np.asarray(A_dict, float)
    A_dict_rot = (Rmat @ A_dict.T).T
    A_dict_rot_scaled = scale * A_dict_rot

    T_dict = _lattice_translations_cart(A_dict_rot_scaled, rng=rng)
    seeds_dict = (placed[None, :, :] + T_dict[:, None, :]).reshape(-1, 3)
    return seeds_dict

def overlap_volume_union_cart(
    A_sys: np.ndarray,
    coms_sys_cart: np.ndarray,
    A_dict: np.ndarray,
    coms_dict_cart: np.ndarray,
    angles_deg: Sequence[float],
    p_ref_cart: np.ndarray,
    rng: int = 1,
    scale: float = 1.0,
) -> float:
    """Overlap volume for union of system seeds and rotated+scaled dict seeds."""
    seeds_sys  = _build_system_seeds_cart(A_sys,  np.asarray(coms_sys_cart, float),  rng=rng)
    seeds_dict = _build_dict_seeds_cart(A_dict, np.asarray(coms_dict_cart, float),
                                        angles_deg, np.asarray(p_ref_cart, float),
                                        rng=rng, scale=scale)
    seeds_all = _unique_rows(np.vstack([seeds_sys, seeds_dict]), atol=1e-8)
    d2 = np.sum((seeds_all - p_ref_cart)**2, axis=1)
    others = seeds_all[d2 > 1e-18]
    verts = _voronoi_cell_vertices_halfspaces_cart(p_ref_cart, others)
    if len(verts) < 4:
        return 0.0
    hull = ConvexHull(verts)
    return float(hull.volume)

def single_set_cell_volume_cart(
    A_cell: np.ndarray,
    coms_cart: np.ndarray,
    p_ref_cart: np.ndarray,
    which: str,
    angles_deg: Optional[Sequence[float]] = None,
    rng: int = 1,
    scale: float = 1.0,
) -> float:
    """
    Volume of Voronoi cell at p_ref_cart for a single seed set.
    which='system' -> use A_cell & coms_cart as-is (scale ignored)
    which='dict'   -> rotate, then uniformly scale by `scale`, then place at p_ref
    """
    if which == 'system':
        seeds = _build_system_seeds_cart(A_cell, np.asarray(coms_cart, float), rng=rng)
    elif which == 'dict':
        assert angles_deg is not None, "angles_deg is required when which='dict'"
        seeds = _build_dict_seeds_cart(A_cell, np.asarray(coms_cart, float),
                                       angles_deg, np.asarray(p_ref_cart, float),
                                       rng=rng, scale=scale)
    else:
        raise ValueError("which must be 'system' or 'dict'.")

    seeds = _unique_rows(seeds, atol=1e-8)
    d2 = np.sum((seeds - p_ref_cart)**2, axis=1)
    others = seeds[d2 > 1e-18]
    verts = _voronoi_cell_vertices_halfspaces_cart(p_ref_cart, others)
    if len(verts) < 4:
        return 0.0
    hull = ConvexHull(verts)
    return float(hull.volume)

def iou_via_union_cart(
    A_sys: np.ndarray,
    coms_sys_cart: np.ndarray,
    A_dict: np.ndarray,
    coms_dict_cart: np.ndarray,
    p_ref_cart: np.ndarray,
    angles_deg: Sequence[float],
    rng: int = 1,
    scale: float = 1.0,
) -> float:
    """IoU = overlap / (V_sys + V_dict - overlap) with rotated+scaled dict seeds."""
    overlap = overlap_volume_union_cart(
        A_sys, coms_sys_cart, A_dict, coms_dict_cart, angles_deg, p_ref_cart, rng=rng, scale=scale
    )
    vol_sys  = single_set_cell_volume_cart(A_sys,  coms_sys_cart,  p_ref_cart, 'system', rng=rng)
    # For a single-set Voronoi cell, rotation does not change the volume; only scale matters.
    vol_dict = single_set_cell_volume_cart(A_dict, coms_dict_cart, p_ref_cart, 'dict',
                                           angles_deg=angles_deg, rng=rng, scale=scale)
    denom = vol_sys + vol_dict - overlap
    if denom <= 0.0:
        return 0.0
    return float(overlap / denom)


# (Optional) per-site optimizer sketch kept for convenience
def optimize_angles_each_system_site(
    A_sys: np.ndarray,
    coms_sys_cart: np.ndarray,
    A_dict: np.ndarray,
    coms_dict_cart: np.ndarray,
    sites_sys_cart: np.ndarray,
    angles_init: Sequence[float] = (0.0, 0.0, 0.0),
    rng: int = 1,
    optimizer: str = "Powell"
) -> List[tuple]:
    """
    For each system lattice point p_ref, align dict[0] to p_ref and optimize angles to maximize IoU.
    Returns: list of (site_index, angles_opt, best_iou).
    """
    from scipy.optimize import minimize
    results = []
    x0 = np.asarray(angles_init, float)
    for idx, p_ref in enumerate(np.asarray(sites_sys_cart, float)):
        def obj(ang):
            return -iou_via_union_cart(A_sys, coms_sys_cart, A_dict, coms_dict_cart, p_ref, ang, rng=rng)
        res = minimize(obj, x0, method=optimizer)
        results.append((idx, res.x, float(-res.fun)))
    return results


# ===========================
# Volume tiling sanity checker
# ===========================

def check_system_voronoi_volume_tiling(
    system: Atoms,
    groups: Optional[List[List[int]]] = None,
    *,
    rel_tol: float = 1e-3,
    abs_tol: float = 1e-6,
    verbose: bool = True,
) -> dict:
    """Check whether the sum of system Voronoi cell volumes equals the supercell volume."""
    from .build_molecule_groups import build_molecule_groups  # local import to avoid cycles

    if groups is None:
        groups = build_molecule_groups(system)

    cells = compute_system_voronoi_cells(system, groups)
    vols = [polyhedron_volume(c) for c in cells]
    sum_vol = float(np.sum(vols))

    cell = system.get_cell().array
    cell_vol = float(abs(np.linalg.det(cell)))

    diff = sum_vol - cell_vol
    rel_err = (abs(diff) / cell_vol) if cell_vol > 0 else np.nan
    ok = (abs(diff) <= abs_tol) or (rel_err <= rel_tol)

    report = {
        "sum_voronoi_volume": sum_vol,
        "cell_volume": cell_vol,
        "abs_diff": float(diff),
        "rel_error": float(rel_err),
        "tolerances": {"abs_tol": abs_tol, "rel_tol": rel_tol},
        "ok": bool(ok),
        "site_volumes": vols,
    }

    if verbose:
        print("\n[Voronoi tiling volume check]")
        print(f"  Sum of Voronoi cell volumes for the system : {sum_vol:.6f} Å^3")
        print(f"  System volume                              : {cell_vol:.6f} Å^3")
        print(f"  Abs diff                                   : {diff:.6e} Å^3")
        print(f"  Rel error                                  : {rel_err:.6e}")
        print(f"  Result                                     : {'PASS' if ok else 'FAIL'} "
              f"(abs_tol={abs_tol:g}, rel_tol={rel_tol:g})\n")

    return report
