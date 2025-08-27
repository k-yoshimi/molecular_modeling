# voronoi_geom.py
# Geometry helpers.
from __future__ import annotations
from typing import List, Sequence
import numpy as np
from scipy.spatial import Voronoi
import trimesh
from ase.atoms import Atoms
from scipy.spatial.transform import Rotation as R

class Polyhedron:
    """Minimal placeholder for a polyhedron shape."""
    def __init__(self, verts: np.ndarray, faces: np.ndarray):
        self.verts = np.asarray(verts, float)  # (V, 3)
        self.faces = np.asarray(faces, int)    # (F, 3 or 4, ...)

def _build_unit_cell_box_mesh(cell: np.ndarray) -> trimesh.Trimesh:
    """Create a convex mesh of the *unit* cell parallelepiped defined by `cell` (3x3, rows=a,b,c)."""
    # corners in fractional coords {0,1}^3
    corners_frac = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=float)
    corners_cart = corners_frac @ cell  # (8,3)
    # convex hull of the 8 corners gives the parallelepiped as a triangular mesh
    hull = trimesh.convex.convex_hull(corners_cart)
    return hull

def _trimesh_to_polyhedron(mesh: trimesh.Trimesh) -> "Polyhedron":
    """Convert a trimesh.Trimesh (triangulated) to Polyhedron."""
    verts = mesh.vertices.view(np.ndarray)
    faces = mesh.faces.view(np.ndarray)   # triangles (F,3)
    return Polyhedron(verts=verts, faces=faces)

def compute_voronoi_cells_from_coms(cell: np.ndarray, coms_cart: np.ndarray) -> List["Polyhedron"]:
    """
    Build Voronoi cells for COM points under periodicity by tiling {-1,0,1}^3 and
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
            # fallback: small cube around COM
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

def rotate_and_place_polyhedron(
    poly: Polyhedron,
    euler_deg: Sequence[float],
    src_center: np.ndarray,
    dst_center: np.ndarray,
) -> Polyhedron:
    """Rotate poly about `src_center` (dictionary COM), then translate to `dst_center` (site COM)."""
    rot = R.from_euler('xyz', np.asarray(euler_deg, float), degrees=True)
    verts = poly.verts
    rotated = rot.apply(verts - src_center) + dst_center
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
    shifted = poly.verts - center
    rotated = rot.apply(shifted) + center
    return Polyhedron(rotated, poly.faces.copy())

def iou_polyhedra(a: Polyhedron, b: Polyhedron) -> float:
    """Intersection-over-Union of two polyhedra volumes.
    TODO: Implement exact boolean ops and volume() function (e.g., trimesh.boolean).
    """
    def diag(p: Polyhedron) -> float:
        d = p.verts.max(0) - p.verts.min(0)
        return float(np.linalg.norm(d))
    da, db = diag(a), diag(b)
    if da <= 0 or db <= 0:
        return 0.0
    return float(min(da, db) / max(da, db))

def wrap_angles_deg(x: np.ndarray) -> np.ndarray:
    """Wrap Euler angles to [-180, 180) for stability."""
    y = (np.asarray(x, float) + 180.0) % 360.0 - 180.0
    y[np.isclose(y, 180.0)] = -180.0
    return y
