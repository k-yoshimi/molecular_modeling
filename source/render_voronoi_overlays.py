# render_voronoi_overlays.py
"""
Render overlays of:
  (1) system single-set Voronoi cell,
  (2) dictionary single-set Voronoi cell (rotated & placed),
  (3) overlap cell obtained from ONE Voronoi built on the union of seeds
      (system seeds ∪ rotated-placed dictionary seeds) around p_ref.

IoU shown in the title is computed via iou_via_union_cart (union-seed method),
consistent with the new optimization pipeline.

Additionally, for each site this script prints and annotates:
  - system single-set cell volume
  - dictionary single-set cell volume (rotated & placed)
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
from ase.atoms import Atoms

# Headless backend for PNG export
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull

from .build_molecule_groups import build_molecule_groups
from .voronoi_geom import (
    Polyhedron,
    compute_system_voronoi_cells,       # system single-set cells (for visualization)
    compute_voronoi_cells_from_coms,    # dict single-set cell (built in dict lattice)
    rotate_and_place_polyhedron,        # rotate about dict[0] and translate to site center
    iou_via_union_cart,                 # IoU via union-seed Voronoi
    single_set_cell_volume_cart,        # <-- volumes consistent with union method
    # Internal helpers to visualize the union-overlap cell
    _build_system_seeds_cart,
    _build_dict_seeds_cart,
    _unique_rows,
    _voronoi_cell_vertices_halfspaces_cart,
)

# ----------------------------- helpers -----------------------------

def _poly_faces_to_verts(poly: Polyhedron) -> List[np.ndarray]:
    """Convert Polyhedron.faces + verts to a list of polygon vertex arrays [(m,3), ...]."""
    faces = []
    for face in np.asarray(poly.faces):
        faces.append(poly.verts[np.asarray(face, int)])
    return faces


def _plot_polyhedron(
    ax,
    poly: Polyhedron,
    *,
    facecolor="C0",
    edgecolor="k",
    alpha=0.35,
    lw=0.5,
    label: Optional[str] = None,
) -> None:
    """Plot a polyhedron to a 3D axis."""
    poly_verts = _poly_faces_to_verts(poly)
    coll = Poly3DCollection(
        poly_verts,
        alpha=alpha,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidths=lw,
        label=label,
    )
    ax.add_collection3d(coll)


def _set_axes_equal_from_points(ax, pts: np.ndarray) -> None:
    """Set 3D axes to equal scale based on a set of points."""
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    centers = (mins + maxs) * 0.5
    extents = (maxs - mins)
    r = float(extents.max()) * 0.55 if extents.max() > 0 else 1.0
    ax.set_xlim(centers[0] - r, centers[0] + r)
    ax.set_ylim(centers[1] - r, centers[1] + r)
    ax.set_zlim(centers[2] - r, centers[2] + r)
    try:
        ax.set_box_aspect([1, 1, 1])
    except Exception:
        pass


def _convex_hull_to_polyhedron(points: np.ndarray) -> Polyhedron:
    """Turn a set of vertices into a triangular Polyhedron via ConvexHull."""
    hull = ConvexHull(points)
    return Polyhedron(verts=hull.points.copy(), faces=hull.simplices.copy())


# ----------------------------- main API -----------------------------

def render_voronoi_overlays(
    system: Atoms,
    *,
    groups: Optional[List[List[int]]] = None,
    angles_deg: Optional[np.ndarray] = None,  # (N,3) final Euler angles in degrees
    sites: Optional[Sequence[int]] = None,    # None = all sites
    out_dir: str | Path = "./voronoi_pngs",
    dpi: int = 220,
    figsize: Tuple[float, float] = (6.0, 5.0),
    azim: float = -60.0,   # camera azimuth
    elev: float = 20.0,    # camera elevation
    rng: int = 1,          # +/- replication range for both lattices (used for union cell & volumes)
    show_overlap: bool = True,  # draw overlap (union) cell
) -> None:
    """
    Render overlays and print volumes per site.

    For each site (p_ref = site's COM (Center of Mass)), we:
      - draw system single-set cell (blue),
      - draw dictionary single-set cell rotated & placed (orange),
      - optionally draw overlap cell from union seeds (green),
      - compute IoU (union method),
      - compute and PRINT volumes: V_sys_single, V_dict_single (both via single_set_cell_volume_cart).

    The title of each PNG includes IoU and both volumes.
    """
    # --- groups & angles auto-wire ---
    if groups is None:
        groups = build_molecule_groups(system)
    N = len(groups)

    # Ensure metadata exists
    if "site_mol_name" not in system.info or "dict_meta" not in system.info:
        raise ValueError(
            "system.info must contain 'site_mol_name' and 'dict_meta'. "
            "Please store them in build_system_from_toml()."
        )
    site_mol_name: List[str] = list(system.info["site_mol_name"])
    dict_meta = system.info["dict_meta"]  # name -> {"lattice_vectors": (3,3), "center_of_mass": (K,3)}

    if angles_deg is None:
        if "optimized_angles_deg" in system.info:
            angles_deg = np.asarray(system.info["optimized_angles_deg"], float)
        else:
            angles_deg = np.zeros((N, 3), dtype=float)
    else:
        angles_deg = np.asarray(angles_deg, float)

    if angles_deg.shape != (N, 3):
        raise ValueError(f"angles_deg must have shape (N,3) with N=len(groups)={N}, got {angles_deg.shape}")

    if sites is None:
        sites = list(range(N))
    else:
        sites = list(map(int, sites))

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- system single-set Voronoi cells computed once (order = groups order) ---
    sys_cells_all = compute_system_voronoi_cells(system, groups)

    # --- Precompute dictionary single-set cell (per molecule name) from its lattice and FIRST COM only (for viz) ---
    dict_cell_cache: dict[str, tuple[Polyhedron, np.ndarray]] = {}
    for name in set(site_mol_name):
        meta = dict_meta[name]
        B_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)
        coms_dict = np.asarray(meta["center_of_mass"], float)
        if coms_dict.ndim != 2 or coms_dict.shape[1] != 3 or coms_dict.shape[0] < 1:
            raise ValueError(f"dict_meta['center_of_mass'] for '{name}' must be (K,3) with K>=1")
        dict_seed_cart = coms_dict[0]  # FIRST COM
        dict_poly = compute_voronoi_cells_from_coms(B_dict, dict_seed_cart[None, :])[0]
        dict_cell_cache[name] = (dict_poly, dict_seed_cart)

    # --- Common data for IoU/volumes ---
    A_sys = system.get_cell().array  # (3,3) rows [a;b;c] in Cartesian Å
    pos = system.get_positions()
    coms_sys_cart_all = np.array([pos[idxs].mean(axis=0) for idxs in groups], dtype=float)  # (N,3)

    # --- Render each selected site ---
    for sid in sites:
        idxs = groups[sid]
        center = system.positions[idxs].mean(axis=0)  # p_ref in system space (Å)

        # (1) system single-set cell (for drawing)
        sys_cell = sys_cells_all[sid]

        # Dictionary data for this site
        mol_name = site_mol_name[sid]
        meta = dict_meta[mol_name]
        A_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)
        coms_dict = np.asarray(meta["center_of_mass"], float)
        dict_poly_master, dict_seed_cart = dict_cell_cache[mol_name]

        # (2) dictionary single-set cell for viz: rotate & place geom
        dict_rot_placed = rotate_and_place_polyhedron(
            dict_poly_master,
            angles_deg[sid],
            src_center=dict_seed_cart,
            dst_center=center,
        )

        # IoU via union-seed method (consistent with optimizer)
        iou_val = iou_via_union_cart(
            A_sys=A_sys,
            coms_sys_cart=coms_sys_cart_all,
            A_dict=A_dict,
            coms_dict_cart=coms_dict,
            p_ref_cart=center,
            angles_deg=angles_deg[sid],
            rng=rng,
        )

        # ---- Volumes (single-set) via halfspace method (consistent with union logic) ----
        V_sys = single_set_cell_volume_cart(
            A_cell=A_sys, coms_cart=coms_sys_cart_all, p_ref_cart=center, which='system', rng=rng
        )
        V_dict = single_set_cell_volume_cart(
            A_cell=A_dict, coms_cart=coms_dict, p_ref_cart=center, which='dict',
            angles_deg=angles_deg[sid], rng=rng
        )

        # Optional: build overlap cell (union-seed Voronoi at p_ref) for visualization
        overlap_poly: Optional[Polyhedron] = None
        if show_overlap:
            seeds_sys = _build_system_seeds_cart(A_sys, coms_sys_cart_all, rng=rng)
            seeds_dict = _build_dict_seeds_cart(A_dict, coms_dict, angles_deg[sid], center, rng=rng)
            seeds_all = _unique_rows(np.vstack([seeds_sys, seeds_dict]), atol=1e-8)
            d2 = np.sum((seeds_all - center) ** 2, axis=1)
            others = seeds_all[d2 > 1e-18]
            verts = _voronoi_cell_vertices_halfspaces_cart(center, others)
            if verts.shape[0] >= 4:
                overlap_poly = _convex_hull_to_polyhedron(verts)

        # Points for scaling
        pts_list = [dict_rot_placed.verts, sys_cell.verts]
        if overlap_poly is not None:
            pts_list.append(overlap_poly.verts)
        all_pts = np.vstack(pts_list)

        # ---- Plot ----
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection="3d")
        ax.view_init(elev=elev, azim=azim)

        # (1) system cell (blue)
        _plot_polyhedron(ax, sys_cell, facecolor="C0", edgecolor="k", alpha=0.25, lw=0.4,
                         label="system Voronoi")
        # (2) dictionary single-set (rotated & placed) (orange)
        _plot_polyhedron(ax, dict_rot_placed, facecolor="C1", edgecolor="k", alpha=0.45, lw=0.4,
                         label="dictionary (rotated)")
        # (3) overlap cell (green)
        if overlap_poly is not None:
            _plot_polyhedron(ax, overlap_poly, facecolor="C2", edgecolor="k", alpha=0.50, lw=0.4,
                             label="overlap")

        _set_axes_equal_from_points(ax, all_pts)
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        ax.set_zlabel("z (Å)")
        ax.set_title(
            f"Site {sid} | Mol: {mol_name} | IoU(union)={iou_val:.3f}"
        )
        ax.legend(loc="upper right", fontsize=9, frameon=True)

        out_path = Path(out_dir) / f"overlay_site_{sid:03d}.png"
        plt.tight_layout()
        plt.savefig(out_path, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)

        # ---- Console print for quick verification ----
        print(f"[render] saved: {out_path}")
        print(
            f"[vol] site {sid:03d} | V_sys = {V_sys:.6f} Å^3 | "
            f"V_dict = {V_dict:.6f} Å^3 | IoU = {iou_val:.6f}"
        )
    print("")
