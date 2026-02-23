from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
from ase.atoms import Atoms

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull

from .build_molecule_groups import build_molecule_groups
from .voronoi_geom import (
    Polyhedron,
    _build_dict_seeds_cart,
    _build_system_seeds_cart,
    _unique_rows,
    _voronoi_cell_vertices_halfspaces_cart,
    compute_system_voronoi_cells,
    compute_voronoi_cells_from_coms,
    iou_via_union_cart,
    rotate_and_place_polyhedron,
    single_set_cell_volume_cart,
)


def _poly_faces_to_verts(poly: Polyhedron) -> List[np.ndarray]:
    """
    Convert face index lists to vertex arrays for plotting.

    Parameters
    ----------
    poly : Polyhedron
        Polyhedron with `verts` and `faces`.

    Returns
    -------
    list[numpy.ndarray]
        Vertices per face suitable for `Poly3DCollection`.
    """
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
    """
    Add a Polyhedron patch to a Matplotlib 3D axis.

    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes
        Target 3D axes.
    poly : Polyhedron
        Polyhedron to render.
    facecolor : str, optional
        Face color.
    edgecolor : str, optional
        Edge color.
    alpha : float, optional
        Transparency in [0,1].
    lw : float, optional
        Edge linewidth.
    label : str, optional
        Legend label.
    """
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
    """
    Set axis limits to equal aspect based on a point cloud.

    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes
        Target 3D axes.
    pts : ndarray
        Point cloud (N,3).
    """
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
    """
    Convert a point cloud to a convex polyhedron via `ConvexHull`.

    Parameters
    ----------
    points : ndarray
        Input points (N,3).

    Returns
    -------
    Polyhedron
        Convex hull polyhedron.
    """
    hull = ConvexHull(points)
    return Polyhedron(verts=hull.points.copy(), faces=hull.simplices.copy())

def render_voronoi_overlays(
    system: Atoms,
    *,
    groups: Optional[List[List[int]]] = None,
    angles_deg: Optional[np.ndarray] = None,
    sites: Optional[Sequence[int]] = None,
    out_dir: str | Path = "./voronoi_pngs",
    dpi: int = 220,
    figsize: Tuple[float, float] = (6.0, 5.0),
    azim: float = -60.0,
    elev: float = 20.0,
    rng: Optional[int] = None,
    show_overlap: bool = True,
) -> None:
    """
    Render Voronoi overlays for selected sites and report volumes.

    Parameters
    ----------
    system : ase.Atoms
        System with `system.info` containing `site_mol_name`, `dict_meta`, and
        optional `optimized_angles_deg` / `optimized_scales`.
    groups : list[list[int]], optional
        Atom index groups per site. If None, computed with `build_molecule_groups`.
    angles_deg : ndarray, optional
        Euler angles (N,3) in degrees; defaults to optimized values or zeros.
    sites : sequence[int], optional
        Site indices to render; defaults to all.
    out_dir : str or Path, optional
        Output directory for PNGs.
    dpi : int, optional
        Figure DPI.
    figsize : tuple[float, float], optional
        Matplotlib figure size in inches.
    azim : float, optional
        Azimuth angle for 3D view.
    elev : float, optional
        Elevation angle for 3D view.
    rng : int, optional
        Random seed for Voronoi sampling; falls back to `system.info["rng"]`.
    show_overlap : bool, optional
        If True, plot overlap polyhedron in orange.
    """
    if rng is None:
        rng_info = system.info.get("rng", {})
        rng = int(rng_info.get("render", rng_info.get("global", 1)))

    if groups is None:
        groups = build_molecule_groups(system)
    N = len(groups)

    if "site_mol_name" not in system.info or "dict_meta" not in system.info:
        raise ValueError("system.info must contain 'site_mol_name' and 'dict_meta'.")
    site_mol_name: List[str] = list(system.info["site_mol_name"])
    dict_meta = system.info["dict_meta"]

    if angles_deg is None:
        if "optimized_angles_deg" in system.info:
            angles_deg = np.asarray(system.info["optimized_angles_deg"], float)
        else:
            angles_deg = np.zeros((N, 3), dtype=float)
    else:
        angles_deg = np.asarray(angles_deg, float)
    if angles_deg.shape != (N, 3):
        raise ValueError(f"angles_deg must have shape (N,3) with N={N}.")

    # Per-site scales (default 1.0)
    if "optimized_scales" in system.info:
        scales = np.asarray(system.info["optimized_scales"], float)
        if scales.shape != (N,):
            scales = np.ones((N,), float)
    else:
        scales = np.ones((N,), float)

    if sites is None:
        sites = list(range(N))
    else:
        sites = list(map(int, sites))

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sys_cells_all = compute_system_voronoi_cells(system, groups)

    # Prepare per-molecule dictionary base poly (scale-independent at this point)
    dict_cell_cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, int]] = {}
    for name in set(site_mol_name):
        meta = dict_meta[name]
        B_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)
        coms_dict = np.asarray(meta["center_of_mass"], float)
        ref_idx   = int(meta.get("voronoi_ref_id", 0))
        dict_seed_cart = coms_dict[ref_idx]
        dict_cell_cache[name] = (B_dict, dict_seed_cart, coms_dict, ref_idx)

    A_sys = system.get_cell().array
    pos = system.get_positions()
    coms_sys_cart_all = np.array([pos[idxs].mean(axis=0) for idxs in groups], dtype=float)

    for sid in sites:
        idxs = groups[sid]
        center = system.positions[idxs].mean(axis=0)
        sys_cell = sys_cells_all[sid]
        s = float(scales[sid])

        mol_name = site_mol_name[sid]
        B_dict, dict_seed_cart, coms_dict, ref_idx = dict_cell_cache[mol_name]
        A_dict = B_dict  # alias

        # Build a scaled dictionary single-cell poly at the origin COM,
        # then rotate and place it to the site center.
        # For the single seed, compute_voronoi_cells_from_coms uses lattice only; scale lattice by s.
        dict_poly_master = compute_voronoi_cells_from_coms(s * B_dict, dict_seed_cart[None, :])[0]
        dict_rot_placed = rotate_and_place_polyhedron(
            dict_poly_master,
            angles_deg[sid],
            src_center=dict_seed_cart,
            dst_center=center,
        )

        iou_val = iou_via_union_cart(
            A_sys=A_sys,
            coms_sys_cart=coms_sys_cart_all,
            A_dict=A_dict,
            coms_dict_cart=coms_dict,
            p_ref_cart=center,
            angles_deg=angles_deg[sid],
            rng=rng,
            scale=s,
            ref_idx=ref_idx,
        )

        V_sys = single_set_cell_volume_cart(
            A_cell=A_sys, coms_cart=coms_sys_cart_all, p_ref_cart=center, which='system', rng=rng
        )
        V_dict = single_set_cell_volume_cart(
            A_cell=A_dict, coms_cart=coms_dict, p_ref_cart=center, which='dict',
            angles_deg=angles_deg[sid], rng=rng, scale=s, ref_idx=ref_idx
        )

        if show_overlap:
            seeds_sys = _build_system_seeds_cart(A_sys, coms_sys_cart_all, rng=rng)
            seeds_dict = _build_dict_seeds_cart(A_dict, coms_dict, angles_deg[sid], center,
                                                rng=rng, scale=s, ref_idx=ref_idx)
            seeds_all = _unique_rows(np.vstack([seeds_sys, seeds_dict]), atol=1e-8)
            d2 = np.sum((seeds_all - center) ** 2, axis=1)
            others = seeds_all[d2 > 1e-18]
            verts = _voronoi_cell_vertices_halfspaces_cart(center, others)
            if verts.shape[0] >= 4:
                overlap_poly = _convex_hull_to_polyhedron(verts)

        pts_list = [dict_rot_placed.verts, sys_cell.verts]
        if overlap_poly is not None:
            pts_list.append(overlap_poly.verts)
        all_pts = np.vstack(pts_list)

        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection="3d")
        ax.view_init(elev=elev, azim=azim)  # type: ignore[attr-defined]

        _plot_polyhedron(ax, sys_cell, facecolor="C0", edgecolor="k", alpha=0.25, lw=0.4,
                         label="system Voronoi")
        _plot_polyhedron(ax, dict_rot_placed, facecolor="C1", edgecolor="k", alpha=0.45, lw=0.4,
                         label=f"dictionary (scaled s={s:.4f})")
        if overlap_poly is not None:
            _plot_polyhedron(ax, overlap_poly, facecolor="C2", edgecolor="k", alpha=0.50, lw=0.4,
                             label="overlap")

        _set_axes_equal_from_points(ax, all_pts)
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        ax.set_zlabel("z (Å)")  # type: ignore[attr-defined]
        ax.set_title(
            f"Site {sid} | Mol: {mol_name} | IoU={iou_val:.3f} | rng={rng}"
        )
        ax.legend(loc="upper right", fontsize=9, frameon=True)

        out_path = Path(out_dir) / f"overlay_site_{sid:03d}.png"
        plt.tight_layout()
        plt.savefig(out_path, bbox_inches="tight", pad_inches=0.3)
        plt.close(fig)

        print(f"[render] saved: {out_path}")
        print(
            f"[vol] site {sid:03d} | V_sys = {V_sys:.6f} Å^3 | "
            f"V_dict(scaled) = {V_dict:.6f} Å^3 | IoU = {iou_val:.6f} | rng = {rng} | s = {s:.6f}"
        )
    print("")
