# render_voronoi_overlays.py
"""
Render overlays of (system-level Voronoi cell) vs (rotated dictionary Voronoi cell) as PNGs.

Dictionary-cell definition (updated):
  - Built PER MOLECULE NAME from the dictionary's lattice vectors and the
    FIRST entry in the dictionary's center_of_mass list.
  - For a given site, the dictionary-cell is rotated about that dictionary COM
    and then translated (placed) to the site's COM before comparing/plotting.

Dependencies:
  - numpy, matplotlib, ase, trimesh (indirectly via voronoi_geom), scipy (installed already)

Data required on `system.info` (populated in build_system_from_toml):
  - system.info["site_mol_name"] : list[str]  # site_id -> molecule name used
  - system.info["dict_meta"]     : dict[name] -> {"lattice_vectors": (3,3), "center_of_mass": (K,3)}
  - optionally system.info["optimized_angles_deg"] : (N,3)  # final Euler angles (deg)
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

from .build_molecule_groups import build_molecule_groups
from .voronoi_geom import (
    Polyhedron,
    compute_system_voronoi_cells,
    compute_voronoi_cells_from_coms,
    rotate_and_place_polyhedron,
    iou_polyhedra,
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
    # Optional: draw vertices
    # ax.scatter(poly.verts[:, 0], poly.verts[:, 1], poly.verts[:, 2], s=6, c=edgecolor)


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
    azim: float = -60.0,  # camera azimuth
    elev: float = 20.0,   # camera elevation
) -> None:
    """
    Render overlays of (system Voronoi cell) vs (dictionary Voronoi cell rotated & placed) as PNGs.

    Parameters
    ----------
    system : Atoms
        Target structure (positions reflect the current orientation).
    groups : list[list[int]] | None
        Per-site atom index groups. If None, inferred from system.arrays['mol_id'].
    angles_deg : (N,3) | None
        Euler angles for each site (deg). If None, uses system.info['optimized_angles_deg']
        if available; otherwise zeros.
    sites : sequence[int] | None
        Subset of site_ids to render. If None, renders all sites.
    out_dir : str | Path
        Output directory to write "overlay_site_XXX.png".
    dpi : int
        Output PNG DPI.
    figsize : (float, float)
        Matplotlib figure size in inches.
    azim, elev : float
        Matplotlib 3D view angles (deg).
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

    # --- system-level Voronoi cells computed once (order = groups order) ---
    sys_cells_all = compute_system_voronoi_cells(system, groups)

    # --- Precompute dictionary-cell (per molecule name) from its lattice and FIRST COM only ---
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

    # --- Render each selected site ---
    for sid in sites:
        idxs = groups[sid]
        center = system.positions[idxs].mean(axis=0)

        # System cell for this site (fixed)
        sys_cell = sys_cells_all[sid]

        # Dictionary-cell for this site's molecule name
        mol_name = site_mol_name[sid]
        dict_poly_master, dict_seed_cart = dict_cell_cache[mol_name]

        # Rotate about dictionary COM and place at site's COM using final Euler angles
        dict_rot_placed = rotate_and_place_polyhedron(
            dict_poly_master,
            angles_deg[sid],
            src_center=dict_seed_cart,
            dst_center=center,
        )

        # IoU for title
        iou_val = iou_polyhedra(dict_rot_placed, sys_cell)

        # Points for scaling
        all_pts = np.vstack([dict_rot_placed.verts, sys_cell.verts])

        # ---- Plot ----
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection="3d")
        ax.view_init(elev=elev, azim=azim)

        # Plot system cell and rotated dictionary cell
        _plot_polyhedron(ax, sys_cell, facecolor="C0", edgecolor="k", alpha=0.25, lw=0.4, label="system Voronoi")
        _plot_polyhedron(ax, dict_rot_placed, facecolor="C1", edgecolor="k", alpha=0.50, lw=0.4, label="dictionary (rot+placed)")

        _set_axes_equal_from_points(ax, all_pts)
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        ax.set_zlabel("z (Å)")
        ax.set_title(f"Site {sid}  |  Mol: {mol_name}  |  IoU = {iou_val:.3f}")
        ax.legend(loc="upper right", fontsize=9, frameon=True)

        out_path = out_dir / f"overlay_site_{sid:03d}.png"
        plt.tight_layout()
        plt.savefig(out_path, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)

        print(f"[render] saved: {out_path}")
    print("")
