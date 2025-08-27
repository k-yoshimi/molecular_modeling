# optimize_voronoi_overlap.py
"""
Two-stage orientation optimization for Voronoi-cell overlap (dictionary-cell = built from
the first COM in the molecule dictionary, using the dictionary's lattice vectors).

Flow:
  - Stage 0: Precompute system-level Voronoi cells once (fixed during optimization).
  - Precompute dictionary-cell (one per molecule name) from dictionary lattice + first COM.
  - For each site: define objective that rotates the dictionary-cell about its dictionary COM
    and places it at the site's COM, then compares to the fixed system Voronoi cell via IoU.

Notes:
    * Geometry parts rely on helpers in voronoi_geom.py.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Optional

import numpy as np
from ase import Atoms

from .voronoi_geom import (
    Polyhedron,
    compute_system_voronoi_cells,
    compute_voronoi_cells_from_coms,   # build Voronoi from given cell + COMs
    rotate_and_place_polyhedron,       # rotate (about dict COM) then place at site COM
    iou_polyhedra,
    wrap_angles_deg,
)
from .build_molecule_groups import build_molecule_groups


def optimize_voronoi_overlap(
    system: Atoms,
    groups: Optional[List[List[int]]] = None,
    x0_euler_deg: Optional[np.ndarray] = None,  # (N, 3) in degrees
    *,
    select_sites: Optional[Sequence[int]] = None,  # optimize only these site_ids
    grid_num: Tuple[int, int, int] = (13, 13, 13),
    angle_bounds: Tuple[float, float] = (-180.0, 180.0),
    log_dir: Optional[str | Path] = None,  # if None, no files are written
) -> Atoms:
    """
    Optimize orientation of dictionary-derived Voronoi cells to maximize overlap with
    system-level Voronoi cells (fixed). Dictionary-cell is defined per molecule name
    from the dictionary lattice vectors and the *first* center_of_mass entry.
    """

    # --- Auto-wire groups when omitted ---
    if groups is None:
        groups = build_molecule_groups(system)  # uses system.arrays['mol_id']
    num_sites = len(groups)

    # --- Auto-wire initial Euler angles when omitted ---
    if x0_euler_deg is None:
        if "site_angles_deg" in system.info:
            x0_euler_deg = np.asarray(system.info["site_angles_deg"], float)  # (N,3)
        else:
            x0_euler_deg = np.zeros((num_sites, 3), dtype=float)
    else:
        x0_euler_deg = np.asarray(x0_euler_deg, float)

    # --- Validate x0 shape ---
    if x0_euler_deg.shape != (num_sites, 3):
        raise ValueError(
            f"x0_euler_deg must have shape (N,3) with N=len(groups). "
            f"Got {x0_euler_deg.shape} and len(groups)={num_sites}."
        )

    # --- Select subset if requested ---
    if select_sites is not None:
        active_site_ids = np.asarray(select_sites, int).tolist()
        active_groups = [groups[i] for i in active_site_ids]
        active_x0 = np.asarray([x0_euler_deg[i] for i in active_site_ids], float)
    else:
        active_site_ids = list(range(num_sites))
        active_groups = groups
        active_x0 = x0_euler_deg

    # --- Ensure dictionary metadata is available ---
    if "site_mol_name" not in system.info or "dict_meta" not in system.info:
        raise ValueError(
            "system.info must contain 'site_mol_name' and 'dict_meta'. "
            "Store them in build_system_from_toml()."
        )
    site_mol_name = list(system.info["site_mol_name"])
    dict_meta = system.info["dict_meta"]  # {name: {'lattice_vectors': (3,3), 'center_of_mass': (K,3)}}

    # --- Stage 0: Precompute system-level Voronoi cells once (fixed) ---
    system_cells = compute_system_voronoi_cells(system, active_groups)
    if len(system_cells) != len(active_groups):
        raise RuntimeError(
            "compute_system_voronoi_cells returned mismatched count: "
            f"{len(system_cells)} vs {len(active_groups)}"
        )

    # --- Precompute dictionary-cell (per molecule name) using dict lattice + first COM only ---
    dict_cell_cache: Dict[str, tuple[Polyhedron, np.ndarray]] = {}
    for name in set(site_mol_name):
        meta = dict_meta[name]
        B_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)  # dictionary lattice
        coms_dict = np.asarray(meta["center_of_mass"], float)              # shape (K,3), K>=1
        if coms_dict.ndim != 2 or coms_dict.shape[1] != 3 or coms_dict.shape[0] < 1:
            raise ValueError(f"dict_meta['center_of_mass'] for '{name}' must be (K,3) with K>=1")
        dict_seed_cart = coms_dict[0]                                      # use *first* COM
        # Build exactly one Voronoi cell in the dictionary space
        dict_poly = compute_voronoi_cells_from_coms(B_dict, dict_seed_cart[None, :])[0]
        dict_cell_cache[name] = (dict_poly, dict_seed_cart)

    # --- Per-site optimization loop (start) ---
    for local_site, (site_id, idxs) in enumerate(zip(active_site_ids, active_groups)):
        # Rotation/placement centers
        center = system.get_positions()[idxs].mean(axis=0)  # site COM in system space
        target = system_cells[local_site]                   # fixed system Voronoi cell
        mol_name = site_mol_name[site_id]                   # which dictionary to use
        dict_poly_master, dict_seed_cart = dict_cell_cache[mol_name]

        # ----- Objective up to here -----
        def objective(angles_deg):
            """Return 1 - IoU (lower is better)."""
            angles_deg = wrap_angles_deg(np.asarray(angles_deg, float))
            dict_rot_placed = rotate_and_place_polyhedron(
                dict_poly_master, angles_deg, src_center=dict_seed_cart, dst_center=center
            )
            iou = iou_polyhedra(dict_rot_placed, target)
            return 1.0 - float(iou)

        # ---------- Stage 1: coarse grid search ----------
        amin, amax = float(angle_bounds[0]), float(angle_bounds[1])
        nx, ny, nz = (int(grid_num[0]), int(grid_num[1]), int(grid_num[2]))
        ax = np.linspace(amin, amax, nx)
        ay = np.linspace(amin, amax, ny)
        az = np.linspace(amin, amax, nz)

        # Start from given initial guess
        best_x = wrap_angles_deg(np.asarray(active_x0[local_site], float))
        best_fx = float(objective(best_x))

        for a in ax:
            for b in ay:
                for c in az:
                    x = np.array([a, b, c], dtype=float)
                    fx = float(objective(x))
                    if fx < best_fx:
                        best_fx, best_x = fx, x

        # ---------- Stage 2: Nelder–Mead refinement ----------
        from scipy.optimize import minimize  # local import to keep module deps light
        res = minimize(
            objective,
            x0=best_x,
            method="Nelder-Mead",
            options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 2000, "maxfev": 10000, "disp": False},
        )
        x_star = wrap_angles_deg(res.x)
        fx_star = float(res.fun)
        if fx_star < best_fx:
            best_x, best_fx = x_star, fx_star

        # ---------- Apply the best rotation to atoms (in-place) ----------
        from scipy.spatial.transform import Rotation as R
        rot = R.from_euler("xyz", best_x, degrees=True)
        pos = system.get_positions()
        sub = pos[idxs]
        pos[idxs] = rot.apply(sub - center) + center
        system.set_positions(pos)

        # Ensure collectors exist (fallback if not created before the loop)
        try:
            final_angles
        except NameError:
            final_angles = np.array(x0_euler_deg, float)
        try:
            best_iou_by_site
        except NameError:
            best_iou_by_site = {sid: np.nan for sid in range(num_sites)}

        # Record final angles and IoU for this site
        final_angles[site_id] = best_x
        best_iou_by_site[site_id] = 1.0 - best_fx

        # Optional: per-site log
        if log_dir is not None:
            outdir = Path(log_dir) / f"site_{site_id}"
            outdir.mkdir(parents=True, exist_ok=True)
            with open(outdir / "res.txt", "w", encoding="utf-8") as f:
                f.write(f"fx = {best_fx}\n")
                for lab, val in zip(["alpha", "beta", "gamma"], best_x):
                    f.write(f"{lab} = {val}\n")

    # --- Store arrays on system.info ---
    system.info["optimized_angles_deg"] = final_angles
    system.info["optimized_iou"] = np.array(
        [best_iou_by_site[sid] for sid in range(num_sites)], dtype=float
    )

    # --- Pretty print summary for *all* lattice sites (even if some were skipped) ---
    print("\n[Optimization summary] Euler angles (deg) per site")
    header = (
        f"{'site':>4}  {'alpha':>9}  {'beta':>9}  {'gamma':>9}  "
        f"{'d_alpha':>9}  {'d_beta':>9}  {'d_gamma':>9}  {'IoU*':>6}"
    )
    print(header)
    print("-" * len(header))
    for sid in range(num_sites):
        init = np.asarray(x0_euler_deg[sid], float)
        final = np.asarray(final_angles[sid], float)
        delta = wrap_angles_deg(final - init)
        iou_star = float(best_iou_by_site[sid])
        print(
            f"{sid:4d}  {final[0]:9.3f}  {final[1]:9.3f}  {final[2]:9.3f}  "
            f"{delta[0]:9.3f}  {delta[1]:9.3f}  {delta[2]:9.3f}  "
            f"{(iou_star if np.isfinite(iou_star) else float('nan')):6.3f}"
        )

    return system
