# optimize_voronoi_overlap.py
"""
Orientation optimization for Voronoi-cell overlap via a single Cartesian Voronoi problem:
for each system site (reference point p_ref), build ONE Voronoi using the union of
  - all system COM (Center of Mass) seeds (replicated by A_sys, +/- rng)
  - all dictionary COM seeds after rotating A_dict and COMs, aligning dict[0] to p_ref,
    then replicating by the rotated lattice (+/- rng)
and evaluate IoU around p_ref.

This replaces the previous method that compared two precomputed meshes via boolean ops.

Flow:
  - Gather system COMs (for ALL groups; even if optimizing only a subset of sites).
  - For each target site:
      * p_ref = site COM (system space).
      * pick the molecule dictionary (by site_mol_name[site]).
      * Objective(angles) = 1 - IoU computed by iou_via_union_cart(...)
      * Coarse grid search -> Nelder–Mead refine.
      * Apply the best rotation to atoms in that site (in-place).
  - Record per-site best angles and IoU into system.info.

Notes:
    * All geometry is Cartesian (Å).
    * Requires voronoi_geom.py providing iou_via_union_cart and wrap_angles_deg.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Optional

import numpy as np
from ase import Atoms

from .voronoi_geom import (
    iou_via_union_cart,  # IoU via union-seed Voronoi at p_ref (Cartesian)
    wrap_angles_deg,     # wrap Euler angles into [-180,180)
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
    rng: int = 1,                                # +/- supercell replication range
    log_dir: Optional[str | Path] = None,        # if None, no files are written
) -> Atoms:
    """
    Optimize dictionary orientation per site to maximize Voronoi-overlap IoU,
    using a single Voronoi problem formed by the union of seed sets (Cartesian).

    Inputs:
        system        : ASE Atoms with system lattice (system.get_cell()) and positions.
        groups        : list of index lists; each group corresponds to one molecule/site.
        x0_euler_deg  : initial Euler angles per site, shape (N,3) in degrees. If None,
                        uses system.info["site_angles_deg"] or zeros.
        select_sites  : subset of site indices to optimize (others kept as-is).
        grid_num      : coarse grid resolution for angles (nx,ny,nz).
        angle_bounds  : per-axis bounds (deg) for grid search.
        rng           : +/- replication range for seeds in both lattices (>=1 recommended).
        log_dir       : directory to write per-site logs.

    Returns:
        system : ASE Atoms with per-site optimal rotation applied in-place.
                 Stores 'optimized_angles_deg' and 'optimized_iou' in system.info.
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

    # --- Select subset if requested (for optimization only) ---
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

    # --- Lattices (Cartesian, rows [a;b;c]) ---
    A_sys = system.get_cell().array  # (3,3)

    # --- COMs for ALL system groups (we need the full set as seeds) ---
    pos = system.get_positions()
    coms_sys_cart_all = np.array([pos[idxs].mean(axis=0) for idxs in groups], dtype=float)  # (N,3)

    # --- Prepare output collectors ---
    final_angles = np.array(x0_euler_deg, float)  # (N,3)
    best_iou_by_site = {sid: np.nan for sid in range(num_sites)}

    # --- Coarse grid spec ---
    amin, amax = float(angle_bounds[0]), float(angle_bounds[1])
    nx, ny, nz = (int(grid_num[0]), int(grid_num[1]), int(grid_num[2]))
    ax = np.linspace(amin, amax, nx)
    ay = np.linspace(amin, amax, ny)
    az = np.linspace(amin, amax, nz)

    # --- Per-site optimization loop ---
    for local_site, (site_id, idxs) in enumerate(zip(active_site_ids, active_groups)):
        # Reference point p_ref = site COM in system space (Å)
        p_ref = pos[idxs].mean(axis=0)

        # Select dictionary for this site
        mol_name = site_mol_name[site_id]
        meta = dict_meta[mol_name]
        A_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)
        coms_dict = np.asarray(meta["center_of_mass"], float)
        if coms_dict.ndim != 2 or coms_dict.shape[1] != 3 or coms_dict.shape[0] < 1:
            raise ValueError(f"dict_meta['center_of_mass'] for '{mol_name}' must be (K,3) with K>=1")

        # Objective: minimize (1 - IoU) using union-seed Voronoi at p_ref
        def objective(angles_deg: Sequence[float]) -> float:
            angles_deg = wrap_angles_deg(np.asarray(angles_deg, float))
            iou = iou_via_union_cart(
                A_sys=A_sys,
                coms_sys_cart=coms_sys_cart_all,
                A_dict=A_dict,
                coms_dict_cart=coms_dict,
                p_ref_cart=p_ref,
                angles_deg=angles_deg,
                rng=rng,
            )
            return 1.0 - float(iou)

        # Start from initial guess for this site
        best_x = wrap_angles_deg(np.asarray(active_x0[local_site], float))
        best_fx = float(objective(best_x))

        # ---------- Stage 1: coarse grid search ----------
        for a in ax:
            for b in ay:
                for c in az:
                    x = np.array([a, b, c], dtype=float)
                    fx = float(objective(x))
                    if fx < best_fx:
                        best_fx, best_x = fx, x

        # ---------- Stage 2: Nelder–Mead refinement ----------
        from scipy.optimize import minimize
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

        # ---------- Apply the best rotation to atoms in this site (in-place) ----------
        from scipy.spatial.transform import Rotation as R
        rot = R.from_euler("xyz", best_x, degrees=True)
        sub = pos[idxs]
        pos[idxs] = rot.apply(sub - p_ref) + p_ref
        system.set_positions(pos)  # update for next sites if geometry affects COMs (usually minor)

        # Record results
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
