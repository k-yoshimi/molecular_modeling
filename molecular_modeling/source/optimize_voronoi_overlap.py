from __future__ import annotations

from pathlib import Path
from typing import List, Sequence, Tuple, Optional

import numpy as np
from ase import Atoms

from .voronoi_geom import (
    iou_via_union_cart,
    single_set_cell_volume_cart,
    wrap_angles_deg,
)
from .build_molecule_groups import build_molecule_groups

def optimize_voronoi_overlap(
    system: Atoms,
    groups: Optional[List[List[int]]] = None,
    x0_euler_deg: Optional[np.ndarray] = None,  # (N, 3) in degrees
    *,
    select_sites: Optional[Sequence[int]] = None,
    grid_num: Tuple[int, int, int] = (13, 13, 13),
    angle_bounds: Tuple[float, float] = (-180.0, 180.0),
    rng: Optional[int] = None,
    log_dir: Optional[str | Path] = None,
) -> Atoms:
    """
    Optimize per-site Euler angles to maximize Voronoi IoU with uniform scaling.

    Parameters
    ----------
    system : ase.Atoms
        Input supercell with `mol_id` array and dictionary metadata in `system.info`.
    groups : list[list[int]], optional
        Atom index groups per site. If None, computed via `build_molecule_groups`.
    x0_euler_deg : ndarray, optional
        Initial Euler angles (N,3) in degrees. Defaults to zeros or `system.info["site_angles_deg"]`.
    select_sites : sequence[int], optional
        Subset of site indices to optimize. Others remain unchanged.
    grid_num : tuple[int, int, int], optional
        Grid counts for coarse search over (alpha, beta, gamma).
    angle_bounds : tuple[float, float], optional
        Inclusive bounds for grid search angles in degrees.
    rng : int, optional
        Random seed for Voronoi sampling; falls back to `system.info["rng"]`.
    log_dir : str or Path, optional
        If provided, writes per-site optimization logs to the directory.

    Returns
    -------
    ase.Atoms
        The same `Atoms` with positions rotated/scaled in place and
        optimization results stored in `system.info` (`optimized_angles_deg`,
        `optimized_scales`, `optimized_iou`).

    Notes
    -----
    Per-site procedure:
      1. Compute system Voronoi cell volume ``V_sys`` at the site COM.
      2. Compute dictionary Voronoi cell volume ``V_dict`` (rotation-invariant).
      3. Set uniform scale ``s = (V_sys / V_dict)^(1/3)`` (or 1 if ``V_dict <= 0``).
      4. Grid search then Nelder–Mead to maximize IoU using rotated+scaled dictionary seeds.
      5. Apply rotation and uniform scaling about the site COM to atom coordinates.
    """
    # Resolve rng
    if rng is None:
        rng_info = system.info.get("rng", {})
        rng = int(rng_info.get("optimize", rng_info.get("global", 1)))

    if groups is None:
        groups = build_molecule_groups(system)
    num_sites = len(groups)

    if x0_euler_deg is None:
        if "site_angles_deg" in system.info:
            x0_euler_deg = np.asarray(system.info["site_angles_deg"], float)
        else:
            x0_euler_deg = np.zeros((num_sites, 3), dtype=float)
    else:
        x0_euler_deg = np.asarray(x0_euler_deg, float)
    if x0_euler_deg.shape != (num_sites, 3):
        raise ValueError(f"x0_euler_deg must have shape (N,3) with N=len(groups)={num_sites}.")

    if select_sites is not None:
        active_site_ids = np.asarray(select_sites, int).tolist()
        active_groups = [groups[i] for i in active_site_ids]
        active_x0 = np.asarray([x0_euler_deg[i] for i in active_site_ids], float)
    else:
        active_site_ids = list(range(num_sites))
        active_groups = groups
        active_x0 = x0_euler_deg

    if "site_mol_name" not in system.info or "dict_meta" not in system.info:
        raise ValueError("system.info must contain 'site_mol_name' and 'dict_meta'.")

    site_mol_name = list(system.info["site_mol_name"])
    dict_meta = system.info["dict_meta"]

    A_sys = system.get_cell().array
    pos = system.get_positions()
    coms_sys_cart_all = np.array([pos[idxs].mean(axis=0) for idxs in groups], dtype=float)

    final_angles = np.array(x0_euler_deg, float)
    best_iou_by_site = {sid: np.nan for sid in range(num_sites)}
    scales_by_site = np.ones((num_sites,), float)

    amin, amax = float(angle_bounds[0]), float(angle_bounds[1])
    nx, ny, nz = (int(grid_num[0]), int(grid_num[1]), int(grid_num[2]))
    ax = np.linspace(amin, amax, nx)
    ay = np.linspace(amin, amax, ny)
    az = np.linspace(amin, amax, nz)

    for local_site, (site_id, idxs) in enumerate(zip(active_site_ids, active_groups)):
        p_ref = pos[idxs].mean(axis=0)
        mol_name = site_mol_name[site_id]
        meta = dict_meta[mol_name]
        A_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)
        coms_dict = np.asarray(meta["center_of_mass"], float)
        ref_idx   = int(meta.get("voronoi_ref_id", 0)) 

        # ---- Step 1-3: scale s by matching volumes ----
        V_sys = single_set_cell_volume_cart(A_sys, coms_sys_cart_all, p_ref, which='system', rng=rng)
        zeros = np.zeros(3, float)
        V_dict = single_set_cell_volume_cart(A_dict, coms_dict, p_ref, which='dict',
                                             angles_deg=zeros, rng=rng, scale=1.0, ref_idx=ref_idx)
        s = float((V_sys / V_dict) ** (1.0 / 3.0)) if V_dict > 0.0 else 1.0
        scales_by_site[site_id] = s

        # ---- Step 4: optimize rotation with scaled dictionary seeds ----
        def objective(angles_deg: Sequence[float]) -> float:
            """
            Objective to minimize: ``1 - IoU`` between system and dictionary Voronoi cells.

            Parameters
            ----------
            angles_deg : sequence[float]
                Euler angles (XYZ, degrees).

            Returns
            -------
            float
                One minus the intersection-over-union.
            """
            ang = wrap_angles_deg(np.asarray(angles_deg, float))
            iou = iou_via_union_cart(
                A_sys=A_sys,
                coms_sys_cart=coms_sys_cart_all,
                A_dict=A_dict,
                coms_dict_cart=coms_dict,
                p_ref_cart=p_ref,
                angles_deg=ang,
                rng=rng,
                scale=s,
                ref_idx=ref_idx,
            )
            return 1.0 - float(iou)

        best_x = wrap_angles_deg(np.asarray(active_x0[local_site], float))
        best_fx = float(objective(best_x))

        # Coarse grid (optional but robust)
        for a in ax:
            for b in ay:
                for c in az:
                    x = np.array([a, b, c], dtype=float)
                    fx = float(objective(x))
                    if fx < best_fx:
                        best_fx, best_x = fx, x

        # Nelder–Mead refine
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

        # ---- Step 5: apply rotation + uniform scaling to site atoms about COM ----
        from scipy.spatial.transform import Rotation as R
        rot = R.from_euler("xyz", best_x, degrees=True)
        sub = pos[idxs]
        # Order doesn't matter for uniform scaling around the same center; use rotate-then-scale.
        pos[idxs] = p_ref + s * (rot.apply(sub - p_ref))
        system.set_positions(pos)

        final_angles[site_id] = best_x
        best_iou_by_site[site_id] = 1.0 - best_fx

        if log_dir is not None:
            outdir = Path(log_dir) / f"site_{site_id}"
            outdir.mkdir(parents=True, exist_ok=True)
            with open(outdir / "res.txt", "w", encoding="utf-8") as f:
                f.write(f"fx = {best_fx}\n")
                for lab, val in zip(["alpha", "beta", "gamma"], best_x):
                    f.write(f"{lab} = {val}\n")
                f.write(f"scale = {s}\n")

    system.info["optimized_angles_deg"] = final_angles
    system.info["optimized_scales"] = scales_by_site
    system.info["optimized_iou"] = np.array(
        [best_iou_by_site[sid] for sid in range(num_sites)], dtype=float
    )

    print("\n[Optimization summary] Euler angles (deg) + scale per site")
    header = (
        f"{'site':>4}  {'alpha':>9}  {'beta':>9}  {'gamma':>9}  "
        f"{'IoU*':>6}  {'scale':>8}  {'rng':>4}"
    )
    print(header)
    print("-" * len(header))
    for sid in range(num_sites):
        iou_star = float(best_iou_by_site[sid])
        s = float(scales_by_site[sid])
        a, b, c = final_angles[sid]
        print(f"{sid:4d}  {a:9.3f}  {b:9.3f}  {c:9.3f}  {iou_star:6.3f}  {s:8.5f}  {rng:4d}")

    return system
