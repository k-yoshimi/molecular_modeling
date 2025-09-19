# pre_scale_by_dict_voronoi.py
from __future__ import annotations
from typing import List, Optional
import numpy as np
from ase import Atoms

from .build_molecule_groups import build_molecule_groups
from .voronoi_geom import single_set_cell_volume_cart

def pre_scale_by_dict_voronoi(
    system: Atoms,
    *,
    groups: Optional[List[List[int]]] = None,
    angles_deg: Optional[np.ndarray] = None,
    rng: Optional[int] = None,
    verbose: bool = True,
) -> Atoms:
    """
    Isotropically scale the *cell volume* to match the sum of dictionary Voronoi
    cell volumes, and move each group's COM accordingly around the origin.

    Behavior:
      - Cell is scaled by s = (V_target / V_current)^(1/3) (isotropic).
      - For each group, COM is mapped as COM' = s * COM (pivot = origin).
      - Intra-molecular displacements (p_i - COM) are preserved (atoms get a pure translation).
      - Intended to run BEFORE per-site rotation optimization.
    """
    # Groups
    if groups is None:
        groups = build_molecule_groups(system)
    N = len(groups)

    # RNG resolve
    if rng is None:
        rng_info = system.info.get("rng", {})
        rng = int(rng_info.get("pre_scale", rng_info.get("global", 1)))

    # Site meta
    if "site_mol_name" not in system.info or "dict_meta" not in system.info:
        raise ValueError(
            "system.info must contain 'site_mol_name' and 'dict_meta'. "
            "Ensure build_system_from_toml() stored them."
        )
    site_mol_name: List[str] = list(system.info["site_mol_name"])
    dict_meta = system.info["dict_meta"]

    # Angles for dictionary placement (default zeros)
    if angles_deg is None:
        angles_deg = np.zeros((N, 3), dtype=float)
    else:
        angles_deg = np.asarray(angles_deg, float)
        if angles_deg.shape != (N, 3):
            raise ValueError(f"angles_deg must have shape (N,3) with N={N}, got {angles_deg.shape}")

    # Current cell & positions
    A_sys = system.get_cell().array.copy()
    V_current = float(abs(np.linalg.det(A_sys)))
    if not np.isfinite(V_current) or V_current <= 0.0:
        raise RuntimeError(f"Current system cell volume is invalid: {V_current}")

    pos = system.get_positions()
    coms_sys_cart_all = np.array([pos[idxs].mean(axis=0) for idxs in groups], dtype=float)

    # Target volume: sum of dictionary Voronoi volumes over sites
    V_target = 0.0
    for sid in range(N):
        mol_name = site_mol_name[sid]
        meta = dict_meta[mol_name]
        A_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)
        coms_dict = np.asarray(meta["center_of_mass"], float)
        if coms_dict.ndim != 2 or coms_dict.shape[1] != 3 or coms_dict.shape[0] < 1:
            raise ValueError(f"dict_meta['center_of_mass'] for '{mol_name}' must be (K,3) with K>=1")

        p_ref = coms_sys_cart_all[sid]
        V_d = single_set_cell_volume_cart(
            A_cell=A_dict,
            coms_cart=coms_dict,
            p_ref_cart=p_ref,
            which="dict",
            angles_deg=angles_deg[sid],
            rng=rng,
        )
        V_target += float(V_d)

    if not np.isfinite(V_target) or V_target <= 0.0:
        raise RuntimeError(f"Computed target dictionary volume is invalid: {V_target}")

    # Isotropic scale factor
    s = (V_target / V_current) ** (1.0 / 3.0)

    # 1) Scale the cell (atoms untouched)
    if abs(s - 1.0) > 1e-15:
        A_new = s * A_sys
        system.set_cell(A_new, scale_atoms=False)
    else:
        A_new = A_sys.copy()

    # 2) Move each group's atoms so that COM -> s * COM (pivot = origin)
    if N > 0 and abs(s - 1.0) > 0.0:
        pos2 = pos.copy()
        for sid, idxs in enumerate(groups):
            COM = pos[idxs].mean(axis=0)         # current COM (Å)
            COM_scaled = s * COM                 # pivot = origin
            delta = COM_scaled - COM             # pure translation
            if delta.any():
                pos2[idxs] = pos2[idxs] + delta
        system.set_positions(pos2)

    # Bookkeeping
    info = dict(system.info.get("pre_scale", {}))
    info.update({
        "rng_used": int(rng),
        "V_target_sum_dict": float(V_target),
        "V_prev": float(V_current),
        "scale_factor": float(s),
        "A_prev": A_sys.copy(),
        "A_new": system.get_cell().array.copy(),
        "angles_used_deg": angles_deg.copy(),
        "pivot": "origin",
        "com_scaled": True,
    })
    system.info["pre_scale"] = info

    if verbose:
        print("[pre_scale] dictionary-sum target volume matching (shape-preserving)")
        print(f"  rng used         = {rng:d}")
        print(f"  V_sys(before)    = {V_current:.9f} Å^3")
        print(f"  V_target(dictΣ)  = {V_target:.9f} Å^3")
        print(f"  scale factor s   = {s:.9f}")
        print(f"  V_sys(after)     = {abs(np.linalg.det(system.get_cell().array)):.9f} Å^3")
        print("  COM scaling      = enabled (pivot=origin)")

    return system
