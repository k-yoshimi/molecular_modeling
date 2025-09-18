from __future__ import annotations
from typing import Optional, Sequence, Tuple, Dict
import numpy as np
from ase.atoms import Atoms
from collections import Counter

from .build_molecule_groups import build_molecule_groups
from .voronoi_geom import single_set_cell_volume_cart

def pre_scale_system_to_sum_dict_voronoi(
    system: Atoms,
    groups: Optional[list[list[int]]] = None,
    *,
    rng: Optional[int] = None,  # <-- was int; now Optional to auto-pick from system.info
    invariant: bool = True,
    angles_deg: Optional[np.ndarray] = None,  # per-site (N,3) if invariant=False
    safety_clamp: Tuple[float, float] = (0.02, 50.0),
    wrap_positions: bool = True,
    verbose: bool = True,
) -> Atoms:
    """
    Uniformly scale the system CELL and SITE COMs, while preserving each atom's
    *relative* position from its site's COM (intramolecular geometry unchanged).

    Target: new system volume equals the sum of dictionary single-set Voronoi volumes.
    """
    # ---- Resolve rng from system.info if not specified ----
    if rng is None:
        rng_info = system.info.get("rng", {})
        rng = int(rng_info.get("pre_scale", rng_info.get("global", 1)))

    if "site_mol_name" not in system.info or "dict_meta" not in system.info:
        raise ValueError("system.info must contain 'site_mol_name' and 'dict_meta' (set in build_system_from_toml).")

    if groups is None:
        groups = build_molecule_groups(system)

    site_mol_name: list[str] = list(system.info["site_mol_name"])
    dict_meta: Dict[str, Dict[str, np.ndarray]] = system.info["dict_meta"]

    A_sys = system.get_cell().array
    V_sys = float(abs(np.linalg.det(A_sys)))
    if V_sys <= 0.0:
        raise ValueError("System cell volume is non-positive; cannot scale.")

    # ---- Build V_target (sum of dictionary Voronoi volumes) ----
    if invariant:
        counts = Counter(site_mol_name)
        per_type_vol: Dict[str, float] = {}
        V_target = 0.0
        for name, cnt in counts.items():
            meta = dict_meta[name]
            A_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)
            coms_dict = np.asarray(meta["center_of_mass"], float)
            if coms_dict.ndim != 2 or coms_dict.shape[1] != 3 or coms_dict.shape[0] < 1:
                raise ValueError(f"dict_meta['center_of_mass'] for '{name}' must be (K,3) with K>=1")
            p_ref = coms_dict[0]
            V_one = single_set_cell_volume_cart(
                A_cell=A_dict, coms_cart=coms_dict, p_ref_cart=p_ref,
                which='dict', angles_deg=(0.0, 0.0, 0.0), rng=rng
            )
            per_type_vol[name] = float(V_one)
            V_target += float(cnt) * float(V_one)
    else:
        if angles_deg is None:
            if "site_angles_deg" in system.info:
                angles_deg = np.asarray(system.info["site_angles_deg"], float)
            else:
                raise ValueError("angles_deg must be provided when invariant=False.")
        angles_deg = np.asarray(angles_deg, float)
        if angles_deg.shape != (len(groups), 3):
            raise ValueError(f"angles_deg must be (N,3) with N=len(groups)={len(groups)}")

        pos = system.get_positions()
        V_target = 0.0
        for sid, idxs in enumerate(groups):
            center = pos[idxs].mean(axis=0)
            name = site_mol_name[sid]
            meta = dict_meta[name]
            A_dict = np.asarray(meta["lattice_vectors"], float).reshape(3, 3)
            coms_dict = np.asarray(meta["center_of_mass"], float)
            if coms_dict.ndim != 2 or coms_dict.shape[1] != 3 or coms_dict.shape[0] < 1:
                raise ValueError(f"dict_meta['center_of_mass'] for '{name}' must be (K,3) with K>=1")
            V_one = single_set_cell_volume_cart(
                A_cell=A_dict, coms_cart=coms_dict, p_ref_cart=center,
                which='dict', angles_deg=angles_deg[sid], rng=rng
            )
            V_target += float(V_one)
        per_type_vol = {}

    if not np.isfinite(V_target) or V_target <= 0.0:
        raise ValueError(f"Computed target dictionary-volume sum is invalid: {V_target}")

    # ---- Compute uniform scale s ----
    s_cubed = V_target / V_sys
    s = float(np.cbrt(s_cubed))
    s_min, s_max = safety_clamp
    if s_min > 0.0 or s_max < float("inf"):
        s = float(np.clip(s, s_min, s_max))

    # ---- Preserve intramolecular geometry ----
    pos0 = system.get_positions()
    if pos0.size == 0:
        return system

    coms0 = np.array([pos0[idxs].mean(axis=0) for idxs in groups], dtype=float)
    rel = np.zeros_like(pos0)
    for sid, idxs in enumerate(groups):
        rel[idxs] = pos0[idxs] - coms0[sid]

    new_cell = s * A_sys
    system.set_cell(new_cell, scale_atoms=False)

    invA0 = np.linalg.inv(A_sys)
    coms_frac = coms0 @ invA0
    coms1 = coms_frac @ new_cell

    pos1 = pos0.copy()
    for sid, idxs in enumerate(groups):
        pos1[idxs] = coms1[sid] + rel[idxs]

    system.set_positions(pos1)
    if system.get_pbc().any():
        try:
            system.wrap(eps=1e-9)
        except Exception:
            pass

    info = {
        "mode": "invariant" if invariant else "per-site",
        "rng": int(rng),
        "V_sys_before": V_sys,
        "V_target_dict_sum": float(V_target),
        "scale_factor": s,
        "safety_clamp": list(safety_clamp),
        "preserve_intramolecular_shape": True,
        "per_type_dict_volume": per_type_vol,
    }
    system.info["pre_scale"] = info

    if verbose:
        V_after = float(abs(np.linalg.det(system.get_cell().array)))
        print("[pre_scale] dictionary-sum target volume matching (shape-preserving)")
        print(f"  rng used         = {rng}")
        print(f"  V_sys(before)    = {V_sys:.6f} Å^3")
        print(f"  V_target(dict)   = {V_target:.6f} Å^3")
        print(f"  s (uniform)      = {s:.6f}  -> new V_sys ≈ {V_after:.6f} Å^3")
        if invariant and per_type_vol:
            for k, v in sorted(per_type_vol.items()):
                print(f"  per-type: {k}: V_dict0 = {v:.6f} Å^3")

    return system
