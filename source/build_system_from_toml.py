import re
import numpy as np
import tomllib
from typing import Dict, List, Tuple, Optional
from scipy.spatial.transform import Rotation as R
from ase import Atoms
from itertools import product
from collections import Counter

def build_system_from_toml(*, verbose: bool = True) -> Atoms:
    """Build an ASE Atoms system from TOML inputs (supercell + substitutions)."""
    # Load input from TOML file
    with open("./input/input.toml", "rb") as f:
        input_data = tomllib.load(f)

    # --- NEW: wire rng settings ---
    # Global default
    rng_global = int(input_data.get("rng", 1))
    # Per-stage overrides (optional)
    rng_pre_scale = int(input_data.get("pre_scale_rng", rng_global))
    rng_optimize  = int(input_data.get("optimize_rng",  rng_global))
    rng_render    = int(input_data.get("render_rng",    rng_global))

    # Load TOML configuration for each molecule
    with open("./data/benzene.toml", "rb") as f:
        config_benzene = tomllib.load(f)
    with open("./data/naphthalene.toml", "rb") as f:
        config_naphthalene = tomllib.load(f)

    # Map molecule names to their configuration dictionaries
    config_map = {
        config_benzene["molecule_name"]: config_benzene,
        config_naphthalene["molecule_name"]: config_naphthalene,
    }
    site_mol_name: List[str] = []

    # Extract parameters
    molecules       = input_data.get("molecules", [])
    center_of_mass  = input_data.get("center_of_mass", [])
    rotations       = input_data.get("rotation", [])
    repeats         = input_data.get("supercell_repeat", [1,1,1])
    raw_cell        = np.array(input_data.get("cell", []), float)
    use_fractional  = input_data.get("use_fractional", False)  # kept for compatibility

    def expand_symbols(tokens):
        """
        Expand tokens like ["C6H6"] or ["C6","H6"] into a flat list of symbols.
        """
        pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
        result = []
        for token in tokens:
            for element, count_str in pattern.findall(token):
                count = int(count_str) if count_str else 1
                result.extend([element] * count)
        return result

    def compute_positions(cfg, com, angles):
        """Compute rotated+translated Cartesian positions for a molecule."""
        pos = np.array(cfg["positions"], float)
        if cfg.get("use_fractional", True):
            lattice = np.array(cfg["lattice_vectors"], float)
            pos = pos.dot(lattice)
        rot = R.from_euler('xyz', angles, degrees=True)
        return rot.apply(pos) + np.array(com, float)

    # Per-cell molecule substitutions
    mol_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"],): s["new_molecule"]
        for s in input_data.get("molecule_substitutions", [])
    }
    # Per-atom substitutions
    atom_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"], s["atom_index"]): s["new_symbol"]
        for s in input_data.get("atom_substitutions", [])
    }
    # Global atom substitutions
    global_atom_subs = {s["old_symbol"]: s["new_symbol"]
                        for s in input_data.get("global_atom_substitutions", [])}
    # Molecule-level (all atoms) substitutions
    molecule_atom_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"],): s["new_symbols"]
        for s in input_data.get("molecule_atom_substitutions", [])
    }

    super_symbols: List[str] = []
    super_positions: List[List[float]] = []
    mol_ids: List[int] = []
    groups: List[List[int]] = []
    site_angles: List[Tuple[float, float, float]] = []

    nx, ny, nz = repeats
    site_counter = 0

    for (i, j, k) in product(range(nx), range(ny), range(nz)):
        shift = i * raw_cell[0] + j * raw_cell[1] + k * raw_cell[2]

        for mol_idx, mol_name in enumerate(molecules):
            site_id = site_counter
            site_counter += 1

            key_mol = (i, j, k, mol_idx)
            use_name = mol_subs.get(key_mol, mol_name)
            cfg = config_map[use_name]
            site_mol_name.append(use_name)

            angles = rotations[mol_idx]  # Euler angles [deg]
            com = center_of_mass[mol_idx]
            pos_arr = compute_positions(cfg, com, angles)
            syms = expand_symbols(cfg["symbols"])

            # substitutions
            syms = [global_atom_subs.get(sym, sym) for sym in syms]
            override_tokens = molecule_atom_subs.get(key_mol)
            if override_tokens is not None:
                expanded = expand_symbols(override_tokens)
                assert len(expanded) == len(syms), (
                    f"Expected {len(syms)} atoms, but got {len(expanded)}"
                )
                syms = expanded

            base_idx = len(super_symbols)
            groups.append([])

            for atom_local_idx, (sym, pos) in enumerate(zip(syms, pos_arr)):
                final_sym = atom_subs.get((i, j, k, mol_idx, atom_local_idx), sym)
                super_symbols.append(final_sym)
                super_positions.append(list(pos + shift))
                atom_global_idx = base_idx + atom_local_idx
                mol_ids.append(site_id)
                groups[site_id].append(atom_global_idx)

            site_angles.append(tuple(angles))

    # Define scaled cell for supercell
    scaled_cell = (raw_cell.T * repeats).T  # (3,3) ndarray

    # Build ASE Atoms object
    system = Atoms(
        symbols=super_symbols,
        positions=super_positions,
        cell=scaled_cell.tolist(),
        pbc=True,
    )

    # ---- Attach arrays / info for downstream grouping and optimization ----
    system.new_array("mol_id", np.array(mol_ids, dtype=int))
    system.info["groups"] = groups
    system.info["site_angles_deg"] = np.array(site_angles, float)  # (num_sites, 3)
    system.info["site_mol_name"] = site_mol_name

    dict_meta = {}
    for name, cfg in config_map.items():
        dict_meta[name] = {
            "lattice_vectors": np.array(cfg["lattice_vectors"], dtype=float),
            "center_of_mass":  np.array(cfg["center_of_mass"],  dtype=float),  # (K,3)
        }
    system.info["dict_meta"] = dict_meta

    # --- NEW: store rng settings so every stage can pick them up ---
    system.info["rng"] = {
        "global":   int(rng_global),
        "pre_scale": int(rng_pre_scale),
        "optimize":  int(rng_optimize),
        "render":    int(rng_render),
    }

    # -------- concise stdout (summary) --------
    if verbose:
        n_atoms = len(super_symbols)
        n_sites = len(groups)
        vol = float(abs(np.linalg.det(scaled_cell)))
        counts = Counter(site_mol_name)
        mol_line = ", ".join(f"{k}×{v}" for k, v in sorted(counts.items()))
        a, b, c = scaled_cell
        print("[build_system_from_toml] Built system.")
        print(f"  atoms={n_atoms}, sites={n_sites}, repeats={tuple(repeats)}, "
              f"cell_volume={vol:.3f} Å^3")
        print(f"  cell(Å): a={a.tolist()}, b={b.tolist()}, c={c.tolist()}")
        print(f"  molecules: {mol_line}")
        print(f"  rng settings: {system.info['rng']}")

    return system
