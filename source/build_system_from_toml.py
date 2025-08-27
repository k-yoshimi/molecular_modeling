import re
import numpy as np
import tomllib
from typing import Dict, List, Tuple
from scipy.spatial.transform import Rotation as R
from ase import Atoms
from itertools import product

def build_system_from_toml() -> Atoms:
    """Build an ASE Atoms system from TOML inputs (supercell + substitutions)."""
   
    # Load input from TOML file
    with open("./input/input.toml", "rb") as f:
        input_data = tomllib.load(f)
    
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
    raw_cell        = np.array(input_data.get("cell", []))
    use_fractional  = input_data.get("use_fractional", False)
    
    # Preprocess shorthand expansion for a given molecule config
    def expand_symbols(tokens):
        """
        Expand a list of chemical-symbol tokens like ["C6","H6"] or ["C6H6"]
        into a flat list of individual symbols: ["C","C",...,"H","H",...].
        """
        # Regex: element = uppercase + optional lowercase, count = digits (optional)
        pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
        result = []
        for token in tokens:
            # findall returns list of (element, count_str) tuples
            for element, count_str in pattern.findall(token):
                count = int(count_str) if count_str else 1
                # repeat each element count times
                result.extend([element] * count)
        return result
    
    # Function to compute positions for a molecule (rotated+translated)
    def compute_positions(cfg, com, angles):
        # Load raw positions from molecule dictionary
        pos = np.array(cfg["positions"])
    
        # If positions are given in fractional coordinates relative to the
        # molecule’s lattice_vectors, convert them to Cartesian.
        if cfg.get("use_fractional", True):
            # lattice_vectors is a list of three vectors [[ax,ay,az],...]
            lattice = np.array(cfg["lattice_vectors"])
            # fractional -> Cartesian: f1*a + f2*b + f3*c
            pos = pos.dot(lattice)
    
        # Apply rotation and then translate by center-of-mass
        rot = R.from_euler('xyz', angles, degrees=True)
        return rot.apply(pos) + np.array(com)
    
    # Build dictionaries for per-cell molecule substitutions
    mol_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"],): s["new_molecule"]
        for s in input_data.get("molecule_substitutions", [])
    }
    
    # Build dictionaries for per-cell, per-atom substitutions
    atom_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"], s["atom_index"]): s["new_symbol"]
        for s in input_data.get("atom_substitutions", [])
    }
    
    # Build dictionary for global atom substitutions
    global_atom_subs = {
        s["old_symbol"]: s["new_symbol"]
        for s in input_data.get("global_atom_substitutions", [])
    }
    
    # Build dictionary for molecule-level (all atoms) substitutions
    molecule_atom_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"],): s["new_symbols"]
        for s in input_data.get("molecule_atom_substitutions", [])
    }
    
    # Build supercell lists  
    super_symbols: List[str] = []
    super_positions: List[List[float]] = []
    mol_ids: List[int] = []              # per-atom molecule (site) id
    groups: List[List[int]] = []         # list of atom-index lists per site
    site_angles: List[Tuple[float, float, float]] = []  # per-site initial Euler angles (deg)

    nx, ny, nz = repeats
    site_counter = 0

    for (i, j, k) in product(range(nx), range(ny), range(nz)):
        shift = i * raw_cell[0] + j * raw_cell[1] + k * raw_cell[2]

        for mol_idx, mol_name in enumerate(molecules):
            # ---- Assign site_id in the *exact* placement order ----
            site_id = site_counter
            site_counter += 1

            # 1) molecule_substitutions
            key_mol = (i, j, k, mol_idx)
            use_name = mol_subs.get(key_mol, mol_name)
            cfg = config_map[use_name]
            site_mol_name.append(use_name)

            # Compute positions and base symbols
            angles = rotations[mol_idx]                  # Euler angles [deg] for this molecule
            com = center_of_mass[mol_idx]
            pos_arr = compute_positions(cfg, com, angles)
            syms = expand_symbols(cfg["symbols"])

            # Apply substitutions (global / molecule-level / per-atom) -- your existing logic
            # 2) global_atom_substitutions
            syms = [global_atom_subs.get(sym, sym) for sym in syms]

            # 3) molecule_atom_substitutions (override all atoms if provided)
            override_tokens = molecule_atom_subs.get(key_mol)
            if override_tokens is not None:
                expanded = expand_symbols(override_tokens)
                assert len(expanded) == len(syms), (
                    f"Expected {len(syms)} atoms, but got {len(expanded)}"
                )
                syms = expanded

            # ---- Start index of this molecule's atoms in the flattened arrays ----
            base_idx = len(super_symbols)
            groups.append([])  # prepare group for this site_id

            # 4) per-atom substitutions + append
            for atom_local_idx, (sym, pos) in enumerate(zip(syms, pos_arr)):
                final_sym = atom_subs.get((i, j, k, mol_idx, atom_local_idx), sym)
                super_symbols.append(final_sym)
                super_positions.append(list(pos + shift))
                # Record per-atom mol_id and group membership
                atom_global_idx = base_idx + atom_local_idx
                mol_ids.append(site_id)
                groups[site_id].append(atom_global_idx)

            # Store per-site initial angles in placement order
            site_angles.append(tuple(angles))

    # Define scaled cell for supercell
    scaled_cell = (raw_cell.T * repeats).T.tolist()

    # Build ASE Atoms object
    system = Atoms(
        symbols=super_symbols,
        positions=super_positions,
        cell=scaled_cell,
        pbc=True,
    )

    # ---- Attach arrays / info for downstream grouping and optimization ----
    # Each entry in mol_ids is an int (site_id) for that atom
    system.new_array("mol_id", np.array(mol_ids, dtype=int))
    # Optional conveniences:
    system.info["groups"] = groups                           # list[list[int]]
    system.info["site_angles_deg"] = np.array(site_angles)    # shape = (num_sites, 3)
    system.info["site_mol_name"] = site_mol_name

    dict_meta = {}
    for name, cfg in config_map.items():
        dict_meta[name] = {
            "lattice_vectors": np.array(cfg["lattice_vectors"], dtype=float),
            "center_of_mass":  np.array(cfg["center_of_mass"],  dtype=float),  # (K,3) を想定
        }
    system.info["dict_meta"] = dict_meta

    # Debug print (optional)
    print("Scaled cell vectors:", scaled_cell)
    print(system)
    return system