# build_system_from_toml.py
import re
import numpy as np
import tomllib
from typing import Dict, List, Tuple, Optional
from scipy.spatial.transform import Rotation as R
from ase import Atoms
from itertools import product
from collections import Counter

def build_system_from_toml(*, verbose: bool = True) -> Atoms:
    """
    Build an ASE Atoms system from TOML inputs (supercell + substitutions).

    Key points in this version:
      - Dictionary TOMLs (e.g., benzene.toml) may include `voronoi_ref_id` (0-based).
        * Atomic positions are first shifted so that they correspond to the
          voronoi_ref_id-th molecule (same cell), then rotated and placed.
        * dict_meta stores voronoi_ref_id for downstream Voronoi/IoU routines.
      - Backward compatible: if `voronoi_ref_id` is absent, ref_id = 0.
    """
    # ----------- Read inputs -----------
    with open("./input/input.toml", "rb") as f:
        input_data = tomllib.load(f)

    # rng: keep only global/optimize/render (compatibility with older keys)
    rng_global   = int(input_data.get("rng", 1))
    rng_optimize = int(input_data.get("optimize_rng",  rng_global))
    rng_render   = int(input_data.get("render_rng",    rng_global))

    # Load dictionary molecule configs (extend here as needed)
    with open("./data/benzene.toml", "rb") as f:
        config_benzene = tomllib.load(f)
    with open("./data/naphthalene.toml", "rb") as f:
        config_naphthalene = tomllib.load(f)

    config_map: Dict[str, dict] = {
        config_benzene["molecule_name"]: config_benzene,
        config_naphthalene["molecule_name"]: config_naphthalene,
    }

    # ----------- Input fields -----------
    molecules:       List[str] = input_data.get("molecules", [])
    center_of_mass:  List[List[float]] = input_data.get("center_of_mass", [])
    rotations:       List[List[float]] = input_data.get("rotation", [])
    repeats:         List[int] = input_data.get("supercell_repeat", [1, 1, 1])
    raw_cell         = np.array(input_data.get("cell", []), float)
    use_fractional   = input_data.get("use_fractional", False)  # kept for compatibility

    # ----------- Helpers -----------
    def expand_symbols(tokens: List[str]) -> List[str]:
        """Expand tokens like ['C6H6', 'Si6O6'] into a flat list of symbols."""
        pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
        result: List[str] = []
        for token in tokens:
            for element, count_str in pattern.findall(token):
                count = int(count_str) if count_str else 1
                result.extend([element] * count)
        return result

    def compute_positions(cfg: dict, com: List[float], angles: List[float]) -> np.ndarray:
        """
        Positions are already for the molecule selected by 'voronoi_ref_id'.
        We must rotate about that molecule's COM, then place at the site COM 'com'.
        """
        pos = np.array(cfg["positions"], float)
        lattice = np.array(cfg["lattice_vectors"], float)
    
        # 1) fractional -> Cartesian if needed
        if cfg.get("use_fractional", True):
            pos = pos @ lattice
    
        # 2) recenter to rotate about the molecule's own COM
        #    (User promises: COM(positions) == center_of_mass[voronoi_ref_id])
        pos_com = pos.mean(axis=0)  # robust even if tiny rounding exists
        pos_centered = pos - pos_com
    
        # 3) rotate about origin (which is now the molecule COM)
        rot = R.from_euler('xyz', angles, degrees=True)
        pos_rot = rot.apply(pos_centered)
    
        # 4) place at the target site COM
        return pos_rot + np.array(com, float)

    # ----------- Substitution tables -----------
    mol_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"],): s["new_molecule"]
        for s in input_data.get("molecule_substitutions", [])
    }
    atom_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"], s["atom_index"]): s["new_symbol"]
        for s in input_data.get("atom_substitutions", [])
    }
    global_atom_subs = {s["old_symbol"]: s["new_symbol"]
                        for s in input_data.get("global_atom_substitutions", [])}
    molecule_atom_subs = {
        tuple(s["cell_index"]) + (s["molecule_index"],): s["new_symbols"]
        for s in input_data.get("molecule_atom_substitutions", [])
    }

    # ----------- Supercell assembly buffers -----------
    super_symbols: List[str] = []
    super_positions: List[List[float]] = []
    mol_ids: List[int] = []
    groups: List[List[int]] = []
    site_angles: List[Tuple[float, float, float]] = []
    site_mol_name: List[str] = []

    nx, ny, nz = map(int, repeats)
    site_counter = 0

    # ----------- Build supercell by tiling + per-site placement -----------
    for (i, j, k) in product(range(nx), range(ny), range(nz)):
        shift = i * raw_cell[0] + j * raw_cell[1] + k * raw_cell[2]
        for mol_idx, mol_name in enumerate(molecules):
            site_id = site_counter
            site_counter += 1

            key_mol = (i, j, k, mol_idx)
            use_name = mol_subs.get(key_mol, mol_name)
            cfg = config_map[use_name]
            site_mol_name.append(use_name)

            # Euler angles & COM for this site
            angles = rotations[mol_idx]
            com = center_of_mass[mol_idx]

            # Build positions in Cartesian, already matched to cfg['voronoi_ref_id']
            pos_arr = compute_positions(cfg, com, angles)

            # Symbols (expand formula tokens if present)
            syms = cfg["symbols"]
            syms = expand_symbols(syms)

            # Global replacements
            syms = [global_atom_subs.get(sym, sym) for sym in syms]

            # Per-molecule override (expanded tokens)
            override_tokens = molecule_atom_subs.get(key_mol)
            if override_tokens is not None:
                expanded = expand_symbols(override_tokens)
                assert len(expanded) == len(syms), (
                    f"Expanded symbols length mismatch at cell={key_mol[:3]} mol={key_mol[3]}: "
                    f"{len(expanded)} != {len(syms)}"
                )
                syms = expanded

            # Append atoms of this site
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

    # ----------- Build ASE Atoms -----------
    scaled_cell = (raw_cell.T * repeats).T
    system = Atoms(
        symbols=super_symbols,
        positions=super_positions,
        cell=scaled_cell.tolist(),
        pbc=True,
    )
    system.new_array("mol_id", np.array(mol_ids, dtype=int))
    system.info["groups"] = groups
    system.info["site_angles_deg"] = np.array(site_angles, float)
    system.info["site_mol_name"] = site_mol_name

    # Pack dictionary meta for downstream Voronoi/IoU:
    #   - lattice_vectors (3x3)
    #   - center_of_mass  (K,3)
    #   - voronoi_ref_id  (int)
    dict_meta: Dict[str, dict] = {}
    for name, cfg in config_map.items():
        dict_meta[name] = {
            "lattice_vectors": np.array(cfg["lattice_vectors"], dtype=float),
            "center_of_mass":  np.array(cfg["center_of_mass"],  dtype=float),
            "voronoi_ref_id":  int(cfg.get("voronoi_ref_id", 0)),
        }
    system.info["dict_meta"] = dict_meta

    # rng info (downstream optimizers/renderers will read this)
    system.info["rng"] = {
        "global":   int(rng_global),
        "optimize": int(rng_optimize),
        "render":   int(rng_render),
    }

    # ----------- Verbose summary -----------
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
