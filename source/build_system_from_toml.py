# build_system_from_toml.py
import re
import tomllib
from collections import Counter
from itertools import product
from typing import Dict, List

import numpy as np
from ase import Atoms
from scipy.spatial.transform import Rotation as R


def build_system_from_toml(*, verbose: bool = True) -> Atoms:
    """
    Build an ASE `Atoms` object from TOML inputs (supercell + substitutions).

    Parameters
    ----------
    verbose : bool, optional
        If True, print a summary of the assembled system.

    Returns
    -------
    ase.Atoms
        Periodic supercell with per-atom `mol_id` and metadata in
        `system.info` (`groups`, `site_mol_name`, `dict_meta`, `rng`).

    Notes
    -----
    - Input (`input/input.toml`)
        - `data_dir` specifies the directory for molecule TOML files (default: `./data`).
        - `use_fractional=true` interprets `center_of_mass` as fractional and converts to
          Cartesian via `cell` (rows a,b,c).
        - Downstream code always uses Cartesian site COMs.
    - Dictionary (`{data_dir}/<name>.toml`)
        - `voronoi_ref_id` selects which molecule COM is the reference (0-based).
        - If `use_fractional=true`, both `positions` and `center_of_mass` are converted to
          Cartesian with that dictionary's `lattice_vectors`.
        - `dict_meta` stores dictionary COMs in Cartesian.
    - Placement
        - For each site: recenter `positions` by its COM, rotate (initially zero), then
          translate to the site COM; lattice tiling handled by `supercell_repeat`.
    """
    # ----------- Read inputs -----------
    with open("./input/input.toml", "rb") as f:
        input_data = tomllib.load(f)

    # rng: keep only global/optimize/render (compatibility with older keys)
    rng_global   = int(input_data.get("rng", 1))
    rng_optimize = int(input_data.get("optimize_rng",  rng_global))
    rng_render   = int(input_data.get("render_rng",    rng_global))

    # Data directory for molecule TOML files (default: ./data)
    data_dir = input_data.get("data_dir", "./data")

    # Load dictionary molecule configs
    config_map: Dict[str, dict] = {}
    for molecule in input_data.get("molecules", []):
        with open(f"{data_dir}/{molecule}.toml", "rb") as f:
            config_map[molecule] = tomllib.load(f)

    # ----------- Input fields -----------
    molecules: List[str] = input_data.get("molecules", [])
    repeats:   List[int] = input_data.get("supercell_repeat", [1, 1, 1])

    raw_cell = np.array(input_data.get("cell", []), float)  # rows a,b,c (Å)
    input_use_fractional = bool(input_data.get("use_fractional", False))

    # site COMs: maybe fractional (if input_use_fractional)
    center_of_mass_in = np.asarray(input_data.get("center_of_mass", []), float)
    if input_use_fractional:
        # fractional -> Cartesian using input cell
        center_of_mass_cart = center_of_mass_in @ raw_cell
    else:
        center_of_mass_cart = center_of_mass_in

    # ----------- Helpers -----------
    def expand_symbols(tokens: List[str]) -> List[str]:
        """
        Expand compact chemical tokens into elemental symbols.

        Parameters
        ----------
        tokens : list[str]
            Strings such as ``"C6H6"`` or ``"Si6O6"``.

        Returns
        -------
        list[str]
            Flattened elemental symbols expanded by count.
        """
        pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
        result: List[str] = []
        for token in tokens:
            for element, count_str in pattern.findall(token):
                count = int(count_str) if count_str else 1
                result.extend([element] * count)
        return result

    def compute_positions(cfg: dict, site_com_cart: np.ndarray, angles_xyz_deg: List[float]) -> np.ndarray:
        """
        Compute rotated and translated Cartesian positions for one site.

        Parameters
        ----------
        cfg : dict
            Molecule dictionary entry with `positions`, `center_of_mass`, `lattice_vectors`,
            and optional `use_fractional` / `voronoi_ref_id`.
        site_com_cart : ndarray
            Target site center-of-mass in Cartesian coordinates.
        angles_xyz_deg : list[float]
            Euler angles (XYZ, degrees) applied about the molecule COM.

        Returns
        -------
        numpy.ndarray
            Rotated and placed Cartesian positions (N_atoms, 3).

        Notes
        -----
        1. If `use_fractional` is True, positions are converted with `lattice_vectors`.
        2. Positions are recentered to their own COM, rotated, then translated to `site_com_cart`.
        """
        pos = np.array(cfg["positions"], float)
        lattice = np.array(cfg["lattice_vectors"], float)

        # 1) fractional -> Cartesian if needed (no wrapping)
        if bool(cfg.get("use_fractional", False)):
            pos = pos @ lattice

        # 2) recenter to molecule's COM (robust even with slight numeric noise)
        pos_com = pos.mean(axis=0)
        pos_centered = pos - pos_com

        # 3) rotate about origin (now molecule COM)
        rot = R.from_euler('xyz', angles_xyz_deg, degrees=True)
        pos_rot = rot.apply(pos_centered)

        # 4) place to site COM (Cartesian)
        return pos_rot + np.asarray(site_com_cart, float)

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
    site_mol_name: List[str] = []

    nx, ny, nz = map(int, repeats)
    site_counter = 0

    # ----------- Build supercell by tiling + per-site placement -----------
    for (i, j, k) in product(range(nx), range(ny), range(nz)):
        # lattice shift for this cell instance
        shift = i * raw_cell[0] + j * raw_cell[1] + k * raw_cell[2]

        for mol_idx, mol_name in enumerate(molecules):
            site_id = site_counter
            site_counter += 1

            key_mol = (i, j, k, mol_idx)
            use_name = mol_subs.get(key_mol, mol_name)
            cfg = config_map[use_name]
            site_mol_name.append(use_name)

            # Euler angles & site COM (ALWAYS Cartesian)
            # No input rotation: place unrotated; optimizer will rotate later.
            angles = [0.0, 0.0, 0.0]
            site_com = center_of_mass_cart[mol_idx]

            # Build positions in Cartesian (already for voronoi_ref_id)
            pos_arr = compute_positions(cfg, site_com, angles)

            # Expand and apply symbol substitutions
            syms = expand_symbols(cfg["symbols"])
            syms = [global_atom_subs.get(sym, sym) for sym in syms]

            # Per-molecule atom override (expanded tokens)
            override_tokens = molecule_atom_subs.get(key_mol)
            if override_tokens is not None:
                expanded = expand_symbols(override_tokens)
                assert len(expanded) == len(syms), (
                    f"Expanded symbols length mismatch at cell={key_mol[:3]} mol={key_mol[3]}: "
                    f"{len(expanded)} != {len(syms)}"
                )
                syms = expanded

            # Append atoms of this site (with per-atom substitutions)
            base_idx = len(super_symbols)
            groups.append([])
            for atom_local_idx, (sym, pos) in enumerate(zip(syms, pos_arr)):
                final_sym = atom_subs.get((i, j, k, mol_idx, atom_local_idx), sym)
                super_symbols.append(final_sym)
                super_positions.append(list(pos + shift))
                atom_global_idx = base_idx + atom_local_idx
                mol_ids.append(site_id)
                groups[site_id].append(atom_global_idx)

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
    system.info["site_mol_name"] = site_mol_name

    # ----------- Pack dictionary meta (ALWAYS Cartesian COMs) -----------
    dict_meta: Dict[str, dict] = {}
    for name, cfg in config_map.items():
        lattice = np.array(cfg["lattice_vectors"], dtype=float)          # rows a,b,c (Å)
        coms    = np.array(cfg["center_of_mass"],  dtype=float)          # (K,3), fractional or Cartesian
        if bool(cfg.get("use_fractional", False)):
            coms_cart = coms @ lattice                                   # fractional -> Cartesian
        else:
            coms_cart = coms
        dict_meta[name] = {
            "lattice_vectors": lattice,                                  # Cartesian basis (Å)
            "center_of_mass":  coms_cart,                                # ALWAYS Cartesian from here on
            "voronoi_ref_id":  int(cfg.get("voronoi_ref_id", 0)),
        }
    system.info["dict_meta"] = dict_meta

    # rng info (downstream readers will use this)
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
        print(f"  rng setting: rng={system.info['rng']['global']}")

    return system
