# Import necessary modules
import re
import numpy as np
import tomllib
from scipy.spatial.transform import Rotation as R
from ase import Atoms
from ase.io import write
from itertools import product

# Load input from TOML file
with open("input/input.toml", "rb") as f:
    input_data = tomllib.load(f)

# Load TOML configuration for each molecule
with open("input/benzene.toml", "rb") as f:
    config_benzene = tomllib.load(f)
with open("input/naphthalene.toml", "rb") as f:
    config_naphthalene = tomllib.load(f)

# Map molecule names to their configuration dictionaries
config_map = {
    config_benzene["molecule_name"]: config_benzene,
    config_naphthalene["molecule_name"]: config_naphthalene,
}

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
    pos = np.array(cfg["positions"])
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
super_symbols = []
super_positions = []
nx, ny, nz = repeats

for (i, j, k) in product(range(nx), range(ny), range(nz)):
    # Compute translation for this cell
    shift = i * raw_cell[0] + j * raw_cell[1] + k * raw_cell[2]

    for mol_idx, mol_name in enumerate(molecules):
        key_mol = (i, j, k, mol_idx)
        # 1) molecule_substitutions
        use_name = mol_subs.get(key_mol, mol_name)
        cfg = config_map[use_name]

        # Compute positions and base symbols
        angles = rotations[mol_idx]
        com = center_of_mass[mol_idx]
        pos_arr = compute_positions(cfg, com, angles)
        syms = expand_symbols(cfg["symbols"])

        # 2) Apply global atom substitutions
        syms = [global_atom_subs.get(sym, sym) for sym in syms]

        # 3) Apply molecule-level (all atoms) substitutions
        override_tokens = molecule_atom_subs.get(key_mol)
        if override_tokens is not None:
            expanded = expand_symbols(override_tokens)
            assert len(expanded) == len(syms), (
                f"Expected {len(syms)} atoms, but got {len(expanded)}"
            )
            syms = expanded

        # 4) Loop over atoms and apply cell/molecule/atom-specific substitutions
        for atom_idx, (sym, pos) in enumerate(zip(syms, pos_arr)):
            final_sym = atom_subs.get((i, j, k, mol_idx, atom_idx), sym)
            super_symbols.append(final_sym)
            super_positions.append(list(pos + shift))


# Define scaled cell for supercell
scaled_cell = (raw_cell.T * repeats).T.tolist()

# Build ASE Atoms object
system = Atoms(
    symbols=super_symbols,
    positions=super_positions,
    cell=scaled_cell,
    pbc=True
)

# Print the summary
print("Scaled cell vectors:", scaled_cell)
print(system)

# Define pseudopotentials for each element
# 1. Get unique chemical species symbols
# get_chemical_symbols() returns a list of strings, so use set() to remove duplicates
unique_species = sorted(set(system.get_chemical_symbols()))
# 2. Dynamically create the pseudopotentials dictionary
# Apply the ".pbe-nc.UPF" suffix to all species
pseudopotentials = {
    elem: f"{elem}.pbe-nc.UPF"
    for elem in unique_species
}

# Set Quantum ESPRESSO namelist parameters
input_data = {
    'control': {
        'calculation': 'scf',        # single-point energy calculation
        'prefix':      'molecule',   # output file prefix
        'outdir':      './tmp',      # working directory
    },
    'system': {
        'ecutwfc':  40,    # plane-wave cutoff energy in Ry
        'ecutrho': 320,    # charge density cutoff in Ry
    },
    'electrons': {
        'conv_thr': 1e-8,  # convergence threshold for electronic SCF
    }
}

# Define k-point grid (Gamma-point only)
kpoints = (1, 1, 1)

# Write Quantum ESPRESSO input file (PWscf)
write('qe.in',
      system,
      format='espresso-in',
      pseudopotentials=pseudopotentials,
      input_data=input_data,
      kpoints=kpoints)

print("Quantum ESPRESSO input file 'qe.in' has been generated.")

# Generate input files for VASP calculations
# --- 1. Write POSCAR file ---
# Use direct (fractional) coordinates
write('POSCAR', system, format='vasp', direct=True)

# --- 2. Define INCAR parameters ---
incar_settings = {
    'PREC':     'Accurate',    # precision level
    'ENCUT':    300,           # plane-wave cutoff energy [eV]
    'EDIFF':    1e-6,          # electronic convergence criterion
    'ISMEAR':   0,             # Gaussian smearing
    'SIGMA':    0.05,          # smearing width [eV]
    'IBRION':   -1,            # no ionic relaxation (static run)
    'NSW':      0,             # number of ionic steps
    'LWAVE':    True,          # write WAVECAR
    'LCHARG':   True,          # write CHGCAR
}

with open('INCAR', 'w') as f:
    for key, val in incar_settings.items():
        f.write(f"{key} = {val}\n")

# --- 3. Write KPOINTS file (Gamma-point only) ---
with open('KPOINTS', 'w') as f:
    f.write("Gamma-point sampling\n")    # comment line
    f.write("0\n")                       # 0 means automatic generation
    f.write("Gamma\n")                   # centering
    f.write("1 1 1\n")                   # k-point grid
    f.write("0 0 0\n")                   # shift

# Generate input files for VESTA
write('POSCAR.vasp', system, format='vasp', direct=True)
write('structure.cif', system, format='cif')