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

# Pattern to split symbol and count, e.g. "C6" -> ("C", "6")
pattern = re.compile(r"^([A-Z][a-z]?)(\d+)$")

# Extract parameters
molecules       = input_data.get("molecules", [])
center_of_mass  = input_data.get("center_of_mass", [])
rotations       = input_data.get("rotation", [])
repeats         = input_data.get("supercell_repeat", [1,1,1])
raw_cell        = np.array(input_data.get("cell", []))

# Preprocess shorthand expansion for a given molecule config
def expand_symbols(symbols):
    result = []
    for spec in symbols:
        m = pattern.match(spec)
        if m:
            result.extend([m.group(1)] * int(m.group(2)))
        else:
            result.append(spec)
    return result

# Function to compute positions for a molecule (rotated+translated)
def compute_positions(cfg, com, angles):
    pos = np.array(cfg["positions"])
    rot = R.from_euler('xyz', angles, degrees=True)
    return rot.apply(pos) + np.array(com)

# Prepare substitution dictionaries keyed by (i,j,k,mol_idx)
mol_subs = {}
for sub in input_data.get("molecule_substitutions", []):
    key = tuple(sub["cell_index"]) + (sub["molecule_index"],)
    mol_subs[key] = sub["new_molecule"]

atom_subs = {}
for sub in input_data.get("atom_substitutions", []):
    key = tuple(sub["cell_index"]) + (sub["molecule_index"], sub["atom_index"])
    atom_subs[key] = sub["new_symbol"]

# Build supercell lists
super_symbols = []
super_positions = []
nx, ny, nz = repeats

for (i, j, k) in product(range(nx), range(ny), range(nz)):
    # compute cell translation
    shift = i*raw_cell[0] + j*raw_cell[1] + k*raw_cell[2]
    # loop over molecules in base cell
    for mol_idx, mol_name in enumerate(molecules):
        # determine actual molecule after substitution
        key_mol = (i, j, k, mol_idx)
        use_name = mol_subs.get(key_mol, mol_name)
        cfg = config_map[use_name]

        # compute base positions and symbols
        angles = rotations[mol_idx]
        com    = center_of_mass[mol_idx]
        pos_arr = compute_positions(cfg, com, angles)
        syms = expand_symbols(cfg["symbols"])

        # append each atom, apply atom-level substitution if any
        for atom_idx, (sym, pos) in enumerate(zip(syms, pos_arr)):
            # global substitution key for atom
            key_atom = (i, j, k, mol_idx, atom_idx)
            final_sym = atom_subs.get(key_atom, sym)
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

# Define pseudopotentials for each element; ensure these UPF files exist
pseudopotentials = {
    'C': 'C.pbe-nc.UPF',
    'H': 'H.pbe-nc.UPF',
    'Si': 'Si.pbe-nc.UPF'
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