# Import necessary modules
import tomllib
from ase import Atoms
from ase.io import write

# Load input
with open("input.toml", "rb") as f:
    input = tomllib.load(f)

# Load TOML configuration from file
with open("benzene.toml", "rb") as f:
    config_benzene = tomllib.load(f)
with open("naphthalene.toml", "rb") as f:
    config_naphthalene = tomllib.load(f)

# Map molecule names to their loaded config dicts
config_map = {
    config_benzene["molecule_name"]: config_benzene,
    config_naphthalene["molecule_name"]: config_naphthalene,
}

# Scale lattice vectors according to supercell_repeat
raw_cell = input["cell"]                # [[a1,a2,a3], [b1,b2,b3], [c1,c2,c3]]
repeats  = input["supercell_repeat"]    # [nx, ny, nz]

# Multiply each cell vector by its corresponding repeat factor
scaled_cell = [
    [coord * rep for coord in vec]
    for vec, rep in zip(raw_cell, repeats)
]

# Prepare lists for all symbols and positions
all_symbols = []
all_positions = []

# Loop over the sequence in input.toml, pairing each molecule with its center-of-mass offset
for mol_name, com in zip(input["molecules"], input["center_of_mass"]):
    # Retrieve the corresponding config
    cfg = config_map.get(mol_name)
    if cfg is None:
        raise ValueError(f"Unknown molecule: {mol_name}")
    
    # Extend symbols and positions with translation
    for symbol, pos in zip(cfg["symbols"], cfg["positions"]):
        # Translate each position by the center-of-mass vector
        translated_pos = [
            pos[0] + com[0],
            pos[1] + com[1],
            pos[2] + com[2],
        ]
        all_symbols.append(symbol)
        all_positions.append(translated_pos)

# Create ASE Atoms with the stretched cell
system = Atoms(
    symbols=all_symbols,
    positions=all_positions,
    cell=scaled_cell,
    pbc=True        # enable periodic boundary conditions
)
# Print the summary
print("Scaled cell vectors:", scaled_cell)
print(system)

# Define pseudopotentials for each element; ensure these UPF files exist
pseudopotentials = {
    'C': 'C.pbe-nc.UPF',
    'H': 'H.pbe-nc.UPF',
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