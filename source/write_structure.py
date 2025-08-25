from ase.atoms import Atoms
from ase.io import write

def write_structure(system: Atoms) -> None:
    """Write multiple simulator-friendly files from an ASE Atoms object."""

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

    return