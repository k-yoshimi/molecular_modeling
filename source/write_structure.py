from pathlib import Path

from ase.atoms import Atoms
from ase.io import write


def write_structure(system: Atoms, *, outdir: str | Path = ".", verbose: bool = True):
    """
    Write simulator-friendly structure files from an ASE `Atoms` object.

    Parameters
    ----------
    system : ase.Atoms
        System to export.
    outdir : str or Path, optional
        Output directory. Created if missing.
    verbose : bool, optional
        If True, print a summary of written files.

    Writes
    ------
    - Quantum ESPRESSO input: `qe.in`
    - VASP inputs: `POSCAR`, `INCAR`, `KPOINTS`, `POSCAR.vasp`
    - CIF: `structure.cif`
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []

    unique_species = sorted(set(system.get_chemical_symbols()))
    pseudopotentials = {elem: f"{elem}.pbe-nc.UPF" for elem in unique_species}

    input_data = {
        'control': {'calculation': 'scf', 'prefix': 'molecule', 'outdir': './tmp'},
        'system':  {'ecutwfc': 40, 'ecutrho': 320},
        'electrons': {'conv_thr': 1e-8},
    }
    kpoints = (1, 1, 1)

    # QE
    qe_in = outdir / "qe.in"
    write(qe_in.as_posix(), system, format='espresso-in',
          pseudopotentials=pseudopotentials, input_data=input_data, kpoints=kpoints)
    written.append(qe_in.as_posix())

    # VASP
    poscar = outdir / "POSCAR"
    write(poscar.as_posix(), system, format='vasp', direct=True)
    written.append(poscar.as_posix())

    incar = outdir / "INCAR"
    with open(incar, "w", encoding="utf-8") as f:
        for k, v in {
            'PREC': 'Accurate', 'ENCUT': 300, 'EDIFF': 1e-6, 'ISMEAR': 0, 'SIGMA': 0.05,
            'IBRION': -1, 'NSW': 0, 'LWAVE': True, 'LCHARG': True,
        }.items():
            f.write(f"{k} = {v}\n")
    written.append(incar.as_posix())

    kpoints_file = outdir / "KPOINTS"
    with open(kpoints_file, "w", encoding="utf-8") as f:
        f.write("Gamma-point sampling\n0\nGamma\n1 1 1\n0 0 0\n")
    written.append(kpoints_file.as_posix())

    # Exports
    poscar_vasp = outdir / "POSCAR.vasp"
    write(poscar_vasp.as_posix(), system, format='vasp', direct=True)
    written.append(poscar_vasp.as_posix())

    cif = outdir / "structure.cif"
    write(cif.as_posix(), system, format='cif')
    written.append(cif.as_posix())

    if verbose:
        names = ", ".join(Path(p).name for p in written)
        print(f"[write_structure] Structure input and export files have been written ({names}).")

    return
