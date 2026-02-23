# source/gui/export_dialog.py
"""
Dialog for exporting structure files.

Allows exporting to Quantum ESPRESSO, VASP, and CIF formats.
"""

import tkinter as tk
from pathlib import Path
from tkinter import filedialog, ttk
from typing import List, Optional

from ase import Atoms
from ase.io import write

from .base_dialog import BaseDialog


class ExportDialog(BaseDialog):
    """
    Dialog for exporting structure files.

    Allows selecting output directory and formats to export.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    atoms : ase.Atoms
        ASE Atoms object to export.
    default_outdir : str, optional
        Default output directory.
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        atoms: Optional[Atoms] = None,
        default_outdir: str = "./structure_out",
    ):
        self.atoms = atoms
        self.default_outdir = default_outdir
        super().__init__(parent, title="Export Structure Files")

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        if self.atoms is None:
            ttk.Label(parent, text="No structure to export.").pack(pady=20)
            return

        # Output directory
        dir_frame = ttk.Frame(parent)
        dir_frame.pack(fill=tk.X, pady=5)

        ttk.Label(dir_frame, text="Output directory:").pack(side=tk.LEFT)

        self.outdir_var = tk.StringVar(value=self.default_outdir)
        self.outdir_entry = ttk.Entry(dir_frame, textvariable=self.outdir_var, width=30)
        self.outdir_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

        browse_btn = ttk.Button(dir_frame, text="Browse...", command=self._browse_dir)
        browse_btn.pack(side=tk.LEFT)

        # Format checkboxes
        ttk.Label(parent, text="Export formats:").pack(anchor=tk.W, pady=(10, 5))

        self.qe_var = tk.BooleanVar(value=True)
        self.vasp_var = tk.BooleanVar(value=True)
        self.cif_var = tk.BooleanVar(value=True)

        ttk.Checkbutton(
            parent,
            text="Quantum ESPRESSO (qe.in)",
            variable=self.qe_var,
        ).pack(anchor=tk.W)

        ttk.Checkbutton(
            parent,
            text="VASP (POSCAR, INCAR, KPOINTS)",
            variable=self.vasp_var,
        ).pack(anchor=tk.W)

        ttk.Checkbutton(
            parent,
            text="CIF (structure.cif)",
            variable=self.cif_var,
        ).pack(anchor=tk.W)

        # Structure info
        n_atoms = len(self.atoms)
        symbols = self.atoms.get_chemical_symbols()
        unique = sorted(set(symbols))
        ttk.Label(
            parent,
            text=f"Structure: {n_atoms} atoms ({', '.join(unique)})",
            foreground="gray",
        ).pack(anchor=tk.W, pady=(10, 0))

    def _browse_dir(self) -> None:
        """Open directory browser."""
        dir_path = filedialog.askdirectory(
            title="Select Output Directory",
            initialdir=self.outdir_var.get() or ".",
        )
        if dir_path:
            self.outdir_var.set(dir_path)

    def validate(self) -> bool:
        """Validate dialog input."""
        if self.atoms is None:
            self.show_error("Error", "No structure to export.")
            return False

        if not any([self.qe_var.get(), self.vasp_var.get(), self.cif_var.get()]):
            self.show_error("Error", "Please select at least one export format.")
            return False

        outdir = self.outdir_var.get().strip()
        if not outdir:
            self.show_error("Error", "Please specify an output directory.")
            return False

        return True

    def apply(self) -> None:
        """Apply dialog result - perform the export."""
        outdir = Path(self.outdir_var.get().strip())
        outdir.mkdir(parents=True, exist_ok=True)

        written: List[str] = []

        try:
            # Quantum ESPRESSO
            if self.qe_var.get():
                qe_path = outdir / "qe.in"
                unique_species = sorted(set(self.atoms.get_chemical_symbols()))
                pseudopotentials = {elem: f"{elem}.pbe-nc.UPF" for elem in unique_species}
                input_data = {
                    'control': {'calculation': 'scf', 'prefix': 'molecule', 'outdir': './tmp'},
                    'system': {'ecutwfc': 40, 'ecutrho': 320},
                    'electrons': {'conv_thr': 1e-8},
                }
                write(
                    qe_path.as_posix(),
                    self.atoms,
                    format='espresso-in',
                    pseudopotentials=pseudopotentials,
                    input_data=input_data,
                    kpoints=(1, 1, 1),
                )
                written.append("qe.in")

            # VASP
            if self.vasp_var.get():
                # POSCAR
                poscar_path = outdir / "POSCAR"
                write(poscar_path.as_posix(), self.atoms, format='vasp', direct=True)
                written.append("POSCAR")

                # INCAR
                incar_path = outdir / "INCAR"
                with open(incar_path, "w", encoding="utf-8") as f:
                    for k, v in {
                        'PREC': 'Accurate',
                        'ENCUT': 300,
                        'EDIFF': 1e-6,
                        'ISMEAR': 0,
                        'SIGMA': 0.05,
                        'IBRION': -1,
                        'NSW': 0,
                        'LWAVE': True,
                        'LCHARG': True,
                    }.items():
                        f.write(f"{k} = {v}\n")
                written.append("INCAR")

                # KPOINTS
                kpoints_path = outdir / "KPOINTS"
                with open(kpoints_path, "w", encoding="utf-8") as f:
                    f.write("Gamma-point sampling\n0\nGamma\n1 1 1\n0 0 0\n")
                written.append("KPOINTS")

                # POSCAR.vasp (alternate)
                poscar_vasp_path = outdir / "POSCAR.vasp"
                write(poscar_vasp_path.as_posix(), self.atoms, format='vasp', direct=True)
                written.append("POSCAR.vasp")

            # CIF
            if self.cif_var.get():
                cif_path = outdir / "structure.cif"
                # Create copy with PBC disabled to prevent atom wrapping
                atoms_copy = self.atoms.copy()
                atoms_copy.set_pbc(False)
                write(cif_path.as_posix(), atoms_copy, format='cif')
                written.append("structure.cif")

            self.result = {
                "outdir": str(outdir),
                "written": written,
                "success": True,
            }

            self.show_info(
                "Export Complete",
                f"Exported {len(written)} file(s) to {outdir}:\n" + "\n".join(written),
            )

        except Exception as e:
            self.result = {
                "outdir": str(outdir),
                "written": written,
                "success": False,
                "error": str(e),
            }
            self.show_error("Export Error", f"Failed to export: {e}")


def show_export_dialog(
    parent: Optional[tk.Widget] = None,
    atoms: Optional[Atoms] = None,
    default_outdir: str = "./structure_out",
) -> Optional[dict]:
    """
    Show the export dialog and return the result.

    Parameters
    ----------
    parent : tk.Widget, optional
        Parent widget.
    atoms : ase.Atoms, optional
        ASE Atoms object to export.
    default_outdir : str, optional
        Default output directory.

    Returns
    -------
    dict or None
        Dialog result with export information, or None if cancelled.
    """
    dialog = ExportDialog(parent, atoms=atoms, default_outdir=default_outdir)
    dialog.wait_window()
    return dialog.result if not dialog.cancelled else None
