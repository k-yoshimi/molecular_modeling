# source/gui/molecule_list_dialog.py
"""
Dialog for displaying molecule/group information in the system.

Shows a list of all molecules in the current system with their properties.
"""

import tkinter as tk
from tkinter import ttk
from typing import Any, Dict, List, Optional

import numpy as np
from ase import Atoms

from .base_dialog import BaseDialog


class MoleculeListDialog(BaseDialog):
    """
    Dialog for displaying molecule groups in the system.

    Shows molecule indices, types, atom counts, and center of mass positions.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    atoms : ase.Atoms
        ASE Atoms object containing molecule information.
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        atoms: Optional[Atoms] = None,
    ):
        self.atoms = atoms
        super().__init__(parent, title="Molecule List")
        # Adjust window size
        self.geometry("600x400")

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        if self.atoms is None:
            ttk.Label(parent, text="No system loaded.").pack(pady=20)
            return

        # Check for required data
        groups = self.atoms.info.get("groups", [])
        site_mol_name = self.atoms.info.get("site_mol_name", [])

        if not groups:
            ttk.Label(parent, text="No molecule groups found.").pack(pady=20)
            return

        # Summary frame
        summary_frame = ttk.Frame(parent)
        summary_frame.pack(fill=tk.X, pady=(0, 10))

        n_molecules = len(groups)
        n_atoms = len(self.atoms)
        ttk.Label(
            summary_frame,
            text=f"Total: {n_molecules} molecules, {n_atoms} atoms",
            font=("TkDefaultFont", 10, "bold"),
        ).pack(anchor=tk.W)

        # Create treeview for molecule list
        columns = ("index", "name", "atoms", "elements", "com_x", "com_y", "com_z")
        self.tree = ttk.Treeview(parent, columns=columns, show="headings", height=15)

        # Define headings
        self.tree.heading("index", text="Index")
        self.tree.heading("name", text="Molecule")
        self.tree.heading("atoms", text="Atoms")
        self.tree.heading("elements", text="Elements")
        self.tree.heading("com_x", text="COM X")
        self.tree.heading("com_y", text="COM Y")
        self.tree.heading("com_z", text="COM Z")

        # Define column widths
        self.tree.column("index", width=50, anchor=tk.CENTER)
        self.tree.column("name", width=100, anchor=tk.W)
        self.tree.column("atoms", width=50, anchor=tk.CENTER)
        self.tree.column("elements", width=120, anchor=tk.W)
        self.tree.column("com_x", width=80, anchor=tk.E)
        self.tree.column("com_y", width=80, anchor=tk.E)
        self.tree.column("com_z", width=80, anchor=tk.E)

        # Add scrollbar
        scrollbar = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)

        # Pack widgets
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Populate tree
        self._populate_tree(groups, site_mol_name)

    def _populate_tree(self, groups: List[List[int]], site_mol_name: List[str]) -> None:
        """Populate the treeview with molecule data."""
        positions = self.atoms.get_positions()
        symbols = self.atoms.get_chemical_symbols()
        n_atoms = len(self.atoms)

        for i, group in enumerate(groups):
            if not group:
                continue

            # Validate group indices
            valid_indices = [idx for idx in group if 0 <= idx < n_atoms]
            if not valid_indices:
                continue

            # Get molecule name
            mol_name = site_mol_name[i] if i < len(site_mol_name) else f"mol_{i}"

            # Get atom count
            atom_count = len(valid_indices)

            # Get unique elements
            mol_symbols = [symbols[idx] for idx in valid_indices]
            unique_elements = sorted(set(mol_symbols))
            elements_str = ", ".join(unique_elements)

            # Calculate center of mass (simple average)
            mol_positions = positions[valid_indices]
            com = mol_positions.mean(axis=0)

            # Insert row
            self.tree.insert(
                "",
                tk.END,
                values=(
                    i,
                    mol_name,
                    atom_count,
                    elements_str,
                    f"{com[0]:.3f}",
                    f"{com[1]:.3f}",
                    f"{com[2]:.3f}",
                ),
            )

    def build_buttons(self, parent: ttk.Frame) -> None:
        """Build buttons (only Close button for this dialog)."""
        close_btn = ttk.Button(parent, text="Close", command=self.on_cancel)
        close_btn.pack(side=tk.RIGHT)


def show_molecule_list_dialog(
    parent: Optional[tk.Widget] = None,
    atoms: Optional[Atoms] = None,
) -> None:
    """
    Show the molecule list dialog.

    Parameters
    ----------
    parent : tk.Widget, optional
        Parent widget.
    atoms : ase.Atoms, optional
        ASE Atoms object containing molecule information.
    """
    dialog = MoleculeListDialog(parent, atoms=atoms)
    dialog.wait_window()
