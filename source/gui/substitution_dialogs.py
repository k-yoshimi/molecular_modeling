# source/gui/substitution_dialogs.py
"""
Dialogs for atom and molecule substitution operations.

Provides four types of substitution:
1. GlobalAtomSubstitutionDialog - Replace all atoms of one element with another
2. MoleculeSubstitutionDialog - Replace a molecule at a specific site
3. MoleculeAtomSubstitutionDialog - Replace atoms in a molecule with a new symbol list
4. AtomSubstitutionDialog - Replace a single atom at a specific position
"""

import tkinter as tk
from tkinter import ttk
from typing import List, Optional, Tuple

from ase import Atoms

from .base_dialog import BaseDialog, LabeledCombobox, LabeledEntry
from .molecule_registry import get_registry

# Common element symbols for quick selection (periodic table order)
COMMON_ELEMENTS = [
    # Period 1-2
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    # Period 3
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    # Period 4
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    # Period 5
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    # Period 6 (selected)
    "Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi",
    # Lanthanides (selected)
    "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
]


class GlobalAtomSubstitutionDialog(BaseDialog):
    """
    Dialog for global element substitution.

    Replaces all atoms of one element with another throughout the system.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    atoms : ase.Atoms
        ASE Atoms object to analyze for available elements.
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        atoms: Optional[Atoms] = None,
    ):
        self.atoms = atoms
        super().__init__(parent, title="Global Element Substitution")

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        # Get current elements from system
        if self.atoms is not None:
            current_elements = sorted(set(self.atoms.get_chemical_symbols()))
        else:
            current_elements = COMMON_ELEMENTS[:10]

        # Info label
        ttk.Label(
            parent,
            text=f"Current elements: {', '.join(current_elements)}",
        ).pack(anchor=tk.W, pady=(0, 10))

        # From element
        self.from_combo = LabeledCombobox(
            parent, "Replace:", current_elements
        )
        self.from_combo.pack(fill=tk.X, pady=2)

        # Arrow label
        ttk.Label(parent, text="    ->").pack(anchor=tk.W)

        # To element
        self.to_combo = LabeledCombobox(
            parent, "With:", COMMON_ELEMENTS
        )
        self.to_combo.pack(fill=tk.X, pady=2)

    def validate(self) -> bool:
        """Validate dialog input."""
        from_elem = self.from_combo.get()
        to_elem = self.to_combo.get()

        if not from_elem or not to_elem:
            self.show_error("Error", "Please select both elements.")
            return False

        if from_elem == to_elem:
            self.show_error("Error", "Source and target elements must be different.")
            return False

        return True

    def apply(self) -> None:
        """Apply dialog result."""
        self.result = {
            "old_symbol": self.from_combo.get(),
            "new_symbol": self.to_combo.get(),
        }


class MoleculeSubstitutionDialog(BaseDialog):
    """
    Dialog for molecule substitution at a specific site.

    Replaces a molecule at a given cell index and molecule index with another molecule.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    atoms : ase.Atoms
        ASE Atoms object containing molecule information.
    selected_mol_index : int, optional
        Pre-selected molecule index (from selected atom's mol_id).
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        atoms: Optional[Atoms] = None,
        selected_mol_index: Optional[int] = None,
    ):
        self.atoms = atoms
        self.selected_mol_index = selected_mol_index
        self.registry = get_registry()
        super().__init__(parent, title="Molecule Substitution")

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        # Cell index
        cell_frame = ttk.Frame(parent)
        cell_frame.pack(fill=tk.X, pady=2)

        ttk.Label(cell_frame, text="Cell index:", width=15).pack(side=tk.LEFT)

        self.cell_i = tk.StringVar(value="0")
        self.cell_j = tk.StringVar(value="0")
        self.cell_k = tk.StringVar(value="0")

        ttk.Label(cell_frame, text="i:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_i, width=5).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Label(cell_frame, text="j:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_j, width=5).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Label(cell_frame, text="k:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_k, width=5).pack(side=tk.LEFT)

        # Molecule index (use selected mol_id if available)
        default_mol_idx = str(self.selected_mol_index) if self.selected_mol_index is not None else "0"
        self.mol_index_entry = LabeledEntry(parent, "Molecule index:", default=default_mol_idx)
        self.mol_index_entry.pack(fill=tk.X, pady=2)

        # New molecule selection
        molecule_names = self.registry.get_molecule_names()
        if not molecule_names:
            molecule_names = ["(no molecules found)"]

        self.new_mol_combo = LabeledCombobox(
            parent, "New molecule:", molecule_names
        )
        self.new_mol_combo.pack(fill=tk.X, pady=2)

    def validate(self) -> bool:
        """Validate dialog input."""
        try:
            int(self.cell_i.get())
            int(self.cell_j.get())
            int(self.cell_k.get())
            int(self.mol_index_entry.get())
        except ValueError:
            self.show_error("Error", "Cell index and molecule index must be integers.")
            return False

        new_mol = self.new_mol_combo.get()
        if not new_mol or new_mol == "(no molecules found)":
            self.show_error("Error", "Please select a new molecule.")
            return False

        return True

    def apply(self) -> None:
        """Apply dialog result."""
        self.result = {
            "cell_index": [
                int(self.cell_i.get()),
                int(self.cell_j.get()),
                int(self.cell_k.get()),
            ],
            "molecule_index": int(self.mol_index_entry.get()),
            "new_molecule": self.new_mol_combo.get(),
        }


class MoleculeAtomSubstitutionDialog(BaseDialog):
    """
    Dialog for replacing atoms in a molecule with a new symbol list.

    Replaces the atomic symbols of a specific molecule with a new set.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    atoms : ase.Atoms
        ASE Atoms object containing molecule information.
    selected_mol_index : int, optional
        Pre-selected molecule index (from selected atom's mol_id).
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        atoms: Optional[Atoms] = None,
        selected_mol_index: Optional[int] = None,
    ):
        self.atoms = atoms
        self.selected_mol_index = selected_mol_index
        super().__init__(parent, title="Molecule Atom Substitution")

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        # Cell index
        cell_frame = ttk.Frame(parent)
        cell_frame.pack(fill=tk.X, pady=2)

        ttk.Label(cell_frame, text="Cell index:", width=15).pack(side=tk.LEFT)

        self.cell_i = tk.StringVar(value="0")
        self.cell_j = tk.StringVar(value="0")
        self.cell_k = tk.StringVar(value="0")

        ttk.Label(cell_frame, text="i:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_i, width=5).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Label(cell_frame, text="j:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_j, width=5).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Label(cell_frame, text="k:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_k, width=5).pack(side=tk.LEFT)

        # Molecule index (use selected mol_id if available)
        default_mol_idx = str(self.selected_mol_index) if self.selected_mol_index is not None else "0"
        self.mol_index_entry = LabeledEntry(parent, "Molecule index:", default=default_mol_idx)
        self.mol_index_entry.pack(fill=tk.X, pady=2)

        # New symbols entry
        ttk.Label(
            parent,
            text="New symbols (comma-separated or compact notation like 'Si6O6'):",
        ).pack(anchor=tk.W, pady=(10, 2))

        self.new_symbols_var = tk.StringVar(value="")
        self.new_symbols_entry = ttk.Entry(parent, textvariable=self.new_symbols_var)
        self.new_symbols_entry.pack(fill=tk.X, pady=2)

        # Help text
        ttk.Label(
            parent,
            text="Example: C6H6 or C,C,C,C,C,C,H,H,H,H,H,H",
            font=("TkDefaultFont", 9),
            foreground="gray",
        ).pack(anchor=tk.W)

    def validate(self) -> bool:
        """Validate dialog input."""
        try:
            int(self.cell_i.get())
            int(self.cell_j.get())
            int(self.cell_k.get())
            int(self.mol_index_entry.get())
        except ValueError:
            self.show_error("Error", "Cell index and molecule index must be integers.")
            return False

        new_symbols = self.new_symbols_var.get().strip()
        if not new_symbols:
            self.show_error("Error", "Please enter new symbols.")
            return False

        return True

    def apply(self) -> None:
        """Apply dialog result."""
        # Parse new symbols
        raw = self.new_symbols_var.get().strip()
        if "," in raw:
            # Comma-separated list
            new_symbols = [s.strip() for s in raw.split(",") if s.strip()]
        else:
            # Compact notation (e.g., "Si6O6")
            new_symbols = [raw]

        self.result = {
            "cell_index": [
                int(self.cell_i.get()),
                int(self.cell_j.get()),
                int(self.cell_k.get()),
            ],
            "molecule_index": int(self.mol_index_entry.get()),
            "new_symbols": new_symbols,
        }


class AtomSubstitutionDialog(BaseDialog):
    """
    Dialog for individual atom substitution.

    Replaces a single atom at a specific position.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    atoms : ase.Atoms
        ASE Atoms object containing molecule information.
    selected_atom : int, optional
        Index of a pre-selected atom.
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        atoms: Optional[Atoms] = None,
        selected_atom: Optional[int] = None,
    ):
        self.atoms = atoms
        self.selected_atom = selected_atom
        super().__init__(parent, title="Individual Atom Substitution")

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        # Cell index
        cell_frame = ttk.Frame(parent)
        cell_frame.pack(fill=tk.X, pady=2)

        ttk.Label(cell_frame, text="Cell index:", width=15).pack(side=tk.LEFT)

        self.cell_i = tk.StringVar(value="0")
        self.cell_j = tk.StringVar(value="0")
        self.cell_k = tk.StringVar(value="0")

        ttk.Label(cell_frame, text="i:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_i, width=5).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Label(cell_frame, text="j:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_j, width=5).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Label(cell_frame, text="k:").pack(side=tk.LEFT)
        ttk.Entry(cell_frame, textvariable=self.cell_k, width=5).pack(side=tk.LEFT)

        # Molecule index
        self.mol_index_entry = LabeledEntry(parent, "Molecule index:", default="0")
        self.mol_index_entry.pack(fill=tk.X, pady=2)

        # Atom index
        default_atom = str(self.selected_atom) if self.selected_atom is not None else "0"
        self.atom_index_entry = LabeledEntry(parent, "Atom index:", default=default_atom)
        self.atom_index_entry.pack(fill=tk.X, pady=2)

        # Show current atom info if available
        if self.atoms is not None and self.selected_atom is not None:
            try:
                current_symbol = self.atoms.get_chemical_symbols()[self.selected_atom]
                pos = self.atoms.get_positions()[self.selected_atom]
                ttk.Label(
                    parent,
                    text=f"Selected: {current_symbol} at ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})",
                    foreground="#66CCFF",  # Light blue, visible on dark backgrounds
                ).pack(anchor=tk.W, pady=(5, 0))
            except IndexError:
                pass

        # New symbol
        self.new_symbol_combo = LabeledCombobox(
            parent, "New symbol:", COMMON_ELEMENTS
        )
        self.new_symbol_combo.pack(fill=tk.X, pady=2)

    def validate(self) -> bool:
        """Validate dialog input."""
        try:
            int(self.cell_i.get())
            int(self.cell_j.get())
            int(self.cell_k.get())
            int(self.mol_index_entry.get())
            int(self.atom_index_entry.get())
        except ValueError:
            self.show_error("Error", "All indices must be integers.")
            return False

        new_symbol = self.new_symbol_combo.get()
        if not new_symbol:
            self.show_error("Error", "Please select a new element.")
            return False

        return True

    def apply(self) -> None:
        """Apply dialog result."""
        self.result = {
            "cell_index": [
                int(self.cell_i.get()),
                int(self.cell_j.get()),
                int(self.cell_k.get()),
            ],
            "molecule_index": int(self.mol_index_entry.get()),
            "atom_index": int(self.atom_index_entry.get()),
            "new_symbol": self.new_symbol_combo.get(),
        }


# Convenience functions for showing dialogs

def show_global_substitution_dialog(
    parent: Optional[tk.Widget] = None,
    atoms: Optional[Atoms] = None,
) -> Optional[dict]:
    """Show global element substitution dialog."""
    dialog = GlobalAtomSubstitutionDialog(parent, atoms=atoms)
    dialog.wait_window()
    return dialog.result if not dialog.cancelled else None


def show_molecule_substitution_dialog(
    parent: Optional[tk.Widget] = None,
    atoms: Optional[Atoms] = None,
    selected_mol_index: Optional[int] = None,
) -> Optional[dict]:
    """Show molecule substitution dialog."""
    dialog = MoleculeSubstitutionDialog(parent, atoms=atoms, selected_mol_index=selected_mol_index)
    dialog.wait_window()
    return dialog.result if not dialog.cancelled else None


def show_molecule_atom_substitution_dialog(
    parent: Optional[tk.Widget] = None,
    atoms: Optional[Atoms] = None,
    selected_mol_index: Optional[int] = None,
) -> Optional[dict]:
    """Show molecule atom substitution dialog."""
    dialog = MoleculeAtomSubstitutionDialog(parent, atoms=atoms, selected_mol_index=selected_mol_index)
    dialog.wait_window()
    return dialog.result if not dialog.cancelled else None


def show_atom_substitution_dialog(
    parent: Optional[tk.Widget] = None,
    atoms: Optional[Atoms] = None,
    selected_atom: Optional[int] = None,
) -> Optional[dict]:
    """Show individual atom substitution dialog."""
    dialog = AtomSubstitutionDialog(parent, atoms=atoms, selected_atom=selected_atom)
    dialog.wait_window()
    return dialog.result if not dialog.cancelled else None
