# source/gui/save_molecule_dialog.py
"""
Dialog for saving selected atoms as a molecule TOML file.

Allows users to select atoms in ASE-GUI, edit the selection,
and write a molecule definition to data/<name>.toml.
"""

import tkinter as tk
from tkinter import ttk, messagebox
from typing import Optional, List, Callable
from pathlib import Path

import numpy as np
from ase import Atoms

from .base_dialog import BaseDialog, LabeledEntry


def _format_toml(
    molecule_name: str,
    symbols: List[str],
    positions: np.ndarray,
    lattice_vectors: np.ndarray,
    num_molecules: int,
    center_of_mass: np.ndarray,
) -> str:
    """
    Format molecule data as TOML string.

    Parameters
    ----------
    molecule_name : str
        Name of the molecule.
    symbols : list[str]
        Element symbols for the molecule atoms.
    positions : numpy.ndarray
        Cartesian positions of the first molecule, shape (N, 3).
    lattice_vectors : numpy.ndarray
        Cell vectors, shape (3, 3).
    num_molecules : int
        Number of molecules in the unit cell.
    center_of_mass : numpy.ndarray
        Center of mass of each molecule, shape (num_molecules, 3).

    Returns
    -------
    str
        TOML-formatted string.
    """
    lines = []
    lines.append(f'molecule_name = "{molecule_name}"')
    lines.append("")

    # symbols
    lines.append("symbols = [")
    for i, sym in enumerate(symbols):
        comma = "," if i < len(symbols) - 1 else ","
        lines.append(f'  "{sym}"{comma}')
    lines.append("]")
    lines.append("")

    lines.append("use_fractional = false")
    lines.append("")

    # positions
    lines.append("# positions of first molecule in the cell")
    lines.append("positions = [")
    for i, pos in enumerate(positions):
        comma = "," if i < len(positions) - 1 else ","
        lines.append(f"  [{pos[0]:12.6f}, {pos[1]:12.6f}, {pos[2]:12.6f}]{comma}")
    lines.append("]")
    lines.append("")

    # lattice_vectors
    lines.append("lattice_vectors = [")
    labels = ["# a-vector", "# b-vector", "# c-vector"]
    for i, vec in enumerate(lattice_vectors):
        comma = "," if i < len(lattice_vectors) - 1 else ","
        lines.append(f"  [{vec[0]:12.6f}, {vec[1]:12.6f}, {vec[2]:12.6f}]{comma}  {labels[i]}")
    lines.append("]")
    lines.append("")

    # num_molecules
    lines.append("# number of molecules in the cell")
    lines.append(f"num_molecules = {num_molecules}")
    lines.append("")

    # center_of_mass
    lines.append("# center of mass of each molecule")
    lines.append("center_of_mass = [")
    for i, com in enumerate(center_of_mass):
        comma = "," if i < len(center_of_mass) - 1 else ","
        lines.append(f"  [{com[0]:12.6f}, {com[1]:12.6f}, {com[2]:12.6f}]{comma}")
    lines.append("]")

    return "\n".join(lines) + "\n"


class SaveMoleculeDialog(BaseDialog):
    """
    Dialog for saving selected atoms as a molecule TOML file.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    atoms : ase.Atoms
        Full atomic structure.
    selected_indices : list[int]
        Indices of selected atoms.
    highlight_callback : callable, optional
        Called with a list of atom indices to highlight in the main GUI.
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        atoms: Atoms,
        selected_indices: List[int],
        highlight_callback: Optional[Callable[[List[int]], None]] = None,
    ):
        self.atoms = atoms
        self.selected_indices = list(selected_indices)
        self._highlight_callback = highlight_callback
        self._original_indices = list(selected_indices)
        super().__init__(parent, title="Save Selection as Molecule")
        self.geometry("500x500")
        # Show initial highlight for all selected atoms
        self._notify_highlight(self.selected_indices)

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        # Molecule name
        self.name_entry = LabeledEntry(parent, "Molecule name:", default="")
        self.name_entry.pack(fill=tk.X, pady=2)
        self.name_entry.entry.focus_set()

        # Selected atoms editor
        sel_frame = ttk.LabelFrame(parent, text="Selected Atoms", padding=5)
        sel_frame.pack(fill=tk.BOTH, expand=True, pady=5)

        # Info label (updated dynamically)
        self.sel_info_label = ttk.Label(sel_frame, text="")
        self.sel_info_label.pack(anchor=tk.W)

        # Atom list (Treeview)
        tree_frame = ttk.Frame(sel_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, pady=(2, 0))

        columns = ("index", "element", "x", "y", "z")
        self.atom_tree = ttk.Treeview(
            tree_frame, columns=columns, show="headings", height=8,
            selectmode="extended",
        )
        self.atom_tree.heading("index", text="Index")
        self.atom_tree.heading("element", text="Elem")
        self.atom_tree.heading("x", text="X")
        self.atom_tree.heading("y", text="Y")
        self.atom_tree.heading("z", text="Z")
        self.atom_tree.column("index", width=50, anchor=tk.CENTER)
        self.atom_tree.column("element", width=50, anchor=tk.CENTER)
        self.atom_tree.column("x", width=90, anchor=tk.E)
        self.atom_tree.column("y", width=90, anchor=tk.E)
        self.atom_tree.column("z", width=90, anchor=tk.E)

        tree_scroll = ttk.Scrollbar(
            tree_frame, orient=tk.VERTICAL, command=self.atom_tree.yview
        )
        self.atom_tree.configure(yscrollcommand=tree_scroll.set)
        self.atom_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Add / Remove buttons
        btn_frame = ttk.Frame(sel_frame)
        btn_frame.pack(fill=tk.X, pady=(4, 0))

        self.add_index_var = tk.StringVar()
        add_entry = ttk.Entry(btn_frame, textvariable=self.add_index_var, width=12)
        add_entry.pack(side=tk.LEFT, padx=(0, 2))

        add_btn = ttk.Button(
            btn_frame, text="Add Index", command=self._on_add_atom
        )
        add_btn.pack(side=tk.LEFT, padx=(0, 8))

        remove_btn = ttk.Button(
            btn_frame, text="Remove Selected", command=self._on_remove_atoms
        )
        remove_btn.pack(side=tk.LEFT)

        self._refresh_atom_tree()

        # Bind Treeview selection to highlight atoms in the main GUI
        self.atom_tree.bind("<<TreeviewSelect>>", self._on_tree_select)

    def _notify_highlight(self, indices: List[int]) -> None:
        """Notify the main GUI to highlight the given atom indices."""
        if self._highlight_callback is not None:
            self._highlight_callback(indices)

    def _on_tree_select(self, event=None) -> None:
        """Highlight atoms corresponding to the Treeview selection in the main GUI."""
        selection = self.atom_tree.selection()
        if selection:
            indices = [
                int(self.atom_tree.item(item, "values")[0])
                for item in selection
            ]
        else:
            # Nothing focused in tree -> show all selected_indices
            indices = list(self.selected_indices)
        self._notify_highlight(indices)

    def _refresh_atom_tree(self) -> None:
        """Refresh the atom list Treeview and info label from self.selected_indices."""
        # Clear tree
        for item in self.atom_tree.get_children():
            self.atom_tree.delete(item)

        symbols = self.atoms.get_chemical_symbols()
        positions = self.atoms.get_positions()

        for idx in self.selected_indices:
            sym = symbols[idx]
            pos = positions[idx]
            self.atom_tree.insert(
                "", tk.END,
                values=(idx, sym, f"{pos[0]:.4f}", f"{pos[1]:.4f}", f"{pos[2]:.4f}"),
            )

        # Update info label
        sel_symbols = [symbols[i] for i in self.selected_indices]
        unique = sorted(set(sel_symbols))
        self.sel_info_label.config(
            text=f"{len(self.selected_indices)} atoms: {', '.join(unique) if unique else '(none)'}"
        )

        # Highlight all selected atoms in the main GUI
        self._notify_highlight(self.selected_indices)

    def _on_add_atom(self) -> None:
        """Add atom indices from the entry field."""
        text = self.add_index_var.get().strip()
        if not text:
            return

        n_atoms = len(self.atoms)
        added = []
        errors = []
        # Support comma/space separated, and ranges like "3-7"
        for token in text.replace(",", " ").split():
            if "-" in token and not token.startswith("-"):
                parts = token.split("-", 1)
                try:
                    start, end = int(parts[0]), int(parts[1])
                    for i in range(start, end + 1):
                        if 0 <= i < n_atoms and i not in self.selected_indices:
                            self.selected_indices.append(i)
                            added.append(i)
                except ValueError:
                    errors.append(token)
            else:
                try:
                    i = int(token)
                    if i < 0 or i >= n_atoms:
                        errors.append(f"{i} (out of range)")
                    elif i not in self.selected_indices:
                        self.selected_indices.append(i)
                        added.append(i)
                except ValueError:
                    errors.append(token)

        if errors:
            self.show_error("Error", f"Invalid indices: {', '.join(errors)}")

        if added:
            self.add_index_var.set("")
            self._refresh_atom_tree()

    def _on_remove_atoms(self) -> None:
        """Remove atoms selected in the Treeview."""
        selection = self.atom_tree.selection()
        if not selection:
            return

        remove_indices = set()
        for item in selection:
            values = self.atom_tree.item(item, "values")
            remove_indices.add(int(values[0]))

        self.selected_indices = [
            i for i in self.selected_indices if i not in remove_indices
        ]
        self._refresh_atom_tree()

    def on_cancel(self) -> None:
        """Handle Cancel: restore original selection highlight."""
        self._notify_highlight(self._original_indices)
        super().on_cancel()

    def validate(self) -> bool:
        """Validate dialog input."""
        if not self.selected_indices:
            self.show_error("Error", "No atoms in the selection.")
            return False

        name = self.name_entry.get().strip()
        if not name:
            self.show_error("Error", "Please enter a molecule name.")
            return False

        if not name.replace("_", "").replace("-", "").isalnum():
            self.show_error(
                "Error",
                "Molecule name must contain only letters, digits, hyphens, and underscores.",
            )
            return False

        # Check if file already exists
        data_dir = Path(__file__).resolve().parent.parent.parent / "data"
        toml_path = data_dir / f"{name}.toml"
        if toml_path.exists():
            if not self.ask_yes_no(
                "File Exists",
                f"{toml_path.name} already exists. Overwrite?",
            ):
                return False

        return True

    def apply(self) -> None:
        """Apply dialog result - write the TOML file."""
        name = self.name_entry.get().strip()

        positions = self.atoms.get_positions()
        symbols = self.atoms.get_chemical_symbols()

        # Molecule positions and symbols
        mol_symbols = [symbols[i] for i in self.selected_indices]
        mol_positions = positions[self.selected_indices]

        # Center of mass (using ASE)
        sub_atoms = self.atoms[self.selected_indices]
        com = sub_atoms.get_center_of_mass()
        com_array = np.array([com])

        # Lattice vectors from the cell
        cell = self.atoms.get_cell()
        lattice = np.array(cell)

        # Positions relative to the molecule's COM
        mol_positions_shifted = mol_positions - com

        # Generate TOML
        toml_str = _format_toml(
            molecule_name=name,
            symbols=mol_symbols,
            positions=mol_positions_shifted,
            lattice_vectors=lattice,
            num_molecules=1,
            center_of_mass=com_array,
        )

        # Write file
        data_dir = Path(__file__).resolve().parent.parent.parent / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        toml_path = data_dir / f"{name}.toml"
        toml_path.write_text(toml_str, encoding="utf-8")

        self.result = {
            "molecule_name": name,
            "toml_path": str(toml_path),
            "num_molecules": 1,
        }


def show_save_molecule_dialog(
    parent: Optional[tk.Widget] = None,
    atoms: Optional[Atoms] = None,
    selected_indices: Optional[List[int]] = None,
    highlight_callback: Optional[Callable[[List[int]], None]] = None,
) -> Optional[dict]:
    """
    Show the save molecule dialog and return the result.

    Parameters
    ----------
    parent : tk.Widget, optional
        Parent widget.
    atoms : ase.Atoms
        Full atomic structure.
    selected_indices : list[int]
        Indices of selected atoms.
    highlight_callback : callable, optional
        Called with a list of atom indices to highlight in the main GUI.

    Returns
    -------
    dict or None
        Dialog result with keys 'molecule_name', 'toml_path', 'num_molecules',
        or None if cancelled.
    """
    if atoms is None or selected_indices is None or len(selected_indices) == 0:
        messagebox.showerror("Error", "No atoms selected.")
        return None

    dialog = SaveMoleculeDialog(
        parent, atoms, selected_indices,
        highlight_callback=highlight_callback,
    )
    dialog.wait_window()
    return dialog.result if not dialog.cancelled else None
