# source/gui/supercell_dialog.py
"""
Dialog for supercell configuration.

Allows setting supercell repeat parameters.
"""

import tkinter as tk
from tkinter import ttk
from typing import Optional, Tuple

import numpy as np
from ase import Atoms

from .base_dialog import BaseDialog, LabeledEntry


class SupercellDialog(BaseDialog):
    """
    Dialog for setting supercell repeat parameters.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    current_repeats : tuple[int, int, int], optional
        Current supercell repeat values.
    atoms : ase.Atoms, optional
        Current atoms object for cell information display.
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        current_repeats: Tuple[int, int, int] = (1, 1, 1),
        atoms: Optional[Atoms] = None,
    ):
        self.current_repeats = current_repeats
        self.atoms = atoms
        super().__init__(parent, title="Set Supercell Repeat")

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        # Current cell info
        if self.atoms is not None:
            cell = self.atoms.get_cell()
            a, b, c = np.linalg.norm(cell, axis=1)
            ttk.Label(
                parent,
                text=f"Current cell: a={a:.2f}, b={b:.2f}, c={c:.2f} (Ang)",
            ).pack(anchor=tk.W, pady=(0, 10))

        # Repeat entries
        ttk.Label(parent, text="Supercell repeat factors:").pack(anchor=tk.W)

        repeat_frame = ttk.Frame(parent)
        repeat_frame.pack(fill=tk.X, pady=5)

        self.nx_var = tk.StringVar(value=str(self.current_repeats[0]))
        self.ny_var = tk.StringVar(value=str(self.current_repeats[1]))
        self.nz_var = tk.StringVar(value=str(self.current_repeats[2]))

        for label, var in [("nx:", self.nx_var), ("ny:", self.ny_var), ("nz:", self.nz_var)]:
            ttk.Label(repeat_frame, text=label).pack(side=tk.LEFT)
            entry = ttk.Entry(repeat_frame, textvariable=var, width=5)
            entry.pack(side=tk.LEFT, padx=(0, 10))

        # Preview
        self.preview_label = ttk.Label(parent, text="")
        self.preview_label.pack(anchor=tk.W, pady=(10, 0))
        self._update_preview()

        # Bind changes
        for var in [self.nx_var, self.ny_var, self.nz_var]:
            var.trace_add("write", lambda *args: self._update_preview())

    def _update_preview(self) -> None:
        """Update the preview label."""
        try:
            nx = int(self.nx_var.get())
            ny = int(self.ny_var.get())
            nz = int(self.nz_var.get())
            total_cells = nx * ny * nz

            if self.atoms is not None:
                n_atoms = len(self.atoms)
                total_atoms = n_atoms * total_cells
                self.preview_label.config(
                    text=f"Result: {total_cells} cells, ~{total_atoms} atoms"
                )
            else:
                self.preview_label.config(text=f"Result: {total_cells} cells")
        except ValueError:
            self.preview_label.config(text="Invalid input")

    def validate(self) -> bool:
        """Validate dialog input."""
        try:
            nx = int(self.nx_var.get())
            ny = int(self.ny_var.get())
            nz = int(self.nz_var.get())

            if nx < 1 or ny < 1 or nz < 1:
                self.show_error("Error", "Repeat factors must be >= 1.")
                return False

            if nx * ny * nz > 1000:
                if not self.ask_yes_no(
                    "Warning",
                    f"This will create {nx * ny * nz} cells. Continue?"
                ):
                    return False

        except ValueError:
            self.show_error("Error", "Repeat factors must be integers.")
            return False

        return True

    def apply(self) -> None:
        """Apply dialog result."""
        self.result = {
            "supercell_repeat": [
                int(self.nx_var.get()),
                int(self.ny_var.get()),
                int(self.nz_var.get()),
            ]
        }


def show_supercell_dialog(
    parent: Optional[tk.Widget] = None,
    current_repeats: Tuple[int, int, int] = (1, 1, 1),
    atoms: Optional[Atoms] = None,
) -> Optional[dict]:
    """
    Show the supercell dialog and return the result.

    Parameters
    ----------
    parent : tk.Widget, optional
        Parent widget.
    current_repeats : tuple[int, int, int], optional
        Current supercell repeat values.
    atoms : ase.Atoms, optional
        Current atoms object for cell information display.

    Returns
    -------
    dict or None
        Dialog result with key 'supercell_repeat', or None if cancelled.
    """
    dialog = SupercellDialog(parent, current_repeats=current_repeats, atoms=atoms)
    dialog.wait_window()
    return dialog.result if not dialog.cancelled else None
