#!/usr/bin/env python3
# molecular_modeling_gui.py
"""
Molecular Modeling GUI - ASE-GUI extension for molecular modeling operations.

Extends ASE-GUI with menus for:
- MoleculeTools: Molecule addition and substitution
- Export: Structure file export (QE, VASP, CIF)

Usage:
    python molecular_modeling_gui.py [structure_file]
"""

import json
import os
import sys
import threading
import numpy as np

from ase.io import read
from ase.gui.gui import GUI
import ase.gui.ui as ui
from ase import Atoms

# i18n support
try:
    from ase.gui.i18n import _
except Exception:
    def _(s): return s

# tkinter imports
try:
    import tkinter as tk
    from tkinter import messagebox, ttk
except Exception:
    tk = None
    messagebox = None
    ttk = None

# Local GUI module imports
from source.gui.molecule_registry import get_registry
from source.gui.add_molecule_dialog import show_add_molecule_dialog
from source.gui.molecule_list_dialog import show_molecule_list_dialog
from source.gui.substitution_dialogs import (
    show_global_substitution_dialog,
    show_molecule_substitution_dialog,
    show_molecule_atom_substitution_dialog,
    show_atom_substitution_dialog,
)
from source.gui.export_dialog import show_export_dialog
from source.gui.save_molecule_dialog import show_save_molecule_dialog
from source.gui.supercell_dialog import show_supercell_dialog
from source.gui.controller import get_controller, reset_controller


def read_cif(filepath: str) -> Atoms:
    """Read a CIF file with fallback cleanup for files that break ASE's parser.

    Parameters
    ----------
    filepath : str
        Path to the CIF file.

    Returns
    -------
    ase.Atoms
        The structure read from the CIF file.
    """
    try:
        return read(filepath)
    except Exception:
        # ASE's CIF parser can fail on files with complex citation loops
        # (e.g. semicolon-delimited multi-line text fields in ICSD CIFs).
        # Retry after stripping non-structural citation data.
        return _read_cif_with_cleanup(filepath)


def _read_cif_with_cleanup(filepath: str) -> Atoms:
    """Read a CIF file after removing citation loops that break ASE's parser."""
    import tempfile
    import os

    with open(filepath) as f:
        lines = f.readlines()

    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        # Skip standalone _citation_* tags
        if line.strip().startswith('_citation'):
            i += 1
            continue
        # Skip loop_ blocks whose headers are all _citation*
        if line.strip() == 'loop_':
            j = i + 1
            while j < len(lines) and lines[j].strip().startswith('_citation'):
                j += 1
            if j > i + 1:
                # Citation loop detected — skip its data rows too
                in_text_field = False
                while j < len(lines):
                    l = lines[j].strip()
                    if in_text_field:
                        if l.startswith(';'):
                            in_text_field = False
                        j += 1
                        continue
                    if l.startswith(';'):
                        in_text_field = True
                        j += 1
                        continue
                    if (l.startswith('_') or l.startswith('loop_')
                            or l.startswith('data_') or l.startswith('#')):
                        break
                    j += 1
                i = j
                continue
        out.append(line)
        i += 1

    fd, tmp = tempfile.mkstemp(suffix='.cif')
    try:
        with os.fdopen(fd, 'w') as f:
            f.writelines(out)
        return read(tmp)
    finally:
        os.unlink(tmp)


class MolecularModelingGUI(GUI):
    """
    Extended ASE-GUI with molecular modeling capabilities.

    Provides menus for:
    - MoleculeTools: Add molecules, list molecules, substitutions
    - Export: Export to various simulation formats
    """

    # Preset window sizes: (label, width, height)
    WINDOW_SIZES = [
        ("Small (450x450)", 450, 450),
        ("Medium (700x700)", 700, 700),
        ("Large (1000x800)", 1000, 800),
        ("Wide (1200x600)", 1200, 600),
    ]

    MAX_UNDO = 50  # Maximum number of undo snapshots
    MAX_RECENT = 10  # Maximum number of recent files
    RECENT_FILE = os.path.join(
        os.path.expanduser("~"), ".molecular_modeling_recent.json"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.controller = get_controller()
        if self.atoms is not None:
            self.controller.set_atoms(self.atoms)
        self._repeat_id = None
        self._undo_stack = []
        self._redo_stack = []
        self._setup_toolbar()
        # Set a reasonable default window size (larger than ASE's 450x450)
        self._resize_window(700, 700)

    # =========================================================================
    # View toolbar
    # =========================================================================

    def _setup_toolbar(self):
        """Add translate/rotate toolbar above the status bar."""
        if tk is None:
            return

        toolbar = tk.Frame(self.window.win)
        toolbar.pack(side=tk.BOTTOM, fill=tk.X, before=self.window.status)

        # --- Translate buttons ---
        tk.Label(toolbar, text="Translate:").pack(side=tk.LEFT, padx=(4, 2))
        translate_buttons = [
            ("\u25C0", -1, 0),   # ◀
            ("\u25B6", +1, 0),   # ▶
            ("\u25B2", 0, +1),   # ▲
            ("\u25BC", 0, -1),   # ▼
        ]
        for label, dx, dy in translate_buttons:
            btn = tk.Button(toolbar, text=label, width=3)
            btn.pack(side=tk.LEFT, padx=1)
            action = lambda dx=dx, dy=dy: self._translate_view(dx, dy)
            btn.bind("<ButtonPress-1>", lambda e, a=action: self._start_repeat(a))
            btn.bind("<ButtonRelease-1>", self._stop_repeat)

        # --- Rotate buttons ---
        tk.Label(toolbar, text="  Rotate:").pack(side=tk.LEFT, padx=(8, 2))
        rotate_buttons = [
            ("X-", 'x', -5),
            ("X+", 'x', +5),
            ("Y-", 'y', -5),
            ("Y+", 'y', +5),
            ("Z-", 'z', -5),
            ("Z+", 'z', +5),
        ]
        for label, axis, angle in rotate_buttons:
            btn = tk.Button(toolbar, text=label, width=3)
            btn.pack(side=tk.LEFT, padx=1)
            action = lambda axis=axis, angle=angle: self._rotate_view(axis, angle)
            btn.bind("<ButtonPress-1>", lambda e, a=action: self._start_repeat(a))
            btn.bind("<ButtonRelease-1>", self._stop_repeat)

        # --- Zoom buttons ---
        tk.Label(toolbar, text="  Zoom:").pack(side=tk.LEFT, padx=(8, 2))
        for label, factor in [("-", 1 / 1.1), ("+", 1.1)]:
            btn = tk.Button(toolbar, text=label, width=3)
            btn.pack(side=tk.LEFT, padx=1)
            action = lambda f=factor: self._zoom_view(f)
            btn.bind("<ButtonPress-1>", lambda e, a=action: self._start_repeat(a))
            btn.bind("<ButtonRelease-1>", self._stop_repeat)

    def _translate_view(self, dx, dy):
        """Translate the view in screen coordinates. dx, dy are -1 or +1."""
        step = 20.0 / self.scale
        self.center -= np.dot(self.axes, np.array([dx * step, -dy * step, 0.0]))
        self.draw()

    def _zoom_view(self, factor):
        """Zoom the view by multiplying scale by factor."""
        self.scale *= factor
        self.draw()

    def _rotate_view(self, axis, angle_deg):
        """Rotate the view around the given axis ('x', 'y', 'z') by angle_deg."""
        angle = np.radians(angle_deg)
        c, s = np.cos(angle), np.sin(angle)
        if axis == 'x':
            R = [[1, 0, 0], [0, c, -s], [0, s, c]]
        elif axis == 'y':
            R = [[c, 0, s], [0, 1, 0], [-s, 0, c]]
        else:  # 'z'
            R = [[c, -s, 0], [s, c, 0], [0, 0, 1]]
        self.axes = np.dot(self.axes, np.array(R))
        self.draw()

    def _start_repeat(self, action):
        """Execute action and schedule repeat every 50ms."""
        action()
        self._repeat_id = self.window.win.after(50, self._start_repeat, action)

    def _stop_repeat(self, event=None):
        """Cancel the repeating action."""
        if self._repeat_id is not None:
            self.window.win.after_cancel(self._repeat_id)
            self._repeat_id = None

    # =========================================================================
    # Enhanced status bar
    # =========================================================================

    def status(self, atoms):
        """Override ASE status to show structure info and molecule details."""
        natoms = len(atoms)
        indices = np.arange(natoms)[self.images.selected[:natoms]]
        n = len(indices)

        if n == 0:
            # Show system overview
            parts = [f"{natoms} atoms"]
            # Element composition
            syms = atoms.get_chemical_symbols()
            unique = sorted(set(syms))
            if unique:
                parts.append("(" + ", ".join(unique) + ")")
            # Molecule groups
            if "mol_id" in atoms.arrays:
                mol_ids = atoms.get_array("mol_id")
                n_groups = len(set(mol_ids))
                parts.append(f"{n_groups} groups")
            # Cell info
            cell = atoms.get_cell()
            lengths = cell.lengths()
            if lengths.any():
                parts.append(
                    f"Cell: {lengths[0]:.1f}x{lengths[1]:.1f}x{lengths[2]:.1f} Å"
                )
            text = " | ".join(parts)
            self.window.update_status_line(text)
            return

        if n == 1:
            # Show atom details + molecule name
            Z = atoms.numbers[indices]
            R = atoms.positions[indices]
            from ase.data import atomic_names, chemical_symbols
            idx = indices[0]
            text = (f" #{idx} {atomic_names[Z[0]]} ({chemical_symbols[Z[0]]}): "
                    f"{R[0][0]:.3f}, {R[0][1]:.3f}, {R[0][2]:.3f} Å")
            # Add mol_id and molecule name
            if "mol_id" in atoms.arrays:
                mol_id = int(atoms.get_array("mol_id")[idx])
                text += f" | mol_id={mol_id}"
                site_names = atoms.info.get("site_mol_name", [])
                if mol_id < len(site_names):
                    text += f" ({site_names[mol_id]})"
            self.window.update_status_line(text)
            return

        from ase.data import chemical_symbols

        if n == 2:
            # Distance with MIC
            Z = atoms.numbers[indices]
            d = atoms.get_distance(int(indices[0]), int(indices[1]), mic=True)
            s0, s1 = chemical_symbols[Z[0]], chemical_symbols[Z[1]]
            text = f" #{indices[0]}({s0})-#{indices[1]}({s1}): {d:.3f} Å"
            # Add molecule info
            text += self._mol_info_for_indices(atoms, indices)
            self.window.update_status_line(text)
            # Log measurement
            self._log_measurement(f"Distance #{indices[0]}({s0})-#{indices[1]}({s1}): {d:.3f} Å")
            return

        if n == 3:
            # Angle with MIC
            Z = atoms.numbers[indices]
            angle = atoms.get_angle(int(indices[0]), int(indices[1]), int(indices[2]), mic=True)
            syms = [chemical_symbols[z] for z in Z]
            text = (f" #{indices[0]}({syms[0]})-#{indices[1]}({syms[1]})"
                    f"-#{indices[2]}({syms[2]}): {angle:.1f}°")
            text += self._mol_info_for_indices(atoms, indices)
            self.window.update_status_line(text)
            self._log_measurement(
                f"Angle #{indices[0]}({syms[0]})-#{indices[1]}({syms[1]})"
                f"-#{indices[2]}({syms[2]}): {angle:.1f}°"
            )
            return

        ordered_indices = [i for i in self.images.selected_ordered
                          if i < len(atoms)]
        if n == 4 and len(ordered_indices) == 4:
            # Dihedral angle with MIC
            Z = atoms.numbers[indices]
            dihedral = atoms.get_dihedral(*ordered_indices, mic=True)
            syms = [chemical_symbols[atoms.numbers[i]] for i in ordered_indices]
            text = (f" Dihedral #{ordered_indices[0]}({syms[0]})"
                    f"→#{ordered_indices[1]}({syms[1]})"
                    f"→#{ordered_indices[2]}({syms[2]})"
                    f"→#{ordered_indices[3]}({syms[3]}): {dihedral:.1f}°")
            self.window.update_status_line(text)
            self._log_measurement(
                f"Dihedral {syms[0]}→{syms[1]}→{syms[2]}→{syms[3]}: {dihedral:.1f}°"
            )
            return

        # Fallback for 5+ atoms
        super().status(atoms)

    def _mol_info_for_indices(self, atoms, indices):
        """Get molecule info string for a set of atom indices."""
        if "mol_id" not in atoms.arrays:
            return ""
        mol_ids = atoms.get_array("mol_id")
        site_names = atoms.info.get("site_mol_name", [])
        mols = set()
        for idx in indices:
            mid = int(mol_ids[idx])
            name = site_names[mid] if mid < len(site_names) else f"mol_{mid}"
            mols.add(name)
        return " | " + ", ".join(sorted(mols))

    def _log_measurement(self, text):
        """Add a measurement to the log."""
        if not hasattr(self, '_measurement_log'):
            self._measurement_log = []
        self._measurement_log.append(text)
        # Keep last 100 measurements
        if len(self._measurement_log) > 100:
            self._measurement_log = self._measurement_log[-100:]

    def show_measurement_log(self, key=None):
        """Show the measurement log dialog."""
        if tk is None:
            return

        log = getattr(self, '_measurement_log', [])

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Measurement Log")
        dialog.resizable(True, True)
        dialog.transient(self.window.win)

        frame = tk.Frame(dialog, padx=10, pady=10)
        frame.pack(fill=tk.BOTH, expand=True)

        text_widget = tk.Text(frame, width=60, height=20, wrap=tk.WORD)
        scrollbar = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=text_widget.yview)
        text_widget.configure(yscrollcommand=scrollbar.set)
        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        if log:
            for i, entry in enumerate(log, 1):
                text_widget.insert(tk.END, f"{i}. {entry}\n")
        else:
            text_widget.insert(tk.END, "No measurements recorded yet.\n"
                              "Select 2 atoms for distance, 3 for angle, "
                              "4 for dihedral angle.")
        text_widget.configure(state=tk.DISABLED)

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(fill=tk.X, pady=5)

        def clear_log():
            self._measurement_log = []
            text_widget.configure(state=tk.NORMAL)
            text_widget.delete("1.0", tk.END)
            text_widget.insert(tk.END, "Log cleared.\n")
            text_widget.configure(state=tk.DISABLED)

        tk.Button(btn_frame, text="Clear", command=clear_log).pack(side=tk.LEFT, padx=10)
        tk.Button(btn_frame, text="Close", command=dialog.destroy).pack(side=tk.RIGHT, padx=10)

    # =========================================================================
    # Selection helpers
    # =========================================================================

    def select_by_element(self, key=None):
        """Open dialog to select atoms by element symbol."""
        if tk is None or self.atoms is None:
            return

        from ase.data import chemical_symbols
        elements = sorted(set(chemical_symbols[z] for z in self.atoms.numbers))
        if not elements:
            return

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Select by Element")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)

        tk.Label(dialog, text="Select elements:", padx=10, pady=5).pack(anchor=tk.W)

        frame = tk.Frame(dialog, padx=10)
        frame.pack(fill=tk.BOTH)

        vars_ = {}
        for sym in elements:
            v = tk.BooleanVar(value=False)
            vars_[sym] = v
            tk.Checkbutton(frame, text=sym, variable=v).pack(anchor=tk.W)

        mode_var = tk.StringVar(value="replace")
        mode_frame = tk.LabelFrame(dialog, text="Mode", padx=10, pady=5)
        mode_frame.pack(fill=tk.X, padx=10, pady=5)
        tk.Radiobutton(mode_frame, text="Replace selection", variable=mode_var,
                       value="replace").pack(anchor=tk.W)
        tk.Radiobutton(mode_frame, text="Add to selection", variable=mode_var,
                       value="add").pack(anchor=tk.W)

        def apply():
            syms = self.atoms.get_chemical_symbols()
            selected = {s for s, v in vars_.items() if v.get()}
            mask = np.array([s in selected for s in syms])
            n = len(self.atoms)
            if mode_var.get() == "replace":
                self.images.selected[:n] = mask
            else:
                self.images.selected[:n] |= mask
            self.draw()
            dialog.destroy()

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        tk.Button(btn_frame, text="Select", command=apply).pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=tk.LEFT, padx=5)

    def select_by_molecule(self, key=None):
        """Open dialog to select atoms by molecule group."""
        if tk is None or self.atoms is None:
            return

        if "mol_id" not in self.atoms.arrays:
            if messagebox:
                messagebox.showinfo("Select by Molecule",
                                    "No molecule groups found (mol_id).")
            return

        mol_ids = self.atoms.get_array("mol_id")
        site_names = self.atoms.info.get("site_mol_name", [])
        groups = self.atoms.info.get("groups", [])
        n_groups = int(mol_ids.max()) + 1

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Select by Molecule")
        dialog.resizable(True, True)
        dialog.transient(self.window.win)

        tk.Label(dialog, text="Select molecule groups:", padx=10, pady=5).pack(anchor=tk.W)

        frame = tk.Frame(dialog, padx=10)
        frame.pack(fill=tk.BOTH, expand=True)

        canvas = tk.Canvas(frame, height=200)
        scroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=canvas.yview)
        inner = tk.Frame(canvas)
        inner.bind("<Configure>",
                    lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=inner, anchor="nw")
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

        vars_ = {}
        for i in range(n_groups):
            name = site_names[i] if i < len(site_names) else f"mol_{i}"
            n_atoms = len(groups[i]) if i < len(groups) else 0
            v = tk.BooleanVar(value=False)
            vars_[i] = v
            tk.Checkbutton(inner, text=f"#{i} {name} ({n_atoms} atoms)",
                          variable=v).pack(anchor=tk.W)

        mode_var = tk.StringVar(value="replace")
        mode_frame = tk.LabelFrame(dialog, text="Mode", padx=10, pady=5)
        mode_frame.pack(fill=tk.X, padx=10, pady=5)
        tk.Radiobutton(mode_frame, text="Replace selection", variable=mode_var,
                       value="replace").pack(anchor=tk.W)
        tk.Radiobutton(mode_frame, text="Add to selection", variable=mode_var,
                       value="add").pack(anchor=tk.W)

        def apply():
            selected_groups = {g for g, v in vars_.items() if v.get()}
            n = len(self.atoms)
            mask = np.array([int(mol_ids[i]) in selected_groups for i in range(n)])
            if mode_var.get() == "replace":
                self.images.selected[:n] = mask
            else:
                self.images.selected[:n] |= mask
            self.draw()
            dialog.destroy()

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        tk.Button(btn_frame, text="Select", command=apply).pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=tk.LEFT, padx=5)

    def select_by_sphere(self, key=None):
        """Open dialog to select atoms within a radius of a selected atom."""
        if tk is None or self.atoms is None:
            return

        n = len(self.atoms)
        indices = np.arange(n)[self.images.selected[:n]]
        if len(indices) == 0:
            if messagebox:
                messagebox.showinfo("Select by Sphere",
                                    "Select one atom as the center first.")
            return

        center_idx = int(indices[0])
        center_pos = self.atoms.positions[center_idx]

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Select by Sphere")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)

        tk.Label(dialog, text=f"Center: atom #{center_idx}", padx=10, pady=5).pack(anchor=tk.W)

        row = tk.Frame(dialog, padx=10)
        row.pack(fill=tk.X)
        tk.Label(row, text="Radius (Å):").pack(side=tk.LEFT)
        r_var = tk.DoubleVar(value=3.0)
        tk.Spinbox(row, from_=0.5, to=50.0, increment=0.5,
                   textvariable=r_var, width=8).pack(side=tk.LEFT, padx=5)

        mode_var = tk.StringVar(value="replace")
        mode_frame = tk.LabelFrame(dialog, text="Mode", padx=10, pady=5)
        mode_frame.pack(fill=tk.X, padx=10, pady=5)
        tk.Radiobutton(mode_frame, text="Replace selection", variable=mode_var,
                       value="replace").pack(anchor=tk.W)
        tk.Radiobutton(mode_frame, text="Add to selection", variable=mode_var,
                       value="add").pack(anchor=tk.W)

        def apply():
            radius = r_var.get()
            dists = np.linalg.norm(self.atoms.positions - center_pos, axis=1)
            mask = dists <= radius
            if mode_var.get() == "replace":
                self.images.selected[:n] = mask
            else:
                self.images.selected[:n] |= mask
            self.draw()
            dialog.destroy()

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        tk.Button(btn_frame, text="Select", command=apply).pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=tk.LEFT, padx=5)

    # =========================================================================
    # Delete atoms/molecules
    # =========================================================================

    def delete_selected_atoms(self, widget=None, data=None):
        """Override ASE's delete to add undo support and update molecule metadata."""
        nselected = sum(self.images.selected[:len(self.atoms)])
        if nselected == 0:
            return
        if messagebox and not messagebox.askyesno(
            "Delete atoms", f"Delete {nselected} selected atom(s)?",
        ):
            return
        self._save_undo_state()
        self._delete_selected_with_metadata()

    def delete_molecule(self, key=None):
        """Delete the entire molecule group of selected atoms."""
        if self.atoms is None:
            return
        n = len(self.atoms)
        indices = np.arange(n)[self.images.selected[:n]]
        if len(indices) == 0:
            if messagebox:
                messagebox.showinfo("Delete Molecule", "Select an atom first.")
            return

        if "mol_id" not in self.atoms.arrays:
            if messagebox:
                messagebox.showinfo("Delete Molecule",
                                    "No molecule groups found (mol_id).")
            return

        mol_ids = self.atoms.get_array("mol_id")
        selected_mols = set(int(mol_ids[i]) for i in indices)
        site_names = self.atoms.info.get("site_mol_name", [])

        names = []
        for mid in sorted(selected_mols):
            name = site_names[mid] if mid < len(site_names) else f"mol_{mid}"
            names.append(f"#{mid} ({name})")

        if messagebox and not messagebox.askyesno(
            "Delete Molecule",
            f"Delete molecule(s): {', '.join(names)}?",
        ):
            return

        self._save_undo_state()
        # Select all atoms in the molecule groups
        mask = np.array([int(mol_ids[i]) in selected_mols for i in range(n)])
        self.images.selected[:n] = mask
        self._delete_selected_with_metadata()

    def _delete_selected_with_metadata(self):
        """Delete selected atoms and rebuild molecule metadata."""
        n = len(self.atoms)
        mask = self.images.selected[:n]
        keep = ~mask

        # Rebuild atoms
        new_atoms = self.atoms[keep]

        # Rebuild molecule metadata if present
        if "mol_id" in self.atoms.arrays:
            old_mol_ids = self.atoms.get_array("mol_id")
            old_groups = self.atoms.info.get("groups", [])
            old_names = self.atoms.info.get("site_mol_name", [])

            # Map old indices to new indices
            old_to_new = {}
            new_idx = 0
            for old_idx in range(n):
                if keep[old_idx]:
                    old_to_new[old_idx] = new_idx
                    new_idx += 1

            # Rebuild groups, removing empty ones
            new_groups = []
            new_names = []
            new_mol_ids = np.zeros(len(new_atoms), dtype=int)
            new_group_id = 0

            for g_idx, group in enumerate(old_groups):
                new_group = [old_to_new[i] for i in group if i in old_to_new]
                if new_group:
                    new_groups.append(new_group)
                    name = old_names[g_idx] if g_idx < len(old_names) else f"mol_{g_idx}"
                    new_names.append(name)
                    for ai in new_group:
                        new_mol_ids[ai] = new_group_id
                    new_group_id += 1

            if "mol_id" in new_atoms.arrays:
                del new_atoms.arrays["mol_id"]
            new_atoms.new_array("mol_id", new_mol_ids)
            new_atoms.info["groups"] = new_groups
            new_atoms.info["site_mol_name"] = new_names

        self.atoms = new_atoms
        self.images[self.frame] = new_atoms
        self.images.selected[:] = False
        self.set_frame(self.frame)
        self.draw()

    # =========================================================================
    # Show/hide atoms
    # =========================================================================

    def toggle_hydrogen(self, key=None):
        """Toggle visibility of hydrogen atoms."""
        if self.atoms is None:
            return
        h_mask = np.array(self.atoms.numbers == 1)
        if not h_mask.any():
            return
        # If any H is visible, hide all; otherwise show all
        n = len(self.atoms)
        visible = self.images.visible[:n]
        if visible[h_mask].any():
            self.images.visible[:n][h_mask] = False
        else:
            self.images.visible[:n][h_mask] = True
        self.draw()

    def toggle_element_visibility(self, key=None):
        """Open dialog to toggle visibility by element."""
        if tk is None or self.atoms is None:
            return

        from ase.data import chemical_symbols
        elements = sorted(set(chemical_symbols[z] for z in self.atoms.numbers))

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Element Visibility")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)

        tk.Label(dialog, text="Toggle element visibility:", padx=10, pady=5).pack(anchor=tk.W)

        n = len(self.atoms)
        visible = self.images.visible[:n]
        syms = self.atoms.get_chemical_symbols()

        frame = tk.Frame(dialog, padx=10)
        frame.pack(fill=tk.BOTH)

        vars_ = {}
        for sym in elements:
            # Check if any atom of this element is visible
            elem_mask = np.array([s == sym for s in syms])
            is_visible = visible[elem_mask].any()
            v = tk.BooleanVar(value=is_visible)
            vars_[sym] = v
            tk.Checkbutton(frame, text=sym, variable=v).pack(anchor=tk.W)

        def apply():
            for sym, v in vars_.items():
                elem_mask = np.array([s == sym for s in syms])
                self.images.visible[:n][elem_mask] = v.get()
            self.draw()
            dialog.destroy()

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        tk.Button(btn_frame, text="Apply", command=apply).pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=tk.LEFT, padx=5)

    # =========================================================================
    # Atom appearance panel
    # =========================================================================

    # Distinct colors for molecule groups (up to 20, then recycle)
    _GROUP_COLORS = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
        '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
    ]

    def show_appearance_panel(self, key=None):
        """Show the atom appearance settings panel."""
        if tk is None:
            return

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Atom Appearance")
        dialog.resizable(True, True)
        dialog.transient(self.window.win)
        dialog.geometry("420x500")

        notebook = ttk.Notebook(dialog)
        notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # --- Tab 1: Color Mode ---
        color_frame = tk.Frame(notebook, padx=10, pady=10)
        notebook.add(color_frame, text="Color Mode")

        color_mode_var = tk.StringVar(value="jmol")
        modes = [
            ("Default (by element)", "jmol"),
            ("By molecule group", "mol_group"),
            ("By tag", "tag"),
        ]
        for label, mode in modes:
            tk.Radiobutton(
                color_frame, text=label, variable=color_mode_var,
                value=mode, command=lambda: self._apply_color_mode(color_mode_var.get()),
            ).pack(anchor=tk.W, pady=2)

        # Custom array modes (if atoms have extra arrays)
        if self.atoms is not None:
            skip = {'numbers', 'positions', 'forces', 'momenta',
                    'initial_charges', 'initial_magmoms', 'tags'}
            for key_name in self.atoms.arrays:
                if key_name not in skip and key_name != 'mol_id':
                    tk.Radiobutton(
                        color_frame, text=f"By '{key_name}'",
                        variable=color_mode_var, value=key_name,
                        command=lambda: self._apply_color_mode(color_mode_var.get()),
                    ).pack(anchor=tk.W, pady=2)

        # --- Tab 2: Element Colors ---
        elem_frame = tk.Frame(notebook, padx=10, pady=10)
        notebook.add(elem_frame, text="Element Colors")

        canvas = tk.Canvas(elem_frame)
        scroll = ttk.Scrollbar(elem_frame, orient=tk.VERTICAL, command=canvas.yview)
        inner = tk.Frame(canvas)

        inner.bind("<Configure>",
                    lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=inner, anchor="nw")
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Get unique elements in current structure
        elements = []
        if self.atoms is not None:
            from ase.data import chemical_symbols
            seen = set()
            for z in self.atoms.numbers:
                sym = chemical_symbols[z]
                if sym not in seen:
                    seen.add(sym)
                    elements.append((z, sym))
            elements.sort(key=lambda x: x[0])

        self._elem_color_labels = {}
        for z, sym in elements:
            row = tk.Frame(inner)
            row.pack(fill=tk.X, pady=1)
            tk.Label(row, text=f"{sym} (Z={z})", width=12, anchor=tk.W).pack(side=tk.LEFT)
            color_hex = self.colors.get(z, '#808080')
            lbl = tk.Label(row, text="  ", bg=color_hex, width=4, relief=tk.SUNKEN)
            lbl.pack(side=tk.LEFT, padx=5)
            self._elem_color_labels[z] = lbl
            btn = tk.Button(
                row, text="Change...",
                command=lambda z=z, lbl=lbl: self._pick_element_color(z, lbl),
            )
            btn.pack(side=tk.LEFT)

        # Reset button
        tk.Button(
            elem_frame, text="Reset to Default",
            command=lambda: self._reset_element_colors(),
        ).pack(pady=5)

        # --- Tab 3: Size ---
        size_frame = tk.Frame(notebook, padx=10, pady=10)
        notebook.add(size_frame, text="Size")

        tk.Label(size_frame, text="Atom Size Scale:").pack(anchor=tk.W, pady=(0, 5))

        current_scale = getattr(self.images, 'atom_scale', 0.89)
        scale_var = tk.DoubleVar(value=current_scale)

        scale_label = tk.Label(size_frame, text=f"{current_scale:.2f}")
        scale_label.pack(anchor=tk.W)

        def update_scale(val):
            v = float(val)
            self.images.atom_scale = v
            scale_label.configure(text=f"{v:.2f}")
            self.draw()

        scale_slider = tk.Scale(
            size_frame, from_=0.1, to=3.0, resolution=0.05,
            orient=tk.HORIZONTAL, variable=scale_var,
            command=update_scale, length=300,
        )
        scale_slider.pack(fill=tk.X, pady=5)

        # Per-element size multipliers
        tk.Label(size_frame, text="\nPer-Element Radius Override:",
                 font=("", 10, "bold")).pack(anchor=tk.W, pady=(10, 5))

        size_inner = tk.Frame(size_frame)
        size_inner.pack(fill=tk.X)

        from ase.data import covalent_radii as default_radii
        for z, sym in elements:
            row = tk.Frame(size_inner)
            row.pack(fill=tk.X, pady=1)
            tk.Label(row, text=f"{sym}:", width=6, anchor=tk.W).pack(side=tk.LEFT)
            current_r = self.images.covalent_radii[z]
            default_r = default_radii[z]
            r_var = tk.DoubleVar(value=current_r)
            tk.Label(row, text=f"(default: {default_r:.2f} Å)", width=18,
                     anchor=tk.W).pack(side=tk.LEFT)
            spin = tk.Spinbox(
                row, from_=0.1, to=5.0, increment=0.05,
                textvariable=r_var, width=8,
                command=lambda z=z, v=r_var: self._set_element_radius(z, v),
            )
            spin.pack(side=tk.LEFT, padx=5)

        # Close button
        tk.Button(dialog, text="Close", command=dialog.destroy).pack(pady=5)

    def _apply_color_mode(self, mode):
        """Apply a color mode to the current view."""
        if mode == "mol_group":
            self._color_by_molecule_group()
        elif mode == "jmol":
            self.colormode = 'jmol'
            self.draw()
        elif mode == "tag":
            self.colormode = 'tag'
            # Set colormode_data for tag mode
            if self.atoms is not None:
                tags = self.atoms.get_tags()
                n_colors = max(len(self._GROUP_COLORS), int(tags.max()) + 1)
                colorscale = (self._GROUP_COLORS * ((n_colors // len(self._GROUP_COLORS)) + 1))[:n_colors]
                self.colormode_data = (colorscale, 0, n_colors - 1)
            self.draw()
        else:
            # Custom array mode
            if self.atoms is not None and mode in self.atoms.arrays:
                data = self.atoms.get_array(mode)
                mn, mx = float(data.min()), float(data.max())
                if mn == mx:
                    mx = mn + 1
                n_colors = len(self._GROUP_COLORS)
                self.colormode = mode
                self.colormode_data = (self._GROUP_COLORS[:n_colors], mn, mx)
                self.draw()

    def _color_by_molecule_group(self):
        """Color atoms by their molecule group (mol_id)."""
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            if messagebox:
                messagebox.showinfo(
                    "Color by Group",
                    "No molecule groups found (mol_id).\n"
                    "Load a structure built from TOML first.",
                )
            return

        mol_ids = self.atoms.get_array("mol_id")
        n_groups = int(mol_ids.max()) + 1
        n_colors = len(self._GROUP_COLORS)
        colorscale = (self._GROUP_COLORS * ((n_groups // n_colors) + 1))[:max(n_groups, 1)]
        self.colormode = 'mol_id'
        self.colormode_data = (colorscale, 0, max(n_groups - 1, 1))
        self.draw()

    def _pick_element_color(self, z, label_widget):
        """Open a color picker for the given element."""
        from tkinter import colorchooser
        current = self.colors.get(z, '#808080')
        result = colorchooser.askcolor(color=current, title=f"Choose color for Z={z}")
        if result and result[1]:
            hex_color = result[1]
            self.colors[z] = hex_color
            label_widget.configure(bg=hex_color)
            if self.colormode == 'jmol':
                self.draw()

    def _reset_element_colors(self):
        """Reset element colors to ASE default (jmol)."""
        from ase.data.colors import jmol_colors
        self.colors = {
            i: ('#{:02X}{:02X}{:02X}'.format(*(int(x * 255) for x in rgb)))
            for i, rgb in enumerate(jmol_colors)
        }
        self.colormode = 'jmol'
        self.draw()
        # Update labels if panel is open
        for z, lbl in getattr(self, '_elem_color_labels', {}).items():
            try:
                lbl.configure(bg=self.colors.get(z, '#808080'))
            except tk.TclError:
                pass

    def _set_element_radius(self, z, var):
        """Set the covalent radius for an element."""
        try:
            r = var.get()
            self.images.covalent_radii[z] = r
            self.draw()
        except (tk.TclError, ValueError):
            pass

    # =========================================================================
    # Recent files
    # =========================================================================

    def _load_recent_files(self):
        """Load the recent files list from disk."""
        try:
            with open(self.RECENT_FILE) as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return []

    def _save_recent_files(self, files):
        """Save the recent files list to disk."""
        try:
            with open(self.RECENT_FILE, "w") as f:
                json.dump(files[:self.MAX_RECENT], f)
        except OSError:
            pass

    def _add_recent_file(self, filepath):
        """Add a file to the recent files list."""
        filepath = os.path.abspath(filepath)
        files = self._load_recent_files()
        # Remove if already present, then prepend
        files = [f for f in files if f != filepath]
        files.insert(0, filepath)
        self._save_recent_files(files)

    def _open_recent_file(self, filepath):
        """Open a recent file."""
        def callback(key=None):
            if not os.path.exists(filepath):
                if messagebox:
                    messagebox.showerror("Error", f"File not found:\n{filepath}")
                return
            try:
                atoms = read_cif(filepath)
            except Exception as e:
                if messagebox:
                    messagebox.showerror("Error", f"Failed to read:\n{e}")
                return
            self._save_undo_state()
            self.images.initialize([atoms])
            self.set_frame()
            self.draw()
            self.controller.set_atoms(self.atoms)
            self._add_recent_file(filepath)
        return callback

    # =========================================================================
    # Undo / Redo
    # =========================================================================

    def _save_undo_state(self):
        """Save the current atoms state to the undo stack."""
        if self.atoms is None:
            return
        snapshot = self.atoms.copy()
        snapshot.info = dict(self.atoms.info)
        # Copy arrays (e.g. mol_id)
        for key in self.atoms.arrays:
            if key not in ('numbers', 'positions'):
                snapshot.set_array(key, self.atoms.get_array(key).copy())
        self._undo_stack.append(snapshot)
        if len(self._undo_stack) > self.MAX_UNDO:
            self._undo_stack.pop(0)
        self._redo_stack.clear()

    def _restore_atoms(self, snapshot):
        """Restore atoms from a snapshot and update the GUI."""
        self.atoms.positions[:] = snapshot.positions
        self.atoms.numbers[:] = snapshot.numbers
        self.atoms.set_cell(snapshot.get_cell())
        self.atoms.set_pbc(snapshot.get_pbc())
        self.atoms.info = dict(snapshot.info)
        # Restore arrays
        for key in snapshot.arrays:
            if key not in ('numbers', 'positions'):
                self.atoms.set_array(key, snapshot.get_array(key).copy())
        # Handle size change (atoms added/removed)
        if len(snapshot) != len(self.atoms):
            self.images.initialize([snapshot.copy()])
            # Copy info and arrays to the new atoms object
            self.atoms = self.images[0]
            self.atoms.info = dict(snapshot.info)
            for key in snapshot.arrays:
                if key not in ('numbers', 'positions'):
                    self.atoms.set_array(key, snapshot.get_array(key).copy())
        self.controller.set_atoms(self.atoms)
        self.images.initialize([self.atoms])
        self.set_frame()
        self.draw()

    def undo(self, key=None):
        """Undo the last operation."""
        if not self._undo_stack:
            if messagebox:
                messagebox.showinfo("Undo", "Nothing to undo.")
            return
        # Save current state to redo stack
        if self.atoms is not None:
            current = self.atoms.copy()
            current.info = dict(self.atoms.info)
            for k in self.atoms.arrays:
                if k not in ('numbers', 'positions'):
                    current.set_array(k, self.atoms.get_array(k).copy())
            self._redo_stack.append(current)
        snapshot = self._undo_stack.pop()
        self._restore_atoms(snapshot)

    def redo(self, key=None):
        """Redo the last undone operation."""
        if not self._redo_stack:
            if messagebox:
                messagebox.showinfo("Redo", "Nothing to redo.")
            return
        # Save current state to undo stack
        if self.atoms is not None:
            current = self.atoms.copy()
            current.info = dict(self.atoms.info)
            for k in self.atoms.arrays:
                if k not in ('numbers', 'positions'):
                    current.set_array(k, self.atoms.get_array(k).copy())
            self._undo_stack.append(current)
        snapshot = self._redo_stack.pop()
        self._restore_atoms(snapshot)

    # =========================================================================
    # Window size control
    # =========================================================================

    def _resize_window(self, width, height):
        """Resize the GUI window to the given width and height in pixels."""
        win = self.window.win
        win.geometry(f"{width}x{height}")
        # Update the internal canvas size tracking so ASE's scale adjusts
        self.window.size[:] = [width, height]

    def _resize_window_preset(self, width, height):
        """Resize to a preset size (used as menu callback)."""
        def callback(key=None):
            self._resize_window(width, height)
        return callback

    def _resize_window_custom(self, key=None):
        """Show a dialog to enter custom window size."""
        if tk is None:
            return

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Window Size")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)
        dialog.grab_set()

        frame = tk.Frame(dialog, padx=10, pady=10)
        frame.pack()

        # Get current window size
        current_w = self.window.win.winfo_width()
        current_h = self.window.win.winfo_height()

        tk.Label(frame, text="Width:").grid(row=0, column=0, sticky=tk.W, pady=2)
        width_var = tk.StringVar(value=str(current_w))
        width_entry = tk.Entry(frame, textvariable=width_var, width=8)
        width_entry.grid(row=0, column=1, padx=4, pady=2)
        tk.Label(frame, text="px").grid(row=0, column=2, sticky=tk.W)

        tk.Label(frame, text="Height:").grid(row=1, column=0, sticky=tk.W, pady=2)
        height_var = tk.StringVar(value=str(current_h))
        height_entry = tk.Entry(frame, textvariable=height_var, width=8)
        height_entry.grid(row=1, column=1, padx=4, pady=2)
        tk.Label(frame, text="px").grid(row=1, column=2, sticky=tk.W)

        def apply_size():
            try:
                w = int(width_var.get())
                h = int(height_var.get())
                if w < 200 or h < 200:
                    if messagebox:
                        messagebox.showerror("Error", "Minimum size is 200x200.",
                                             parent=dialog)
                    return
                self._resize_window(w, h)
                dialog.destroy()
            except ValueError:
                if messagebox:
                    messagebox.showerror("Error", "Please enter valid integers.",
                                         parent=dialog)

        btn_frame = tk.Frame(frame)
        btn_frame.grid(row=2, column=0, columnspan=3, pady=(8, 0))
        tk.Button(btn_frame, text="OK", width=8, command=apply_size).pack(
            side=tk.LEFT, padx=4)
        tk.Button(btn_frame, text="Cancel", width=8,
                  command=dialog.destroy).pack(side=tk.LEFT, padx=4)

        width_entry.focus_set()
        width_entry.select_range(0, tk.END)
        dialog.bind('<Return>', lambda e: apply_size())
        dialog.bind('<Escape>', lambda e: dialog.destroy())

    # =========================================================================
    # MoleculeTools methods
    # =========================================================================

    # Supported file formats for import
    _IMPORT_FILETYPES = [
        ("All supported", "*.cif *.xyz *.in *.POSCAR *.vasp *.pdb *.gen"),
        ("CIF files", "*.cif"),
        ("XYZ files", "*.xyz"),
        ("QE input files", "*.in *.pwi"),
        ("VASP files", "*.POSCAR *.vasp"),
        ("PDB files", "*.pdb"),
        ("DFTB+ GEN files", "*.gen"),
        ("All files", "*.*"),
    ]

    def load_cif_dialog(self, key=None):
        """Load a CIF file via file dialog and display the structure."""
        from tkinter import filedialog

        root = self.window.win if hasattr(self, "window") else None
        filename = filedialog.askopenfilename(
            parent=root,
            title="Load CIF File",
            filetypes=[("CIF files", "*.cif"), ("All files", "*.*")],
        )
        if not filename:
            return

        try:
            atoms = read_cif(filename)
        except Exception as e:
            if messagebox:
                messagebox.showerror(
                    "Error", f"Failed to read CIF file:\n{e}"
                )
            return

        self._save_undo_state()
        self._add_recent_file(filename)
        self.images.initialize([atoms])
        self.set_frame()
        self.draw()
        self.controller.set_atoms(self.atoms)

    def load_structure_dialog(self, key=None):
        """Load a structure file in any ASE-supported format."""
        from tkinter import filedialog

        root = self.window.win if hasattr(self, "window") else None
        filename = filedialog.askopenfilename(
            parent=root,
            title="Load Structure File",
            filetypes=self._IMPORT_FILETYPES,
        )
        if not filename:
            return

        try:
            # Try CIF-specific reader first for .cif files
            if filename.lower().endswith('.cif'):
                atoms = read_cif(filename)
            else:
                atoms = read(filename)
        except Exception as e:
            if messagebox:
                messagebox.showerror(
                    "Error", f"Failed to read file:\n{e}"
                )
            return

        self._save_undo_state()
        self._add_recent_file(filename)
        self.images.initialize([atoms])
        self.set_frame()
        self.draw()
        self.controller.set_atoms(self.atoms)

        if messagebox:
            fmt = os.path.splitext(filename)[1]
            messagebox.showinfo(
                "Load Structure",
                f"Loaded {len(atoms)} atoms from {os.path.basename(filename)} ({fmt})",
            )

    def add_molecule_dialog(self, key=None):
        """Show add molecule dialog."""
        root = self.window.win if hasattr(self, "window") else None
        cell = self.atoms.get_cell() if self.atoms is not None else None

        result = show_add_molecule_dialog(parent=root, cell=cell)
        if result:
            self._save_undo_state()
            # Ensure controller has current atoms reference
            self.controller.set_atoms(self.atoms)

            add_mode = result.get("add_mode", "single")

            if add_mode == "unitcell":
                # Add all molecules from unit cell
                positions = result.get("positions", [])
                lattice = result.get("lattice")
                molecule_name = result["molecule_name"]

                # Set cell from molecule's lattice if current cell is default
                if lattice and self.atoms is not None:
                    current_cell = self.atoms.get_cell()
                    if np.allclose(current_cell, np.diag([20, 20, 20])):
                        self.atoms.set_cell(lattice)

                success_count = 0
                for pos in positions:
                    success = self.controller.add_molecule(
                        molecule_name=molecule_name,
                        position=pos,
                    )
                    if success:
                        success_count += 1

                if success_count > 0:
                    self.images.initialize([self.atoms])
                    self.set_frame()
                    self.draw()
                    if messagebox:
                        messagebox.showinfo(
                            "Success",
                            f"Added {success_count} {molecule_name} molecule(s)"
                        )
                else:
                    if messagebox:
                        messagebox.showerror("Error", "Failed to add molecules.")
            else:
                # Single molecule mode
                success = self.controller.add_molecule(
                    molecule_name=result["molecule_name"],
                    position=result["position"],
                )
                if success:
                    self.images.initialize([self.atoms])
                    self.set_frame()
                    self.draw()
                    if messagebox:
                        messagebox.showinfo(
                            "Success",
                            f"Added {result['molecule_name']} at {result['position']}"
                        )
                else:
                    if messagebox:
                        messagebox.showerror("Error", "Failed to add molecule.")

    def list_molecules_dialog(self, key=None):
        """Show molecule list dialog."""
        root = self.window.win if hasattr(self, "window") else None
        show_molecule_list_dialog(parent=root, atoms=self.atoms)

    def save_molecule_dialog(self, key=None):
        """Show save selection as molecule dialog."""
        if self.atoms is None or len(self.atoms) == 0:
            if messagebox:
                messagebox.showerror("Error", "No structure loaded.")
            return

        # Get selected atom indices
        selected_indices = []
        if hasattr(self, 'images') and self.images.selected.any():
            selected_indices = list(np.where(self.images.selected)[0])

        if not selected_indices:
            if messagebox:
                messagebox.showerror("Error", "No atoms selected. Select atoms first.")
            return

        root = self.window.win if hasattr(self, "window") else None

        def highlight_atoms(indices):
            """Update ASE GUI selection highlight."""
            self.images.selected[:] = False
            for i in indices:
                if 0 <= i < len(self.images.selected):
                    self.images.selected[i] = True
            self.draw()

        result = show_save_molecule_dialog(
            parent=root,
            atoms=self.atoms,
            selected_indices=selected_indices,
            highlight_callback=highlight_atoms,
        )
        if result:
            if messagebox:
                messagebox.showinfo(
                    "Success",
                    f"Saved molecule '{result['molecule_name']}' "
                    f"({result['num_molecules']} molecule(s)) to\n"
                    f"{result['toml_path']}",
                )

    def global_substitution_dialog(self, key=None):
        """Show global element substitution dialog."""
        root = self.window.win if hasattr(self, "window") else None
        result = show_global_substitution_dialog(parent=root, atoms=self.atoms)
        if result:
            self._save_undo_state()
            self.controller.set_atoms(self.atoms)
            count = self.controller.global_atom_substitution(
                old_symbol=result["old_symbol"],
                new_symbol=result["new_symbol"],
            )
            self.set_frame()
            self.draw()
            if messagebox:
                messagebox.showinfo(
                    "Substitution",
                    f"Replaced {count} {result['old_symbol']} -> {result['new_symbol']}"
                )

    def molecule_substitution_dialog(self, key=None):
        """Show molecule substitution dialog."""
        root = self.window.win if hasattr(self, "window") else None

        # Get mol_id from selected atom if any
        selected_mol_index = None
        if hasattr(self, 'images') and self.images.selected.any():
            selected_indices = np.where(self.images.selected)[0]
            if len(selected_indices) >= 1 and self.atoms is not None:
                if "mol_id" in self.atoms.arrays:
                    mol_ids = self.atoms.get_array("mol_id")
                    selected_mol_index = int(mol_ids[selected_indices[0]])

        result = show_molecule_substitution_dialog(
            parent=root, atoms=self.atoms, selected_mol_index=selected_mol_index
        )
        if result:
            self._save_undo_state()
            # Use molecule_index directly (cell_index is for TOML-based systems)
            mol_idx = result["molecule_index"]
            new_molecule = result["new_molecule"]

            self.controller.set_atoms(self.atoms)
            success = self.controller.substitute_molecule(
                molecule_index=mol_idx,
                new_molecule_name=new_molecule,
            )

            if success:
                # Get the new atoms object from controller and reinitialize
                new_atoms = self.controller.get_atoms()
                self.images.initialize([new_atoms])
                self.set_frame()
                self.draw()
                if messagebox:
                    messagebox.showinfo(
                        "Success",
                        f"Molecule {mol_idx} replaced with {new_molecule}"
                    )
            else:
                if messagebox:
                    messagebox.showerror(
                        "Error",
                        f"Failed to substitute molecule {mol_idx} with {new_molecule}"
                    )

    def molecule_atom_substitution_dialog(self, key=None):
        """Show molecule atom substitution dialog."""
        root = self.window.win if hasattr(self, "window") else None

        # Get mol_id from selected atom if any
        selected_mol_index = None
        if hasattr(self, 'images') and self.images.selected.any():
            selected_indices = np.where(self.images.selected)[0]
            if len(selected_indices) >= 1 and self.atoms is not None:
                if "mol_id" in self.atoms.arrays:
                    mol_ids = self.atoms.get_array("mol_id")
                    # Use the mol_id of the first selected atom
                    selected_mol_index = int(mol_ids[selected_indices[0]])

        result = show_molecule_atom_substitution_dialog(
            parent=root, atoms=self.atoms, selected_mol_index=selected_mol_index
        )
        if result:
            self._save_undo_state()
            mol_idx = result["molecule_index"]
            self.controller.set_atoms(self.atoms)
            success, error_msg = self.controller.substitute_molecule_atoms(
                molecule_index=mol_idx,
                new_symbols=result["new_symbols"],
            )
            if success:
                self.set_frame()
                self.draw()
                if messagebox:
                    messagebox.showinfo("Success", "Molecule atoms substituted.")
            else:
                if messagebox:
                    messagebox.showerror("Error", error_msg)

    def atom_substitution_dialog(self, key=None):
        """Show individual atom substitution dialog."""
        root = self.window.win if hasattr(self, "window") else None

        # Get selected atom if any
        selected = None
        if hasattr(self, 'images') and self.images.selected.any():
            selected_indices = np.where(self.images.selected)[0]
            if len(selected_indices) == 1:
                selected = int(selected_indices[0])

        result = show_atom_substitution_dialog(
            parent=root,
            atoms=self.atoms,
            selected_atom=selected,
        )
        if result:
            self._save_undo_state()
            # For simplicity, use atom_index directly
            atom_idx = result["atom_index"]
            self.controller.set_atoms(self.atoms)
            success = self.controller.substitute_atom(
                atom_index=atom_idx,
                new_symbol=result["new_symbol"],
            )
            if success:
                self.set_frame()
                self.draw()
                if messagebox:
                    messagebox.showinfo("Success", f"Atom {atom_idx} -> {result['new_symbol']}")
            else:
                if messagebox:
                    messagebox.showerror("Error", "Failed to substitute atom.")

    def supercell_dialog(self, key=None):
        """Show supercell builder dialog and apply."""
        if self.atoms is None or len(self.atoms) == 0:
            if messagebox:
                messagebox.showerror("Error", "No structure loaded.")
            return

        root = self.window.win if hasattr(self, "window") else None
        result = show_supercell_dialog(
            parent=root, current_repeats=(1, 1, 1), atoms=self.atoms,
        )
        if result:
            repeat = result["supercell_repeat"]
            self._save_undo_state()
            from ase.build import make_supercell
            P = np.diag(repeat)
            new_atoms = make_supercell(self.atoms, P)
            # Replicate mol_id array if present
            if "mol_id" in self.atoms.arrays:
                old_ids = self.atoms.get_array("mol_id")
                n_orig = len(self.atoms)
                n_cells = repeat[0] * repeat[1] * repeat[2]
                max_id = old_ids.max() + 1 if len(old_ids) > 0 else 0
                new_ids = []
                for cell_i in range(n_cells):
                    new_ids.extend(old_ids + cell_i * max_id)
                if len(new_ids) == len(new_atoms):
                    new_atoms.set_array("mol_id", np.array(new_ids, dtype=int))
            self.images.initialize([new_atoms])
            self.set_frame()
            self.draw()
            self.controller.set_atoms(self.atoms)
            if messagebox:
                messagebox.showinfo(
                    "Supercell",
                    f"Built {repeat[0]}x{repeat[1]}x{repeat[2]} supercell\n"
                    f"({len(new_atoms)} atoms)",
                )

    # =========================================================================
    # Export methods
    # =========================================================================

    def export_all_dialog(self, key=None):
        """Show export dialog for all formats."""
        root = self.window.win if hasattr(self, "window") else None
        result = show_export_dialog(parent=root, atoms=self.atoms)
        # Dialog handles its own success/error messages

    def export_qe(self, key=None):
        """Export to Quantum ESPRESSO format."""
        self._export_single_format("qe")

    def export_vasp(self, key=None):
        """Export to VASP format."""
        self._export_single_format("vasp")

    def export_cif(self, key=None):
        """Export to CIF format."""
        self._export_single_format("cif")

    def _export_single_format(self, fmt: str):
        """Export to a single format."""
        from tkinter import filedialog
        from pathlib import Path
        from ase.io import write

        if self.atoms is None:
            if messagebox:
                messagebox.showerror("Error", "No structure to export.")
            return

        root = self.window.win if hasattr(self, "window") else None

        if fmt == "qe":
            filename = filedialog.asksaveasfilename(
                parent=root,
                title="Save Quantum ESPRESSO Input",
                defaultextension=".in",
                filetypes=[("QE Input", "*.in"), ("All files", "*.*")],
            )
            if filename:
                unique_species = sorted(set(self.atoms.get_chemical_symbols()))
                pseudopotentials = {elem: f"{elem}.pbe-nc.UPF" for elem in unique_species}
                input_data = {
                    'control': {'calculation': 'scf', 'prefix': 'molecule', 'outdir': './tmp'},
                    'system': {'ecutwfc': 40, 'ecutrho': 320},
                    'electrons': {'conv_thr': 1e-8},
                }
                write(filename, self.atoms, format='espresso-in',
                      pseudopotentials=pseudopotentials, input_data=input_data, kpoints=(1, 1, 1))
                if messagebox:
                    messagebox.showinfo("Export", f"Saved: {filename}")

        elif fmt == "vasp":
            filename = filedialog.asksaveasfilename(
                parent=root,
                title="Save VASP POSCAR",
                defaultextension="",
                initialfile="POSCAR",
                filetypes=[("POSCAR", "POSCAR*"), ("All files", "*.*")],
            )
            if filename:
                write(filename, self.atoms, format='vasp', direct=True)
                if messagebox:
                    messagebox.showinfo("Export", f"Saved: {filename}")

        elif fmt == "cif":
            filename = filedialog.asksaveasfilename(
                parent=root,
                title="Save CIF",
                defaultextension=".cif",
                filetypes=[("CIF files", "*.cif"), ("All files", "*.*")],
            )
            if filename:
                # Create a copy to avoid modifying the original
                atoms_copy = self.atoms.copy()
                # Disable PBC to prevent atom wrapping when re-reading
                # This preserves molecular structures that span cell boundaries
                atoms_copy.set_pbc(False)
                write(filename, atoms_copy, format='cif')
                if messagebox:
                    messagebox.showinfo("Export", f"Saved: {filename}")

    def show_voronoi_results(self, key=None):
        """Show Voronoi optimization results in a dialog."""
        if tk is None:
            return
        if self.atoms is None:
            if messagebox:
                messagebox.showerror("Error", "No structure loaded.")
            return

        iou = self.atoms.info.get("optimized_iou")
        angles = self.atoms.info.get("optimized_angles_deg")
        scales = self.atoms.info.get("optimized_scales")

        if iou is None:
            if messagebox:
                messagebox.showinfo(
                    "Voronoi Results",
                    "No optimization results found.\n"
                    "Run Pipeline → Run Optimize (Voronoi) first.",
                )
            return

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Voronoi Optimization Results")
        dialog.resizable(True, True)
        dialog.transient(self.window.win)

        frame = tk.Frame(dialog, padx=10, pady=10)
        frame.pack(fill=tk.BOTH, expand=True)

        # Header
        site_names = self.atoms.info.get("site_mol_name", [])

        # Create treeview table
        columns = ("site", "molecule", "iou", "alpha", "beta", "gamma", "scale")
        tree = ttk.Treeview(frame, columns=columns, show="headings", height=10)
        tree.heading("site", text="Site")
        tree.heading("molecule", text="Molecule")
        tree.heading("iou", text="IoU")
        tree.heading("alpha", text="α (°)")
        tree.heading("beta", text="β (°)")
        tree.heading("gamma", text="γ (°)")
        tree.heading("scale", text="Scale")

        tree.column("site", width=50, anchor=tk.CENTER)
        tree.column("molecule", width=100, anchor=tk.W)
        tree.column("iou", width=80, anchor=tk.CENTER)
        tree.column("alpha", width=80, anchor=tk.CENTER)
        tree.column("beta", width=80, anchor=tk.CENTER)
        tree.column("gamma", width=80, anchor=tk.CENTER)
        tree.column("scale", width=80, anchor=tk.CENTER)

        for i in range(len(iou)):
            mol_name = site_names[i] if i < len(site_names) else "?"
            ang = angles[i] if angles is not None else [0, 0, 0]
            sc = scales[i] if scales is not None else 1.0
            tree.insert("", tk.END, values=(
                i, mol_name, f"{iou[i]:.4f}",
                f"{ang[0]:.1f}", f"{ang[1]:.1f}", f"{ang[2]:.1f}",
                f"{sc:.4f}",
            ))

        scrollbar = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=tree.yview)
        tree.configure(yscrollcommand=scrollbar.set)
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Summary
        avg_iou = np.mean(iou)
        summary = tk.Label(
            dialog,
            text=f"Average IoU: {avg_iou:.4f}  |  {len(iou)} sites",
            anchor=tk.W, padx=10, pady=5,
        )
        summary.pack(fill=tk.X)

        tk.Button(dialog, text="Close", command=dialog.destroy).pack(pady=5)

    # =========================================================================
    # Project save/load
    # =========================================================================

    def save_project(self, key=None):
        """Save project files (input.toml + data/*.toml) to a chosen directory."""
        if tk is None:
            return
        from tkinter import filedialog
        import shutil

        src_dir = os.getcwd()
        input_src = os.path.join(src_dir, "input", "input.toml")
        data_src = os.path.join(src_dir, "data")

        if not os.path.isfile(input_src):
            if messagebox:
                messagebox.showerror(
                    "Error",
                    "No input/input.toml found in the working directory.\n"
                    "Save is only available when a TOML project exists.",
                )
            return

        dest_dir = filedialog.askdirectory(
            title="Save Project To...",
            parent=self.window.win,
        )
        if not dest_dir:
            return

        try:
            # Create input/ and data/ directories
            dest_input = os.path.join(dest_dir, "input")
            dest_data = os.path.join(dest_dir, "data")
            os.makedirs(dest_input, exist_ok=True)
            os.makedirs(dest_data, exist_ok=True)

            # Copy input.toml
            shutil.copy2(input_src, os.path.join(dest_input, "input.toml"))

            # Copy data/*.toml files
            if os.path.isdir(data_src):
                for fname in os.listdir(data_src):
                    if fname.endswith(".toml"):
                        shutil.copy2(
                            os.path.join(data_src, fname),
                            os.path.join(dest_data, fname),
                        )

            if messagebox:
                messagebox.showinfo(
                    "Save Project",
                    f"Project saved to:\n{dest_dir}",
                )
        except Exception as e:
            if messagebox:
                messagebox.showerror("Error", f"Failed to save project:\n{e}")

    def load_project(self, key=None):
        """Load a project from a directory containing input/input.toml + data/*.toml."""
        if tk is None:
            return
        from tkinter import filedialog

        proj_dir = filedialog.askdirectory(
            title="Select Project Directory...",
            parent=self.window.win,
        )
        if not proj_dir:
            return

        input_file = os.path.join(proj_dir, "input", "input.toml")
        if not os.path.isfile(input_file):
            if messagebox:
                messagebox.showerror(
                    "Error",
                    f"No input/input.toml found in:\n{proj_dir}",
                )
            return

        try:
            # Temporarily change to project directory and build
            original_dir = os.getcwd()
            os.chdir(proj_dir)
            try:
                from source.build_system_from_toml import build_system_from_toml
                new_atoms = build_system_from_toml(verbose=True)
            finally:
                os.chdir(original_dir)

            self._save_undo_state()
            self.atoms = new_atoms
            self.images[self.frame] = new_atoms
            self.set_frame(self.frame)
            self.draw()

            if messagebox:
                messagebox.showinfo(
                    "Load Project",
                    f"Loaded project from:\n{proj_dir}\n"
                    f"{len(new_atoms)} atoms, "
                    f"{len(new_atoms.info.get('groups', []))} molecule sites",
                )
        except Exception as e:
            if messagebox:
                messagebox.showerror("Error", f"Failed to load project:\n{e}")

    # =========================================================================
    # Pipeline execution
    # =========================================================================

    def _run_pipeline_in_thread(self, task_name, task_func):
        """Run a pipeline task in a background thread with a progress dialog."""
        if tk is None:
            return

        if self.atoms is None or len(self.atoms) == 0:
            if messagebox:
                messagebox.showerror("Error", "No structure loaded.")
            return

        # Check mol_id exists
        if "mol_id" not in self.atoms.arrays:
            if messagebox:
                messagebox.showerror(
                    "Error",
                    "No molecule groups found (mol_id).\n"
                    "Load a structure built from TOML first.",
                )
            return

        # Create progress dialog
        progress = tk.Toplevel(self.window.win)
        progress.title(task_name)
        progress.resizable(False, False)
        progress.transient(self.window.win)
        progress.grab_set()
        progress.protocol("WM_DELETE_WINDOW", lambda: None)  # Prevent closing

        frame = tk.Frame(progress, padx=20, pady=20)
        frame.pack()
        status_label = tk.Label(frame, text=f"Running {task_name}...", width=40)
        status_label.pack(pady=(0, 10))

        result = {"error": None, "message": None}

        def run():
            try:
                msg = task_func()
                result["message"] = msg
            except Exception as e:
                result["error"] = str(e)

        def check_thread():
            if thread.is_alive():
                self.window.win.after(100, check_thread)
            else:
                progress.destroy()
                if result["error"]:
                    if messagebox:
                        messagebox.showerror("Error", result["error"])
                else:
                    # Update the view with potentially modified atoms
                    self.images.initialize([self.atoms])
                    self.set_frame()
                    self.draw()
                    if messagebox and result["message"]:
                        messagebox.showinfo("Complete", result["message"])

        thread = threading.Thread(target=run, daemon=True)
        thread.start()
        self.window.win.after(100, check_thread)

    def run_optimize(self, key=None):
        """Run Voronoi optimization on the current structure."""
        from source.build_molecule_groups import build_molecule_groups
        from source.optimize_voronoi_overlap import optimize_voronoi_overlap

        self._save_undo_state()

        def task():
            groups = build_molecule_groups(self.atoms)
            rngs = self.atoms.info.get("rng", {})
            rng_opt = int(rngs.get("optimize", rngs.get("global", 1)))
            self.atoms = optimize_voronoi_overlap(
                self.atoms, groups=groups, rng=rng_opt,
            )
            # Format results
            iou = self.atoms.info.get("optimized_iou", [])
            angles = self.atoms.info.get("optimized_angles_deg", [])
            lines = ["Voronoi optimization complete.\n"]
            for i, (val, ang) in enumerate(zip(iou, angles)):
                lines.append(
                    f"  Site {i}: IoU={val:.4f}, "
                    f"angles=({ang[0]:.1f}, {ang[1]:.1f}, {ang[2]:.1f})°"
                )
            return "\n".join(lines)

        self._run_pipeline_in_thread("Optimize Voronoi", task)

    def run_render(self, key=None):
        """Render Voronoi overlay visualizations."""
        from source.build_molecule_groups import build_molecule_groups
        from source.render_voronoi_overlays import render_voronoi_overlays
        from tkinter import filedialog

        root = self.window.win if hasattr(self, "window") else None
        out_dir = filedialog.askdirectory(
            parent=root,
            title="Select output directory for Voronoi PNGs",
            initialdir="./voronoi_pngs",
        )
        if not out_dir:
            return

        def task():
            groups = build_molecule_groups(self.atoms)
            rngs = self.atoms.info.get("rng", {})
            rng_ren = int(rngs.get("render", rngs.get("global", 1)))
            render_voronoi_overlays(
                self.atoms, groups=groups, out_dir=out_dir, rng=rng_ren,
            )
            return f"Voronoi PNGs saved to:\n{out_dir}"

        self._run_pipeline_in_thread("Render Voronoi", task)

    def run_export(self, key=None):
        """Export structure files (QE, VASP, CIF)."""
        from source.write_structure import write_structure
        from tkinter import filedialog

        root = self.window.win if hasattr(self, "window") else None
        out_dir = filedialog.askdirectory(
            parent=root,
            title="Select output directory for structure files",
            initialdir="./structure_out",
        )
        if not out_dir:
            return

        def task():
            write_structure(self.atoms, outdir=out_dir, verbose=True)
            return f"Structure files saved to:\n{out_dir}"

        self._run_pipeline_in_thread("Export Structures", task)

    def run_full_pipeline(self, key=None):
        """Run the full pipeline: Optimize → Render → Export."""
        from source.build_molecule_groups import build_molecule_groups
        from source.optimize_voronoi_overlap import optimize_voronoi_overlap
        from source.voronoi_geom import check_system_voronoi_volume_tiling
        from source.render_voronoi_overlays import render_voronoi_overlays
        from source.write_structure import write_structure
        from tkinter import filedialog

        root = self.window.win if hasattr(self, "window") else None
        out_dir = filedialog.askdirectory(
            parent=root,
            title="Select output directory",
            initialdir=".",
        )
        if not out_dir:
            return

        self._save_undo_state()

        def task():
            import os
            groups = build_molecule_groups(self.atoms)
            rngs = self.atoms.info.get("rng", {})
            rng_opt = int(rngs.get("optimize", rngs.get("global", 1)))
            rng_ren = int(rngs.get("render", rngs.get("global", 1)))

            # Optimize
            self.atoms = optimize_voronoi_overlap(
                self.atoms, groups=groups, rng=rng_opt,
            )

            # Check tiling
            check_system_voronoi_volume_tiling(
                self.atoms, groups, rel_tol=1e-3, abs_tol=1e-6, verbose=True,
            )

            # Render
            png_dir = os.path.join(out_dir, "voronoi_pngs")
            os.makedirs(png_dir, exist_ok=True)
            render_voronoi_overlays(
                self.atoms, groups=groups, out_dir=png_dir, rng=rng_ren,
            )

            # Export
            struct_dir = os.path.join(out_dir, "structure_out")
            os.makedirs(struct_dir, exist_ok=True)
            write_structure(self.atoms, outdir=struct_dir, verbose=True)

            # Format results
            iou = self.atoms.info.get("optimized_iou", [])
            lines = ["Full pipeline complete.\n"]
            for i, val in enumerate(iou):
                lines.append(f"  Site {i}: IoU={val:.4f}")
            lines.append(f"\nOutputs saved to:\n{out_dir}")
            return "\n".join(lines)

        self._run_pipeline_in_thread("Full Pipeline", task)

    # =========================================================================
    # Menu construction
    # =========================================================================

    def get_menu_data(self):
        """Build extended menu with MoleculeTools and Export."""
        menu = super().get_menu_data()
        M = ui.MenuItem

        # MoleculeTools menu
        menu.append((
            _('MoleculeTools'),
            [
                M(_('Undo'), self.undo, 'Ctrl+Z'),
                M(_('Redo'), self.redo, 'Ctrl+Y'),
                M(_('Load Structure...'),
                  self.load_structure_dialog, 'Ctrl+L'),
                M(_('Load CIF...'),
                  self.load_cif_dialog),
                M(_('Add Molecule...'),
                  self.add_molecule_dialog, 'Ctrl+J'),
                M(_('List Molecules/Groups'),
                  self.list_molecules_dialog),
                M(_('Save Selection as Molecule...'),
                  self.save_molecule_dialog),
                M(_('Global Element Substitution...'),
                  self.global_substitution_dialog),
                M(_('Molecule Substitution...'),
                  self.molecule_substitution_dialog),
                M(_('Molecule Atom Substitution...'),
                  self.molecule_atom_substitution_dialog),
                M(_('Individual Atom Substitution...'),
                  self.atom_substitution_dialog),
                M(_('Build Supercell...'),
                  self.supercell_dialog),
                M('---'),
                M(_('Save Project...'),
                  self.save_project, 'Ctrl+K'),
                M(_('Load Project...'),
                  self.load_project, 'Ctrl+T'),
                M('---'),
                M(_('Measurement Log...'),
                  self.show_measurement_log),
            ]
        ))

        # Recent Files menu
        recent = self._load_recent_files()
        if recent:
            recent_items = []
            for filepath in recent:
                label = os.path.basename(filepath)
                recent_items.append(M(_(label), self._open_recent_file(filepath)))
            menu.append((_('Recent Files'), recent_items))

        # Export menu
        menu.append((
            _('Export'),
            [
                M(_('Export All (QE, VASP, CIF)...'),
                  self.export_all_dialog, 'Ctrl+D'),
                M(_('Quantum ESPRESSO...'),
                  self.export_qe),
                M(_('VASP...'),
                  self.export_vasp),
                M(_('CIF...'),
                  self.export_cif),
            ]
        ))

        # Selection menu
        menu.append((
            _('Selection'),
            [
                M(_('Select by Element...'), self.select_by_element),
                M(_('Select by Molecule...'), self.select_by_molecule),
                M(_('Select by Sphere...'), self.select_by_sphere),
                M('---'),
                M(_('Delete Selected Atoms'), self.delete_selected_atoms),
                M(_('Delete Molecule...'), self.delete_molecule),
            ]
        ))

        # Pipeline menu
        menu.append((
            _('Pipeline'),
            [
                M(_('Run Optimize (Voronoi)'),
                  self.run_optimize, 'Ctrl+P'),
                M(_('Run Render (Voronoi PNGs)...'),
                  self.run_render),
                M(_('Run Export (QE, VASP, CIF)...'),
                  self.run_export),
                M(_('Run All (Optimize → Render → Export)...'),
                  self.run_full_pipeline),
                M(_('Show Optimization Results...'),
                  self.show_voronoi_results, 'Ctrl+H'),
            ]
        ))

        # View Settings menu
        view_items = [
            M(_('Atom Appearance...'), self.show_appearance_panel),
            M(_('Color by Molecule Group'), self._color_by_molecule_group),
            M(_('Reset Colors'), self._reset_element_colors),
            M('---'),
            M(_('Toggle Hydrogen'), self.toggle_hydrogen),
            M(_('Element Visibility...'), self.toggle_element_visibility),
            M('---'),
        ]
        for label, w, h in self.WINDOW_SIZES:
            view_items.append(M(_(label), self._resize_window_preset(w, h)))
        view_items.append(M(_('Custom Size...'), self._resize_window_custom))
        menu.append((_('View Settings'), view_items))

        return menu


def launch_gui(atoms: Atoms) -> Atoms:
    """
    Launch the molecular modeling GUI with a given atoms object.

    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object to visualize and edit.

    Returns
    -------
    ase.Atoms
        The modified atoms object after GUI is closed.

    Example
    -------
    >>> from source.build_system_from_toml import build_system_from_toml
    >>> from molecular_modeling_gui import launch_gui
    >>> system = build_system_from_toml(verbose=True)
    >>> modified_system = launch_gui(system)
    >>> # Use modified_system for further processing
    """
    # Reset controller to fresh state and set atoms
    controller = reset_controller()
    controller.set_atoms(atoms)

    # Launch GUI
    gui = MolecularModelingGUI(images=[atoms])
    gui.run()

    # Return the modified atoms after GUI closes
    return gui.atoms


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        # Create empty structure with reasonable default cell
        atoms = Atoms()
        atoms.set_cell([20, 20, 20])  # 20x20x20 Angstrom box
        atoms.set_pbc(True)
        # Initialize empty groups info
        atoms.info["groups"] = []
        atoms.info["site_mol_name"] = []
    else:
        atoms = read(sys.argv[1])

    launch_gui(atoms)


if __name__ == "__main__":
    main()
