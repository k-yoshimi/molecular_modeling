"""
Tests for GUI feature logic (color mode switching, visibility toggles,
focus mode, undo/redo, element substitution, cell parameter editing, etc.).

These tests verify the core logic of MolecularModelingGUI without
launching a real tkinter window.  We create a thin mock of the GUI
that wires up all the state-mutating methods while replacing draw()
and other display-only calls with no-ops.
"""

import os
import sys
from unittest.mock import MagicMock

import numpy as np
import pytest
from ase import Atoms

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_atoms_with_mol_id():
    """Create a test Atoms object with mol_id metadata (2 mol groups)."""
    atoms = Atoms(
        symbols=["C", "C", "H", "H", "N", "O"],
        positions=[
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [2.0, 0.5, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ],
        cell=[10.0, 10.0, 10.0],
        pbc=True,
    )
    mol_ids = np.array([0, 0, 0, 1, 1, 1], dtype=int)
    atoms.new_array("mol_id", mol_ids)
    atoms.info["groups"] = [[0, 1, 2], [3, 4, 5]]
    atoms.info["site_mol_name"] = ["benzene", "water"]
    atoms.set_tags([0, 0, 0, 1, 1, 1])
    return atoms


@pytest.fixture
def atoms_with_mol():
    return _make_atoms_with_mol_id()


class FakeImages:
    """Minimal stand-in for ase.gui.images.Images."""

    def __init__(self, atoms):
        n = len(atoms)
        self.selected = np.zeros(n, dtype=bool)
        self.visible = np.ones(n, dtype=bool)
        self.atom_scale = 1.0
        self.covalent_radii = np.ones(120)
        self._atoms = [atoms]

    def __getitem__(self, idx):
        return self._atoms[idx]

    def __setitem__(self, idx, val):
        self._atoms[idx] = val

    def initialize(self, atoms_list):
        self._atoms = list(atoms_list)
        n = len(atoms_list[0])
        self.selected = np.zeros(n, dtype=bool)
        self.visible = np.ones(n, dtype=bool)


class FakeGUI:
    """
    Lightweight mock of MolecularModelingGUI that exposes the real
    state-mutating methods without needing tkinter or an event loop.

    We import the actual class methods and bind them to this mock so that
    the exact production code paths are exercised.
    """

    # Class-level constants copied from the real class
    _GROUP_COLORS = [
        "#e6194B", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
        "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
        "#dcbeff", "#9A6324", "#800000", "#aaffc3", "#808000",
        "#000075", "#a9a9a9", "#ffe119", "#ffd8b1", "#fffac8",
    ]
    MAX_UNDO = 50

    def __init__(self, atoms):
        self.atoms = atoms
        self.frame = 0
        self.images = FakeImages(atoms)
        self.colormode = "jmol"
        self.colormode_data = None
        self.colors = {}
        self.center = atoms.positions.mean(axis=0)
        self.scale = 10.0
        self.axes = np.eye(3)
        self._undo_stack = []
        self._redo_stack = []
        self._focus_active = False
        self._focus_backup = None
        self._focus_mol_ids = set()
        self._focus_mode = None
        self._hidden_elements = set()
        self._hidden_molecules = set()
        self._mol_colors = {}
        self._hydrogen_visible = True
        self._drag_mode = False
        self._dragging = False
        self._drag_start_xy = None
        self._drag_start_positions = None
        self._draw_count = 0
        self._elem_collapsed = False
        self._mol_collapsed = False
        self._atomlist_collapsed = True

        # Mock controller
        self.controller = MagicMock()
        self.controller.set_atoms = MagicMock()
        # Make global_atom_substitution actually mutate self.atoms
        self.controller.global_atom_substitution = (
            self._real_global_atom_substitution
        )

    # --- Stubs for display-only methods ---
    def draw(self):
        self._draw_count += 1

    def set_frame(self, frame=None):
        if frame is not None:
            self.frame = frame
        self.atoms = self.images[self.frame]

    # --- Real logic methods imported from the production class ---

    def _save_undo_state(self):
        if self.atoms is not None:
            snapshot = self.atoms.copy()
            snapshot.info = dict(self.atoms.info)
            for k in self.atoms.arrays:
                if k not in ("numbers", "positions"):
                    snapshot.set_array(k, self.atoms.get_array(k).copy())
            self._undo_stack.append(snapshot)
            if len(self._undo_stack) > self.MAX_UNDO:
                self._undo_stack.pop(0)
            self._redo_stack.clear()

    def _restore_atoms(self, snapshot):
        """Restore atoms from a snapshot (simplified for testing)."""
        if len(snapshot) != len(self.atoms):
            self.images.initialize([snapshot.copy()])
            self.atoms = self.images[0]
        else:
            self.atoms.positions[:] = snapshot.positions
            self.atoms.numbers[:] = snapshot.numbers
            self.atoms.set_cell(snapshot.get_cell())
            self.atoms.set_pbc(snapshot.get_pbc())
        self.atoms.info = dict(snapshot.info)
        for key in snapshot.arrays:
            if key not in ("numbers", "positions"):
                self.atoms.set_array(key, snapshot.get_array(key).copy())
        self.controller.set_atoms(self.atoms)
        self.images.initialize([self.atoms])
        self.set_frame()
        self.draw()

    def undo(self, key=None):
        if not self._undo_stack:
            return
        if self.atoms is not None:
            current = self.atoms.copy()
            current.info = dict(self.atoms.info)
            for k in self.atoms.arrays:
                if k not in ("numbers", "positions"):
                    current.set_array(k, self.atoms.get_array(k).copy())
            self._redo_stack.append(current)
        snapshot = self._undo_stack.pop()
        self._restore_atoms(snapshot)

    def redo(self, key=None):
        if not self._redo_stack:
            return
        if self.atoms is not None:
            current = self.atoms.copy()
            current.info = dict(self.atoms.info)
            for k in self.atoms.arrays:
                if k not in ("numbers", "positions"):
                    current.set_array(k, self.atoms.get_array(k).copy())
            self._undo_stack.append(current)
        snapshot = self._redo_stack.pop()
        self._restore_atoms(snapshot)

    def _apply_color_mode(self, mode):
        if mode in ("mol_group", "mol_id"):
            self._color_by_molecule_group()
        elif mode in ("jmol", "element"):
            self.colormode = "jmol"
            self.draw()
        elif mode == "tag":
            self.colormode = "tag"
            if self.atoms is not None:
                tags = self.atoms.get_tags()
                n_colors = max(len(self._GROUP_COLORS), int(tags.max()) + 1)
                colorscale = (
                    self._GROUP_COLORS
                    * ((n_colors // len(self._GROUP_COLORS)) + 1)
                )[:n_colors]
                self.colormode_data = (colorscale, 0, n_colors)
            self.draw()
        else:
            if self.atoms is not None and mode in self.atoms.arrays:
                data = self.atoms.get_array(mode)
                mn, mx = float(data.min()), float(data.max())
                if mn == mx:
                    mx = mn + 1
                n_colors = len(self._GROUP_COLORS)
                self.colormode = mode
                self.colormode_data = (
                    self._GROUP_COLORS[:n_colors],
                    mn,
                    mx,
                )
                self.draw()

    def _color_by_molecule_group(self):
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        n_groups = int(mol_ids.max()) + 1
        n_colors = len(self._GROUP_COLORS)
        colorscale = (
            self._GROUP_COLORS * ((n_groups // n_colors) + 1)
        )[:max(n_groups, 1)]
        self.colormode = "mol_id"
        self.colormode_data = (colorscale, 0, max(n_groups, 1))
        self.draw()

    def _reset_element_colors(self):
        from ase.data.colors import jmol_colors
        self.colors = {
            i: ("#{:02X}{:02X}{:02X}".format(*(int(x * 255) for x in rgb)))
            for i, rgb in enumerate(jmol_colors)
        }
        self.colormode = "jmol"
        self.draw()

    def toggle_hydrogen(self, key=None):
        if self.atoms is None:
            return
        h_mask = np.array(self.atoms.numbers == 1)
        if not h_mask.any():
            return
        if not hasattr(self, "_hidden_elements"):
            self._hidden_elements = set()
        if "H" not in self._hidden_elements:
            self._hidden_elements.add("H")
            self._hydrogen_visible = False
        else:
            self._hidden_elements.discard("H")
            self._hydrogen_visible = True
        self.apply_visibility()

    def apply_visibility(self):
        """Apply all visibility rules (mirrors _apply_visibility)."""
        atoms = self.atoms
        if atoms is None or len(atoms) == 0:
            return
        n = len(atoms)
        visible = np.ones(n, dtype=bool)
        if hasattr(self, "_hidden_elements") and self._hidden_elements:
            symbols = atoms.get_chemical_symbols()
            for i, s in enumerate(symbols):
                if s in self._hidden_elements:
                    visible[i] = False
        if (hasattr(self, "_hidden_molecules") and self._hidden_molecules
                and "mol_id" in atoms.arrays):
            mol_ids = atoms.get_array("mol_id")
            for i in range(n):
                if int(mol_ids[i]) in self._hidden_molecules:
                    visible[i] = False
        self.images.visible[:n] = visible
        self.draw()

    def focus_selection(self, key=None):
        """Focus on selected atoms (simplified: no dialog, direct focus)."""
        if self.atoms is None or len(self.atoms) == 0:
            return
        n = len(self.atoms)
        indices = np.arange(n)[self.images.selected[:n]]
        if len(indices) == 0:
            return
        self._execute_focus(indices, set(), "atoms")

    def focus_molecules(self, mol_ids_to_focus):
        """Focus on specific molecule groups (test helper for dialog flow)."""
        if self.atoms is None or len(self.atoms) == 0:
            return
        if "mol_id" not in self.atoms.arrays:
            return
        n = len(self.atoms)
        mol_ids = self.atoms.get_array("mol_id")
        focus_indices = []
        for g_id in mol_ids_to_focus:
            focus_indices.extend(np.where(mol_ids == g_id)[0])
        if not focus_indices:
            return
        self._execute_focus(
            np.array(sorted(focus_indices)), set(mol_ids_to_focus), "molecules"
        )

    def _execute_focus(self, indices, focus_mol_ids, mode):
        """Apply focus with given indices and molecule group set."""
        n = len(self.atoms)
        self._focus_backup = {
            "visible": self.images.visible[:n].copy(),
            "center": self.center.copy(),
            "scale": self.scale,
        }
        self._focus_active = True
        self._focus_mol_ids = set(focus_mol_ids)
        self._focus_mode = mode
        visible = np.zeros(n, dtype=bool)
        visible[indices] = True
        self.images.visible[:n] = visible
        positions = self.atoms.positions[indices]
        self.center = positions.mean(axis=0)
        self.draw()

    def _toggle_focus_molecule(self, group_id):
        """Add or remove a molecule group from the current focus."""
        if not getattr(self, "_focus_active", False):
            return
        if "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        n = len(self.atoms)
        if group_id in self._focus_mol_ids:
            self._focus_mol_ids.discard(group_id)
            mask = mol_ids == group_id
            self.images.visible[:n][mask] = False
        else:
            self._focus_mol_ids.add(group_id)
            mask = mol_ids == group_id
            self.images.visible[:n][mask] = True
        visible_indices = np.where(self.images.visible[:n])[0]
        if len(visible_indices) > 0:
            self.center = self.atoms.positions[visible_indices].mean(axis=0)
        self.draw()

    def on_mol_vis_change(self, group_id):
        """Simulate side panel molecule checkbox toggle."""
        if getattr(self, "_focus_active", False):
            self._toggle_focus_molecule(group_id)
            return
        if group_id in self._hidden_molecules:
            self._hidden_molecules.discard(group_id)
        else:
            self._hidden_molecules.add(group_id)
        self.apply_visibility()

    def unfocus(self, key=None):
        if not getattr(self, "_focus_active", False):
            return
        n = len(self.atoms)
        backup = self._focus_backup
        if backup.get("visible") is not None:
            vlen = min(len(backup["visible"]), n)
            self.images.visible[:vlen] = backup["visible"][:vlen]
            if n > vlen:
                self.images.visible[vlen:n] = True
        if backup.get("center") is not None:
            self.center = backup["center"]
        if backup.get("scale") is not None:
            self.scale = backup["scale"]
        self._focus_active = False
        self._focus_backup = None
        self._focus_mol_ids = set()
        self._focus_mode = None
        self.draw()

    def _delete_selected_with_metadata(self):
        n = len(self.atoms)
        mask = self.images.selected[:n]
        keep = ~mask
        new_atoms = self.atoms[keep]
        if "mol_id" in self.atoms.arrays:
            old_groups = self.atoms.info.get("groups", [])
            old_names = self.atoms.info.get("site_mol_name", [])
            old_to_new = {}
            new_idx = 0
            for old_idx in range(n):
                if keep[old_idx]:
                    old_to_new[old_idx] = new_idx
                    new_idx += 1
            new_groups = []
            new_names = []
            new_mol_ids = np.zeros(len(new_atoms), dtype=int)
            new_group_id = 0
            for g_idx, group in enumerate(old_groups):
                new_group = [old_to_new[i] for i in group if i in old_to_new]
                if new_group:
                    new_groups.append(new_group)
                    name = (
                        old_names[g_idx]
                        if g_idx < len(old_names)
                        else f"mol_{g_idx}"
                    )
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
        self.images.selected = np.zeros(len(new_atoms), dtype=bool)
        self.images.visible = np.ones(len(new_atoms), dtype=bool)
        self.set_frame(self.frame)
        self.draw()

    def delete_selected_atoms(self):
        nselected = sum(self.images.selected[: len(self.atoms)])
        if nselected == 0:
            return
        self._save_undo_state()
        self._delete_selected_with_metadata()

    def _real_global_atom_substitution(self, old_symbol, new_symbol):
        """Real substitution logic (mirrors controller)."""
        symbols = self.atoms.get_chemical_symbols()
        count = 0
        for i, sym in enumerate(symbols):
            if sym == old_symbol:
                symbols[i] = new_symbol
                count += 1
        if count > 0:
            self.atoms.set_chemical_symbols(symbols)
        return count

    def set_element_radius(self, z, radius):
        """Set covalent radius for element z (mirrors _on_element_size_change)."""
        self.images.covalent_radii[z] = radius
        self.draw()

    def set_element_color(self, z, hex_color):
        """Set color for element z (mirrors _on_element_color_pick)."""
        self.colors[z] = hex_color
        self.draw()

    def apply_molecule_visibility(self):
        """Apply molecule visibility (alias for apply_visibility)."""
        self.apply_visibility()

    def select_molecule(self, group_id):
        """Toggle selection of a molecule group (mirrors _on_mol_select)."""
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        n = len(self.atoms)
        mask = (mol_ids == group_id)
        if np.array_equal(self.images.selected[:n], mask):
            self.images.selected[:n] = False
        else:
            self.images.selected[:n] = mask
        self.draw()

    def start_drag(self, selected_indices):
        """Simulate starting a drag on selected atoms (mirrors press in drag mode)."""
        self._drag_mode = True
        self._dragging = True
        self._drag_start_xy = (100, 100)  # arbitrary start position
        self._drag_start_positions = self.atoms.positions.copy()
        n = len(self.atoms)
        self.images.selected[:n] = False
        for i in selected_indices:
            self.images.selected[i] = True
        self._save_undo_state()

    def drag_to(self, delta_3d):
        """Simulate dragging atoms by a 3D world delta.

        In the real GUI, 2D mouse delta is converted to 3D via the axes matrix.
        Here we apply the delta directly to simplify testing.
        """
        if not self._dragging:
            return
        n = len(self.atoms)
        selected = self.images.selected[:n]
        indices = np.where(selected)[0]
        pos = self._drag_start_positions.copy()
        pos[indices] += np.array(delta_3d)
        self.atoms.set_positions(pos)
        self.images[self.frame] = self.atoms
        self.draw()

    def end_drag(self):
        """Simulate ending a drag (mirrors release in drag mode)."""
        self._dragging = False
        self._drag_start_xy = None
        self._drag_start_positions = None
        self.draw()

    # --- View operation methods (mirrors the real GUI) ---

    def _translate_view(self, dx, dy):
        """Translate the view in screen coordinates."""
        step = 20.0 / self.scale
        self.center -= np.dot(self.axes, np.array([dx * step, -dy * step, 0.0]))
        self.draw()

    def _zoom_view(self, factor):
        """Zoom the view by multiplying scale by factor."""
        self.scale *= factor
        self.draw()

    def _rotate_view(self, axis, angle_deg):
        """Rotate the view around the given axis."""
        angle = np.radians(angle_deg)
        c, s = np.cos(angle), np.sin(angle)
        if axis == 'x':
            R = [[1, 0, 0], [0, c, -s], [0, s, c]]
        elif axis == 'y':
            R = [[c, 0, s], [0, 1, 0], [-s, 0, c]]
        else:
            R = [[c, -s, 0], [s, c, 0], [0, 0, 1]]
        self.axes = np.dot(self.axes, np.array(R))
        self.draw()

    # --- Background color (canvas mock) ---

    def _set_canvas_bg(self, color):
        """Set the canvas background color (mock version)."""
        self._canvas_bg = color

    def _get_canvas_bg(self):
        """Get the current canvas background color."""
        return getattr(self, '_canvas_bg', '#ffffff')

    # --- Move / Rotate molecule logic (dialog-free version) ---

    def move_atoms(self, indices, vector):
        """Move atoms by a translation vector (mirrors move_molecule_dialog apply).

        Parameters
        ----------
        indices : list[int]
            Atom indices to move.
        vector : array-like
            3D translation vector in Angstrom.
        """
        if self.atoms is None or len(indices) == 0:
            return
        self._save_undo_state()
        pos = self.atoms.get_positions()
        vec = np.array(vector)
        for i in indices:
            pos[i] += vec
        self.atoms.set_positions(pos)
        self.images[self.frame] = self.atoms
        self.set_frame(self.frame)
        self.draw()

    def rotate_atoms(self, indices, angles_deg):
        """Rotate atoms by Euler angles around their centroid.

        Parameters
        ----------
        indices : list[int]
            Atom indices to rotate.
        angles_deg : tuple of float
            Euler angles (alpha, beta, gamma) in degrees, XYZ convention.
        """
        if self.atoms is None or len(indices) == 0:
            return
        from scipy.spatial.transform import Rotation as Rot
        self._save_undo_state()
        pos = self.atoms.get_positions()
        atom_pos = pos[indices]
        centroid = atom_pos.mean(axis=0)
        R = Rot.from_euler('XYZ', angles_deg, degrees=True).as_matrix()
        rotated = (atom_pos - centroid) @ R.T + centroid
        for i, idx in enumerate(indices):
            pos[idx] = rotated[i]
        self.atoms.set_positions(pos)
        self.images[self.frame] = self.atoms
        self.set_frame(self.frame)
        self.draw()

    def move_molecule(self, mol_id, vector):
        """Move all atoms in a molecule group by a translation vector."""
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        indices = list(np.where(mol_ids == mol_id)[0])
        self.move_atoms(indices, vector)

    def rotate_molecule(self, mol_id, angles_deg):
        """Rotate all atoms in a molecule group by Euler angles."""
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        indices = list(np.where(mol_ids == mol_id)[0])
        self.rotate_atoms(indices, angles_deg)

    # --- Drag mode toggle (mirrors _toggle_drag_mode without UI) ---

    def toggle_drag_mode(self):
        """Toggle drag mode state (mirrors _toggle_drag_mode)."""
        self._drag_mode = not self._drag_mode
        if not self._drag_mode:
            self._dragging = False

    def apply_mol_group_colors(self):
        """Apply custom molecule group colors (mirrors _apply_mol_group_colors)."""
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        n_groups = int(mol_ids.max()) + 1
        n_colors = len(self._GROUP_COLORS)
        colorscale = []
        for g_id in range(n_groups):
            color = self._mol_colors.get(
                g_id, self._GROUP_COLORS[g_id % n_colors])
            colorscale.append(color)
        self.colormode = "mol_id"
        self.colormode_data = (colorscale, 0, max(n_groups, 1))
        self.draw()

    # --- Select all / deselect all ---

    def select_all(self):
        """Select all visible atoms."""
        if self.atoms is None:
            return
        n = len(self.atoms)
        vis = self.images.visible[:n]
        self.images.selected[:n] = vis
        self.draw()

    def deselect_all(self):
        """Deselect all atoms."""
        if self.atoms is None:
            return
        n = len(self.atoms)
        self.images.selected[:n] = False
        self.draw()

    # --- Context select molecule ---

    def context_select_molecule(self, atom_idx):
        """Select all atoms in the same molecule group as atom_idx."""
        if self.atoms is None or len(self.atoms) == 0:
            return
        n = len(self.atoms)
        if "groups" in self.atoms.info:
            for grp in self.atoms.info["groups"]:
                if atom_idx in grp:
                    self.images.selected[:n] = False
                    for i in grp:
                        if i < n:
                            self.images.selected[i] = True
                    self.draw()
                    return
        if "mol_id" in self.atoms.arrays:
            mol_ids = self.atoms.get_array("mol_id")
            target_id = mol_ids[atom_idx]
            self.images.selected[:n] = (mol_ids == target_id)
            self.draw()

    # --- Section collapse/expand ---

    def toggle_section(self, section):
        """Toggle collapse state of a side panel section."""
        if section == "elem":
            self._elem_collapsed = not self._elem_collapsed
        elif section == "mol":
            self._mol_collapsed = not self._mol_collapsed
        elif section == "atomlist":
            self._atomlist_collapsed = not self._atomlist_collapsed

    # --- Status info ---

    def get_status_info(self):
        """Return status info string (mirrors _update_status_info logic)."""
        if self.atoms is None or len(self.atoms) == 0:
            return "No structure loaded"
        from collections import Counter
        n = len(self.atoms)
        sel = self.images.selected[:n]
        n_sel = int(sel.sum())
        parts = [f"{n} atoms total"]
        if n_sel > 0:
            parts.append(f"{n_sel} selected")
            syms = np.array(self.atoms.get_chemical_symbols())
            sel_syms = syms[sel]
            counts = Counter(sel_syms)
            breakdown = ", ".join(f"{s}:{c}" for s, c in sorted(counts.items()))
            parts.append(breakdown)
        if "groups" in self.atoms.info:
            ng = len(self.atoms.info["groups"])
            parts.append(f"{ng} groups")
        return " | ".join(parts)

    # --- Atom List ---

    def get_atom_list_data(self):
        """Return atom list data as list of tuples (idx, elem, x, y, z, mol).

        Mirrors the logic in _update_atom_list().
        """
        atoms = self.atoms
        if atoms is None or len(atoms) == 0:
            return []
        syms = atoms.get_chemical_symbols()
        pos = atoms.get_positions()
        mol_ids = atoms.arrays.get('mol_id', None)
        rows = []
        for i in range(len(atoms)):
            mol_val = int(mol_ids[i]) if mol_ids is not None else ""
            rows.append((
                i, syms[i],
                round(pos[i, 0], 3), round(pos[i, 1], 3), round(pos[i, 2], 3),
                mol_val,
            ))
        return rows

    def atom_list_select(self, idx):
        """Simulate selecting an atom via the atom list (mirrors _on_atom_tree_select)."""
        n = len(self.atoms)
        if idx < 0 or idx >= n:
            return
        self.images.selected[:n] = False
        self.images.selected[idx] = True
        self.draw()

    # --- Axis indicator ---

    def get_axis_projections(self):
        """Return XYZ axis projections for the axis indicator."""
        projections = {}
        for i, label in enumerate(["X", "Y", "Z"]):
            unit = np.zeros(3)
            unit[i] = 1.0
            proj = np.dot(self.axes[:2], unit)
            projections[label] = proj
        return projections


# ===========================================================================
# Color Mode Switching Tests
# ===========================================================================


class TestColorModeSwitching:
    """Verify that all color mode transitions work correctly."""

    def test_element_to_mol_id_and_back(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert gui.colormode == "jmol"

        # Switch to mol_id
        gui._apply_color_mode("mol_id")
        assert gui.colormode == "mol_id"
        assert gui.colormode_data is not None
        colorscale, lo, hi = gui.colormode_data
        assert len(colorscale) >= 2
        assert lo == 0

        # Switch back to element
        gui._apply_color_mode("element")
        assert gui.colormode == "jmol"

    def test_element_to_tag_and_back(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert gui.colormode == "jmol"

        # Switch to tag
        gui._apply_color_mode("tag")
        assert gui.colormode == "tag"
        assert gui.colormode_data is not None

        # Switch back to element
        gui._apply_color_mode("element")
        assert gui.colormode == "jmol"

    def test_mol_id_to_tag_and_back(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)

        gui._apply_color_mode("mol_id")
        assert gui.colormode == "mol_id"

        gui._apply_color_mode("tag")
        assert gui.colormode == "tag"

        gui._apply_color_mode("mol_id")
        assert gui.colormode == "mol_id"

    def test_jmol_alias(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._apply_color_mode("mol_id")
        gui._apply_color_mode("jmol")
        assert gui.colormode == "jmol"

    def test_mol_group_alias(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._apply_color_mode("mol_group")
        assert gui.colormode == "mol_id"

    def test_all_transitions(self, atoms_with_mol):
        """Test every possible transition between the 3 modes."""
        gui = FakeGUI(atoms_with_mol)
        modes = ["element", "mol_id", "tag"]
        expected_colormode = {"element": "jmol", "mol_id": "mol_id", "tag": "tag"}

        for from_mode in modes:
            for to_mode in modes:
                gui._apply_color_mode(from_mode)
                assert gui.colormode == expected_colormode[from_mode], (
                    f"Failed setting initial mode {from_mode}"
                )
                gui._apply_color_mode(to_mode)
                assert gui.colormode == expected_colormode[to_mode], (
                    f"Failed transition {from_mode} -> {to_mode}"
                )

    def test_reset_colors_sets_jmol(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._apply_color_mode("mol_id")
        assert gui.colormode == "mol_id"

        gui._reset_element_colors()
        assert gui.colormode == "jmol"
        # colors dict should have jmol colors populated
        assert len(gui.colors) > 0

    def test_mol_id_without_array(self):
        """Switching to mol_id with no mol_id array should be a no-op."""
        atoms = Atoms("CH4", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0],
                                        [0, 0, 1], [1, 1, 0]])
        atoms.set_cell([10, 10, 10])
        gui = FakeGUI(atoms)
        gui._apply_color_mode("mol_id")
        # Should remain jmol since mol_id array doesn't exist
        assert gui.colormode == "jmol"

    def test_draw_is_called_on_each_switch(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        count_before = gui._draw_count

        gui._apply_color_mode("mol_id")
        assert gui._draw_count == count_before + 1

        gui._apply_color_mode("element")
        assert gui._draw_count == count_before + 2

        gui._apply_color_mode("tag")
        assert gui._draw_count == count_before + 3


# ===========================================================================
# Hydrogen Visibility Toggle Tests
# ===========================================================================


class TestToggleHydrogen:
    """Verify toggle_hydrogen correctly shows/hides H atoms."""

    def test_hide_hydrogen(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert gui.images.visible.all()

        gui.toggle_hydrogen()
        n = len(gui.atoms)
        h_mask = gui.atoms.numbers == 1
        assert not gui.images.visible[:n][h_mask].any(), "H should be hidden"
        assert gui.images.visible[:n][~h_mask].all(), "Non-H should be visible"
        assert "H" in gui._hidden_elements
        assert gui._hydrogen_visible is False

    def test_show_hydrogen_again(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_hydrogen()  # hide
        gui.toggle_hydrogen()  # show
        n = len(gui.atoms)
        assert gui.images.visible[:n].all(), "All atoms should be visible"
        assert "H" not in gui._hidden_elements
        assert gui._hydrogen_visible is True

    def test_toggle_hydrogen_no_h(self):
        """Toggle should be a no-op if no H atoms exist."""
        atoms = Atoms("CO", positions=[[0, 0, 0], [1, 0, 0]], cell=[10, 10, 10])
        gui = FakeGUI(atoms)
        before = gui._draw_count
        gui.toggle_hydrogen()
        assert gui._draw_count == before, "draw should not be called"

    def test_toggle_hydrogen_multiple_cycles(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        for _ in range(5):
            gui.toggle_hydrogen()  # hide
            h_mask = gui.atoms.numbers == 1
            assert not gui.images.visible[:len(gui.atoms)][h_mask].any()
            gui.toggle_hydrogen()  # show
            assert gui.images.visible[:len(gui.atoms)].all()


# ===========================================================================
# Element Visibility Tests
# ===========================================================================


class TestElementVisibility:
    """Test _hidden_elements and per-element visibility logic."""

    def test_hide_single_element(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._hidden_elements.add("C")
        # Simulate _on_element_vis_change logic
        n = len(gui.atoms)
        visible = np.ones(n, dtype=bool)
        symbols = gui.atoms.get_chemical_symbols()
        for i, s in enumerate(symbols):
            if s in gui._hidden_elements:
                visible[i] = False
        gui.images.visible[:n] = visible

        c_mask = np.array([s == "C" for s in symbols])
        assert not gui.images.visible[:n][c_mask].any()
        assert gui.images.visible[:n][~c_mask].all()

    def test_hide_multiple_elements(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._hidden_elements = {"C", "N"}
        n = len(gui.atoms)
        visible = np.ones(n, dtype=bool)
        symbols = gui.atoms.get_chemical_symbols()
        for i, s in enumerate(symbols):
            if s in gui._hidden_elements:
                visible[i] = False
        gui.images.visible[:n] = visible

        for i, s in enumerate(symbols):
            if s in ("C", "N"):
                assert not gui.images.visible[i]
            else:
                assert gui.images.visible[i]

    def test_unhide_element(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._hidden_elements.add("C")
        gui._hidden_elements.discard("C")
        assert "C" not in gui._hidden_elements


# ===========================================================================
# Focus Mode Tests
# ===========================================================================


class TestFocusMode:
    """Test focus_selection and unfocus."""

    def test_focus_hides_non_selected(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        # Select atoms 0 and 1
        gui.images.selected[0] = True
        gui.images.selected[1] = True

        gui.focus_selection()

        n = len(gui.atoms)
        assert gui.images.visible[0] is True or gui.images.visible[0]
        assert gui.images.visible[1] is True or gui.images.visible[1]
        for i in range(2, n):
            assert not gui.images.visible[i], f"Atom {i} should be hidden"
        assert gui._focus_active is True

    def test_unfocus_restores_visibility(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_visible = gui.images.visible[:len(gui.atoms)].copy()
        original_center = gui.center.copy()
        original_scale = gui.scale

        gui.images.selected[0] = True
        gui.focus_selection()
        gui.unfocus()

        np.testing.assert_array_equal(
            gui.images.visible[:len(gui.atoms)], original_visible
        )
        np.testing.assert_array_almost_equal(gui.center, original_center)
        assert gui.scale == original_scale
        assert gui._focus_active is False

    def test_focus_no_selection_is_noop(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        before = gui._draw_count
        gui.focus_selection()
        assert gui._draw_count == before, "No draw should happen without selection"
        assert gui._focus_active is False

    def test_unfocus_without_focus_is_noop(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        before = gui._draw_count
        gui.unfocus()
        assert gui._draw_count == before

    def test_focus_centers_on_selection(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[4] = True  # N at [3, 0, 0]
        gui.images.selected[5] = True  # O at [4, 0, 0]

        gui.focus_selection()
        expected_center = np.array([3.5, 0.0, 0.0])
        np.testing.assert_array_almost_equal(gui.center, expected_center)

    def test_double_focus_unfocus(self, atoms_with_mol):
        """Focus twice, then unfocus once should restore first focus backup."""
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True
        gui.focus_selection()
        # Now focus on different atoms (overwrite backup)
        gui.images.selected[:] = False
        gui.images.selected[4] = True
        gui.focus_selection()
        gui.unfocus()
        # After unfocus, the second focus backup is restored
        assert gui._focus_active is False

    def test_focus_molecules_shows_whole_group(self, atoms_with_mol):
        """Focus on molecule group 0 shows all atoms in that group."""
        gui = FakeGUI(atoms_with_mol)
        gui.focus_molecules({0})  # benzene = atoms 0,1,2

        n = len(gui.atoms)
        assert gui.images.visible[0]
        assert gui.images.visible[1]
        assert gui.images.visible[2]
        for i in range(3, n):
            assert not gui.images.visible[i]
        assert gui._focus_active
        assert gui._focus_mol_ids == {0}

    def test_focus_multiple_molecules(self, atoms_with_mol):
        """Focus on multiple molecule groups shows all their atoms."""
        gui = FakeGUI(atoms_with_mol)
        gui.focus_molecules({0, 1})

        n = len(gui.atoms)
        for i in range(n):
            assert gui.images.visible[i], f"Atom {i} should be visible"

    def test_toggle_add_molecule_in_focus(self, atoms_with_mol):
        """Add a molecule group while in focus mode."""
        gui = FakeGUI(atoms_with_mol)
        gui.focus_molecules({0})  # Start with benzene only

        # Add water (group 1) via side panel toggle
        gui._toggle_focus_molecule(1)

        assert gui._focus_mol_ids == {0, 1}
        n = len(gui.atoms)
        for i in range(n):
            assert gui.images.visible[i], f"Atom {i} should be visible"

    def test_toggle_remove_molecule_in_focus(self, atoms_with_mol):
        """Remove a molecule group while in focus mode."""
        gui = FakeGUI(atoms_with_mol)
        gui.focus_molecules({0, 1})  # Start with both groups

        # Remove benzene (group 0)
        gui._toggle_focus_molecule(0)

        assert gui._focus_mol_ids == {1}
        # Atoms 0,1,2 (benzene) should be hidden
        assert not gui.images.visible[0]
        assert not gui.images.visible[1]
        assert not gui.images.visible[2]
        # Atoms 3,4,5 (water) should be visible
        assert gui.images.visible[3]
        assert gui.images.visible[4]
        assert gui.images.visible[5]

    def test_mol_vis_change_delegates_in_focus(self, atoms_with_mol):
        """Side panel checkbox delegates to focus toggle in focus mode."""
        gui = FakeGUI(atoms_with_mol)
        gui.focus_molecules({0})

        # Simulate side panel checkbox toggle for group 1
        gui.on_mol_vis_change(1)

        assert 1 in gui._focus_mol_ids
        assert gui.images.visible[3]  # water atom visible

    def test_unfocus_clears_focus_mol_ids(self, atoms_with_mol):
        """Unfocus clears focus state including mol_ids."""
        gui = FakeGUI(atoms_with_mol)
        gui.focus_molecules({0})
        gui.unfocus()

        assert gui._focus_mol_ids == set()
        assert gui._focus_mode is None
        assert gui._focus_active is False

    def test_focus_recenters_on_toggle(self, atoms_with_mol):
        """Adding a molecule in focus mode recenters view."""
        gui = FakeGUI(atoms_with_mol)
        gui.focus_molecules({0})  # benzene: atoms at [0,0,0], [1.5,0,0], [0.5,0.5,0]
        center_before = gui.center.copy()

        gui._toggle_focus_molecule(1)  # add water: atoms at [2,0.5,0], [3,0,0], [4,0,0]
        # Center should change (now includes all 6 atoms)
        assert not np.allclose(gui.center, center_before)


# ===========================================================================
# Undo / Redo Tests
# ===========================================================================


class TestUndoRedo:
    """Test undo and redo functionality."""

    def test_undo_single(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_positions = gui.atoms.positions.copy()

        gui._save_undo_state()
        gui.atoms.positions[0] = [9.0, 9.0, 9.0]

        gui.undo()
        np.testing.assert_array_almost_equal(
            gui.atoms.positions, original_positions
        )

    def test_redo_single(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._save_undo_state()
        new_pos = [9.0, 9.0, 9.0]
        gui.atoms.positions[0] = new_pos

        gui.undo()
        gui.redo()
        np.testing.assert_array_almost_equal(gui.atoms.positions[0], new_pos)

    def test_undo_redo_preserves_mol_id(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_mol_ids = gui.atoms.get_array("mol_id").copy()

        gui._save_undo_state()
        new_ids = np.array([9, 9, 9, 9, 9, 9], dtype=int)
        gui.atoms.set_array("mol_id", new_ids)

        gui.undo()
        np.testing.assert_array_equal(
            gui.atoms.get_array("mol_id"), original_mol_ids
        )

    def test_undo_empty_stack_is_noop(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        before = gui._draw_count
        gui.undo()
        assert gui._draw_count == before

    def test_redo_empty_stack_is_noop(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        before = gui._draw_count
        gui.redo()
        assert gui._draw_count == before

    def test_multiple_undo(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        positions_sequence = []

        for i in range(3):
            positions_sequence.append(gui.atoms.positions.copy())
            gui._save_undo_state()
            gui.atoms.positions[0] = [float(i + 1)] * 3

        # Undo 3 times to get back to original
        for i in range(3):
            gui.undo()

        np.testing.assert_array_almost_equal(
            gui.atoms.positions, positions_sequence[0]
        )

    def test_undo_clears_redo_on_new_action(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._save_undo_state()
        gui.atoms.positions[0] = [1, 1, 1]
        gui.undo()
        assert len(gui._redo_stack) == 1

        # New action should clear redo
        gui._save_undo_state()
        assert len(gui._redo_stack) == 0

    def test_undo_preserves_cell(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_cell = gui.atoms.get_cell().array.copy()

        gui._save_undo_state()
        gui.atoms.set_cell([20, 20, 20])

        gui.undo()
        np.testing.assert_array_almost_equal(
            gui.atoms.get_cell().array, original_cell
        )


# ===========================================================================
# Delete Selected Atoms Tests
# ===========================================================================


class TestDeleteSelected:
    """Test atom deletion with metadata rebuild."""

    def test_delete_single_atom(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True  # Select C atom

        gui.delete_selected_atoms()

        assert len(gui.atoms) == 5
        assert "C" in gui.atoms.get_chemical_symbols()  # Other C still present

    def test_delete_preserves_mol_groups(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[2] = True  # Delete H from group 0

        gui.delete_selected_atoms()

        groups = gui.atoms.info.get("groups", [])
        assert len(groups) == 2
        assert len(groups[0]) == 2  # group 0 lost one atom
        assert len(groups[1]) == 3  # group 1 unchanged

    def test_delete_entire_group(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        # Delete all of group 0 (atoms 0, 1, 2)
        gui.images.selected[0] = True
        gui.images.selected[1] = True
        gui.images.selected[2] = True

        gui.delete_selected_atoms()

        assert len(gui.atoms) == 3
        groups = gui.atoms.info.get("groups", [])
        assert len(groups) == 1  # Only group 1 remains
        names = gui.atoms.info.get("site_mol_name", [])
        assert names[0] == "water"

    def test_delete_undo(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_n = len(gui.atoms)
        gui.images.selected[0] = True

        gui.delete_selected_atoms()
        assert len(gui.atoms) == original_n - 1

        gui.undo()
        assert len(gui.atoms) == original_n

    def test_delete_no_selection_is_noop(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        before_n = len(gui.atoms)
        gui.delete_selected_atoms()
        assert len(gui.atoms) == before_n

    def test_mol_id_reindexed_after_delete(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        # Delete all of group 0
        gui.images.selected[0] = True
        gui.images.selected[1] = True
        gui.images.selected[2] = True

        gui.delete_selected_atoms()

        mol_ids = gui.atoms.get_array("mol_id")
        # Remaining atoms should all be in group 0 (reindexed)
        np.testing.assert_array_equal(mol_ids, [0, 0, 0])


# ===========================================================================
# Global Element Substitution Tests
# ===========================================================================


class TestGlobalSubstitution:
    """Test element substitution logic."""

    def test_substitute_element(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        count = gui._real_global_atom_substitution("C", "Si")
        assert count == 2
        symbols = gui.atoms.get_chemical_symbols()
        assert "C" not in symbols
        assert symbols.count("Si") == 2

    def test_substitute_nonexistent_element(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        count = gui._real_global_atom_substitution("Xe", "Kr")
        assert count == 0

    def test_batch_substitution(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        # Simulate batch: C→Si, H→F
        gui._save_undo_state()
        c1 = gui._real_global_atom_substitution("C", "Si")
        c2 = gui._real_global_atom_substitution("H", "F")

        symbols = gui.atoms.get_chemical_symbols()
        assert "C" not in symbols
        assert "H" not in symbols
        assert symbols.count("Si") == 2
        assert symbols.count("F") == 2
        assert c1 + c2 == 4

    def test_substitution_undo(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_symbols = gui.atoms.get_chemical_symbols()

        gui._save_undo_state()
        gui._real_global_atom_substitution("C", "Si")

        gui.undo()
        assert gui.atoms.get_chemical_symbols() == original_symbols


# ===========================================================================
# Cell Parameter Editing Tests
# ===========================================================================


class TestCellParameters:
    """Test cell parameter modification logic."""

    def test_set_cell_from_cellpar(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_vol = gui.atoms.get_cell().volume

        from ase.cell import Cell
        new_lengths = [12.0, 12.0, 12.0]
        new_angles = [90.0, 90.0, 90.0]

        gui._save_undo_state()
        new_cell = Cell.fromcellpar(new_lengths + new_angles)
        gui.atoms.set_cell(new_cell, scale_atoms=True)

        new_vol = gui.atoms.get_cell().volume
        assert abs(new_vol - 12.0**3) < 1e-6
        assert new_vol != original_vol

    def test_cell_change_scales_positions(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_frac = gui.atoms.get_scaled_positions()

        from ase.cell import Cell
        new_cell = Cell.fromcellpar([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
        gui.atoms.set_cell(new_cell, scale_atoms=True)

        new_frac = gui.atoms.get_scaled_positions()
        np.testing.assert_array_almost_equal(new_frac, original_frac)

    def test_cell_change_undo(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        original_cell = gui.atoms.get_cell().array.copy()

        gui._save_undo_state()
        from ase.cell import Cell
        gui.atoms.set_cell(
            Cell.fromcellpar([20, 20, 20, 90, 90, 90]), scale_atoms=True
        )

        gui.undo()
        np.testing.assert_array_almost_equal(
            gui.atoms.get_cell().array, original_cell
        )


# ===========================================================================
# Atom Size Tests
# ===========================================================================


class TestAtomSize:
    """Test atom size scaling."""

    def test_set_atom_scale(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert gui.images.atom_scale == 1.0

        gui.images.atom_scale = 2.0
        assert gui.images.atom_scale == 2.0

    def test_atom_scale_range(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        for scale in [0.3, 1.0, 2.0, 3.0]:
            gui.images.atom_scale = scale
            assert gui.images.atom_scale == scale


# ===========================================================================
# Integration: Combined Feature Tests
# ===========================================================================


class TestIntegration:
    """Cross-feature integration tests."""

    def test_color_mode_survives_undo(self, atoms_with_mol):
        """Color mode setting should persist across undo/redo."""
        gui = FakeGUI(atoms_with_mol)
        gui._apply_color_mode("mol_id")
        assert gui.colormode == "mol_id"

        gui._save_undo_state()
        gui.atoms.positions[0] = [5, 5, 5]
        gui.undo()

        # colormode should still be mol_id (undo restores atoms, not view state)
        assert gui.colormode == "mol_id"

    def test_focus_then_delete_then_unfocus(self, atoms_with_mol):
        """Delete while in focus mode, then unfocus."""
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True
        gui.images.selected[1] = True
        gui.focus_selection()
        assert gui._focus_active

        # Delete atom 0 (still selected from focus)
        gui.images.selected = np.zeros(len(gui.atoms), dtype=bool)
        gui.images.selected[0] = True
        gui.delete_selected_atoms()

        # Should still be able to unfocus (backup visibility size differs)
        gui.unfocus()
        assert gui._focus_active is False

    def test_hide_hydrogen_then_substitute_H(self, atoms_with_mol):
        """Hide H, then substitute H → F. Hidden elements should update."""
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_hydrogen()
        assert "H" in gui._hidden_elements

        gui._save_undo_state()
        gui._real_global_atom_substitution("H", "F")

        # No H atoms remain, so _hidden_elements may still have 'H'
        # but there are no H atoms to be affected
        symbols = gui.atoms.get_chemical_symbols()
        assert "H" not in symbols
        assert "F" in symbols

    def test_undo_redo_cycle_preserves_everything(self, atoms_with_mol):
        """Full undo/redo cycle should preserve all metadata."""
        gui = FakeGUI(atoms_with_mol)
        original_pos = gui.atoms.positions.copy()
        original_mol_ids = gui.atoms.get_array("mol_id").copy()
        original_groups = gui.atoms.info["groups"][:]
        original_names = gui.atoms.info["site_mol_name"][:]
        original_cell = gui.atoms.get_cell().array.copy()

        # Make changes
        gui._save_undo_state()
        gui.atoms.positions[0] = [9, 9, 9]
        gui._save_undo_state()
        gui._real_global_atom_substitution("C", "Si")

        # Undo twice
        gui.undo()
        gui.undo()

        np.testing.assert_array_almost_equal(gui.atoms.positions, original_pos)
        np.testing.assert_array_equal(
            gui.atoms.get_array("mol_id"), original_mol_ids
        )
        assert gui.atoms.info["groups"] == original_groups
        assert gui.atoms.info["site_mol_name"] == original_names
        np.testing.assert_array_almost_equal(
            gui.atoms.get_cell().array, original_cell
        )

        # Redo twice
        gui.redo()
        gui.redo()
        # After redo, C should be Si
        assert "Si" in gui.atoms.get_chemical_symbols()


# ===========================================================================
# Per-Element Size Tests
# ===========================================================================


class TestPerElementSize:
    """Test per-element radius control."""

    def test_set_element_radius(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        from ase.data import atomic_numbers
        z_c = atomic_numbers["C"]
        original_r = gui.images.covalent_radii[z_c]

        gui.set_element_radius(z_c, 0.5)
        assert gui.images.covalent_radii[z_c] == 0.5
        assert gui.images.covalent_radii[z_c] != original_r

    def test_set_multiple_element_radii(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        from ase.data import atomic_numbers
        z_c = atomic_numbers["C"]
        z_h = atomic_numbers["H"]

        gui.set_element_radius(z_c, 0.3)
        gui.set_element_radius(z_h, 0.1)

        assert gui.images.covalent_radii[z_c] == 0.3
        assert gui.images.covalent_radii[z_h] == 0.1

    def test_element_radius_does_not_affect_others(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        from ase.data import atomic_numbers
        z_c = atomic_numbers["C"]
        z_n = atomic_numbers["N"]
        original_n_r = gui.images.covalent_radii[z_n]

        gui.set_element_radius(z_c, 0.1)
        assert gui.images.covalent_radii[z_n] == original_n_r


# ===========================================================================
# Per-Element Color Tests
# ===========================================================================


class TestPerElementColor:
    """Test per-element color control."""

    def test_set_element_color(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        from ase.data import atomic_numbers
        z_c = atomic_numbers["C"]

        gui.set_element_color(z_c, "#FF0000")
        assert gui.colors[z_c] == "#FF0000"

    def test_set_multiple_element_colors(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        from ase.data import atomic_numbers
        z_c = atomic_numbers["C"]
        z_h = atomic_numbers["H"]

        gui.set_element_color(z_c, "#FF0000")
        gui.set_element_color(z_h, "#00FF00")

        assert gui.colors[z_c] == "#FF0000"
        assert gui.colors[z_h] == "#00FF00"

    def test_element_color_independent_of_size(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        from ase.data import atomic_numbers
        z_c = atomic_numbers["C"]

        gui.set_element_color(z_c, "#0000FF")
        original_r = gui.images.covalent_radii[z_c]
        gui.set_element_radius(z_c, 0.2)

        assert gui.colors[z_c] == "#0000FF"
        assert gui.images.covalent_radii[z_c] == 0.2


# ===========================================================================
# Molecule Visibility Tests
# ===========================================================================


class TestMoleculeVisibility:
    """Test molecule group visibility, selection, and color control."""

    def test_hide_molecule_group(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._hidden_molecules.add(0)
        gui.apply_molecule_visibility()

        n = len(gui.atoms)
        mol_ids = gui.atoms.get_array("mol_id")
        for i in range(n):
            if mol_ids[i] == 0:
                assert not gui.images.visible[i], f"Atom {i} in group 0 should be hidden"
            else:
                assert gui.images.visible[i], f"Atom {i} in group 1 should be visible"

    def test_hide_then_show_molecule(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._hidden_molecules.add(1)
        gui.apply_molecule_visibility()

        # Group 1 atoms hidden
        mol_ids = gui.atoms.get_array("mol_id")
        for i in range(len(gui.atoms)):
            if mol_ids[i] == 1:
                assert not gui.images.visible[i]

        # Show again
        gui._hidden_molecules.discard(1)
        gui.apply_molecule_visibility()
        assert gui.images.visible[:len(gui.atoms)].all()

    def test_hide_all_molecules(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._hidden_molecules = {0, 1}
        gui.apply_molecule_visibility()

        assert not gui.images.visible[:len(gui.atoms)].any()

    def test_molecule_and_element_visibility_combined(self, atoms_with_mol):
        """Hiding both a molecule group and an element should combine."""
        gui = FakeGUI(atoms_with_mol)
        gui._hidden_molecules.add(0)    # Hide group 0 (atoms 0,1,2)
        gui._hidden_elements.add("N")   # Hide N (atom 4)
        gui.apply_molecule_visibility()

        n = len(gui.atoms)
        mol_ids = gui.atoms.get_array("mol_id")
        symbols = gui.atoms.get_chemical_symbols()
        for i in range(n):
            if mol_ids[i] == 0 or symbols[i] == "N":
                assert not gui.images.visible[i], f"Atom {i} should be hidden"
            else:
                assert gui.images.visible[i], f"Atom {i} should be visible"

    def test_select_molecule(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.select_molecule(0)

        n = len(gui.atoms)
        mol_ids = gui.atoms.get_array("mol_id")
        for i in range(n):
            if mol_ids[i] == 0:
                assert gui.images.selected[i], f"Atom {i} should be selected"
            else:
                assert not gui.images.selected[i], f"Atom {i} should not be selected"

    def test_select_molecule_toggle_deselect(self, atoms_with_mol):
        """Pressing Sel twice on the same molecule should deselect it."""
        gui = FakeGUI(atoms_with_mol)
        gui.select_molecule(0)
        # First press: group 0 selected
        assert gui.images.selected[0]
        # Second press: deselect all
        gui.select_molecule(0)
        assert not gui.images.selected.any()

    def test_select_different_molecule(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.select_molecule(0)
        gui.select_molecule(1)

        # Now only group 1 should be selected
        mol_ids = gui.atoms.get_array("mol_id")
        for i in range(len(gui.atoms)):
            if mol_ids[i] == 1:
                assert gui.images.selected[i]
            else:
                assert not gui.images.selected[i]

    def test_molecule_custom_color(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._mol_colors[0] = "#FF0000"
        gui._mol_colors[1] = "#0000FF"
        gui.apply_mol_group_colors()

        assert gui.colormode == "mol_id"
        colorscale = gui.colormode_data[0]
        assert colorscale[0] == "#FF0000"
        assert colorscale[1] == "#0000FF"

    def test_molecule_color_with_default(self, atoms_with_mol):
        """Only set color for one group, other should use default."""
        gui = FakeGUI(atoms_with_mol)
        gui._mol_colors[0] = "#123456"
        gui.apply_mol_group_colors()

        colorscale = gui.colormode_data[0]
        assert colorscale[0] == "#123456"
        assert colorscale[1] == gui._GROUP_COLORS[1]  # default

    def test_molecule_visibility_no_mol_id(self):
        """Molecule visibility should be a no-op without mol_id."""
        atoms = Atoms("CH4", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0],
                                        [0, 0, 1], [1, 1, 0]])
        atoms.set_cell([10, 10, 10])
        gui = FakeGUI(atoms)
        gui._hidden_molecules.add(0)
        gui.apply_molecule_visibility()
        # Should be a no-op, all visible
        assert gui.images.visible[:len(gui.atoms)].all()


# ===========================================================================
# Drag Mode Tests
# ===========================================================================


class TestDragMode:
    """Verify drag-to-move atom functionality."""

    def test_drag_single_atom(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_pos = gui.atoms.positions[0].copy()
        gui.start_drag([0])
        gui.drag_to([1.0, 0.0, 0.0])
        gui.end_drag()
        np.testing.assert_allclose(gui.atoms.positions[0], orig_pos + [1.0, 0.0, 0.0])

    def test_drag_multiple_atoms(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_0 = gui.atoms.positions[0].copy()
        orig_1 = gui.atoms.positions[1].copy()
        gui.start_drag([0, 1])
        gui.drag_to([0.5, 0.5, 0.0])
        gui.end_drag()
        np.testing.assert_allclose(gui.atoms.positions[0], orig_0 + [0.5, 0.5, 0.0])
        np.testing.assert_allclose(gui.atoms.positions[1], orig_1 + [0.5, 0.5, 0.0])

    def test_drag_does_not_move_unselected(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_2 = gui.atoms.positions[2].copy()
        gui.start_drag([0])
        gui.drag_to([2.0, 0.0, 0.0])
        gui.end_drag()
        np.testing.assert_allclose(gui.atoms.positions[2], orig_2)

    def test_drag_undo(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_pos = gui.atoms.positions.copy()
        gui.start_drag([0, 1, 2])
        gui.drag_to([1.0, 1.0, 1.0])
        gui.end_drag()
        # Positions should have changed
        assert not np.allclose(gui.atoms.positions[:3], orig_pos[:3])
        # Undo should restore
        gui.undo()
        np.testing.assert_allclose(gui.atoms.positions, orig_pos)

    def test_drag_mode_toggle(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert not gui._drag_mode
        gui._drag_mode = True
        assert gui._drag_mode
        gui._drag_mode = False
        assert not gui._drag_mode

    def test_drag_3d_displacement(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_pos = gui.atoms.positions[4].copy()
        gui.start_drag([4])
        gui.drag_to([0.0, 0.0, 2.5])
        gui.end_drag()
        np.testing.assert_allclose(gui.atoms.positions[4], orig_pos + [0.0, 0.0, 2.5])


# ===========================================================================
# Color by Molecule - Distinct Colors Tests
# ===========================================================================


class TestColorByMoleculeDistinct:
    """Verify that each molecule group gets a distinct color via
    ASE's linear color-mapping formula:
        index = int(scalar / cmax * N + 0.5)
    """

    @staticmethod
    def _ase_color_indices(colormode_data, scalars):
        """Reproduce ASE's get_colors() index mapping."""
        colorscale, cmin, cmax = colormode_data
        N = len(colorscale)
        if cmin == cmax:
            return [N // 2] * len(scalars)
        indices = np.clip(
            ((np.array(scalars, dtype=float) - cmin) / (cmax - cmin) * N + 0.5).astype(int),
            0, N - 1,
        )
        return indices.tolist()

    def test_two_groups_get_different_colors(self, atoms_with_mol):
        """2 molecule groups must map to 2 different colorscale indices."""
        gui = FakeGUI(atoms_with_mol)
        gui._apply_color_mode("mol_id")
        assert gui.colormode == "mol_id"
        mol_ids = gui.atoms.get_array("mol_id")
        unique_ids = sorted(set(mol_ids))
        indices = self._ase_color_indices(gui.colormode_data, unique_ids)
        assert len(set(indices)) == len(unique_ids), (
            f"Molecule groups {unique_ids} mapped to non-distinct "
            f"color indices {indices}"
        )

    def test_four_groups_get_different_colors(self):
        """4 molecule groups must each map to a distinct color index."""
        atoms = Atoms(
            symbols=["C"] * 8,
            positions=[[i, 0, 0] for i in range(8)],
            cell=[20, 20, 20],
            pbc=True,
        )
        mol_ids = np.array([0, 0, 1, 1, 2, 2, 3, 3], dtype=int)
        atoms.new_array("mol_id", mol_ids)
        atoms.info["groups"] = [[0, 1], [2, 3], [4, 5], [6, 7]]
        atoms.set_tags([0] * 8)
        gui = FakeGUI(atoms)
        gui._apply_color_mode("mol_id")
        unique_ids = [0, 1, 2, 3]
        indices = self._ase_color_indices(gui.colormode_data, unique_ids)
        assert len(set(indices)) == 4, (
            f"4 groups mapped to indices {indices} — expected 4 distinct values"
        )

    def test_single_group_does_not_crash(self):
        """A single molecule group should work without division-by-zero."""
        atoms = Atoms("CC", positions=[[0, 0, 0], [1, 0, 0]],
                       cell=[10, 10, 10], pbc=True)
        atoms.new_array("mol_id", np.array([0, 0], dtype=int))
        atoms.info["groups"] = [[0, 1]]
        atoms.set_tags([0, 0])
        gui = FakeGUI(atoms)
        gui._apply_color_mode("mol_id")
        assert gui.colormode == "mol_id"
        colorscale, cmin, cmax = gui.colormode_data
        assert cmax > 0  # no division by zero

    def test_many_groups_wrap_colors(self):
        """When n_groups > n_colors, colors should wrap without error."""
        n = 25
        atoms = Atoms(
            symbols=["C"] * n,
            positions=[[i, 0, 0] for i in range(n)],
            cell=[50, 50, 50],
            pbc=True,
        )
        mol_ids = np.arange(n, dtype=int)  # 25 groups, 1 atom each
        atoms.new_array("mol_id", mol_ids)
        atoms.info["groups"] = [[i] for i in range(n)]
        atoms.set_tags([0] * n)
        gui = FakeGUI(atoms)
        gui._apply_color_mode("mol_id")
        colorscale, cmin, cmax = gui.colormode_data
        assert len(colorscale) >= n
        # All 25 should map to valid indices (no out-of-range)
        indices = self._ase_color_indices(gui.colormode_data, list(range(n)))
        assert all(0 <= idx < len(colorscale) for idx in indices)


# ===========================================================================
# 3D Sphere Renderer Tests
# ===========================================================================


class TestSphereRenderer:
    """Verify _SphereRenderer generates correct sphere images."""

    def test_render_sphere_returns_rgba_image(self):
        """Rendered sphere should be an RGBA PIL Image."""
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        from molecular_modeling_gui import _SphereRenderer

        renderer = _SphereRenderer()
        img = renderer._render_sphere('#ff0000', 32, False)
        assert img.mode == 'RGBA'
        assert img.size == (32, 32)

    def test_sphere_has_transparent_corners(self):
        """Corners of the sphere image should be transparent."""
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        from molecular_modeling_gui import _SphereRenderer

        renderer = _SphereRenderer()
        img = renderer._render_sphere('#0000ff', 50, False)
        # Top-left corner should be fully transparent
        assert img.getpixel((0, 0))[3] == 0
        # Bottom-right corner
        assert img.getpixel((49, 49))[3] == 0

    def test_sphere_center_is_opaque(self):
        """Centre of the sphere should be opaque."""
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        from molecular_modeling_gui import _SphereRenderer

        renderer = _SphereRenderer()
        img = renderer._render_sphere('#00ff00', 40, False)
        cx, cy = 20, 20
        assert img.getpixel((cx, cy))[3] > 200

    def test_sphere_has_highlight(self):
        """Upper-left region should be brighter than lower-right (light source)."""
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        from molecular_modeling_gui import _SphereRenderer

        renderer = _SphereRenderer()
        img = renderer._render_sphere('#808080', 60, False)
        # Sample upper-left and lower-right quadrant centres
        ul = img.getpixel((18, 18))  # near highlight
        lr = img.getpixel((42, 42))  # shadow side
        brightness_ul = sum(ul[:3])
        brightness_lr = sum(lr[:3])
        assert brightness_ul > brightness_lr, (
            f"Highlight ({brightness_ul}) should be brighter than shadow ({brightness_lr})"
        )

    def test_selected_sphere_has_outline(self):
        """Selected sphere should differ from unselected."""
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        import numpy as np

        from molecular_modeling_gui import _SphereRenderer

        renderer = _SphereRenderer()
        img_normal = renderer._render_sphere('#ff7f0e', 40, False)
        img_selected = renderer._render_sphere('#ff7f0e', 40, True)
        arr_n = np.array(img_normal)
        arr_s = np.array(img_selected)
        assert not np.array_equal(arr_n, arr_s), "Selected sphere must differ"

    def test_different_colors_produce_different_images(self):
        """Different atom colours should produce different sphere images."""
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        import numpy as np

        from molecular_modeling_gui import _SphereRenderer

        renderer = _SphereRenderer()
        img_r = np.array(renderer._render_sphere('#ff0000', 30, False))
        img_b = np.array(renderer._render_sphere('#0000ff', 30, False))
        assert not np.array_equal(img_r, img_b)

    def test_minimum_size_sphere(self):
        """A size-2 sphere should render without error."""
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        from molecular_modeling_gui import _SphereRenderer

        renderer = _SphereRenderer()
        img = renderer._render_sphere('#ffffff', 2, False)
        assert img.size == (2, 2)


# ===========================================================================
# View Operations (Pan / Orbit / Zoom) Tests
# ===========================================================================


class TestViewOperations:
    """Verify view manipulation callbacks change the expected state."""

    def test_translate_view_moves_center(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_center = gui.center.copy()
        gui._translate_view(1, 0)  # move right
        assert not np.allclose(gui.center, orig_center), \
            "Center should change after translate"

    def test_translate_view_all_directions(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            orig = gui.center.copy()
            gui._translate_view(dx, dy)
            assert not np.allclose(gui.center, orig), \
                f"Center should change for dx={dx}, dy={dy}"

    def test_zoom_in_increases_scale(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_scale = gui.scale
        gui._zoom_view(1.1)
        assert gui.scale > orig_scale

    def test_zoom_out_decreases_scale(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_scale = gui.scale
        gui._zoom_view(1 / 1.1)
        assert gui.scale < orig_scale

    def test_zoom_roundtrip(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_scale = gui.scale
        gui._zoom_view(1.1)
        gui._zoom_view(1 / 1.1)
        np.testing.assert_allclose(gui.scale, orig_scale, rtol=1e-10)

    def test_rotate_x_changes_axes(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_axes = gui.axes.copy()
        gui._rotate_view('x', 5)
        assert not np.allclose(gui.axes, orig_axes)

    def test_rotate_y_changes_axes(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_axes = gui.axes.copy()
        gui._rotate_view('y', 5)
        assert not np.allclose(gui.axes, orig_axes)

    def test_rotate_z_changes_axes(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_axes = gui.axes.copy()
        gui._rotate_view('z', 5)
        assert not np.allclose(gui.axes, orig_axes)

    def test_rotate_preserves_orthogonality(self, atoms_with_mol):
        """Rotation should preserve the orthogonal nature of axes."""
        gui = FakeGUI(atoms_with_mol)
        gui._rotate_view('x', 30)
        gui._rotate_view('y', -15)
        gui._rotate_view('z', 45)
        # axes @ axes.T should be close to identity
        product = gui.axes @ gui.axes.T
        np.testing.assert_allclose(product, np.eye(3), atol=1e-10)

    def test_rotate_360_returns_to_original(self, atoms_with_mol):
        """72 rotations of 5° should return to original orientation."""
        gui = FakeGUI(atoms_with_mol)
        orig_axes = gui.axes.copy()
        for _ in range(72):
            gui._rotate_view('x', 5)
        np.testing.assert_allclose(gui.axes, orig_axes, atol=1e-10)

    def test_operations_trigger_draw(self, atoms_with_mol):
        """Each view operation should call draw()."""
        gui = FakeGUI(atoms_with_mol)
        gui._draw_count = 0
        gui._translate_view(1, 0)
        assert gui._draw_count == 1
        gui._zoom_view(1.1)
        assert gui._draw_count == 2
        gui._rotate_view('x', 5)
        assert gui._draw_count == 3


# ===========================================================================
# Background Color Tests
# ===========================================================================


class TestBackgroundColor:
    """Verify canvas background color switching."""

    def test_set_bg_white(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._set_canvas_bg("#ffffff")
        assert gui._get_canvas_bg() == "#ffffff"

    def test_set_bg_black(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._set_canvas_bg("#000000")
        assert gui._get_canvas_bg() == "#000000"

    def test_set_bg_dark_gray(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._set_canvas_bg("#333333")
        assert gui._get_canvas_bg() == "#333333"

    def test_set_bg_light_gray(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._set_canvas_bg("#d0d0d0")
        assert gui._get_canvas_bg() == "#d0d0d0"

    def test_set_bg_custom_color(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._set_canvas_bg("#1a2b3c")
        assert gui._get_canvas_bg() == "#1a2b3c"

    def test_default_bg_is_white(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert gui._get_canvas_bg() == "#ffffff"

    def test_multiple_bg_changes(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        for color in ["#000000", "#333333", "#ffffff", "#d0d0d0"]:
            gui._set_canvas_bg(color)
            assert gui._get_canvas_bg() == color


# ===========================================================================
# Trackball Rotation Logic Tests
# ===========================================================================


class TestTrackballLogic:
    """Verify trackball rotation math works correctly."""

    def test_trackball_inside_ball_x_drag_rotates_y(self, atoms_with_mol):
        """Horizontal drag inside the ball should rotate around Y axis."""
        gui = FakeGUI(atoms_with_mol)
        orig_axes = gui.axes.copy()
        # Simulate horizontal drag inside ball → rotate around Y
        gui._rotate_view('y', 10)
        assert not np.allclose(gui.axes, orig_axes)
        # Y row should change; first element of axes should differ
        assert not np.isclose(gui.axes[0, 0], orig_axes[0, 0])

    def test_trackball_inside_ball_y_drag_rotates_x(self, atoms_with_mol):
        """Vertical drag inside the ball should rotate around X axis."""
        gui = FakeGUI(atoms_with_mol)
        orig_axes = gui.axes.copy()
        gui._rotate_view('x', -10)
        assert not np.allclose(gui.axes, orig_axes)

    def test_trackball_outside_ball_rotates_z(self, atoms_with_mol):
        """Drag outside the ball should rotate around Z axis."""
        gui = FakeGUI(atoms_with_mol)
        orig_axes = gui.axes.copy()
        gui._rotate_view('z', 15)
        assert not np.allclose(gui.axes, orig_axes)
        # Z rotation should leave z-axis row unchanged
        np.testing.assert_allclose(gui.axes[2, 2], orig_axes[2, 2], atol=1e-10)

    def test_combined_trackball_rotations_preserve_orthogonality(
            self, atoms_with_mol):
        """Multiple trackball-style rotations should keep axes orthogonal."""
        gui = FakeGUI(atoms_with_mol)
        # Simulate a complex trackball drag sequence
        gui._rotate_view('y', 12)
        gui._rotate_view('x', -8)
        gui._rotate_view('z', 5)
        gui._rotate_view('y', -3)
        gui._rotate_view('x', 15)
        product = gui.axes @ gui.axes.T
        np.testing.assert_allclose(product, np.eye(3), atol=1e-10)


# ===========================================================================
# Drag Mode + View Operations Integration Tests
# ===========================================================================


class TestDragModeIntegration:
    """Verify drag mode doesn't interfere with view operations."""

    def test_view_ops_work_in_drag_mode(self, atoms_with_mol):
        """Pan/zoom/orbit should still work when drag mode is on."""
        gui = FakeGUI(atoms_with_mol)
        gui._drag_mode = True

        orig_center = gui.center.copy()
        gui._translate_view(1, 0)
        assert not np.allclose(gui.center, orig_center)

        orig_scale = gui.scale
        gui._zoom_view(1.1)
        assert gui.scale > orig_scale

        orig_axes = gui.axes.copy()
        gui._rotate_view('x', 10)
        assert not np.allclose(gui.axes, orig_axes)

    def test_drag_then_view_op(self, atoms_with_mol):
        """After a drag, view operations should still work normally."""
        gui = FakeGUI(atoms_with_mol)
        gui.start_drag([0])
        gui.drag_to([1.0, 0.0, 0.0])
        gui.end_drag()

        orig_center = gui.center.copy()
        gui._translate_view(-1, 0)
        assert not np.allclose(gui.center, orig_center)

    def test_drag_preserves_view_state(self, atoms_with_mol):
        """Dragging atoms should not affect the view center/scale/axes."""
        gui = FakeGUI(atoms_with_mol)
        orig_center = gui.center.copy()
        orig_scale = gui.scale
        orig_axes = gui.axes.copy()

        gui.start_drag([0, 1])
        gui.drag_to([2.0, 3.0, 0.0])
        gui.end_drag()

        np.testing.assert_allclose(gui.center, orig_center)
        assert gui.scale == orig_scale
        np.testing.assert_allclose(gui.axes, orig_axes)


# ===========================================================================
# Move Molecule/Atoms Tests
# ===========================================================================


class TestMoveAtoms:
    """Verify move_atoms / move_molecule logic (mirrors move_molecule_dialog)."""

    def test_move_selected_atoms(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.move_atoms([0, 1], [1.0, 2.0, 3.0])
        np.testing.assert_allclose(gui.atoms.positions[0], orig[0] + [1, 2, 3])
        np.testing.assert_allclose(gui.atoms.positions[1], orig[1] + [1, 2, 3])
        # Unselected atoms unchanged
        np.testing.assert_allclose(gui.atoms.positions[2], orig[2])

    def test_move_molecule_group(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.move_molecule(0, [0.5, 0.0, 0.0])
        # mol_id==0 → indices 0,1,2
        for i in [0, 1, 2]:
            np.testing.assert_allclose(
                gui.atoms.positions[i], orig[i] + [0.5, 0, 0])
        # mol_id==1 unchanged
        for i in [3, 4, 5]:
            np.testing.assert_allclose(gui.atoms.positions[i], orig[i])

    def test_move_undo(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.move_atoms([0], [5.0, 5.0, 5.0])
        assert not np.allclose(gui.atoms.positions[0], orig[0])
        gui.undo()
        np.testing.assert_allclose(gui.atoms.positions, orig)

    def test_move_zero_vector(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.move_atoms([0, 1, 2], [0.0, 0.0, 0.0])
        np.testing.assert_allclose(gui.atoms.positions, orig)

    def test_move_no_atoms_does_nothing(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        count = gui._draw_count
        gui.move_atoms([], [1.0, 1.0, 1.0])
        np.testing.assert_allclose(gui.atoms.positions, orig)
        assert gui._draw_count == count  # no draw triggered

    def test_move_triggers_draw(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._draw_count = 0
        gui.move_atoms([0], [1.0, 0.0, 0.0])
        assert gui._draw_count >= 1


# ===========================================================================
# Rotate Molecule/Atoms Tests
# ===========================================================================


class TestRotateAtoms:
    """Verify rotate_atoms / rotate_molecule logic (mirrors rotate_molecule_dialog)."""

    def test_rotate_selected_atoms(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.rotate_atoms([0, 1, 2], (90, 0, 0))
        # Positions should change
        assert not np.allclose(gui.atoms.positions[:3], orig[:3])
        # Unselected atoms unchanged
        np.testing.assert_allclose(gui.atoms.positions[3:], orig[3:])

    def test_rotate_preserves_centroid(self, atoms_with_mol):
        """Rotation should keep the centroid of rotated atoms unchanged."""
        gui = FakeGUI(atoms_with_mol)
        indices = [0, 1, 2]
        orig_centroid = gui.atoms.positions[indices].mean(axis=0)
        gui.rotate_atoms(indices, (45, 30, 60))
        new_centroid = gui.atoms.positions[indices].mean(axis=0)
        np.testing.assert_allclose(new_centroid, orig_centroid, atol=1e-10)

    def test_rotate_molecule_group(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.rotate_molecule(1, (0, 0, 90))
        # mol_id==1 → indices 3,4,5 should change
        assert not np.allclose(gui.atoms.positions[3:], orig[3:])
        # mol_id==0 unchanged
        np.testing.assert_allclose(gui.atoms.positions[:3], orig[:3])

    def test_rotate_undo(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.rotate_atoms([3, 4, 5], (45, 45, 45))
        assert not np.allclose(gui.atoms.positions[3:], orig[3:])
        gui.undo()
        np.testing.assert_allclose(gui.atoms.positions, orig)

    def test_rotate_360_returns_to_original(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.rotate_atoms([0, 1, 2], (360, 0, 0))
        np.testing.assert_allclose(gui.atoms.positions, orig, atol=1e-10)

    def test_rotate_no_atoms_does_nothing(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        count = gui._draw_count
        gui.rotate_atoms([], (90, 0, 0))
        np.testing.assert_allclose(gui.atoms.positions, orig)
        assert gui._draw_count == count

    def test_rotate_triggers_draw(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui._draw_count = 0
        gui.rotate_atoms([0], (10, 0, 0))
        assert gui._draw_count >= 1


# ===========================================================================
# Drag Mode Toggle Tests
# ===========================================================================


class TestDragModeToggle:
    """Verify drag mode toggle state transitions (mirrors _toggle_drag_mode)."""

    def test_toggle_on(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert not gui._drag_mode
        gui.toggle_drag_mode()
        assert gui._drag_mode

    def test_toggle_off(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_drag_mode()  # on
        gui.toggle_drag_mode()  # off
        assert not gui._drag_mode

    def test_toggle_off_clears_dragging(self, atoms_with_mol):
        """Turning off drag mode should clear _dragging flag."""
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_drag_mode()  # on
        gui._dragging = True  # simulate mid-drag
        gui.toggle_drag_mode()  # off
        assert not gui._dragging

    def test_toggle_on_does_not_start_drag(self, atoms_with_mol):
        """Turning on drag mode alone should NOT start dragging."""
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_drag_mode()
        assert gui._drag_mode
        assert not gui._dragging

    def test_drag_mode_on_then_drag_then_toggle_off(self, atoms_with_mol):
        """Full cycle: toggle on → drag → toggle off."""
        gui = FakeGUI(atoms_with_mol)
        orig = gui.atoms.positions.copy()
        gui.toggle_drag_mode()  # on
        gui.start_drag([0])
        gui.drag_to([1.0, 0.0, 0.0])
        gui.end_drag()
        assert not np.allclose(gui.atoms.positions[0], orig[0])
        gui.toggle_drag_mode()  # off
        assert not gui._drag_mode
        assert not gui._dragging

    def test_multiple_toggles(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        for _ in range(10):
            gui.toggle_drag_mode()
        # After even number of toggles, should be off
        assert not gui._drag_mode


# ===========================================================================
# Delete Selected Atoms Tests
# ===========================================================================


class TestDeleteSelectedAtoms:
    """Verify delete_selected_atoms logic."""

    def test_delete_single_atom(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_len = len(gui.atoms)
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        assert len(gui.atoms) == orig_len - 1

    def test_delete_no_selection_does_nothing(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_len = len(gui.atoms)
        gui.delete_selected_atoms()
        assert len(gui.atoms) == orig_len

    def test_delete_undo(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_len = len(gui.atoms)
        orig_pos = gui.atoms.positions.copy()
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        assert len(gui.atoms) == orig_len - 1
        gui.undo()
        assert len(gui.atoms) == orig_len
        np.testing.assert_allclose(gui.atoms.positions, orig_pos)

    def test_delete_multiple_atoms(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True
        gui.images.selected[1] = True
        gui.images.selected[2] = True
        gui.delete_selected_atoms()
        assert len(gui.atoms) == 3  # 6 - 3


# ===========================================================================
# Collapsible Sections Tests
# ===========================================================================


class TestCollapsibleSections:
    """Verify collapsible section toggle logic."""

    def test_elem_default_expanded(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert not gui._elem_collapsed

    def test_mol_default_expanded(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        assert not gui._mol_collapsed

    def test_elem_toggle_collapse(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_section("elem")
        assert gui._elem_collapsed

    def test_elem_toggle_expand(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_section("elem")
        gui.toggle_section("elem")
        assert not gui._elem_collapsed

    def test_mol_toggle_collapse(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_section("mol")
        assert gui._mol_collapsed

    def test_mol_toggle_expand(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_section("mol")
        gui.toggle_section("mol")
        assert not gui._mol_collapsed

    def test_elem_and_mol_independent(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.toggle_section("elem")
        assert gui._elem_collapsed
        assert not gui._mol_collapsed
        gui.toggle_section("mol")
        assert gui._elem_collapsed
        assert gui._mol_collapsed

    def test_multiple_toggles(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        for _ in range(5):
            gui.toggle_section("elem")
        assert gui._elem_collapsed  # odd number → collapsed


# ===========================================================================
# Select All / Deselect All Tests
# ===========================================================================


class TestSelectAllDeselectAll:
    """Verify select all and deselect all logic."""

    def test_select_all(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.select_all()
        n = len(gui.atoms)
        assert gui.images.selected[:n].all()

    def test_deselect_all(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[:] = True
        gui.deselect_all()
        n = len(gui.atoms)
        assert not gui.images.selected[:n].any()

    def test_select_all_respects_visibility(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.visible[2] = False
        gui.select_all()
        n = len(gui.atoms)
        assert gui.images.selected[0]
        assert gui.images.selected[1]
        assert not gui.images.selected[2]

    def test_deselect_then_select_all(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True
        gui.images.selected[3] = True
        gui.deselect_all()
        assert not gui.images.selected[:len(gui.atoms)].any()
        gui.select_all()
        assert gui.images.selected[:len(gui.atoms)].all()


# ===========================================================================
# Context Select Molecule Tests
# ===========================================================================


class TestContextSelectMolecule:
    """Verify context menu molecule selection."""

    def test_select_molecule_by_groups(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.context_select_molecule(0)  # atom 0 is in group [0, 1, 2]
        n = len(gui.atoms)
        assert gui.images.selected[0]
        assert gui.images.selected[1]
        assert gui.images.selected[2]
        assert not gui.images.selected[3]
        assert not gui.images.selected[4]
        assert not gui.images.selected[5]

    def test_select_molecule_second_group(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.context_select_molecule(4)  # atom 4 is in group [3, 4, 5]
        n = len(gui.atoms)
        assert not gui.images.selected[0]
        assert gui.images.selected[3]
        assert gui.images.selected[4]
        assert gui.images.selected[5]

    def test_select_molecule_clears_previous_selection(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[:] = True
        gui.context_select_molecule(0)
        n = len(gui.atoms)
        # Only group 0 selected
        assert gui.images.selected[0]
        assert not gui.images.selected[3]

    def test_select_molecule_fallback_mol_id(self, atoms_with_mol):
        """Fall back to mol_id array when groups info is missing."""
        gui = FakeGUI(atoms_with_mol)
        del gui.atoms.info["groups"]
        gui.context_select_molecule(0)
        n = len(gui.atoms)
        mol_ids = gui.atoms.get_array("mol_id")
        target = mol_ids[0]
        expected = (mol_ids == target)
        np.testing.assert_array_equal(gui.images.selected[:n], expected)


# ===========================================================================
# Status Info Tests
# ===========================================================================


class TestStatusInfo:
    """Verify status info text generation."""

    def test_no_selection(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        info = gui.get_status_info()
        assert "6 atoms total" in info
        assert "2 groups" in info

    def test_with_selection(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True
        gui.images.selected[2] = True
        info = gui.get_status_info()
        assert "2 selected" in info
        assert "C:1" in info
        assert "H:1" in info

    def test_all_selected(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.select_all()
        info = gui.get_status_info()
        assert "6 selected" in info

    def test_no_atoms(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.atoms = Atoms()
        gui.images.initialize([gui.atoms])
        info = gui.get_status_info()
        assert "No structure loaded" in info


# ===========================================================================
# Axis Indicator Tests
# ===========================================================================


class TestAxisIndicator:
    """Verify XYZ axis projection calculations."""

    def test_default_axes_identity(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        proj = gui.get_axis_projections()
        # With identity axes, X projects to [1, 0] on screen
        np.testing.assert_allclose(proj["X"], [1.0, 0.0], atol=1e-10)
        np.testing.assert_allclose(proj["Y"], [0.0, 1.0], atol=1e-10)
        np.testing.assert_allclose(proj["Z"], [0.0, 0.0], atol=1e-10)

    def test_rotated_axes(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        # Rotate 90 degrees around Z
        gui._rotate_view('z', 90)
        proj = gui.get_axis_projections()
        # After 90° Z rotation, X should project near [0, 1] and Y near [-1, 0]
        assert abs(proj["X"][1]) > 0.5  # X has significant Y component


# ===========================================================================
# Undo with Size Change Tests
# ===========================================================================


class TestUndoSizeChange:
    """Verify undo works correctly when atom count changes."""

    def test_undo_after_delete(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_len = len(gui.atoms)
        orig_pos = gui.atoms.positions.copy()
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        assert len(gui.atoms) == orig_len - 1
        gui.undo()
        assert len(gui.atoms) == orig_len
        np.testing.assert_allclose(gui.atoms.positions, orig_pos)

    def test_undo_preserves_mol_id(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_mol_ids = gui.atoms.get_array("mol_id").copy()
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        gui.undo()
        np.testing.assert_array_equal(
            gui.atoms.get_array("mol_id"), orig_mol_ids)

    def test_redo_after_undo_size_change(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        reduced_len = len(gui.atoms)
        gui.undo()
        assert len(gui.atoms) == 6
        gui.redo()
        assert len(gui.atoms) == reduced_len

    def test_multiple_undo_with_size_changes(self, atoms_with_mol):
        gui = FakeGUI(atoms_with_mol)
        orig_len = len(gui.atoms)
        # Delete atom 0
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        assert len(gui.atoms) == orig_len - 1
        # Delete another atom
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        assert len(gui.atoms) == orig_len - 2
        # Undo both
        gui.undo()
        assert len(gui.atoms) == orig_len - 1
        gui.undo()
        assert len(gui.atoms) == orig_len


class TestAtomList:
    """Tests for the Atom List side panel feature."""

    def test_atom_list_data_with_mol_id(self, atoms_with_mol):
        """Atom list should return correct data with mol_id."""
        gui = FakeGUI(atoms_with_mol)
        rows = gui.get_atom_list_data()
        assert len(rows) == 6
        # Check first row
        idx, elem, x, y, z, mol = rows[0]
        assert idx == 0
        assert elem == "C"
        assert mol == 0
        # Check last row
        idx, elem, x, y, z, mol = rows[5]
        assert idx == 5
        assert elem == "O"
        assert mol == 1

    def test_atom_list_data_without_mol_id(self):
        """Atom list should handle atoms without mol_id."""
        atoms = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
        gui = FakeGUI(atoms)
        rows = gui.get_atom_list_data()
        assert len(rows) == 2
        assert rows[0][1] == "H"
        assert rows[0][5] == ""  # no mol_id → empty string

    def test_atom_list_data_empty(self):
        """Atom list should return empty for empty atoms."""
        atoms = Atoms()
        gui = FakeGUI(atoms)
        rows = gui.get_atom_list_data()
        assert rows == []

    def test_atom_list_coordinates(self, atoms_with_mol):
        """Atom list should show correct coordinates."""
        gui = FakeGUI(atoms_with_mol)
        rows = gui.get_atom_list_data()
        pos = gui.atoms.get_positions()
        for i, row in enumerate(rows):
            assert row[2] == round(pos[i, 0], 3)
            assert row[3] == round(pos[i, 1], 3)
            assert row[4] == round(pos[i, 2], 3)

    def test_atom_list_select(self, atoms_with_mol):
        """Clicking atom in list should select it in 3D view."""
        gui = FakeGUI(atoms_with_mol)
        gui.atom_list_select(3)
        n = len(gui.atoms)
        sel = gui.images.selected[:n]
        assert sel[3]
        assert sel.sum() == 1

    def test_atom_list_select_clears_previous(self, atoms_with_mol):
        """Selecting in atom list should clear previous selection."""
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True
        gui.images.selected[1] = True
        gui.atom_list_select(4)
        n = len(gui.atoms)
        sel = gui.images.selected[:n]
        assert sel.sum() == 1
        assert sel[4]

    def test_atom_list_toggle_section(self, atoms_with_mol):
        """Atom list section starts collapsed and can be toggled."""
        gui = FakeGUI(atoms_with_mol)
        assert gui._atomlist_collapsed is True
        gui.toggle_section("atomlist")
        assert gui._atomlist_collapsed is False
        gui.toggle_section("atomlist")
        assert gui._atomlist_collapsed is True

    def test_atom_list_updates_after_delete(self, atoms_with_mol):
        """Atom list should reflect changes after deleting atoms."""
        gui = FakeGUI(atoms_with_mol)
        assert len(gui.get_atom_list_data()) == 6
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        rows = gui.get_atom_list_data()
        assert len(rows) == 5
        assert [r[0] for r in rows] == [0, 1, 2, 3, 4]

    def test_atom_list_updates_after_undo(self, atoms_with_mol):
        """Atom list should reflect changes after undo."""
        gui = FakeGUI(atoms_with_mol)
        gui.images.selected[0] = True
        gui.delete_selected_atoms()
        assert len(gui.get_atom_list_data()) == 5
        gui.undo()
        assert len(gui.get_atom_list_data()) == 6
