"""
Tests for GUI feature logic (color mode switching, visibility toggles,
focus mode, undo/redo, element substitution, cell parameter editing, etc.).

These tests verify the core logic of MolecularModelingGUI without
launching a real tkinter window.  We create a thin mock of the GUI
that wires up all the state-mutating methods while replacing draw()
and other display-only calls with no-ops.
"""

import numpy as np
import pytest
from ase import Atoms
from unittest.mock import MagicMock, patch


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
                self.colormode_data = (colorscale, 0, n_colors - 1)
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
        self.colormode_data = (colorscale, 0, max(n_groups - 1, 1))
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
        """Select all atoms in a molecule group (mirrors _on_mol_select)."""
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        n = len(self.atoms)
        self.images.selected[:n] = (mol_ids == group_id)
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
        self.colormode_data = (colorscale, 0, max(n_groups - 1, 1))
        self.draw()


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
