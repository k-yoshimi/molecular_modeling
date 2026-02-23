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
from source.gui.theme import COLORS, FONTS, SPACING, configure_theme


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


class _SphereRenderer:
    """Generate and cache 3D-shaded sphere images using Pillow.

    Each sphere is rendered with a radial gradient simulating a light source
    from the upper-left, producing a realistic 3D appearance.  Images are
    cached by (hex_color, diameter, selected) so that identical spheres are
    only rendered once.
    """

    def __init__(self):
        self._cache = {}  # (color, diameter, selected) -> PhotoImage

    def clear_cache(self):
        self._cache.clear()

    def _hex_to_rgb(self, hex_color):
        h = hex_color.lstrip('#')
        if len(h) == 3:
            h = ''.join(c * 2 for c in h)
        return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))

    def get_sphere(self, hex_color, diameter, selected, tk_root):
        """Return a tk.PhotoImage of a shaded sphere.

        Parameters
        ----------
        hex_color : str   e.g. '#ff7f0e'
        diameter  : int   pixel diameter (>= 2)
        selected  : bool  draw selection highlight ring
        tk_root   : tk widget for PhotoImage master
        """
        key = (hex_color, diameter, selected)
        img = self._cache.get(key)
        if img is not None:
            return img

        from PIL import Image, ImageDraw

        size = max(int(diameter), 2)
        # Render at 1x (fast enough with per-pixel math via numpy)
        img_pil = self._render_sphere(hex_color, size, selected)
        # Convert to PhotoImage
        from PIL import ImageTk
        photo = ImageTk.PhotoImage(img_pil, master=tk_root)
        self._cache[key] = photo
        return photo

    def _render_sphere(self, hex_color, size, selected):
        from PIL import Image
        import numpy as np

        r, g, b = self._hex_to_rgb(hex_color)
        s = size
        # Create coordinate grids centered on the sphere
        y_coords, x_coords = np.mgrid[0:s, 0:s]
        cx, cy = (s - 1) / 2.0, (s - 1) / 2.0
        radius = s / 2.0

        # Distance from center, normalised to [0, 1]
        dx = (x_coords - cx) / radius
        dy = (y_coords - cy) / radius
        dist2 = dx * dx + dy * dy

        # Mask: inside the sphere
        inside = dist2 <= 1.0

        # Surface normal z-component (hemisphere)
        nz = np.sqrt(np.clip(1.0 - dist2, 0, 1))

        # Light direction (upper-left, towards viewer)
        lx, ly, lz = -0.4, -0.5, 0.75
        ln = np.sqrt(lx * lx + ly * ly + lz * lz)
        lx, ly, lz = lx / ln, ly / ln, lz / ln

        # Diffuse lighting (Lambertian)
        diffuse = np.clip(dx * lx + dy * ly + nz * lz, 0, 1)

        # Specular highlight (Phong-like)
        # Reflection of light vector about surface normal
        rx = lx - 2 * (dx * lx + dy * ly + nz * lz) * dx
        ry = ly - 2 * (dx * lx + dy * ly + nz * lz) * dy
        rz = lz - 2 * (dx * lx + dy * ly + nz * lz) * nz
        # View direction is (0, 0, 1)
        spec = np.clip(rz, 0, 1) ** 30  # shininess exponent

        # Ambient + diffuse + specular
        ambient = 0.25
        intensity = np.clip(ambient + 0.65 * diffuse + 0.6 * spec, 0, 1)

        # Apply to base colour
        out_r = np.clip(r * intensity + 255 * 0.6 * spec, 0, 255).astype(np.uint8)
        out_g = np.clip(g * intensity + 255 * 0.6 * spec, 0, 255).astype(np.uint8)
        out_b = np.clip(b * intensity + 255 * 0.6 * spec, 0, 255).astype(np.uint8)

        # Alpha channel: antialiased circle edge
        # Compute pixel distance from circle boundary
        pixel_dist = (np.sqrt(dist2) - 1.0) * radius  # in pixels
        # Smooth falloff over ~0.8 pixel at the edge
        alpha = np.clip(1.0 - pixel_dist * 1.5, 0, 1)
        alpha = (alpha * 255).astype(np.uint8)

        # Compose RGBA image
        rgba = np.zeros((s, s, 4), dtype=np.uint8)
        rgba[:, :, 0] = out_r
        rgba[:, :, 1] = out_g
        rgba[:, :, 2] = out_b
        rgba[:, :, 3] = alpha

        img = Image.fromarray(rgba)

        # Draw selection ring
        if selected:
            from PIL import ImageDraw
            draw = ImageDraw.Draw(img)
            margin = max(1, s // 20)
            lw = max(2, s // 15)
            bbox = [margin, margin, s - 1 - margin, s - 1 - margin]
            draw.ellipse(bbox, outline=(0, 69, 0, 255), width=lw)

        return img


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
        configure_theme(self.window.win)
        self.controller = get_controller()
        if self.atoms is not None:
            self.controller.set_atoms(self.atoms)
        self._repeat_id = None
        self._undo_stack = []
        self._redo_stack = []
        self._drag_mode = False
        self._dragging = False
        self._drag_start_xy = None
        self._drag_start_positions = None
        self._setup_3d_rendering()
        self._setup_file_panel()
        self._setup_toolbar()
        self._setup_bottom_panel()
        self._setup_side_panel()
        self._setup_context_menu()
        self._setup_keyboard_shortcuts()
        # Set a reasonable default window size (larger than ASE's 450x450)
        self._resize_window(700, 700)

    # =========================================================================
    # 3D sphere rendering
    # =========================================================================

    def _setup_3d_rendering(self):
        """Install the 3D sphere renderer, replacing ASE's flat circles."""
        if tk is None:
            return
        self._sphere_renderer = _SphereRenderer()
        self._3d_mode = True  # enabled by default
        # Keep references to PhotoImages so they aren't garbage-collected
        self._sphere_images = []
        # Save original circle method
        self._original_circle = self.window.circle
        # Replace the window's circle method
        self.window.circle = self._draw_sphere

    def _draw_sphere(self, color, selected, *bbox):
        """Draw a 3D-shaded sphere image instead of a flat circle."""
        if not self._3d_mode:
            self._original_circle(color, selected, *bbox)
            return
        x1, y1, x2, y2 = (int(v) for v in bbox)
        diameter = max(x2 - x1, y2 - y1, 2)
        try:
            photo = self._sphere_renderer.get_sphere(
                color, diameter, selected, self.window.win)
            self._sphere_images.append(photo)
            cx = (x1 + x2) // 2
            cy = (y1 + y2) // 2
            self.window.canvas.create_image(cx, cy, image=photo)
        except Exception:
            # Fallback to flat circle on any error
            self._original_circle(color, selected, *bbox)

    def set_atoms(self, atoms):
        """Override to extend bonds to atom centers when atoms are hidden."""
        super().set_atoms(atoms)
        if not self.showing_bonds():
            return
        n = len(atoms)
        visible = self.images.visible
        nbonds = len(self.X_bonds)
        if nbonds == 0:
            return
        # Recompute bonds so hidden atoms get bonds drawn to their centers
        bonds = self.get_bonds(atoms)
        if len(bonds) == 0:
            return
        ncell = len(self.X_cell)
        P = atoms.positions
        cell = atoms.cell
        Af = self.images.repeat[:, np.newaxis] * cell
        r = 0.65 * self.get_covalent_radii()
        for i, bond in enumerate(bonds):
            a0, a1 = bond[0], bond[1]
            hidden0 = not visible[a0] if a0 < len(visible) else False
            hidden1 = not visible[a1] if a1 < len(visible) else False
            if not hidden0 and not hidden1:
                continue
            # Recompute this bond's start/end
            pa = P[a0]
            pb = P[a1] + np.dot(bond[2:], Af)
            b = pb - pa
            d = np.sqrt((b ** 2).sum())
            if d < 1e-10:
                continue
            x0 = 0.0 if hidden0 else r[a0] / d
            x1 = 0.0 if hidden1 else r[a1] / d
            start = pa + b * x0
            seg = b * (1.0 - x0 - x1)
            if bond[2:].any():
                seg *= 0.5
            self.X_bonds[i] = start
            self.B[ncell + i] = start + seg

    def draw(self, status=True):
        """Override draw to clear sphere image references before redrawing."""
        self._sphere_images = []
        super().draw(status)
        self._update_status_info()

    def _toggle_3d_mode(self):
        """Toggle between 3D sphere and flat circle rendering."""
        self._3d_mode = self._3d_mode_var.get()
        self.draw()

    # =========================================================================
    # Status info, axis indicator, context menu, keyboard shortcuts, screenshot
    # =========================================================================

    def _update_status_info(self):
        """Update the status info label in the side panel."""
        if not hasattr(self, '_status_info_label'):
            return
        if self.atoms is None or len(self.atoms) == 0:
            self._status_info_label.config(text="No structure loaded")
            return
        n = len(self.atoms)
        sel = self.images.selected[:n]
        n_sel = int(sel.sum())
        parts = [f"{n} atoms total"]
        if n_sel > 0:
            parts.append(f"{n_sel} selected")
            # Show element breakdown of selection
            syms = np.array(self.atoms.get_chemical_symbols())
            sel_syms = syms[sel]
            from collections import Counter
            counts = Counter(sel_syms)
            breakdown = ", ".join(f"{s}:{c}" for s, c in sorted(counts.items()))
            parts.append(breakdown)
        # Molecule group count
        if "groups" in self.atoms.info:
            ng = len(self.atoms.info["groups"])
            parts.append(f"{ng} groups")
        self._status_info_label.config(text=" | ".join(parts))

    def _setup_context_menu(self):
        """Create the right-click context menu."""
        if tk is None:
            return
        self._context_menu = tk.Menu(self.window.win, tearoff=0)
        self._context_menu.add_command(
            label="Select Molecule Group", command=self._context_select_molecule)
        self._context_menu.add_separator()
        self._context_menu.add_command(
            label="Move...", command=self.move_molecule_dialog)
        self._context_menu.add_command(
            label="Rotate...", command=self.rotate_molecule_dialog)
        self._context_menu.add_separator()
        self._context_menu.add_command(
            label="Delete Selected", command=self.delete_selected_atoms)
        self._context_menu.add_command(
            label="Delete Molecule", command=self.delete_molecule)

        # Bind right-click
        canvas = self.window.canvas
        # macOS uses Button-2 for right-click in some configurations
        canvas.bind("<Button-3>", self._show_context_menu)
        canvas.bind("<Button-2>", self._show_context_menu)

    def _show_context_menu(self, event):
        """Show the right-click context menu at the cursor position."""
        if not hasattr(self, '_context_menu'):
            return
        try:
            self._context_menu.tk_popup(event.x_root, event.y_root)
        finally:
            self._context_menu.grab_release()

    def _context_select_molecule(self):
        """Select all atoms in the molecule group of the most recently clicked atom."""
        if self.atoms is None or len(self.atoms) == 0:
            return
        n = len(self.atoms)
        sel = self.images.selected[:n]
        if not sel.any():
            return
        # Use the first selected atom to find its molecule group
        idx = int(np.where(sel)[0][0])
        if "groups" in self.atoms.info:
            for grp in self.atoms.info["groups"]:
                if idx in grp:
                    self.images.selected[:n] = False
                    for i in grp:
                        if i < n:
                            self.images.selected[i] = True
                    self.draw()
                    return
        # Fallback: use mol_id array
        if "mol_id" in self.atoms.arrays:
            mol_ids = self.atoms.get_array("mol_id")
            target_id = mol_ids[idx]
            self.images.selected[:n] = (mol_ids == target_id)
            self.draw()

    def _setup_keyboard_shortcuts(self):
        """Bind additional keyboard shortcuts."""
        if tk is None:
            return
        win = self.window.win
        win.bind("<Delete>", lambda e: self.delete_selected_atoms())
        win.bind("<BackSpace>", lambda e: self.delete_selected_atoms())
        win.bind("<Escape>", self._deselect_all)
        # Ctrl+A — select all
        win.bind("<Control-a>", self._select_all)
        win.bind("<Control-A>", self._select_all)
        # Ctrl+Shift+S — screenshot
        win.bind("<Control-Shift-S>", lambda e: self.save_screenshot())
        win.bind("<Control-Shift-s>", lambda e: self.save_screenshot())

    def _select_all(self, event=None):
        """Select all visible atoms."""
        if self.atoms is None:
            return
        n = len(self.atoms)
        vis = self.images.visible[:n]
        self.images.selected[:n] = vis
        self.draw()

    def _deselect_all(self, event=None):
        """Deselect all atoms."""
        if self.atoms is None:
            return
        n = len(self.atoms)
        self.images.selected[:n] = False
        self.draw()

    def save_screenshot(self, key=None):
        """Save current view as a PNG file."""
        if tk is None:
            return
        from tkinter import filedialog
        filepath = filedialog.asksaveasfilename(
            parent=self.window.win,
            title="Save Screenshot",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")],
            initialfile="screenshot.png",
        )
        if not filepath:
            return
        try:
            # Use PostScript capture and convert to PNG
            canvas = self.window.canvas
            # Try to use Pillow for direct canvas capture
            try:
                from PIL import Image
                import io
                ps = canvas.postscript(colormode='color')
                # Try ghostscript conversion
                try:
                    img = Image.open(io.BytesIO(ps.encode('utf-8')))
                    img.save(filepath, 'PNG')
                except Exception:
                    # Fallback: save as EPS and note it
                    eps_path = filepath.rsplit('.', 1)[0] + '.eps'
                    with open(eps_path, 'w') as f:
                        f.write(ps)
                    if messagebox:
                        messagebox.showinfo(
                            "Screenshot",
                            f"Saved as EPS (Pillow/Ghostscript not available):\n{eps_path}")
                    return
            except ImportError:
                # No Pillow — save as PostScript
                ps = canvas.postscript(colormode='color')
                eps_path = filepath.rsplit('.', 1)[0] + '.eps'
                with open(eps_path, 'w') as f:
                    f.write(ps)
                if messagebox:
                    messagebox.showinfo(
                        "Screenshot",
                        f"Saved as EPS (install Pillow for PNG):\n{eps_path}")
                return
            if messagebox:
                messagebox.showinfo("Screenshot", f"Saved:\n{filepath}")
        except Exception as e:
            if messagebox:
                messagebox.showerror("Screenshot Error", str(e))

    # =========================================================================
    # File panel (top)
    # =========================================================================

    def _setup_file_panel(self):
        """Add a file toolbar at the top of the window with icon buttons."""
        if tk is None:
            return

        file_bar = tk.Frame(self.window.win, relief=tk.FLAT, borderwidth=0)
        file_bar.pack(side=tk.TOP, fill=tk.X)

        # --- Input ---
        input_frame = ttk.LabelFrame(file_bar, text="Input", padding=(4, 2))
        input_frame.pack(side=tk.LEFT, padx=2, pady=2)

        _input_icons = [
            ("\U0001F4C2", "Load Structure File (Ctrl+L)\nCIF/XYZ/QE/VASP/PDB/GEN\n(atoms only, no molecule metadata)", self.load_structure_dialog),
            ("\U0001F4C1", "Load TOML Project (Ctrl+T)\ninput.toml + data/*.toml\n(with molecule groups & metadata)", self.load_project),
        ]
        for icon, tip, cmd in _input_icons:
            btn = tk.Button(input_frame, text=icon, command=cmd,
                            font=FONTS["icon_large"], width=2, relief=tk.FLAT)
            btn.pack(side=tk.LEFT, padx=2)
            self._add_tooltip(btn, tip)

        # --- Output ---
        output_frame = ttk.LabelFrame(file_bar, text="Output", padding=(4, 2))
        output_frame.pack(side=tk.LEFT, padx=2, pady=2)

        _output_icons = [
            ("\U0001F4BE", "Save TOML Project (Ctrl+K)\nExport input.toml + data/*.toml", self.save_project),
            ("\U0001F4E4", "Export Structure (Ctrl+D)\nOutput as QE / VASP / CIF for calculation", self.export_all_dialog),
            ("\u2B07",     "Save Selection as Molecule TOML\nSave selected atoms to data/*.toml for reuse", self.save_molecule_dialog),
            ("\U0001F4F7", "Screenshot (Ctrl+Shift+S)\nSave current view as PNG", self.save_screenshot),
        ]
        for icon, tip, cmd in _output_icons:
            btn = tk.Button(output_frame, text=icon, command=cmd,
                            font=FONTS["icon_large"], width=2, relief=tk.FLAT)
            btn.pack(side=tk.LEFT, padx=2)
            self._add_tooltip(btn, tip)

        # --- Undo / Redo ---
        edit_frame = ttk.LabelFrame(file_bar, text="Edit", padding=(4, 2))
        edit_frame.pack(side=tk.LEFT, padx=2, pady=2)

        _edit_icons = [
            ("\u21A9", "Undo (Ctrl+Z)", lambda: self.undo()),
            ("\u21AA", "Redo (Ctrl+Y)", lambda: self.redo()),
        ]
        for icon, tip, cmd in _edit_icons:
            btn = tk.Button(edit_frame, text=icon, command=cmd,
                            font=FONTS["icon_large"], width=2, relief=tk.FLAT)
            btn.pack(side=tk.LEFT, padx=2)
            self._add_tooltip(btn, tip)

        # --- Pipeline ---
        pipe_frame = ttk.LabelFrame(file_bar, text="Pipeline", padding=(4, 2))
        pipe_frame.pack(side=tk.LEFT, padx=2, pady=2)

        _pipe_icons = [
            ("\u2699",     "Run Optimize (Ctrl+P)\nVoronoi optimization + auto render", self.run_optimize),
            ("\u25B6",     "Run All\nOptimize \u2192 Render \u2192 Export", self.run_full_pipeline),
            ("\U0001F4CA", "Show Voronoi Results (Ctrl+H)",   self.show_voronoi_results),
            ("\U0001F5BC", "View Voronoi PNGs\nBrowse saved overlay images", lambda: self._view_voronoi_pngs()),
        ]
        for icon, tip, cmd in _pipe_icons:
            btn = tk.Button(pipe_frame, text=icon, command=cmd,
                            font=FONTS["icon_large"], width=2, relief=tk.FLAT)
            btn.pack(side=tk.LEFT, padx=2)
            self._add_tooltip(btn, tip)

    # =========================================================================
    # Bottom panel (Molecule Tools + Structure Tools)
    # =========================================================================

    def _setup_bottom_panel(self):
        """Add a bottom panel with Molecule Tools and Structure Tools icons."""
        if tk is None:
            return

        bottom_bar = tk.Frame(self.window.win, relief=tk.FLAT, borderwidth=0)
        bottom_bar.pack(side=tk.BOTTOM, fill=tk.X,
                        before=self.window.status)

        # --- Molecule Tools ---
        mol_frame = ttk.LabelFrame(bottom_bar, text="Molecule Tools",
                                   padding=(4, 2))
        mol_frame.pack(side=tk.LEFT, padx=2, pady=2)

        _mol_icons = [
            ("E\u20DD",  "Select by Element",   self.select_by_element),
            ("M\u20DD",  "Select by Molecule",  self.select_by_molecule),
            ("\u25CE",    "Select by Sphere",    self.select_by_sphere),
            ("\u2725",    "Move Molecule/Atoms", self.move_molecule_dialog),
            ("\u21BB",    "Rotate Molecule/Atoms", self.rotate_molecule_dialog),
            ("\u2716",    "Delete Selected Atoms", self.delete_selected_atoms),
            ("\u2298",    "Delete Molecule",     self.delete_molecule),
        ]
        for icon, tip, cmd in _mol_icons:
            btn = tk.Button(mol_frame, text=icon, command=cmd,
                            font=FONTS["icon_medium"], width=2, relief=tk.FLAT)
            btn.pack(side=tk.LEFT, padx=1)
            self._add_tooltip(btn, tip)

        ttk.Separator(mol_frame, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=4, pady=2)

        # Drag mode toggle
        self._drag_btn = tk.Button(
            mol_frame, text="\u270B", font=FONTS["icon_medium"], width=2,
            command=self._toggle_drag_mode, relief=tk.RAISED)
        self._drag_btn.pack(side=tk.LEFT, padx=1)
        self._add_tooltip(self._drag_btn, "Drag mode: move atoms by dragging")

        # --- Structure Tools ---
        struct_frame = ttk.LabelFrame(bottom_bar, text="Structure Tools",
                                      padding=(4, 2))
        struct_frame.pack(side=tk.LEFT, padx=2, pady=2)

        _struct_icons = [
            ("\u2795",    "Add Molecule (Ctrl+J)",                        self.add_molecule_dialog),
            ("\u2630",    "List Molecules/Groups",                        self.list_molecules_dialog),
            ("\U0001F504", "Molecule Substitution\n(Replace whole molecule)", self.molecule_substitution_dialog),
            ("\u2747",    "Atom Substitution\n(Replace atoms within molecule)", self.molecule_atom_substitution_dialog),
            ("\u27F3",    "Global Element Substitution\n(Change element across all)", self.global_substitution_dialog),
            ("\u25A6",    "Build Supercell",                              self.supercell_dialog),
            ("\u2B1C",    "Cell Parameters",                              self.cell_parameter_dialog),
        ]
        for icon, tip, cmd in _struct_icons:
            btn = tk.Button(struct_frame, text=icon, command=cmd,
                            font=FONTS["icon_medium"], width=2, relief=tk.FLAT)
            btn.pack(side=tk.LEFT, padx=1)
            self._add_tooltip(btn, tip)

    # =========================================================================
    # View toolbar
    # =========================================================================

    def _setup_toolbar(self):
        """Add a toolbar with View Operations above the status bar.

        Contains: Pan (d-pad) | Orbit (trackball + X/Y/Z) | Zoom | BG
        """
        if tk is None:
            return

        toolbar = tk.Frame(self.window.win)
        toolbar.pack(side=tk.BOTTOM, fill=tk.X, before=self.window.status)

        # ── View Operations ─────────────────────────────────
        view_frame = ttk.LabelFrame(toolbar, text="View Operations",
                                    padding=(4, 2))
        view_frame.pack(fill=tk.X, padx=2, pady=2)

        # --- Pan (directional pad layout) ---
        tk.Label(view_frame, text="Pan:").pack(side=tk.LEFT, padx=(0, 2))
        tpad = tk.Frame(view_frame)
        tpad.pack(side=tk.LEFT, padx=(0, 4))

        # Top row: up
        tpad_top = tk.Frame(tpad)
        tpad_top.pack()
        btn_up = tk.Button(tpad_top, text="\u25B2", width=2)
        btn_up.pack()
        action_up = lambda: self._translate_view(0, +1)
        btn_up.bind("<ButtonPress-1>",
                     lambda e, a=action_up: self._start_repeat(a))
        btn_up.bind("<ButtonRelease-1>", self._stop_repeat)
        self._add_tooltip(btn_up, "Pan view up")

        # Middle row: left, right
        tpad_mid = tk.Frame(tpad)
        tpad_mid.pack()
        btn_left = tk.Button(tpad_mid, text="\u25C0", width=2)
        btn_left.pack(side=tk.LEFT)
        action_left = lambda: self._translate_view(-1, 0)
        btn_left.bind("<ButtonPress-1>",
                       lambda e, a=action_left: self._start_repeat(a))
        btn_left.bind("<ButtonRelease-1>", self._stop_repeat)
        self._add_tooltip(btn_left, "Pan view left")

        btn_right = tk.Button(tpad_mid, text="\u25B6", width=2)
        btn_right.pack(side=tk.LEFT)
        action_right = lambda: self._translate_view(+1, 0)
        btn_right.bind("<ButtonPress-1>",
                        lambda e, a=action_right: self._start_repeat(a))
        btn_right.bind("<ButtonRelease-1>", self._stop_repeat)
        self._add_tooltip(btn_right, "Pan view right")

        # Bottom row: down
        tpad_bot = tk.Frame(tpad)
        tpad_bot.pack()
        btn_down = tk.Button(tpad_bot, text="\u25BC", width=2)
        btn_down.pack()
        action_down = lambda: self._translate_view(0, -1)
        btn_down.bind("<ButtonPress-1>",
                       lambda e, a=action_down: self._start_repeat(a))
        btn_down.bind("<ButtonRelease-1>", self._stop_repeat)
        self._add_tooltip(btn_down, "Pan view down")

        # --- Separator ---
        ttk.Separator(view_frame, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=4, pady=2)

        # --- Orbit: Trackball + axis buttons ---
        tk.Label(view_frame, text="Orbit:").pack(side=tk.LEFT, padx=(0, 2))
        self._setup_trackball(view_frame)

        ttk.Separator(view_frame, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=4, pady=2)

        # Axis buttons (for precise single-axis rotation)
        rpad = tk.Frame(view_frame)
        rpad.pack(side=tk.LEFT, padx=(0, 4))

        for axis_label, axis in [("X", 'x'), ("Y", 'y'), ("Z", 'z')]:
            rrow = tk.Frame(rpad)
            rrow.pack()
            btn_neg = tk.Button(rrow, text=f"{axis_label}\u2212", width=2)
            btn_neg.pack(side=tk.LEFT, padx=1)
            action_neg = lambda a=axis: self._rotate_view(a, -5)
            btn_neg.bind("<ButtonPress-1>",
                          lambda e, a=action_neg: self._start_repeat(a))
            btn_neg.bind("<ButtonRelease-1>", self._stop_repeat)
            self._add_tooltip(btn_neg, f"Orbit -5\u00b0 around {axis_label}")

            btn_pos = tk.Button(rrow, text=f"{axis_label}+", width=2)
            btn_pos.pack(side=tk.LEFT, padx=1)
            action_pos = lambda a=axis: self._rotate_view(a, +5)
            btn_pos.bind("<ButtonPress-1>",
                          lambda e, a=action_pos: self._start_repeat(a))
            btn_pos.bind("<ButtonRelease-1>", self._stop_repeat)
            self._add_tooltip(btn_pos, f"Orbit +5\u00b0 around {axis_label}")

        # --- Separator ---
        ttk.Separator(view_frame, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=4, pady=2)

        # --- Zoom buttons ---
        tk.Label(view_frame, text="Zoom:").pack(side=tk.LEFT, padx=(0, 2))
        for label, factor, tip in [("\u2212", 1 / 1.1, "Zoom out"),
                                    ("+", 1.1, "Zoom in")]:
            btn = tk.Button(view_frame, text=label, width=2)
            btn.pack(side=tk.LEFT, padx=1)
            action = lambda f=factor: self._zoom_view(f)
            btn.bind("<ButtonPress-1>",
                      lambda e, a=action: self._start_repeat(a))
            btn.bind("<ButtonRelease-1>", self._stop_repeat)
            self._add_tooltip(btn, tip)

        # --- Separator ---
        ttk.Separator(view_frame, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=4, pady=2)

        # --- Background color presets ---
        tk.Label(view_frame, text="BG:").pack(side=tk.LEFT, padx=(0, 2))
        bg_presets = [
            ("#ffffff", "White", "W"),
            ("#000000", "Black", "B"),
            ("#333333", "Dark gray", "D"),
            ("#d0d0d0", "Light gray", "L"),
        ]
        for color, tip, label in bg_presets:
            btn = tk.Label(view_frame, text=label, width=2,
                           bg=color, fg="#888" if color != "#000000" else "#777",
                           relief=tk.RAISED, borderwidth=1, cursor="hand2")
            btn.pack(side=tk.LEFT, padx=1)
            btn.bind("<Button-1>",
                     lambda e, c=color: self._set_canvas_bg(c))
            self._add_tooltip(btn, f"Background: {tip}")

        # Custom color picker
        btn_custom = tk.Button(view_frame, text="\u2026", width=2,
                               command=self._pick_canvas_bg)
        btn_custom.pack(side=tk.LEFT, padx=1)
        self._add_tooltip(btn_custom, "Custom background color")

    def _set_canvas_bg(self, color):
        """Set the main canvas background color."""
        try:
            self.window.canvas.configure(bg=color)
            self.window.bg = color
        except Exception:
            pass

    def _pick_canvas_bg(self):
        """Open a color picker for the canvas background."""
        from tkinter import colorchooser
        current = self.window.bg
        result = colorchooser.askcolor(color=current,
                                       title="Canvas Background Color")
        if result and result[1]:
            self._set_canvas_bg(result[1])

    def _toggle_drag_mode_menu(self, key=None):
        """Menu callback for toggling drag mode."""
        self._toggle_drag_mode()

    def _toggle_drag_mode(self):
        """Toggle drag mode for moving atoms with mouse drag."""
        self._drag_mode = not self._drag_mode
        if self._drag_mode:
            self._drag_btn.configure(relief=tk.SUNKEN, bg=COLORS["bg_active"])
            try:
                self.window.canvas.configure(cursor="fleur")
            except Exception:
                pass
        else:
            # Use platform-independent default background color
            if not hasattr(self, '_default_btn_bg'):
                self._default_btn_bg = self._drag_btn.cget("background")
            self._drag_btn.configure(
                relief=tk.RAISED, bg=self._default_btn_bg)
            self._dragging = False
            try:
                self.window.canvas.configure(cursor="")
            except Exception:
                pass

    # =========================================================================
    # Drag mode: move atoms by mouse dragging
    # =========================================================================

    def press(self, event):
        """Override ASE press to handle drag mode.

        Always calls super().press() first so that ASE's internal state
        (self.button, self.xy, self.t0, etc.) is set correctly.  Then,
        if drag mode is active and the click hits a selected atom, we
        start a drag operation.
        """
        super().press(event)
        if self._drag_mode and event.button == self.b1:
            n = len(self.atoms) if self.atoms is not None else 0
            if n > 0 and hasattr(self, 'P') and self.P is not None:
                selected = self.images.selected[:n]
                if selected.any():
                    click = np.array([event.x, event.y])
                    d = self.P[:n] - click
                    r = self.get_covalent_radii()[:n]
                    hit = (d**2).sum(1) < (self.scale * r)**2
                    if (hit & selected).any():
                        self._dragging = True
                        self._drag_start_xy = (event.x, event.y)
                        self._drag_start_positions = self.atoms.positions.copy()
                        self._save_undo_state()

    def move(self, event):
        """Override ASE move to handle atom dragging."""
        if self._drag_mode and self._dragging:
            dx = event.x - self._drag_start_xy[0]
            dy = event.y - self._drag_start_xy[1]
            delta = np.dot(self.axes, np.array([dx, -dy, 0.0])) / self.scale

            n = len(self.atoms)
            selected = self.images.selected[:n]
            indices = np.where(selected)[0]

            pos = self._drag_start_positions.copy()
            pos[indices] += delta
            self.atoms.set_positions(pos)
            self.set_frame(self.frame)
            self.draw(status=False)
            return
        super().move(event)

    def release(self, event):
        """Override ASE release to add Shift+click toggle selection."""
        if self._drag_mode and self._dragging:
            self._dragging = False
            self._drag_start_xy = None
            self._drag_start_positions = None
            self.controller.set_atoms(self.atoms)
            self.draw()
            return

        # Shift+click: toggle-add selection (add if unselected, remove if selected)
        if (event.button == self.b1
                and hasattr(event, 'modifier') and event.modifier == 'shift'
                and event.time < self.t0 + 200):
            selected = self.images.selected
            selected_ordered = self.images.selected_ordered
            d = self.P - self.xy
            r = self.get_covalent_radii()
            hit = np.less((d ** 2).sum(1), (self.scale * r) ** 2)
            for a in self.indices[::-1]:
                if a < len(self.atoms) and hit[a]:
                    # Toggle: deselect if already selected, select if not
                    selected[a] = not selected[a]
                    if selected[a]:
                        selected_ordered += [a]
                    elif a in selected_ordered:
                        selected_ordered = [x for x in selected_ordered if x != a]
                    break
            self.images.selected_ordered = selected_ordered
            self.draw()
            return

        super().release(event)

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
    # Trackball widget for orbit rotation
    # =========================================================================

    def _setup_trackball(self, parent):
        """Create a virtual trackball canvas for interactive orbit rotation."""
        size = 70
        self._tb_size = size
        self._tb_radius = size // 2 - 4  # leave margin for border
        self._tb_center = size // 2

        tb_canvas = tk.Canvas(parent, width=size, height=size,
                              bg=COLORS["tb_bg"], highlightthickness=1,
                              highlightbackground=COLORS["tb_border"],
                              cursor="fleur")
        tb_canvas.pack(side=tk.LEFT, padx=2)
        self._tb_canvas = tb_canvas

        # Draw the trackball circle and crosshair
        self._tb_draw_idle()

        # Bind mouse events
        tb_canvas.bind("<ButtonPress-1>", self._tb_press)
        tb_canvas.bind("<B1-Motion>", self._tb_motion)
        tb_canvas.bind("<ButtonRelease-1>", self._tb_release)

        self._tb_last_xy = None
        self._add_tooltip(tb_canvas, "Drag to orbit (inside=XY, outside=Z)")

    def _tb_draw_idle(self):
        """Draw the trackball in its idle state."""
        c = self._tb_canvas
        c.delete("all")
        cx = cy = self._tb_center
        r = self._tb_radius

        # Outer circle (the "ball")
        c.create_oval(cx - r, cy - r, cx + r, cy + r,
                      outline=COLORS["tb_ball_outline"], width=1,
                      fill=COLORS["tb_ball_fill"])
        # Inner shading circle for 3D look
        sr = int(r * 0.85)
        c.create_oval(cx - sr - 2, cy - sr - 2, cx + sr - 2, cy + sr - 2,
                      outline="", fill=COLORS["tb_inner_fill"])
        # Crosshair
        c.create_line(cx - r + 4, cy, cx + r - 4, cy,
                      fill=COLORS["tb_crosshair"], dash=(2, 2))
        c.create_line(cx, cy - r + 4, cx, cy + r - 4,
                      fill=COLORS["tb_crosshair"], dash=(2, 2))
        # Center dot
        c.create_oval(cx - 2, cy - 2, cx + 2, cy + 2,
                      fill=COLORS["tb_center_dot"], outline="")

    def _tb_press(self, event):
        """Start trackball rotation."""
        self._tb_last_xy = (event.x, event.y)
        self._tb_axes0 = self.axes.copy()

    def _tb_motion(self, event):
        """Handle trackball drag — apply rotation to the view."""
        if self._tb_last_xy is None:
            return

        cx = cy = self._tb_center
        r = self._tb_radius

        # Previous and current positions relative to center
        px, py = self._tb_last_xy[0] - cx, self._tb_last_xy[1] - cy
        nx, ny = event.x - cx, event.y - cy

        # Distance from center (normalised to radius)
        prev_dist = np.sqrt(px * px + py * py) / r
        curr_dist = np.sqrt(nx * nx + ny * ny) / r

        # Sensitivity: degrees per pixel
        sensitivity = 0.8

        dx = (event.x - self._tb_last_xy[0]) * sensitivity
        dy = (event.y - self._tb_last_xy[1]) * sensitivity

        if prev_dist > 1.0 or curr_dist > 1.0:
            # Outside the ball → Z-axis rotation
            # Use the angle between previous and current positions
            angle_prev = np.arctan2(py, px)
            angle_curr = np.arctan2(ny, nx)
            dz = np.degrees(angle_curr - angle_prev)
            self._rotate_view('z', dz)
        else:
            # Inside the ball → X/Y rotation
            # Horizontal drag → rotate around Y, vertical drag → rotate around X
            if abs(dx) > 0.1:
                self._rotate_view('y', dx)
            if abs(dy) > 0.1:
                self._rotate_view('x', -dy)

        self._tb_last_xy = (event.x, event.y)

        # Update visual feedback — show drag indicator
        self._tb_draw_active(event.x, event.y)

    def _tb_release(self, event):
        """End trackball rotation."""
        self._tb_last_xy = None
        self._tb_draw_idle()

    def _tb_draw_active(self, mx, my):
        """Draw the trackball with an active drag indicator."""
        c = self._tb_canvas
        c.delete("all")
        cx = cy = self._tb_center
        r = self._tb_radius

        # Distance from center
        dx, dy = mx - cx, my - cy
        dist = np.sqrt(dx * dx + dy * dy)
        outside = dist > r

        # Outer circle
        ring_color = COLORS["tb_active_z"] if outside else COLORS["tb_active_xy"]
        c.create_oval(cx - r, cy - r, cx + r, cy + r,
                      outline=ring_color, width=2, fill=COLORS["tb_ball_fill"])

        # Inner shading
        sr = int(r * 0.85)
        c.create_oval(cx - sr - 2, cy - sr - 2, cx + sr - 2, cy + sr - 2,
                      outline="", fill=COLORS["tb_inner_fill"])
        # Crosshair
        c.create_line(cx - r + 4, cy, cx + r - 4, cy,
                      fill=COLORS["tb_crosshair"], dash=(2, 2))
        c.create_line(cx, cy - r + 4, cx, cy + r - 4,
                      fill=COLORS["tb_crosshair"], dash=(2, 2))

        # Draw line from center to mouse position (clamped to canvas)
        tx = max(0, min(mx, self._tb_size))
        ty = max(0, min(my, self._tb_size))
        c.create_line(cx, cy, tx, ty, fill=ring_color, width=1)

        # Drag point
        c.create_oval(tx - 3, ty - 3, tx + 3, ty + 3,
                      fill=ring_color, outline="")
        # Center dot
        c.create_oval(cx - 2, cy - 2, cx + 2, cy + 2,
                      fill=COLORS["tb_center_dot"], outline="")

    @staticmethod
    def _add_tooltip(widget, text):
        """Add a hover tooltip to a widget. Returns a setter to update text."""
        tip = None
        text_ref = [text]

        def show(event):
            nonlocal tip
            if tip is not None:
                return
            x = widget.winfo_rootx() + widget.winfo_width() // 2
            y = widget.winfo_rooty() + widget.winfo_height() + 4
            tip = tk.Toplevel(widget)
            tip.wm_overrideredirect(True)
            tip.wm_geometry(f"+{x}+{y}")
            label = tk.Label(
                tip, text=text_ref[0],
                background=COLORS["tooltip_bg"],
                foreground=COLORS["tooltip_fg"],
                relief=tk.SOLID, borderwidth=1,
                font=FONTS["tooltip"],
                padx=SPACING["tooltip_padx"],
                pady=SPACING["tooltip_pady"],
            )
            label.pack()

        def hide(event):
            nonlocal tip
            if tip is not None:
                tip.destroy()
                tip = None

        def set_text(new_text):
            text_ref[0] = new_text

        widget.bind("<Enter>", show)
        widget.bind("<Leave>", hide)
        return set_text

    @staticmethod
    def _create_scrollable_frame(parent, height=None):
        """Create a scrollable frame with canvas and scrollbar.

        Returns (container, inner_frame) where inner_frame is the frame
        to pack widgets into. The container should be packed into the parent.
        """
        container = tk.Frame(parent)
        canvas = tk.Canvas(container, highlightthickness=0)
        if height is not None:
            canvas.configure(height=height)
        scrollbar = ttk.Scrollbar(
            container, orient=tk.VERTICAL, command=canvas.yview)
        inner = tk.Frame(canvas)
        inner.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=inner, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        return container, inner

    # =========================================================================
    # Side panel
    # =========================================================================

    def _setup_side_panel(self):
        """Add a side panel on the right with quick-access view controls."""
        if tk is None:
            return

        # Repack canvas to make room for side panel on the right
        canvas = self.window.canvas
        canvas.pack_forget()

        # Side panel with separator border
        self._side_panel_border = tk.Frame(
            self.window.win, width=252)
        self._side_panel_border.pack(side=tk.RIGHT, fill=tk.Y)
        self._side_panel_border.pack_propagate(False)
        ttk.Separator(self._side_panel_border, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y)
        self._side_panel = tk.Frame(self._side_panel_border, width=250,
                                    relief=tk.FLAT, borderwidth=0)
        self._side_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self._side_panel.pack_propagate(False)
        self._side_panel_visible = True

        # Repack canvas to fill remaining space
        canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # --- Scrollable wrapper (vertical + horizontal) for the entire side panel ---
        sp_canvas = tk.Canvas(self._side_panel, highlightthickness=0)
        sp_vscroll = ttk.Scrollbar(self._side_panel, orient=tk.VERTICAL,
                                    command=sp_canvas.yview)
        sp_hscroll = ttk.Scrollbar(self._side_panel, orient=tk.HORIZONTAL,
                                    command=sp_canvas.xview)
        panel = tk.Frame(sp_canvas)
        self._side_panel_inner = panel
        self._side_panel_canvas = sp_canvas

        panel.bind(
            "<Configure>",
            lambda e: sp_canvas.configure(scrollregion=sp_canvas.bbox("all")))
        sp_canvas.create_window((0, 0), window=panel, anchor="nw")
        sp_canvas.configure(yscrollcommand=sp_vscroll.set,
                            xscrollcommand=sp_hscroll.set)

        sp_hscroll.pack(side=tk.BOTTOM, fill=tk.X)
        sp_vscroll.pack(side=tk.RIGHT, fill=tk.Y)
        sp_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Mouse wheel scrolling on the side panel
        def _on_mousewheel(event):
            try:
                if not sp_canvas.winfo_exists():
                    return
            except tk.TclError:
                return
            # macOS uses event.delta directly; Linux uses Button-4/5
            if event.num == 4:
                sp_canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                sp_canvas.yview_scroll(1, "units")
            else:
                # macOS / Windows: delta is ±120 (Windows) or ±1 (macOS)
                delta = event.delta
                if abs(delta) > 10:
                    delta = delta // 120
                sp_canvas.yview_scroll(-delta, "units")

        self._sp_mousewheel_handler = _on_mousewheel

        # Bind to the canvas and propagate to all children
        def _bind_wheel(widget):
            widget.bind("<MouseWheel>", _on_mousewheel)       # macOS/Windows
            widget.bind("<Button-4>", _on_mousewheel)          # Linux scroll up
            widget.bind("<Button-5>", _on_mousewheel)          # Linux scroll down

        _bind_wheel(sp_canvas)
        _bind_wheel(panel)
        self._sp_bind_wheel = _bind_wheel

        # Rebind wheel to new children when panel is reconfigured
        def _bind_children(event):
            for child in panel.winfo_children():
                _bind_wheel_recursive(child)

        def _bind_wheel_recursive(widget):
            _bind_wheel(widget)
            for child in widget.winfo_children():
                _bind_wheel_recursive(child)

        panel.bind("<Map>", _bind_children)
        self._sp_bind_wheel_recursive = _bind_wheel_recursive

        # --- Title ---
        ttk.Label(panel, text="View Controls",
                  style="SectionHeading.TLabel"
                  ).pack(fill=tk.X, padx=4, pady=(6, 4))

        ttk.Separator(panel, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=2)

        # --- Focus buttons (stacked vertically, fixed width for scrolling) ---
        _btn_width = 28  # fixed character width for buttons
        self._focus_btn = ttk.Button(panel, text="Focus Selection (Ctrl+F)",
                                     width=_btn_width,
                                     command=self.focus_selection)
        self._focus_btn.pack(anchor=tk.W, padx=8, pady=(2, 1))
        self._add_tooltip(self._focus_btn, "Show only selected atoms")
        self._unfocus_btn = ttk.Button(panel, text="Restore All Atoms",
                                       width=_btn_width,
                                       command=self.unfocus, state=tk.DISABLED)
        self._unfocus_btn.pack(anchor=tk.W, padx=8, pady=(1, 2))
        self._add_tooltip(self._unfocus_btn, "Exit focus mode (Ctrl+Shift+F)")

        ttk.Separator(panel, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=2)

        # --- Color mode ---
        ttk.Label(panel, text="Color Mode:").pack(anchor=tk.W, padx=8, pady=(2, 0))
        self._color_mode_var = tk.StringVar(value="element")
        modes = [("By Element", "element"), ("By Molecule", "mol_id"), ("By Tag", "tag")]
        for label, val in modes:
            ttk.Radiobutton(panel, text=label, variable=self._color_mode_var,
                            value=val, command=self._on_color_mode_change
                            ).pack(anchor=tk.W, padx=16)

        ttk.Separator(panel, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=4)

        # --- Show bonds checkbox ---
        show_bonds = False
        try:
            show_bonds = bool(self.window['toggle-show-bonds'])
        except Exception:
            pass
        self._show_bonds_var = tk.BooleanVar(value=show_bonds)
        ttk.Checkbutton(panel, text="Show Bonds (Ctrl+B)",
                         variable=self._show_bonds_var,
                         command=self._on_toggle_bonds
                         ).pack(anchor=tk.W, padx=8, pady=(2, 0))

        # --- 3D Spheres checkbox ---
        self._3d_mode_var = tk.BooleanVar(value=self._3d_mode)
        ttk.Checkbutton(panel, text="3D Spheres",
                         variable=self._3d_mode_var,
                         command=self._toggle_3d_mode
                         ).pack(anchor=tk.W, padx=8, pady=(2, 0))

        ttk.Separator(panel, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=4)

        # --- Element visibility & size (collapsible) ---
        self._elem_collapsed = False
        elem_header = tk.Frame(panel)
        elem_header.pack(fill=tk.X, padx=4, pady=(2, 0))
        self._elem_toggle_label = tk.Label(
            elem_header, text="\u25BC Element Controls:",
            cursor="hand2", font=FONTS["section_heading"])
        self._elem_toggle_label.pack(side=tk.LEFT, padx=4)
        self._elem_toggle_label.bind("<Button-1>", lambda e: self._toggle_section("elem"))

        self._elem_content_frame = tk.Frame(panel)
        self._elem_content_frame.pack(fill=tk.X, padx=4, pady=2)

        self._elem_vis_frame = tk.Frame(self._elem_content_frame)
        self._elem_vis_frame.pack(fill=tk.X)

        self._elem_vis_vars = {}
        self._elem_size_vars = {}
        self._update_element_visibility_panel()

        # Atom size (all atoms) — inside element content
        tk.Label(self._elem_content_frame, text="Size: All Atoms").pack(
            anchor=tk.W, padx=4, pady=(2, 0))
        self._atom_size_var = tk.DoubleVar(value=1.0)
        size_scale = tk.Scale(self._elem_content_frame, from_=0.3, to=3.0,
                              resolution=0.1, orient=tk.HORIZONTAL,
                              variable=self._atom_size_var,
                              command=self._on_atom_size_change, length=200)
        size_scale.pack(padx=4, pady=2)

        ttk.Separator(panel, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=2)

        # --- Molecule visibility (collapsible) ---
        self._mol_collapsed = False
        mol_header = tk.Frame(panel)
        mol_header.pack(fill=tk.X, padx=4, pady=(2, 0))
        self._mol_toggle_label = tk.Label(
            mol_header, text="\u25BC Molecule Controls:",
            cursor="hand2", font=FONTS["section_heading"])
        self._mol_toggle_label.pack(side=tk.LEFT, padx=4)
        self._mol_toggle_label.bind("<Button-1>", lambda e: self._toggle_section("mol"))

        self._mol_content_frame = tk.Frame(panel)
        self._mol_content_frame.pack(fill=tk.X, padx=4, pady=2)

        self._mol_vis_frame = tk.Frame(self._mol_content_frame)
        self._mol_vis_frame.pack(fill=tk.X)

        self._mol_vis_vars = {}
        self._mol_color_btns = {}
        self._mol_colors = {}
        self._update_molecule_visibility_panel()

        ttk.Separator(panel, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=2)

        # --- Atom List (collapsible, collapsed by default) ---
        self._atomlist_collapsed = True
        atomlist_header = tk.Frame(panel)
        atomlist_header.pack(fill=tk.X, padx=4, pady=(2, 0))
        self._atomlist_toggle_label = tk.Label(
            atomlist_header, text="\u25B6 Atom List:",
            cursor="hand2", font=FONTS["section_heading"])
        self._atomlist_toggle_label.pack(side=tk.LEFT, padx=4)
        self._atomlist_toggle_label.bind(
            "<Button-1>", lambda e: self._toggle_section("atomlist"))

        self._atomlist_content_frame = tk.Frame(panel)
        # Start collapsed — do NOT pack content frame

        # Treeview for atom table
        cols = ("idx", "elem", "x", "y", "z", "mol")
        self._atom_tree = ttk.Treeview(
            self._atomlist_content_frame, columns=cols, show="headings",
            height=10, selectmode="browse")
        self._atom_tree.heading("idx", text="#")
        self._atom_tree.heading("elem", text="Elem")
        self._atom_tree.heading("x", text="X")
        self._atom_tree.heading("y", text="Y")
        self._atom_tree.heading("z", text="Z")
        self._atom_tree.heading("mol", text="Mol")
        self._atom_tree.column("idx", width=30, anchor=tk.E)
        self._atom_tree.column("elem", width=35, anchor=tk.CENTER)
        self._atom_tree.column("x", width=55, anchor=tk.E)
        self._atom_tree.column("y", width=55, anchor=tk.E)
        self._atom_tree.column("z", width=55, anchor=tk.E)
        self._atom_tree.column("mol", width=30, anchor=tk.E)

        tree_scroll = ttk.Scrollbar(
            self._atomlist_content_frame, orient=tk.VERTICAL,
            command=self._atom_tree.yview)
        self._atom_tree.configure(yscrollcommand=tree_scroll.set)
        self._atom_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Click row → select that atom in the 3D view
        self._atom_tree.bind("<<TreeviewSelect>>", self._on_atom_tree_select)

        ttk.Separator(panel, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=4, pady=2)

        # --- Status info label (selected atoms / measurement) ---
        self._status_info_label = tk.Label(
            panel, text="",
            font=FONTS["body"],
            wraplength=230, justify=tk.LEFT, anchor=tk.W)
        self._status_info_label.pack(fill=tk.X, padx=8, pady=(2, 4))


    def _on_color_mode_change(self):
        """Handle color mode change from side panel."""
        mode = self._color_mode_var.get()
        self._apply_color_mode(mode)

    def _on_toggle_bonds(self):
        """Toggle bond display from side panel checkbox."""
        show = self._show_bonds_var.get()
        try:
            self.window['toggle-show-bonds'] = show
        except Exception:
            pass
        self.set_frame()

    def toggle_show_bonds(self, key=None):
        """Toggle bond display and sync side panel checkbox."""
        super().toggle_show_bonds(key)
        if hasattr(self, '_show_bonds_var'):
            try:
                state = bool(self.window['toggle-show-bonds'])
            except Exception:
                state = False
            self._show_bonds_var.set(state)

    def _on_atom_size_change(self, value=None):
        """Handle atom size change from side panel."""
        scale = self._atom_size_var.get()
        self.images.atom_scale = scale
        self.draw()

    def _toggle_section(self, section):
        """Toggle collapse/expand of a side panel section."""
        if section == "elem":
            self._elem_collapsed = not self._elem_collapsed
            if self._elem_collapsed:
                self._elem_content_frame.pack_forget()
                self._elem_toggle_label.config(text="\u25B6 Element Controls:")
            else:
                # Re-pack after the header's parent
                self._elem_content_frame.pack(
                    fill=tk.X, padx=4, pady=2,
                    after=self._elem_toggle_label.master)
                self._elem_toggle_label.config(text="\u25BC Element Controls:")
        elif section == "mol":
            self._mol_collapsed = not self._mol_collapsed
            if self._mol_collapsed:
                self._mol_content_frame.pack_forget()
                self._mol_toggle_label.config(text="\u25B6 Molecule Controls:")
            else:
                self._mol_content_frame.pack(
                    fill=tk.X, padx=4, pady=2,
                    after=self._mol_toggle_label.master)
                self._mol_toggle_label.config(text="\u25BC Molecule Controls:")
        elif section == "atomlist":
            self._atomlist_collapsed = not self._atomlist_collapsed
            if self._atomlist_collapsed:
                self._atomlist_content_frame.pack_forget()
                self._atomlist_toggle_label.config(text="\u25B6 Atom List:")
            else:
                self._atomlist_content_frame.pack(
                    fill=tk.X, padx=4, pady=2,
                    after=self._atomlist_toggle_label.master)
                self._atomlist_toggle_label.config(text="\u25BC Atom List:")
                self._update_atom_list()

    def _update_atom_list(self):
        """Populate the atom list Treeview with current atom data."""
        if not hasattr(self, '_atom_tree'):
            return
        tree = self._atom_tree
        # Clear existing rows
        for item in tree.get_children():
            tree.delete(item)
        atoms = self.atoms
        if atoms is None or len(atoms) == 0:
            return
        syms = atoms.get_chemical_symbols()
        pos = atoms.get_positions()
        mol_ids = atoms.arrays.get('mol_id', None)
        for i in range(len(atoms)):
            mol_val = int(mol_ids[i]) if mol_ids is not None else ""
            tree.insert("", tk.END, iid=str(i), values=(
                i, syms[i],
                f"{pos[i, 0]:.3f}", f"{pos[i, 1]:.3f}", f"{pos[i, 2]:.3f}",
                mol_val,
            ))

    def _on_atom_tree_select(self, event=None):
        """When user clicks a row in the atom list, select that atom in 3D view."""
        if not hasattr(self, '_atom_tree'):
            return
        sel = self._atom_tree.selection()
        if not sel:
            return
        try:
            idx = int(sel[0])
        except (ValueError, IndexError):
            return
        n = len(self.atoms)
        if idx < 0 or idx >= n:
            return
        # Select only this atom
        self.images.selected[:n] = False
        self.images.selected[idx] = True
        self.draw()

    def _update_element_visibility_panel(self):
        """Rebuild the element visibility checkboxes, color buttons, and size sliders."""
        if not hasattr(self, '_elem_vis_frame'):
            return
        for w in self._elem_vis_frame.winfo_children():
            w.destroy()
        self._elem_vis_vars.clear()
        self._elem_size_vars.clear()
        if not hasattr(self, '_elem_color_btns'):
            self._elem_color_btns = {}
        self._elem_color_btns.clear()

        atoms = self.atoms
        if atoms is None or len(atoms) == 0:
            return

        from ase.data import atomic_numbers, covalent_radii as default_radii
        from ase.data.colors import jmol_colors

        symbols = sorted(set(atoms.get_chemical_symbols()))
        for sym in symbols:
            row = tk.Frame(self._elem_vis_frame)
            row.pack(fill=tk.X, pady=1)

            # Visibility checkbox
            vis_var = tk.BooleanVar(value=True)
            if hasattr(self, '_hidden_elements') and sym in self._hidden_elements:
                vis_var.set(False)
            self._elem_vis_vars[sym] = vis_var
            cb = ttk.Checkbutton(row, text=f"{sym:>2s}", variable=vis_var,
                                 command=lambda s=sym: self._on_element_vis_change(s))
            cb.pack(side=tk.LEFT, padx=(0, 2))

            # Color swatch (Label for reliable bg color on all platforms)
            z = atomic_numbers[sym]
            current_color = self._get_element_color(z)
            color_lbl = tk.Label(row, text="  ", width=2, bg=current_color,
                                 relief=tk.RAISED, borderwidth=1, cursor="hand2")
            color_lbl.pack(side=tk.LEFT, padx=1)
            color_lbl.bind("<Button-1>",
                           lambda e, z_=z, s=sym: self._on_element_color_pick(z_, s))
            self._elem_color_btns[sym] = color_lbl
            color_tip_setter = self._add_tooltip(color_lbl, f"{sym}: {current_color}")
            if not hasattr(self, '_elem_tip_setters'):
                self._elem_tip_setters = {}
            self._elem_tip_setters[sym] = color_tip_setter

            # Size slider
            default_r = default_radii[z]
            current_r = self.images.covalent_radii[z]
            size_var = tk.DoubleVar(value=round(current_r, 2))
            self._elem_size_vars[sym] = (z, size_var)

            slider = tk.Scale(
                row, from_=0.05, to=default_r * 3, resolution=0.05,
                orient=tk.HORIZONTAL, variable=size_var, length=100,
                showvalue=False,
                command=lambda val, z_=z, s=sym: self._on_element_size_change(z_, val, s),
            )
            slider.pack(side=tk.LEFT, padx=1)
            size_tip_setter = self._add_tooltip(slider, f"{sym} size ({current_r:.2f} Å)")
            if not hasattr(self, '_elem_size_tip_setters'):
                self._elem_size_tip_setters = {}
            self._elem_size_tip_setters[sym] = size_tip_setter

        # Bind mouse wheel to new element widgets
        if hasattr(self, '_sp_bind_wheel_recursive'):
            self._sp_bind_wheel_recursive(self._elem_vis_frame)

        # Also update molecule panel if it exists
        if hasattr(self, '_mol_vis_frame'):
            self._update_molecule_visibility_panel()

    def _get_element_color(self, z):
        """Get the current display color for element with atomic number z."""
        from ase.data.colors import jmol_colors
        return self.colors.get(
            z,
            '#{:02X}{:02X}{:02X}'.format(*(int(x * 255) for x in jmol_colors[z])),
        )

    def _update_color_swatch(self, widget, hex_color, tip_setter=None, tip_text=None):
        """Update a color swatch widget's background and tooltip."""
        try:
            widget.configure(bg=hex_color)
        except Exception:
            pass
        if tip_setter and tip_text:
            tip_setter(tip_text)

    def _on_element_color_pick(self, z, symbol):
        """Open color picker for an element from the side panel."""
        from tkinter import colorchooser
        current = self._get_element_color(z)
        result = colorchooser.askcolor(color=current, title=f"Color for {symbol}")
        if result and result[1]:
            hex_color = result[1]
            self.colors[z] = hex_color
            self._update_color_swatch(
                self._elem_color_btns.get(symbol),
                hex_color,
                getattr(self, '_elem_tip_setters', {}).get(symbol),
                f"{symbol}: {hex_color}",
            )
            self.draw()

    def _on_element_size_change(self, z, value=None, symbol=None):
        """Handle per-element size change from side panel slider."""
        r = float(value) if value is not None else 0.5
        self.images.covalent_radii[z] = r
        if symbol and hasattr(self, '_elem_size_tip_setters'):
            tip_setter = self._elem_size_tip_setters.get(symbol)
            if tip_setter:
                tip_setter(f"{symbol} size ({r:.2f} Å)")
        self.draw()

    def _on_element_vis_change(self, symbol):
        """Handle element visibility toggle from side panel."""
        if not hasattr(self, '_hidden_elements'):
            self._hidden_elements = set()

        var = self._elem_vis_vars.get(symbol)
        if var is None:
            return

        if var.get():
            self._hidden_elements.discard(symbol)
        else:
            self._hidden_elements.add(symbol)

        self._apply_visibility()

    def _update_molecule_visibility_panel(self):
        """Rebuild the molecule visibility controls based on current mol_id."""
        if not hasattr(self, '_mol_vis_frame'):
            return
        for w in self._mol_vis_frame.winfo_children():
            w.destroy()
        self._mol_vis_vars.clear()
        self._mol_color_btns.clear()

        atoms = self.atoms
        if atoms is None or len(atoms) == 0:
            return
        if "mol_id" not in atoms.arrays:
            tk.Label(self._mol_vis_frame, text="(no mol_id)", fg="gray"
                     ).pack(anchor=tk.W)
            return

        mol_ids = atoms.get_array("mol_id")
        site_names = atoms.info.get("site_mol_name", [])
        n_groups = int(mol_ids.max()) + 1
        n_colors = len(self._GROUP_COLORS)

        if not hasattr(self, '_hidden_molecules'):
            self._hidden_molecules = set()

        for g_id in range(n_groups):
            name = site_names[g_id] if g_id < len(site_names) else f"mol_{g_id}"
            n_atoms = int((mol_ids == g_id).sum())

            row = tk.Frame(self._mol_vis_frame)
            row.pack(fill=tk.X, pady=1)

            # Visibility checkbox
            vis_var = tk.BooleanVar(value=(g_id not in self._hidden_molecules))
            self._mol_vis_vars[g_id] = vis_var
            cb = ttk.Checkbutton(row, text=f"{name}({n_atoms})", variable=vis_var,
                                 command=lambda gid=g_id: self._on_mol_vis_change(gid))
            cb.pack(side=tk.LEFT)

            # Color swatch (Label for reliable bg color on all platforms)
            color = self._mol_colors.get(g_id, self._GROUP_COLORS[g_id % n_colors])
            color_lbl = tk.Label(row, text="  ", width=2, bg=color,
                                 relief=tk.RAISED, borderwidth=1, cursor="hand2")
            color_lbl.pack(side=tk.LEFT, padx=1)
            color_lbl.bind("<Button-1>",
                           lambda e, gid=g_id, n=name: self._on_mol_color_pick(gid, n))
            self._mol_color_btns[g_id] = color_lbl
            tip_setter = self._add_tooltip(color_lbl, f"{name}: {color}")
            if not hasattr(self, '_mol_tip_setters'):
                self._mol_tip_setters = {}
            self._mol_tip_setters[g_id] = tip_setter

            # Select button
            sel_btn = ttk.Button(row, text="Select", width=5,
                                 command=lambda gid=g_id: self._on_mol_select(gid))
            sel_btn.pack(side=tk.LEFT, padx=1)
            self._add_tooltip(sel_btn, f"Select {name}")

        # Bind mouse wheel to new molecule widgets
        if hasattr(self, '_sp_bind_wheel_recursive'):
            self._sp_bind_wheel_recursive(self._mol_vis_frame)

    def _on_mol_vis_change(self, group_id):
        """Toggle molecule group visibility (or focus inclusion in focus mode)."""
        # In focus mode, delegate to focus toggle
        if getattr(self, '_focus_active', False):
            self._toggle_focus_molecule(group_id)
            return

        if not hasattr(self, '_hidden_molecules'):
            self._hidden_molecules = set()

        var = self._mol_vis_vars.get(group_id)
        if var is None:
            return

        if var.get():
            self._hidden_molecules.discard(group_id)
        else:
            self._hidden_molecules.add(group_id)

        self._apply_visibility()

    def _apply_visibility(self):
        """Apply all visibility rules (element + molecule hiding).

        Always starts from all-visible and re-applies hiding rules,
        so that unhiding elements/molecules works correctly.
        """
        atoms = self.atoms
        if atoms is None or len(atoms) == 0:
            return
        n = len(atoms)
        visible = np.ones(n, dtype=bool)

        # Apply element visibility
        if hasattr(self, '_hidden_elements') and self._hidden_elements:
            symbols = atoms.get_chemical_symbols()
            for i, s in enumerate(symbols):
                if s in self._hidden_elements:
                    visible[i] = False

        # Apply molecule visibility
        if (hasattr(self, '_hidden_molecules') and self._hidden_molecules
                and "mol_id" in atoms.arrays):
            mol_ids = atoms.get_array("mol_id")
            for i in range(n):
                if int(mol_ids[i]) in self._hidden_molecules:
                    visible[i] = False

        self.images.visible[:n] = visible
        self.draw()

    # Keep backward compatibility alias
    _apply_molecule_visibility = _apply_visibility

    def _on_mol_color_pick(self, group_id, name):
        """Open color picker for a molecule group."""
        from tkinter import colorchooser
        n_colors = len(self._GROUP_COLORS)
        current = self._mol_colors.get(
            group_id, self._GROUP_COLORS[group_id % n_colors])
        result = colorchooser.askcolor(color=current, title=f"Color for {name}")
        if result and result[1]:
            hex_color = result[1]
            self._mol_colors[group_id] = hex_color
            self._update_color_swatch(
                self._mol_color_btns.get(group_id),
                hex_color,
                getattr(self, '_mol_tip_setters', {}).get(group_id),
                f"{name}: {hex_color}",
            )
            # Apply color if in mol_id color mode
            if self.colormode == 'mol_id':
                self._apply_mol_group_colors()

    def _apply_mol_group_colors(self):
        """Apply custom molecule group colors."""
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
        self.colormode = 'mol_id'
        self.colormode_data = (colorscale, 0, max(n_groups - 1, 1))
        self.draw()

    def _on_mol_select(self, group_id):
        """Toggle selection of all atoms in a molecule group."""
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        n = len(self.atoms)
        mask = (mol_ids == group_id)
        # If already exactly this group is selected, deselect all
        if np.array_equal(self.images.selected[:n], mask):
            self.images.selected[:n] = False
        else:
            self.images.selected[:n] = mask
        self.draw()

    def _highlight_molecule_group(self, mol_idx):
        """Select all atoms in a molecule group and switch to mol_id coloring.

        Called after Move/Rotate with Scope='Whole molecule group' so the
        user can see which molecule was affected.
        """
        if self.atoms is None or "mol_id" not in self.atoms.arrays:
            return
        mol_ids = self.atoms.get_array("mol_id")
        n = len(self.atoms)
        self.images.selected[:n] = (mol_ids == mol_idx)
        self._color_by_molecule_group()
        # Sync side panel color mode radio button if it exists
        if hasattr(self, '_color_mode_var'):
            self._color_mode_var.set("mol_id")

    def toggle_side_panel(self, key=None):
        """Toggle the side panel visibility."""
        if not hasattr(self, '_side_panel_border'):
            return
        if self._side_panel_visible:
            self._side_panel_border.pack_forget()
            self._side_panel_visible = False
        else:
            # Repack: unpack canvas, pack side panel border, repack canvas
            canvas = self.window.canvas
            canvas.pack_forget()
            self._side_panel_border.pack(side=tk.RIGHT, fill=tk.Y)
            canvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            self._side_panel_visible = True
            self._update_element_visibility_panel()
            self._update_molecule_visibility_panel()
            if not getattr(self, '_atomlist_collapsed', True):
                self._update_atom_list()

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

        ttk.Button(btn_frame, text="Clear", command=clear_log).pack(side=tk.LEFT, padx=10)
        ttk.Button(btn_frame, text="Close", command=dialog.destroy).pack(side=tk.RIGHT, padx=10)

    # =========================================================================
    # Selection helpers
    # =========================================================================

    def select_by_element(self, key=None):
        """Open dialog to select atoms by element symbol.

        Checkboxes toggle in real time — the main viewer highlights
        matching atoms immediately.
        """
        if tk is None or self.atoms is None:
            return

        from ase.data import chemical_symbols
        n = len(self.atoms)
        syms = self.atoms.get_chemical_symbols()
        elements = sorted(set(chemical_symbols[z] for z in self.atoms.numbers))
        if not elements:
            return

        saved_selection = self.images.selected[:n].copy()
        # Count atoms per element for display
        elem_counts = {}
        for s in syms:
            elem_counts[s] = elem_counts.get(s, 0) + 1

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Select by Element")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)

        tk.Label(dialog, text="Select elements:", padx=10, pady=5
                 ).pack(anchor=tk.W)

        frame = tk.Frame(dialog, padx=10)
        frame.pack(fill=tk.BOTH)

        vars_ = {}
        for sym in elements:
            v = tk.BooleanVar(value=False)
            vars_[sym] = v
            cnt = elem_counts.get(sym, 0)
            ttk.Checkbutton(frame, text=f"{sym}  ({cnt} atoms)",
                            variable=v,
                            command=lambda: _preview()
                            ).pack(anchor=tk.W)

        # --- Info label ---
        info_var = tk.StringVar(value="0 atoms selected")
        tk.Label(dialog, textvariable=info_var, padx=10,
                 font=FONTS["fixed"]).pack(anchor=tk.W, pady=(4, 0))

        mode_var = tk.StringVar(value="replace")
        mode_frame = ttk.LabelFrame(dialog, text="Mode", padding=(10, 5))
        mode_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Radiobutton(mode_frame, text="Replace selection",
                        variable=mode_var, value="replace",
                        command=lambda: _preview()).pack(anchor=tk.W)
        ttk.Radiobutton(mode_frame, text="Add to selection",
                        variable=mode_var, value="add",
                        command=lambda: _preview()).pack(anchor=tk.W)

        def _preview():
            chosen = {s for s, v in vars_.items() if v.get()}
            mask = np.array([s in chosen for s in syms])
            if mode_var.get() == "replace":
                self.images.selected[:n] = mask
            else:
                self.images.selected[:n] = saved_selection | mask
            self.draw()
            info_var.set(f"{int(self.images.selected[:n].sum())} atoms selected")

        def apply():
            dialog.destroy()

        def cancel():
            self.images.selected[:n] = saved_selection
            self.draw()
            dialog.destroy()

        dialog.protocol("WM_DELETE_WINDOW", cancel)

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="Apply", command=apply
                   ).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=cancel
                   ).pack(side=tk.LEFT, padx=5)

    def select_by_molecule(self, key=None):
        """Open dialog to select atoms by molecule group.

        Checkboxes toggle in real time — the main viewer highlights
        matching molecule groups immediately.  Each row shows the
        molecule colour swatch with a click-to-change colour picker.
        """
        if tk is None or self.atoms is None:
            return

        if "mol_id" not in self.atoms.arrays:
            if messagebox:
                messagebox.showinfo("Select by Molecule",
                                    "No molecule groups found (mol_id).")
            return

        n = len(self.atoms)
        mol_ids = self.atoms.get_array("mol_id")
        site_names = self.atoms.info.get("site_mol_name", [])
        groups = self.atoms.info.get("groups", [])
        n_groups = int(mol_ids.max()) + 1
        n_colors = len(self._GROUP_COLORS)
        saved_selection = self.images.selected[:n].copy()

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Select by Molecule")
        dialog.resizable(True, True)
        dialog.transient(self.window.win)

        tk.Label(dialog, text="Select molecule groups:", padx=10, pady=5
                 ).pack(anchor=tk.W)

        frame, inner = self._create_scrollable_frame(dialog, height=200)
        frame.pack(fill=tk.BOTH, expand=True, padx=10)

        vars_ = {}
        color_labels = {}
        for i in range(n_groups):
            name = site_names[i] if i < len(site_names) else f"mol_{i}"
            n_atoms = len(groups[i]) if i < len(groups) else 0
            color = self._mol_colors.get(
                i, self._GROUP_COLORS[i % n_colors])

            row = tk.Frame(inner)
            row.pack(fill=tk.X, pady=1)

            v = tk.BooleanVar(value=False)
            vars_[i] = v

            # Color swatch
            clbl = tk.Label(row, text="  ", width=2, bg=color,
                            relief=tk.RAISED, borderwidth=1, cursor="hand2")
            clbl.pack(side=tk.LEFT, padx=(0, 4))
            color_labels[i] = clbl

            # Click swatch to pick colour
            def _pick_color(event, gid=i, nm=name, lbl=clbl):
                from tkinter import colorchooser
                cur = self._mol_colors.get(
                    gid, self._GROUP_COLORS[gid % n_colors])
                result = colorchooser.askcolor(
                    color=cur, title=f"Color for {nm}",
                    parent=dialog)
                if result and result[1]:
                    hex_c = result[1]
                    self._mol_colors[gid] = hex_c
                    lbl.configure(bg=hex_c)
                    # Sync side panel swatch
                    if gid in self._mol_color_btns:
                        self._mol_color_btns[gid].configure(bg=hex_c)
                    if self.colormode == 'mol_id':
                        self._apply_mol_group_colors()

            clbl.bind("<Button-1>", _pick_color)
            self._add_tooltip(clbl, f"{name}: click to change color")

            # Checkbox
            ttk.Checkbutton(row, text=f"#{i} {name} ({n_atoms} atoms)",
                            variable=v,
                            command=lambda: _preview()
                            ).pack(side=tk.LEFT)

        # --- Info label ---
        info_var = tk.StringVar(value="0 atoms selected")
        tk.Label(dialog, textvariable=info_var, padx=10,
                 font=FONTS["fixed"]).pack(anchor=tk.W, pady=(4, 0))

        mode_var = tk.StringVar(value="replace")
        mode_frame = ttk.LabelFrame(dialog, text="Mode", padding=(10, 5))
        mode_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Radiobutton(mode_frame, text="Replace selection",
                        variable=mode_var, value="replace",
                        command=lambda: _preview()).pack(anchor=tk.W)
        ttk.Radiobutton(mode_frame, text="Add to selection",
                        variable=mode_var, value="add",
                        command=lambda: _preview()).pack(anchor=tk.W)

        def _preview():
            chosen = {g for g, v in vars_.items() if v.get()}
            mask = np.array([int(mol_ids[i]) in chosen for i in range(n)])
            if mode_var.get() == "replace":
                self.images.selected[:n] = mask
            else:
                self.images.selected[:n] = saved_selection | mask
            self.draw()
            count = int(self.images.selected[:n].sum())
            n_mol = len(chosen)
            info_var.set(
                f"{count} atoms selected ({n_mol} molecule group(s))")

        def apply():
            dialog.destroy()

        def cancel():
            self.images.selected[:n] = saved_selection
            self.draw()
            dialog.destroy()

        dialog.protocol("WM_DELETE_WINDOW", cancel)

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="Apply", command=apply
                   ).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=cancel
                   ).pack(side=tk.LEFT, padx=5)

    def select_by_sphere(self, key=None):
        """Open dialog to select atoms within a radius of a selected atom.

        Shows a real-time preview: as the radius changes, the main viewer
        highlights the atoms inside the sphere and an info label shows the
        count and element breakdown.
        """
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
        saved_selection = self.images.selected[:n].copy()
        symbols = self.atoms.get_chemical_symbols()
        dists = np.linalg.norm(self.atoms.positions - center_pos, axis=1)

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Select by Sphere")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)

        center_sym = symbols[center_idx]
        tk.Label(dialog,
                 text=f"Center: atom #{center_idx} ({center_sym})",
                 padx=10, pady=5).pack(anchor=tk.W)

        # --- Radius control ---
        r_frame = tk.Frame(dialog, padx=10)
        r_frame.pack(fill=tk.X)
        tk.Label(r_frame, text="Radius (\u00c5):").pack(side=tk.LEFT)
        r_var = tk.DoubleVar(value=3.0)
        r_scale = tk.Scale(r_frame, from_=0.5, to=20.0, resolution=0.1,
                           orient=tk.HORIZONTAL, variable=r_var,
                           command=lambda _: _preview(), length=200)
        r_scale.pack(side=tk.LEFT, padx=4)
        r_entry = tk.Entry(r_frame, textvariable=r_var, width=6,
                           justify=tk.RIGHT)
        r_entry.pack(side=tk.LEFT, padx=2)
        r_entry.bind("<Return>", lambda e: _preview())
        r_entry.bind("<FocusOut>", lambda e: _preview())

        # --- Info label ---
        info_var = tk.StringVar(value="")
        info_label = tk.Label(dialog, textvariable=info_var, justify=tk.LEFT,
                              padx=10, font=FONTS["fixed"])
        info_label.pack(anchor=tk.W, padx=10, pady=(4, 0))

        # --- Mode ---
        mode_var = tk.StringVar(value="replace")
        mode_frame = ttk.LabelFrame(dialog, text="Mode", padding=(10, 5))
        mode_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Radiobutton(mode_frame, text="Replace selection",
                        variable=mode_var, value="replace").pack(anchor=tk.W)
        ttk.Radiobutton(mode_frame, text="Add to selection",
                        variable=mode_var, value="add").pack(anchor=tk.W)

        def _preview():
            radius = r_var.get()
            mask = dists <= radius
            # Show preview in main viewer
            self.images.selected[:n] = mask
            self.draw()
            # Update info
            count = int(mask.sum())
            elem_counts = {}
            for i, m in enumerate(mask):
                if m:
                    s = symbols[i]
                    elem_counts[s] = elem_counts.get(s, 0) + 1
            parts = [f"{s}: {c}" for s, c in
                     sorted(elem_counts.items(), key=lambda x: -x[1])]
            info_var.set(
                f"{count} atoms in sphere (r={radius:.1f} \u00c5)\n"
                + "  ".join(parts)
            )

        def apply():
            radius = r_var.get()
            mask = dists <= radius
            if mode_var.get() == "replace":
                self.images.selected[:n] = mask
            else:
                self.images.selected[:n] = saved_selection | mask
            self.draw()
            dialog.destroy()

        def cancel():
            self.images.selected[:n] = saved_selection
            self.draw()
            dialog.destroy()

        dialog.protocol("WM_DELETE_WINDOW", cancel)

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="Select", command=apply
                   ).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=cancel
                   ).pack(side=tk.LEFT, padx=5)

        # Initial preview
        _preview()

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

        self.images.initialize([new_atoms])
        self.images.selected[:] = False
        self.set_frame(self.frame)
        self.draw()

    # =========================================================================
    # Move/rotate molecule
    # =========================================================================

    def _get_target_indices(self, scope, mol_idx, fallback_indices):
        """Return atom indices based on scope selection."""
        if scope == "molecule" and mol_idx is not None:
            groups = self.atoms.info.get("groups", [])
            if mol_idx < len(groups):
                return groups[mol_idx]
        return list(fallback_indices)

    def move_molecule_dialog(self, key=None):
        """Open interactive dialog to translate selected molecule/atoms.

        Provides both direct entry fields and +/- buttons for each axis.
        Changes are previewed in real time; Apply confirms, Cancel reverts.
        """
        if tk is None or self.atoms is None:
            return

        n = len(self.atoms)
        indices = np.arange(n)[self.images.selected[:n]]
        if len(indices) == 0:
            if messagebox:
                messagebox.showinfo("Move Molecule", "Select atoms first.")
            return

        mol_idx = None
        if "mol_id" in self.atoms.arrays:
            mol_ids = self.atoms.get_array("mol_id")
            mol_idx = int(mol_ids[indices[0]])

        # Save initial state for preview / cancel / scope switch
        saved_positions = self.atoms.positions.copy()
        saved_selection = self.images.selected[:n].copy()

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Move Molecule/Atoms")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)

        # --- Scope ---
        scope_var = tk.StringVar(value="selected")

        def _on_scope_change(*_args):
            if scope_var.get() == "molecule" and mol_idx is not None:
                self._highlight_molecule_group(mol_idx)
            else:
                self.images.selected[:n] = saved_selection
                self._restore_original_coloring()
            self.draw()

        scope_var.trace_add("write", _on_scope_change)

        scope_frame = ttk.LabelFrame(dialog, text="Scope", padding=(10, 4))
        scope_frame.pack(fill=tk.X, padx=10, pady=(8, 4))
        if mol_idx is not None:
            ttk.Radiobutton(scope_frame, text="Whole molecule group",
                            variable=scope_var, value="molecule"
                            ).pack(anchor=tk.W)
        ttk.Radiobutton(scope_frame, text="Selected atoms only",
                        variable=scope_var, value="selected"
                        ).pack(anchor=tk.W)

        # --- Translation controls ---
        ctrl_frame = ttk.LabelFrame(dialog, text="Translation (Å)",
                                    padding=(10, 4))
        ctrl_frame.pack(fill=tk.X, padx=10, pady=4)

        vars_ = []
        step_var = tk.DoubleVar(value=0.1)

        def nudge(axis_idx, sign):
            """Increment axis value by ±step and preview."""
            step = step_var.get()
            old = vars_[axis_idx].get()
            vars_[axis_idx].set(round(old + sign * step, 6))
            preview()

        for i, label in enumerate(("X:", "Y:", "Z:")):
            row = tk.Frame(ctrl_frame)
            row.pack(fill=tk.X, pady=1)
            tk.Label(row, text=label, width=2).pack(side=tk.LEFT)
            tk.Button(row, text="\u2212", width=2,
                      command=lambda ax=i: nudge(ax, -1),
                      repeatdelay=400, repeatinterval=50
                      ).pack(side=tk.LEFT, padx=1)
            v = tk.DoubleVar(value=0.0)
            vars_.append(v)
            entry = tk.Entry(row, textvariable=v, width=10, justify=tk.RIGHT)
            entry.pack(side=tk.LEFT, padx=2)
            entry.bind("<Return>", lambda e: preview())
            entry.bind("<FocusOut>", lambda e: preview())
            tk.Button(row, text="+", width=2,
                      command=lambda ax=i: nudge(ax, +1),
                      repeatdelay=400, repeatinterval=50
                      ).pack(side=tk.LEFT, padx=1)

        # Step size
        step_frame = tk.Frame(ctrl_frame)
        step_frame.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(step_frame, text="Step:").pack(side=tk.LEFT)
        for val in (0.01, 0.1, 0.5, 1.0):
            ttk.Radiobutton(step_frame, text=str(val), variable=step_var,
                            value=val).pack(side=tk.LEFT)
        ttk.Label(step_frame, text="Å").pack(side=tk.LEFT)

        # --- Preview logic ---
        def preview():
            atom_indices = self._get_target_indices(
                scope_var.get(), mol_idx, indices)
            vec = np.array([v.get() for v in vars_])
            pos = saved_positions.copy()
            for idx in atom_indices:
                pos[idx] += vec
            self.atoms.set_positions(pos)
            self.set_frame(self.frame)
            self.draw()

        def apply():
            self._save_undo_state_with(saved_positions)
            is_molecule = scope_var.get() == "molecule" and mol_idx is not None
            if is_molecule:
                self._highlight_molecule_group(mol_idx)
            self.draw()
            dialog.destroy()

        def cancel():
            self.atoms.set_positions(saved_positions)
            self.images.selected[:n] = saved_selection
            self._restore_original_coloring()
            self.set_frame(self.frame)
            self.draw()
            dialog.destroy()

        def reset():
            for v in vars_:
                v.set(0.0)
            self.atoms.set_positions(saved_positions)
            self.set_frame(self.frame)
            self.draw()

        dialog.protocol("WM_DELETE_WINDOW", cancel)

        # --- Buttons ---
        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=8)
        ttk.Button(btn_frame, text="Apply", width=8,
                   command=apply).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Reset", width=8,
                   command=reset).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Cancel", width=8,
                   command=cancel).pack(side=tk.LEFT, padx=4)

    def rotate_molecule_dialog(self, key=None):
        """Open interactive dialog to rotate selected molecule/atoms.

        Provides both direct entry fields and +/- buttons for each Euler
        angle.  Changes are previewed in real time.
        """
        if tk is None or self.atoms is None:
            return

        n = len(self.atoms)
        indices = np.arange(n)[self.images.selected[:n]]
        if len(indices) == 0:
            if messagebox:
                messagebox.showinfo("Rotate Molecule", "Select atoms first.")
            return

        mol_idx = None
        if "mol_id" in self.atoms.arrays:
            mol_ids = self.atoms.get_array("mol_id")
            mol_idx = int(mol_ids[indices[0]])

        saved_positions = self.atoms.positions.copy()
        saved_selection = self.images.selected[:n].copy()

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Rotate Molecule/Atoms")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)

        # --- Scope ---
        scope_var = tk.StringVar(value="selected")

        def _on_scope_change(*_args):
            if scope_var.get() == "molecule" and mol_idx is not None:
                self._highlight_molecule_group(mol_idx)
            else:
                self.images.selected[:n] = saved_selection
                self._restore_original_coloring()
            self.draw()

        scope_var.trace_add("write", _on_scope_change)

        scope_frame = ttk.LabelFrame(dialog, text="Scope", padding=(10, 4))
        scope_frame.pack(fill=tk.X, padx=10, pady=(8, 4))
        if mol_idx is not None:
            ttk.Radiobutton(scope_frame, text="Whole molecule group",
                            variable=scope_var, value="molecule"
                            ).pack(anchor=tk.W)
        ttk.Radiobutton(scope_frame, text="Selected atoms only",
                        variable=scope_var, value="selected"
                        ).pack(anchor=tk.W)

        # --- Rotation controls ---
        ctrl_frame = ttk.LabelFrame(dialog, text="Euler Angles (deg, XYZ)",
                                    padding=(10, 4))
        ctrl_frame.pack(fill=tk.X, padx=10, pady=4)

        # Trackball + numeric controls side by side
        ctrl_inner = tk.Frame(ctrl_frame)
        ctrl_inner.pack(fill=tk.X)

        # --- Trackball for interactive rotation ---
        tb_frame = tk.Frame(ctrl_inner)
        tb_frame.pack(side=tk.LEFT, padx=(0, 8))

        tb_size = 90
        tb_radius = tb_size // 2 - 4
        tb_center = tb_size // 2
        rot_tb_last = [None]  # mutable for closure

        tb_canvas = tk.Canvas(tb_frame, width=tb_size, height=tb_size,
                              bg=COLORS["tb_bg"], highlightthickness=1,
                              highlightbackground=COLORS["tb_border"],
                              cursor="fleur")
        tb_canvas.pack()
        tk.Label(tb_frame, text="Drag to rotate",
                 font=FONTS["icon_small"]).pack()

        def _rot_tb_draw_idle():
            tb_canvas.delete("all")
            cx = cy = tb_center
            r = tb_radius
            tb_canvas.create_oval(cx - r, cy - r, cx + r, cy + r,
                                  outline=COLORS["tb_ball_outline"], width=1,
                                  fill=COLORS["tb_ball_fill"])
            sr = int(r * 0.85)
            tb_canvas.create_oval(cx - sr - 2, cy - sr - 2,
                                  cx + sr - 2, cy + sr - 2,
                                  outline="", fill=COLORS["tb_inner_fill"])
            tb_canvas.create_line(cx - r + 4, cy, cx + r - 4, cy,
                                  fill=COLORS["tb_crosshair"], dash=(2, 2))
            tb_canvas.create_line(cx, cy - r + 4, cx, cy + r - 4,
                                  fill=COLORS["tb_crosshair"], dash=(2, 2))
            tb_canvas.create_oval(cx - 2, cy - 2, cx + 2, cy + 2,
                                  fill=COLORS["tb_center_dot"], outline="")

        def _rot_tb_draw_active(mx, my):
            tb_canvas.delete("all")
            cx = cy = tb_center
            r = tb_radius
            dx, dy = mx - cx, my - cy
            dist = np.sqrt(dx * dx + dy * dy)
            outside = dist > r
            ring_color = (COLORS["tb_active_z"] if outside
                          else COLORS["tb_active_xy"])
            tb_canvas.create_oval(cx - r, cy - r, cx + r, cy + r,
                                  outline=ring_color, width=2,
                                  fill=COLORS["tb_ball_fill"])
            sr = int(r * 0.85)
            tb_canvas.create_oval(cx - sr - 2, cy - sr - 2,
                                  cx + sr - 2, cy + sr - 2,
                                  outline="", fill=COLORS["tb_inner_fill"])
            tb_canvas.create_line(cx - r + 4, cy, cx + r - 4, cy,
                                  fill=COLORS["tb_crosshair"], dash=(2, 2))
            tb_canvas.create_line(cx, cy - r + 4, cx, cy + r - 4,
                                  fill=COLORS["tb_crosshair"], dash=(2, 2))
            tx = max(0, min(mx, tb_size))
            ty = max(0, min(my, tb_size))
            tb_canvas.create_line(cx, cy, tx, ty, fill=ring_color, width=1)
            tb_canvas.create_oval(tx - 3, ty - 3, tx + 3, ty + 3,
                                  fill=ring_color, outline="")
            tb_canvas.create_oval(cx - 2, cy - 2, cx + 2, cy + 2,
                                  fill=COLORS["tb_center_dot"], outline="")

        def _rot_tb_press(event):
            rot_tb_last[0] = (event.x, event.y)

        def _rot_tb_motion(event):
            if rot_tb_last[0] is None:
                return
            cx = cy = tb_center
            r = tb_radius
            px, py = rot_tb_last[0][0] - cx, rot_tb_last[0][1] - cy
            nx, ny = event.x - cx, event.y - cy
            prev_dist = np.sqrt(px * px + py * py) / r
            curr_dist = np.sqrt(nx * nx + ny * ny) / r
            sensitivity = 0.8
            ddx = (event.x - rot_tb_last[0][0]) * sensitivity
            ddy = (event.y - rot_tb_last[0][1]) * sensitivity
            if prev_dist > 1.0 or curr_dist > 1.0:
                angle_prev = np.arctan2(py, px)
                angle_curr = np.arctan2(ny, nx)
                dz = np.degrees(angle_curr - angle_prev)
                vars_[2].set(round(vars_[2].get() + dz, 2))
            else:
                if abs(ddy) > 0.1:
                    vars_[0].set(round(vars_[0].get() - ddy, 2))
                if abs(ddx) > 0.1:
                    vars_[1].set(round(vars_[1].get() + ddx, 2))
            rot_tb_last[0] = (event.x, event.y)
            _rot_tb_draw_active(event.x, event.y)
            preview()

        def _rot_tb_release(event):
            rot_tb_last[0] = None
            _rot_tb_draw_idle()

        tb_canvas.bind("<ButtonPress-1>", _rot_tb_press)
        tb_canvas.bind("<B1-Motion>", _rot_tb_motion)
        tb_canvas.bind("<ButtonRelease-1>", _rot_tb_release)
        _rot_tb_draw_idle()

        # --- Numeric controls ---
        numeric_frame = tk.Frame(ctrl_inner)
        numeric_frame.pack(side=tk.LEFT, fill=tk.X, expand=True)

        vars_ = []
        step_var = tk.DoubleVar(value=5.0)

        def nudge(axis_idx, sign):
            step = step_var.get()
            old = vars_[axis_idx].get()
            vars_[axis_idx].set(round(old + sign * step, 4))
            preview()

        for i, label in enumerate(("\u03b1 (X):", "\u03b2 (Y):", "\u03b3 (Z):")):
            row = tk.Frame(numeric_frame)
            row.pack(fill=tk.X, pady=1)
            tk.Label(row, text=label, width=6).pack(side=tk.LEFT)
            tk.Button(row, text="\u2212", width=2,
                      command=lambda ax=i: nudge(ax, -1),
                      repeatdelay=400, repeatinterval=50
                      ).pack(side=tk.LEFT, padx=1)
            v = tk.DoubleVar(value=0.0)
            vars_.append(v)
            entry = tk.Entry(row, textvariable=v, width=10, justify=tk.RIGHT)
            entry.pack(side=tk.LEFT, padx=2)
            entry.bind("<Return>", lambda e: preview())
            entry.bind("<FocusOut>", lambda e: preview())
            tk.Button(row, text="+", width=2,
                      command=lambda ax=i: nudge(ax, +1),
                      repeatdelay=400, repeatinterval=50
                      ).pack(side=tk.LEFT, padx=1)

        # Step size
        step_frame = tk.Frame(numeric_frame)
        step_frame.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(step_frame, text="Step:").pack(side=tk.LEFT)
        for val in (1, 5, 15, 45, 90):
            ttk.Radiobutton(step_frame, text=f"{val}\u00b0", variable=step_var,
                            value=float(val)).pack(side=tk.LEFT)

        # --- Preview logic ---
        def preview():
            from scipy.spatial.transform import Rotation as Rot
            atom_indices = self._get_target_indices(
                scope_var.get(), mol_idx, indices)
            a, b, g = [v.get() for v in vars_]

            pos = saved_positions.copy()
            target_pos = pos[atom_indices]
            com = target_pos.mean(axis=0)

            rot = Rot.from_euler('xyz', [a, b, g], degrees=True)
            rotated = rot.apply(target_pos - com) + com

            for j, ai in enumerate(atom_indices):
                pos[ai] = rotated[j]
            self.atoms.set_positions(pos)
            self.set_frame(self.frame)
            self.draw()

        def apply():
            self._save_undo_state_with(saved_positions)
            is_molecule = scope_var.get() == "molecule" and mol_idx is not None
            if is_molecule:
                self._highlight_molecule_group(mol_idx)
            self.draw()
            dialog.destroy()

        def cancel():
            self.atoms.set_positions(saved_positions)
            self.images.selected[:n] = saved_selection
            self._restore_original_coloring()
            self.set_frame(self.frame)
            self.draw()
            dialog.destroy()

        def reset():
            for v in vars_:
                v.set(0.0)
            self.atoms.set_positions(saved_positions)
            self.set_frame(self.frame)
            self.draw()

        dialog.protocol("WM_DELETE_WINDOW", cancel)

        # --- Buttons ---
        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=8)
        ttk.Button(btn_frame, text="Apply", width=8,
                   command=apply).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Reset", width=8,
                   command=reset).pack(side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Cancel", width=8,
                   command=cancel).pack(side=tk.LEFT, padx=4)

    # =========================================================================
    # Focus mode
    # =========================================================================

    def focus_selection(self, key=None):
        """Focus on selected atoms/molecules with auto-toggle preview.

        The preview window is non-modal.  Clicking atoms in the main view
        automatically toggles them: focused atoms are removed, unfocused
        atoms are added.  The main view highlights the current focus set.
        """
        if self.atoms is None or len(self.atoms) == 0:
            return

        n = len(self.atoms)
        selected = self.images.selected[:n]
        indices = np.arange(n)[selected]
        if len(indices) == 0:
            if messagebox:
                messagebox.showinfo("Focus", "Select atoms first.")
            return

        has_mol_id = "mol_id" in self.atoms.arrays
        if has_mol_id:
            mol_ids = self.atoms.get_array("mol_id")
            site_names = self.atoms.info.get("site_mol_name", [])
            n_groups = int(mol_ids.max()) + 1
            selected_groups = set(int(mol_ids[i]) for i in indices)
        else:
            mol_ids = None
            n_groups = 0
            selected_groups = set()

        # If no molecules or only one group, skip dialog
        if n_groups <= 1:
            self._execute_focus(indices, set(), "atoms")
            return

        # Save backup for unfocus (including original selection)
        backup = {
            'visible': self.images.visible[:n].copy(),
            'selected': self.images.selected[:n].copy(),
            'center': self.center.copy() if hasattr(self, 'center') else None,
            'scale': getattr(self, 'scale', None),
        }

        # --- Focus set: the set of atom indices to show ---
        # Initialize with exactly the selected atoms (not expanded to molecules)
        focus_set = set(int(i) for i in indices)

        # Atom colors for preview
        from ase.data.colors import jmol_colors
        from ase.data import covalent_radii as default_radii
        positions = self.atoms.positions
        numbers = self.atoms.numbers

        # --- Build preview window (non-modal) ---
        dialog = tk.Toplevel(self.window.win)
        dialog.title("Focus Selection — Preview")

        # Matplotlib 3D preview (top)
        import matplotlib
        matplotlib.use('TkAgg')
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        fig = Figure(figsize=(5, 4), dpi=90)
        ax = fig.add_subplot(111, projection='3d')
        canvas_widget = FigureCanvasTkAgg(fig, master=dialog)
        canvas_widget.get_tk_widget().pack(fill=tk.BOTH, expand=True,
                                           padx=5, pady=(5, 0))

        # --- Controls (bottom) ---
        ctrl_frame = tk.Frame(dialog)
        ctrl_frame.pack(fill=tk.X, padx=5, pady=5)

        # Status label
        status_label = tk.Label(ctrl_frame, text="", anchor=tk.W)
        status_label.pack(fill=tk.X, padx=5, pady=(0, 4))

        # Unit selector (Atom / Molecule)
        unit_frame = ttk.LabelFrame(ctrl_frame, text="Unit", padding=(8, 4))
        unit_frame.pack(fill=tk.X, padx=5, pady=(0, 5))

        unit_var = tk.StringVar(value="atom")
        ttk.Radiobutton(unit_frame, text="Atom", variable=unit_var,
                         value="atom").pack(side=tk.LEFT)
        ttk.Radiobutton(unit_frame, text="Molecule", variable=unit_var,
                         value="molecule").pack(side=tk.LEFT, padx=(10, 0))

        # Molecule checkboxes
        mol_frame = ttk.LabelFrame(ctrl_frame, text="Molecules", padding=(8, 4))
        mol_frame.pack(fill=tk.X, padx=5, pady=(0, 5))

        n_colors = len(self._GROUP_COLORS)
        mol_vars = {}
        for g_id in range(n_groups):
            name = site_names[g_id] if g_id < len(site_names) else f"mol_{g_id}"
            n_atoms = int((mol_ids == g_id).sum())
            group_atoms = set(int(j) for j in np.where(mol_ids == g_id)[0])
            var = tk.BooleanVar(value=(group_atoms.issubset(focus_set)
                                       and len(group_atoms) > 0))
            mol_vars[g_id] = var
            row = tk.Frame(mol_frame)
            row.pack(fill=tk.X, anchor=tk.W)
            color = self._mol_colors.get(g_id, self._GROUP_COLORS[g_id % n_colors])
            tk.Label(row, text="  ", bg=color, width=2,
                     relief=tk.RAISED, borderwidth=1).pack(side=tk.LEFT, padx=(0, 4))
            ttk.Checkbutton(row, text=f"{name} ({n_atoms} atoms)",
                            variable=var,
                            command=lambda gid=g_id: on_mol_checkbox(gid)
                            ).pack(side=tk.LEFT)

        # Help text
        tk.Label(ctrl_frame,
                 text="Click atoms in main view to toggle focus\n"
                      "(focused \u2192 remove, unfocused \u2192 add)",
                 font=FONTS["body"],
                 justify=tk.LEFT).pack(padx=5, pady=(0, 4))

        def sync_main_view():
            """Highlight focus set as selected atoms in main view."""
            sel = np.zeros(len(self.images.selected), dtype=bool)
            for i in focus_set:
                if 0 <= i < n:
                    sel[i] = True
            self.images.selected[:] = sel

        def on_mol_checkbox(g_id):
            """Handle molecule checkbox toggle."""
            var = mol_vars[g_id]
            group_atoms = set(int(i) for i in np.where(mol_ids == g_id)[0])
            if var.get():
                focus_set.update(group_atoms)
            else:
                focus_set.difference_update(group_atoms)
            sync_main_view()
            self.draw()
            draw_preview()

        def sync_mol_checkboxes():
            """Update molecule checkboxes to reflect focus_set."""
            for g_id, var in mol_vars.items():
                group_atoms = set(int(i) for i in np.where(mol_ids == g_id)[0])
                var.set(group_atoms.issubset(focus_set) and len(group_atoms) > 0)

        def on_click_detected():
            """Handle a click in the main view: auto-toggle.

            ASE-GUI click patterns vs. expected (focus_set as selected):
            - Regular click on atom A: cur = {A} (clears others, selects A)
            - Ctrl+click on atom A: toggles A in cur (keeps others)
            - Shift+click on atom A: adds A to cur (keeps others)
            - Spurious (window focus loss etc.): cur = {} (everything cleared)

            We detect the clicked atom and toggle it in/out of focus_set.
            """
            cur = self.images.selected[:n]
            # Build expected selection from focus_set
            expected = np.zeros(n, dtype=bool)
            for i in focus_set:
                if 0 <= i < n:
                    expected[i] = True

            newly_on = np.where(cur & ~expected)[0]
            newly_off = np.where(~cur & expected)[0]
            total_diff = len(newly_on) + len(newly_off)

            if total_diff == 0:
                return

            clicked_atom = None

            if total_diff == 1:
                # Ctrl+click or Shift+click: exactly one atom toggled
                if len(newly_on) == 1:
                    clicked_atom = int(newly_on[0])
                else:
                    clicked_atom = int(newly_off[0])
            elif len(newly_on) == 1:
                # Regular click on unfocused atom: A added, focus atoms cleared
                clicked_atom = int(newly_on[0])
            elif len(newly_on) == 0 and int(cur.sum()) == 1:
                # Regular click on focused atom: only that atom stays selected
                clicked_atom = int(np.where(cur)[0][0])

            if clicked_atom is None:
                # Can't determine clicked atom (e.g. all cleared) → re-sync
                sync_main_view()
                return

            unit = unit_var.get()

            if unit == "molecule" and mol_ids is not None:
                # Expand to whole molecule group
                g_id = int(mol_ids[clicked_atom])
                clicked = set(int(j) for j in np.where(mol_ids == g_id)[0])
            else:
                clicked = {clicked_atom}

            # Auto-toggle: if clicked atom was in focus, remove; else add
            if clicked_atom in focus_set:
                focus_set.difference_update(clicked)
            else:
                focus_set.update(clicked)

            sync_mol_checkboxes()
            sync_main_view()
            draw_preview()

        def draw_preview():
            """Redraw the 3D preview."""
            ax.clear()
            visible = np.zeros(n, dtype=bool)
            for i in focus_set:
                if 0 <= i < n:
                    visible[i] = True
            vis_idx = np.where(visible)[0]

            status_label.configure(text=f"Focus: {len(vis_idx)} / {n} atoms")

            if len(vis_idx) == 0:
                ax.set_title("No atoms in focus")
                canvas_widget.draw_idle()
                return

            # Visible atoms (colored)
            vis_pos = positions[vis_idx]
            vis_colors = [jmol_colors[numbers[i]] for i in vis_idx]
            vis_sizes = [default_radii[numbers[i]] * 120 for i in vis_idx]
            ax.scatter(vis_pos[:, 0], vis_pos[:, 1], vis_pos[:, 2],
                       c=vis_colors, s=vis_sizes, alpha=1.0,
                       edgecolors='black', linewidths=0.3, depthshade=True)

            # Hidden atoms (faded)
            hid_idx = np.where(~visible)[0]
            if len(hid_idx) > 0:
                hid_pos = positions[hid_idx]
                hid_sizes = [default_radii[numbers[i]] * 60 for i in hid_idx]
                ax.scatter(hid_pos[:, 0], hid_pos[:, 1], hid_pos[:, 2],
                           c='lightgray', s=hid_sizes, alpha=0.25,
                           edgecolors='none', depthshade=True)

            # Auto-scale
            mid = positions.mean(axis=0)
            max_range = (positions.max(axis=0) - positions.min(axis=0)).max() / 2
            if max_range < 0.1:
                max_range = 5.0
            ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
            ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
            ax.set_zlim(mid[2] - max_range, mid[2] + max_range)
            ax.set_xlabel('X (Å)')
            ax.set_ylabel('Y (Å)')
            ax.set_zlabel('Z (Å)')
            ax.set_title(f"{len(vis_idx)} atoms shown")
            canvas_widget.draw_idle()

        # --- Hook into draw() to detect selection changes ---
        original_draw = self.draw
        drawing = [False]

        def hooked_draw():
            if drawing[0]:
                original_draw()
                return
            try:
                if not dialog.winfo_exists():
                    original_draw()
                    return
            except tk.TclError:
                original_draw()
                return
            cur = self.images.selected[:n].copy()
            # Compare with expected focus set selection
            expected = np.zeros(n, dtype=bool)
            for i in focus_set:
                if 0 <= i < n:
                    expected[i] = True
            if not np.array_equal(cur, expected):
                drawing[0] = True
                try:
                    on_click_detected()
                finally:
                    drawing[0] = False
            # Always render with the (possibly corrected) selection
            original_draw()

        self.draw = hooked_draw

        # Highlight focus set in main view and draw initial preview
        sync_main_view()
        original_draw()
        draw_preview()

        def on_ok():
            self.draw = original_draw
            checked_groups = {gid for gid, v in mol_vars.items() if v.get()}
            dialog.destroy()
            focus_indices = np.array(sorted(focus_set))
            if len(focus_indices) == 0:
                # Restore original selection
                self.images.selected[:n] = backup['selected']
                self.draw()
                return
            self._execute_focus(
                focus_indices, checked_groups, "molecules", backup=backup)

        def on_cancel():
            self.draw = original_draw
            # Restore original selection
            self.images.selected[:n] = backup['selected']
            self.draw()
            dialog.destroy()

        # OK / Cancel buttons
        btn_frame = tk.Frame(ctrl_frame)
        btn_frame.pack(pady=(5, 0))
        ttk.Button(btn_frame, text="OK", width=8, command=on_ok).pack(
            side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", width=8,
                   command=on_cancel).pack(side=tk.LEFT, padx=5)

        dialog.protocol("WM_DELETE_WINDOW", on_cancel)
        dialog.geometry("500x700+%d+%d" % (
            self.window.win.winfo_x() + 50, self.window.win.winfo_y() + 50))

    def _execute_focus(self, indices, focus_mol_ids, mode, backup=None):
        """Apply focus with given indices and molecule group set.

        Args:
            indices: atom indices to make visible
            focus_mol_ids: set of mol_id values included in focus
            mode: "atoms" or "molecules"
            backup: pre-saved backup dict (from preview); if None, create one
        """
        n = len(self.atoms)

        # Save pre-focus state (use provided backup if available)
        if backup is not None:
            self._focus_backup = backup
        else:
            self._focus_backup = {
                'visible': self.images.visible[:n].copy(),
                'center': self.center.copy() if hasattr(self, 'center') else None,
                'scale': getattr(self, 'scale', None),
            }
        self._focus_active = True
        self._focus_mol_ids = set(focus_mol_ids)
        self._focus_mode = mode

        # Update button states
        if hasattr(self, '_focus_btn'):
            self._focus_btn.configure(state=tk.DISABLED)
            # ttk buttons don't support bg; use state only
            self._unfocus_btn.configure(state=tk.NORMAL)

        # Set visibility
        visible = np.zeros(n, dtype=bool)
        visible[indices] = True
        self.images.visible[:n] = visible

        # Update side panel molecule checkboxes to reflect focus state
        self._sync_mol_panel_to_focus()

        # Center and zoom on visible atoms
        self._center_and_zoom(indices)

    def _center_and_zoom(self, indices):
        """Center view and auto-zoom on given atom indices."""
        positions = self.atoms.positions[indices]
        self.center = positions.mean(axis=0)

        if len(indices) > 1:
            span = positions.max(axis=0) - positions.min(axis=0)
            max_span = max(span.max(), 1.0)
            canvas_size = min(self.window.canvas.winfo_width(),
                              self.window.canvas.winfo_height())
            if canvas_size > 0:
                self.scale = canvas_size / (max_span * 3.0)

        self.draw()

    def _sync_mol_panel_to_focus(self):
        """Update side panel molecule checkboxes to match focus state."""
        if not hasattr(self, '_mol_vis_vars'):
            return
        for g_id, var in self._mol_vis_vars.items():
            var.set(g_id in self._focus_mol_ids)

    def _toggle_focus_molecule(self, group_id):
        """Add or remove a molecule group from the current focus."""
        if not getattr(self, '_focus_active', False):
            return
        if "mol_id" not in self.atoms.arrays:
            return

        mol_ids = self.atoms.get_array("mol_id")
        n = len(self.atoms)

        if group_id in self._focus_mol_ids:
            # Remove from focus
            self._focus_mol_ids.discard(group_id)
            mask = (mol_ids == group_id)
            self.images.visible[:n][mask] = False
        else:
            # Add to focus
            self._focus_mol_ids.add(group_id)
            mask = (mol_ids == group_id)
            self.images.visible[:n][mask] = True

        # Re-center on all visible atoms
        visible_indices = np.where(self.images.visible[:n])[0]
        if len(visible_indices) > 0:
            self._center_and_zoom(visible_indices)
        else:
            self.draw()

    def unfocus(self, key=None):
        """Exit focus mode: restore all atom visibility."""
        if not getattr(self, '_focus_active', False):
            return

        n = len(self.atoms)
        backup = self._focus_backup

        # Restore visibility
        if backup.get('visible') is not None:
            vlen = min(len(backup['visible']), n)
            self.images.visible[:vlen] = backup['visible'][:vlen]
            if n > vlen:
                self.images.visible[vlen:n] = True

        # Restore selection
        if backup.get('selected') is not None:
            slen = min(len(backup['selected']), n)
            self.images.selected[:slen] = backup['selected'][:slen]

        # Restore view
        if backup.get('center') is not None:
            self.center = backup['center']
        if backup.get('scale') is not None:
            self.scale = backup['scale']

        self._focus_active = False
        self._focus_backup = None
        self._focus_mol_ids = set()
        self._focus_mode = None

        # Reset button states
        if hasattr(self, '_unfocus_btn'):
            self._unfocus_btn.configure(state=tk.DISABLED)
            self._focus_btn.configure(state=tk.NORMAL)

        # Restore side panel molecule checkboxes to normal visibility state
        if hasattr(self, '_mol_vis_vars') and hasattr(self, '_hidden_molecules'):
            for g_id, var in self._mol_vis_vars.items():
                var.set(g_id not in self._hidden_molecules)

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
        if not hasattr(self, '_hidden_elements'):
            self._hidden_elements = set()
        # Toggle: if H currently visible, hide; otherwise show
        if 'H' not in self._hidden_elements:
            self._hidden_elements.add('H')
            self._hydrogen_visible = False
        else:
            self._hidden_elements.discard('H')
            self._hydrogen_visible = True
        # Sync Element Controls checkbox
        if hasattr(self, '_elem_vis_vars') and 'H' in self._elem_vis_vars:
            self._elem_vis_vars['H'].set(self._hydrogen_visible)
        self._apply_visibility()

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
            ttk.Checkbutton(frame, text=sym, variable=v).pack(anchor=tk.W)

        def apply():
            for sym, v in vars_.items():
                elem_mask = np.array([s == sym for s in syms])
                self.images.visible[:n][elem_mask] = v.get()
            self.draw()
            dialog.destroy()

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="Apply", command=apply).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=dialog.destroy).pack(side=tk.LEFT, padx=5)

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
            ttk.Radiobutton(
                color_frame, text=label, variable=color_mode_var,
                value=mode, command=lambda: self._apply_color_mode(color_mode_var.get()),
            ).pack(anchor=tk.W, pady=2)

        # Custom array modes (if atoms have extra arrays)
        if self.atoms is not None:
            skip = {'numbers', 'positions', 'forces', 'momenta',
                    'initial_charges', 'initial_magmoms', 'tags'}
            for key_name in self.atoms.arrays:
                if key_name not in skip and key_name != 'mol_id':
                    ttk.Radiobutton(
                        color_frame, text=f"By '{key_name}'",
                        variable=color_mode_var, value=key_name,
                        command=lambda: self._apply_color_mode(color_mode_var.get()),
                    ).pack(anchor=tk.W, pady=2)

        # --- Tab 2: Element Colors ---
        elem_frame = tk.Frame(notebook, padx=10, pady=10)
        notebook.add(elem_frame, text="Element Colors")

        scroll_container, inner = self._create_scrollable_frame(elem_frame)
        scroll_container.pack(fill=tk.BOTH, expand=True)

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
            btn = ttk.Button(
                row, text="Change...",
                command=lambda z=z, lbl=lbl: self._pick_element_color(z, lbl),
            )
            btn.pack(side=tk.LEFT)

        # Reset button
        ttk.Button(
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
                 font=FONTS["section_heading"]).pack(anchor=tk.W, pady=(10, 5))

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
        ttk.Button(dialog, text="Close", command=dialog.destroy).pack(pady=5)

    def _apply_color_mode(self, mode):
        """Apply a color mode to the current view."""
        if mode in ("mol_group", "mol_id"):
            self._color_by_molecule_group()
        elif mode in ("jmol", "element"):
            self.colormode = 'jmol'
            self.draw()
        elif mode == "tag":
            self.colormode = 'tag'
            # Set colormode_data for tag mode
            if self.atoms is not None:
                tags = self.atoms.get_tags()
                n_colors = max(len(self._GROUP_COLORS), int(tags.max()) + 1)
                colorscale = (self._GROUP_COLORS * ((n_colors // len(self._GROUP_COLORS)) + 1))[:n_colors]
                self.colormode_data = (colorscale, 0, n_colors)
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

    def _restore_original_coloring(self):
        """Restore color mode to element (jmol) and sync side panel."""
        self.colormode = 'jmol'
        if hasattr(self, '_color_mode_var'):
            self._color_mode_var.set("element")

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
        # Use cmax=n_groups so that ASE's linear mapping formula
        # int(k/cmax * N + 0.5) = int(k + 0.5) = k for each integer mol_id k,
        # producing a distinct color for each molecule group.
        self.colormode_data = (colorscale, 0, max(n_groups, 1))
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

    @staticmethod
    def _create_snapshot(atoms):
        """Create a deep copy snapshot of an Atoms object with all metadata."""
        snapshot = atoms.copy()
        snapshot.info = dict(atoms.info)
        for key in atoms.arrays:
            if key not in ('numbers', 'positions'):
                snapshot.set_array(key, atoms.get_array(key).copy())
        return snapshot

    def _refresh_view(self, update_panels=False):
        """Reinitialize images and redraw the view.

        Call after modifying self.atoms to update the display.
        If update_panels is True, also refresh side panel controls.
        """
        self.images.initialize([self.atoms])
        self.set_frame()
        if update_panels:
            self._update_element_visibility_panel()
        if not getattr(self, '_atomlist_collapsed', True):
            self._update_atom_list()
        self.draw()

    def _save_undo_state(self):
        """Save the current atoms state to the undo stack."""
        if self.atoms is None:
            return
        self._undo_stack.append(self._create_snapshot(self.atoms))
        if len(self._undo_stack) > self.MAX_UNDO:
            self._undo_stack.pop(0)
        self._redo_stack.clear()

    def _save_undo_state_with(self, original_positions):
        """Save undo state using original positions (before preview changes).

        Used by interactive dialogs where the atoms have already been
        modified for preview but we want the undo snapshot to hold the
        pre-dialog state.
        """
        if self.atoms is None:
            return
        snapshot = self._create_snapshot(self.atoms)
        snapshot.set_positions(original_positions)
        self._undo_stack.append(snapshot)
        if len(self._undo_stack) > self.MAX_UNDO:
            self._undo_stack.pop(0)
        self._redo_stack.clear()

    def _restore_atoms(self, snapshot):
        """Restore atoms from a snapshot and update the GUI."""
        if len(snapshot) != len(self.atoms):
            # Size changed — must reinitialize with a full copy
            restored = snapshot.copy()
            restored.info = dict(snapshot.info)
            for key in snapshot.arrays:
                if key not in ('numbers', 'positions'):
                    restored.set_array(key, snapshot.get_array(key).copy())
            self.images.initialize([restored])
            self.set_frame()
        else:
            self.atoms.positions[:] = snapshot.positions
            self.atoms.numbers[:] = snapshot.numbers
            self.atoms.set_cell(snapshot.get_cell())
            self.atoms.set_pbc(snapshot.get_pbc())
            self.atoms.info = dict(snapshot.info)
            for key in snapshot.arrays:
                if key not in ('numbers', 'positions'):
                    self.atoms.set_array(key, snapshot.get_array(key).copy())
        self.controller.set_atoms(self.atoms)
        self._refresh_view(update_panels=True)

    def undo(self, key=None):
        """Undo the last operation."""
        if not self._undo_stack:
            if messagebox:
                messagebox.showinfo("Undo", "Nothing to undo.")
            return
        if self.atoms is not None:
            self._redo_stack.append(self._create_snapshot(self.atoms))
        self._restore_atoms(self._undo_stack.pop())

    def redo(self, key=None):
        """Redo the last undone operation."""
        if not self._redo_stack:
            if messagebox:
                messagebox.showinfo("Redo", "Nothing to redo.")
            return
        if self.atoms is not None:
            self._undo_stack.append(self._create_snapshot(self.atoms))
        self._restore_atoms(self._redo_stack.pop())

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
        ttk.Button(btn_frame, text="OK", width=8, command=apply_size).pack(
            side=tk.LEFT, padx=4)
        ttk.Button(btn_frame, text="Cancel", width=8,
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
        self.controller.set_atoms(self.atoms)
        self._refresh_view(update_panels=True)

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
        self.controller.set_atoms(self.atoms)
        self._refresh_view(update_panels=True)

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

    def batch_substitution_dialog(self, key=None):
        """Show batch substitution dialog for multiple element replacements at once."""
        if tk is None or self.atoms is None or len(self.atoms) == 0:
            return

        elements = sorted(set(self.atoms.get_chemical_symbols()))

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Batch Element Substitution")
        dialog.resizable(False, True)
        dialog.transient(self.window.win)

        tk.Label(dialog, text="Define element substitution rules:",
                 font=FONTS["section_heading"], padx=10, pady=8).pack(anchor=tk.W)

        # Scrollable frame for rules
        rules_frame = tk.Frame(dialog)
        rules_frame.pack(fill=tk.BOTH, expand=True, padx=10)

        rule_widgets = []

        def add_rule(old="", new=""):
            row = tk.Frame(rules_frame)
            row.pack(fill=tk.X, pady=2)

            old_var = tk.StringVar(value=old)
            new_var = tk.StringVar(value=new)

            old_combo = ttk.Combobox(row, textvariable=old_var, values=elements,
                                     width=8, state="readonly")
            old_combo.pack(side=tk.LEFT, padx=2)
            tk.Label(row, text="→").pack(side=tk.LEFT, padx=4)
            new_entry = tk.Entry(row, textvariable=new_var, width=8)
            new_entry.pack(side=tk.LEFT, padx=2)

            def remove():
                row.destroy()
                rule_widgets.remove((old_var, new_var, row))

            ttk.Button(row, text="\u00d7", command=remove, width=2).pack(side=tk.LEFT, padx=4)
            rule_widgets.append((old_var, new_var, row))

        # Add initial empty rule
        add_rule()

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="+ Add Rule", command=add_rule).pack(side=tk.LEFT, padx=5)

        ttk.Separator(dialog, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=10, pady=5)

        def apply_batch():
            rules = []
            for old_var, new_var, _ in rule_widgets:
                old_sym = old_var.get().strip()
                new_sym = new_var.get().strip()
                if old_sym and new_sym:
                    rules.append((old_sym, new_sym))

            if not rules:
                if messagebox:
                    messagebox.showwarning("Batch Substitution",
                                           "No valid rules defined.", parent=dialog)
                return

            self._save_undo_state()
            self.controller.set_atoms(self.atoms)
            total = 0
            results = []
            for old_sym, new_sym in rules:
                count = self.controller.global_atom_substitution(
                    old_symbol=old_sym, new_symbol=new_sym)
                total += count
                results.append(f"{old_sym} → {new_sym}: {count} atoms")

            self._refresh_view(update_panels=True)

            if messagebox:
                msg = "\n".join(results) + f"\n\nTotal: {total} atoms replaced"
                messagebox.showinfo("Batch Substitution", msg, parent=dialog)
            dialog.destroy()

        action_frame = tk.Frame(dialog)
        action_frame.pack(pady=8)
        ttk.Button(action_frame, text="Apply All", command=apply_batch, width=10
                   ).pack(side=tk.LEFT, padx=5)
        ttk.Button(action_frame, text="Cancel", command=dialog.destroy, width=10
                   ).pack(side=tk.LEFT, padx=5)

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
    # Cell parameter editor
    # =========================================================================

    def cell_parameter_dialog(self, key=None):
        """Show cell parameter display/editor dialog."""
        if tk is None or self.atoms is None:
            return

        cell = self.atoms.get_cell()
        lengths = cell.lengths()  # a, b, c
        angles = cell.angles()    # alpha, beta, gamma

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Cell Parameters")
        dialog.resizable(False, False)
        dialog.transient(self.window.win)
        dialog.grab_set()
        dialog.lift()
        dialog.focus_force()

        # --- Lattice lengths ---
        tk.Label(dialog, text="Lattice lengths (\u00c5):", font=FONTS["section_heading"]
                 ).grid(row=0, column=0, columnspan=4, sticky=tk.W, padx=10, pady=(10, 2))

        len_vars = []
        for col, (label, val) in enumerate(zip(("a:", "b:", "c:"), lengths)):
            tk.Label(dialog, text=label).grid(row=1, column=col * 2, sticky=tk.E, padx=5)
            v = tk.DoubleVar(value=round(val, 6))
            len_vars.append(v)
            e = tk.Entry(dialog, textvariable=v, width=12)
            e.grid(row=1, column=col * 2 + 1, padx=2, pady=2)

        # --- Lattice angles ---
        tk.Label(dialog, text="Lattice angles (\u00b0):", font=FONTS["section_heading"]
                 ).grid(row=2, column=0, columnspan=4, sticky=tk.W, padx=10, pady=(10, 2))

        ang_vars = []
        for col, (label, val) in enumerate(zip(("α:", "β:", "γ:"), angles)):
            tk.Label(dialog, text=label).grid(row=3, column=col * 2, sticky=tk.E, padx=5)
            v = tk.DoubleVar(value=round(val, 4))
            ang_vars.append(v)
            e = tk.Entry(dialog, textvariable=v, width=12)
            e.grid(row=3, column=col * 2 + 1, padx=2, pady=2)

        # --- Lattice vectors (read-only) ---
        tk.Label(dialog, text="Lattice vectors:", font=FONTS["section_heading"]
                 ).grid(row=4, column=0, columnspan=4, sticky=tk.W, padx=10, pady=(10, 2))

        vec_text = tk.Text(dialog, height=3, width=50, state=tk.NORMAL, font=FONTS["fixed_small"])
        vec_text.grid(row=5, column=0, columnspan=6, padx=10, pady=2)
        for i in range(3):
            vec_text.insert(tk.END,
                            f"  v{i+1} = ({cell[i][0]:10.6f}, {cell[i][1]:10.6f}, {cell[i][2]:10.6f})\n")
        vec_text.config(state=tk.DISABLED)

        # --- Volume ---
        vol = cell.volume
        tk.Label(dialog, text=f"Volume: {vol:.4f} Å³"
                 ).grid(row=6, column=0, columnspan=4, sticky=tk.W, padx=10, pady=4)

        # --- Buttons ---
        btn_frame = tk.Frame(dialog)
        btn_frame.grid(row=7, column=0, columnspan=6, pady=10)

        def apply_changes():
            try:
                new_lengths = [v.get() for v in len_vars]
                new_angles = [v.get() for v in ang_vars]
            except tk.TclError:
                if messagebox:
                    messagebox.showerror("Error", "Invalid numeric input.", parent=dialog)
                return

            # Validate
            for val in new_lengths:
                if val <= 0:
                    if messagebox:
                        messagebox.showerror("Error", "Lengths must be positive.", parent=dialog)
                    return
            for val in new_angles:
                if val <= 0 or val >= 180:
                    if messagebox:
                        messagebox.showerror("Error", "Angles must be between 0° and 180°.",
                                             parent=dialog)
                    return

            self._save_undo_state()
            from ase.cell import Cell
            new_cell = Cell.fromcellpar(new_lengths + new_angles)
            self.atoms.set_cell(new_cell, scale_atoms=True)
            self.images.initialize([self.atoms])
            self.set_frame()
            self.draw()
            self.controller.set_atoms(self.atoms)

            # Update vector display
            c = self.atoms.get_cell()
            vec_text.config(state=tk.NORMAL)
            vec_text.delete("1.0", tk.END)
            for i in range(3):
                vec_text.insert(tk.END,
                                f"  v{i+1} = ({c[i][0]:10.6f}, {c[i][1]:10.6f}, {c[i][2]:10.6f})\n")
            vec_text.config(state=tk.DISABLED)

            # Update volume
            dialog.children  # not needed, update label below

            if messagebox:
                messagebox.showinfo("Cell Parameters",
                                    f"Cell updated.\nVolume: {c.volume:.4f} Å³",
                                    parent=dialog)

        ttk.Button(btn_frame, text="Apply", command=apply_changes, width=10
                   ).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Close", command=dialog.destroy, width=10
                   ).pack(side=tk.LEFT, padx=5)

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

        # Export CSV button
        def export_csv():
            from tkinter import filedialog as fd
            path = fd.asksaveasfilename(
                parent=dialog, title="Export CSV",
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            )
            if not path:
                return
            with open(path, "w") as f:
                f.write("Site,Molecule,IoU,Alpha,Beta,Gamma,Scale\n")
                for i in range(len(iou)):
                    mol_name = site_names[i] if i < len(site_names) else "?"
                    ang = angles[i] if angles is not None else [0, 0, 0]
                    sc = scales[i] if scales is not None else 1.0
                    f.write(f"{i},{mol_name},{iou[i]:.4f},"
                            f"{ang[0]:.1f},{ang[1]:.1f},{ang[2]:.1f},{sc:.4f}\n")
            if messagebox:
                messagebox.showinfo("Export", f"Saved: {path}", parent=dialog)

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=5)
        ttk.Button(btn_frame, text="Export CSV", command=export_csv).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="View PNGs", command=self._view_voronoi_pngs).pack(
            side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Close", command=dialog.destroy).pack(side=tk.LEFT, padx=5)

    def _view_voronoi_pngs(self):
        """Open a viewer for Voronoi overlay PNGs."""
        if tk is None:
            return

        # Look for PNG files in common locations
        png_dir = None
        for candidate in [
            os.path.join(os.getcwd(), "source", "voronoi_pngs"),
            os.path.join(os.getcwd(), "voronoi_pngs"),
        ]:
            if os.path.isdir(candidate):
                png_dir = candidate
                break

        if png_dir is None:
            from tkinter import filedialog
            png_dir = filedialog.askdirectory(
                parent=self.window.win,
                title="Select Voronoi PNGs directory",
            )
            if not png_dir:
                return

        pngs = sorted([f for f in os.listdir(png_dir) if f.lower().endswith('.png')])
        if not pngs:
            if messagebox:
                messagebox.showinfo("Voronoi PNGs", "No PNG files found in the directory.")
            return

        dialog = tk.Toplevel(self.window.win)
        dialog.title("Voronoi Overlay PNGs")
        dialog.resizable(True, True)
        dialog.transient(self.window.win)

        # Navigation
        idx_var = tk.IntVar(value=0)
        label_var = tk.StringVar(value=f"1 / {len(pngs)}: {pngs[0]}")

        tk.Label(dialog, textvariable=label_var, font=FONTS["section_heading"]
                 ).pack(pady=4)

        # Image display
        try:
            from PIL import Image, ImageTk
            has_pil = True
        except ImportError:
            has_pil = False

        if has_pil:
            canvas = tk.Canvas(dialog, width=600, height=500)
            canvas.pack(fill=tk.BOTH, expand=True)
            photo_ref = [None]

            def show_image(idx):
                img = Image.open(os.path.join(png_dir, pngs[idx]))
                # Fit to canvas
                cw, ch = 600, 500
                img.thumbnail((cw, ch), Image.LANCZOS)
                photo = ImageTk.PhotoImage(img)
                canvas.delete("all")
                canvas.create_image(cw // 2, ch // 2, image=photo, anchor=tk.CENTER)
                photo_ref[0] = photo  # keep reference
                label_var.set(f"{idx + 1} / {len(pngs)}: {pngs[idx]}")
        else:
            tk.Label(dialog, text="(Install Pillow for image preview)",
                     fg="gray").pack()

            def show_image(idx):
                label_var.set(f"{idx + 1} / {len(pngs)}: {pngs[idx]}")

        def prev_img():
            i = idx_var.get()
            if i > 0:
                idx_var.set(i - 1)
                show_image(i - 1)

        def next_img():
            i = idx_var.get()
            if i < len(pngs) - 1:
                idx_var.set(i + 1)
                show_image(i + 1)

        nav_frame = tk.Frame(dialog)
        nav_frame.pack(pady=5)
        ttk.Button(nav_frame, text="< Prev", command=prev_img, width=8).pack(side=tk.LEFT, padx=5)
        ttk.Button(nav_frame, text="Next >", command=next_img, width=8).pack(side=tk.LEFT, padx=5)
        ttk.Button(nav_frame, text="Close", command=dialog.destroy, width=8).pack(
            side=tk.LEFT, padx=5)

        show_image(0)

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
                from molecular_modeling.build_system_from_toml import build_system_from_toml
                new_atoms = build_system_from_toml(verbose=True)
            finally:
                os.chdir(original_dir)

            self._save_undo_state()
            self.images.initialize([new_atoms])
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

    def _run_pipeline_in_thread(self, task_name, task_func, on_complete=None):
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
                    if on_complete:
                        on_complete()

        thread = threading.Thread(target=run, daemon=True)
        thread.start()
        self.window.win.after(100, check_thread)

    def run_optimize(self, key=None):
        """Run Voronoi optimization with automatic Before/After rendering."""
        from source.build_molecule_groups import build_molecule_groups
        from source.optimize_voronoi_overlap import optimize_voronoi_overlap
        from source.render_voronoi_overlays import render_voronoi_overlays

        self._save_undo_state()

        # Output directory
        out_dir = os.path.join(os.getcwd(), "voronoi_pngs")
        before_dir = os.path.join(out_dir, "before")
        after_dir = os.path.join(out_dir, "after")

        def task():
            groups = build_molecule_groups(self.atoms)
            rngs = self.atoms.info.get("rng", {})
            rng_opt = int(rngs.get("optimize", rngs.get("global", 1)))
            rng_ren = int(rngs.get("render", rngs.get("global", 1)))
            N = len(groups)

            # --- Before: render with initial angles (zeros) ---
            before_angles = np.zeros((N, 3), dtype=float)
            if "optimized_angles_deg" in self.atoms.info:
                before_angles = np.asarray(
                    self.atoms.info["optimized_angles_deg"], float)
            render_voronoi_overlays(
                self.atoms, groups=groups, angles_deg=before_angles,
                out_dir=before_dir, rng=rng_ren,
            )

            # --- Optimize ---
            optimize_voronoi_overlap(
                self.atoms, groups=groups, rng=rng_opt,
            )

            # --- After: render with optimized angles ---
            render_voronoi_overlays(
                self.atoms, groups=groups, out_dir=after_dir, rng=rng_ren,
            )

            # --- Generate Before/After comparison images ---
            self._render_comparison_pngs(before_dir, after_dir, out_dir)

            # Format results
            iou = self.atoms.info.get("optimized_iou", [])
            angles = self.atoms.info.get("optimized_angles_deg", [])
            lines = ["Voronoi optimization complete.\n"]
            for i, (val, ang) in enumerate(zip(iou, angles)):
                lines.append(
                    f"  Site {i}: IoU={val:.4f}, "
                    f"angles=({ang[0]:.1f}, {ang[1]:.1f}, {ang[2]:.1f})\u00b0"
                )
            lines.append(f"\nPNGs saved to: {out_dir}")
            return "\n".join(lines)

        def _on_complete():
            """Show Voronoi PNG viewer after optimization."""
            self._view_voronoi_pngs(png_dir=out_dir)

        self._run_pipeline_in_thread("Optimize Voronoi", task,
                                     on_complete=_on_complete)

    @staticmethod
    def _render_comparison_pngs(before_dir, after_dir, out_dir):
        """Generate side-by-side Before/After comparison PNGs."""
        try:
            from PIL import Image
        except ImportError:
            return  # Pillow not available, skip comparison

        before_pngs = sorted(
            f for f in os.listdir(before_dir) if f.lower().endswith('.png'))
        after_pngs = sorted(
            f for f in os.listdir(after_dir) if f.lower().endswith('.png'))

        for bname, aname in zip(before_pngs, after_pngs):
            bimg = Image.open(os.path.join(before_dir, bname))
            aimg = Image.open(os.path.join(after_dir, aname))

            # Resize to same height
            h = max(bimg.height, aimg.height)
            if bimg.height != h:
                bimg = bimg.resize(
                    (int(bimg.width * h / bimg.height), h), Image.LANCZOS)
            if aimg.height != h:
                aimg = aimg.resize(
                    (int(aimg.width * h / aimg.height), h), Image.LANCZOS)

            # Create side-by-side image with labels
            gap = 20
            total_w = bimg.width + gap + aimg.width
            label_h = 30
            comp = Image.new("RGB", (total_w, h + label_h), "white")

            # Draw labels
            from PIL import ImageDraw
            draw = ImageDraw.Draw(comp)
            draw.text((bimg.width // 2 - 30, 4), "Before", fill="black")
            draw.text((bimg.width + gap + aimg.width // 2 - 20, 4),
                      "After", fill="black")

            comp.paste(bimg, (0, label_h))
            comp.paste(aimg, (bimg.width + gap, label_h))

            site_id = bname.replace("overlay_site_", "").replace(".png", "")
            comp_path = os.path.join(out_dir, f"comparison_site_{site_id}.png")
            comp.save(comp_path)
            bimg.close()
            aimg.close()
            comp.close()

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

        png_dir = os.path.join(out_dir, "voronoi_pngs")
        before_dir = os.path.join(png_dir, "before")
        after_dir = os.path.join(png_dir, "after")

        def task():
            import os as _os
            groups = build_molecule_groups(self.atoms)
            rngs = self.atoms.info.get("rng", {})
            rng_opt = int(rngs.get("optimize", rngs.get("global", 1)))
            rng_ren = int(rngs.get("render", rngs.get("global", 1)))
            N = len(groups)

            # Before render
            before_angles = np.zeros((N, 3), dtype=float)
            if "optimized_angles_deg" in self.atoms.info:
                before_angles = np.asarray(
                    self.atoms.info["optimized_angles_deg"], float)
            render_voronoi_overlays(
                self.atoms, groups=groups, angles_deg=before_angles,
                out_dir=before_dir, rng=rng_ren,
            )

            # Optimize (modifies atoms in place)
            optimize_voronoi_overlap(
                self.atoms, groups=groups, rng=rng_opt,
            )

            # Check tiling
            check_system_voronoi_volume_tiling(
                self.atoms, groups, rel_tol=1e-3, abs_tol=1e-6, verbose=True,
            )

            # After render
            render_voronoi_overlays(
                self.atoms, groups=groups, out_dir=after_dir, rng=rng_ren,
            )

            # Comparison PNGs
            self._render_comparison_pngs(before_dir, after_dir, png_dir)

            # Export
            struct_dir = _os.path.join(out_dir, "structure_out")
            _os.makedirs(struct_dir, exist_ok=True)
            write_structure(self.atoms, outdir=struct_dir, verbose=True)

            # Format results
            iou = self.atoms.info.get("optimized_iou", [])
            lines = ["Full pipeline complete.\n"]
            for i, val in enumerate(iou):
                lines.append(f"  Site {i}: IoU={val:.4f}")
            lines.append(f"\nOutputs saved to:\n{out_dir}")
            return "\n".join(lines)

        def _on_complete():
            self._view_voronoi_pngs(png_dir=png_dir)

        self._run_pipeline_in_thread("Full Pipeline", task,
                                     on_complete=_on_complete)

    # =========================================================================
    # Help menu actions
    # =========================================================================

    def show_keyboard_shortcuts(self, key=None):
        """Show a dialog with keyboard shortcuts."""
        if tk is None:
            return
        dialog = tk.Toplevel(self.window.win)
        dialog.title("Keyboard Shortcuts")
        dialog.transient(self.window.win)

        shortcuts = [
            ("Ctrl+Z", "Undo"),
            ("Ctrl+Y", "Redo"),
            ("Ctrl+A", "Select all atoms"),
            ("Escape", "Deselect all"),
            ("Delete / Backspace", "Delete selected atoms"),
            ("", ""),
            ("Ctrl+L", "Load structure file"),
            ("Ctrl+J", "Add molecule"),
            ("Ctrl+D", "Export to all formats"),
            ("Ctrl+P", "Run Voronoi optimization"),
            ("Ctrl+H", "Show optimization results"),
            ("Ctrl+K", "Save project"),
            ("Ctrl+T", "Load project"),
            ("Ctrl+B", "Show/hide bonds"),
            ("Ctrl+F", "Focus on selection"),
            ("Ctrl+Shift+F", "Exit focus mode"),
            ("Ctrl+Shift+S", "Save screenshot"),
            ("", ""),
            ("Left click", "Select atom"),
            ("Shift+click", "Add to selection"),
            ("Ctrl+click", "Toggle selection"),
            ("Right-click", "Context menu"),
            ("Left drag", "Rotate view"),
            ("Right drag", "Translate view"),
            ("Mouse wheel", "Zoom"),
        ]

        frame = tk.Frame(dialog, padx=16, pady=12)
        frame.pack(fill=tk.BOTH, expand=True)

        tk.Label(frame, text="Keyboard Shortcuts",
                 font=("TkDefaultFont", 12, "bold")
                 ).pack(anchor=tk.W, pady=(0, 8))

        for key_str, desc in shortcuts:
            if not key_str:
                ttk.Separator(frame, orient=tk.HORIZONTAL).pack(
                    fill=tk.X, pady=4)
                continue
            row = tk.Frame(frame)
            row.pack(fill=tk.X, pady=1)
            tk.Label(row, text=key_str, font=FONTS["fixed"],
                     width=16, anchor=tk.W).pack(side=tk.LEFT)
            tk.Label(row, text=desc, anchor=tk.W).pack(side=tk.LEFT, padx=(8, 0))

        tk.Label(frame, text="\nStandard ASE-GUI shortcuts (Ctrl+O, Ctrl+S, etc.) "
                 "are also available.",
                 fg="gray", wraplength=350, justify=tk.LEFT).pack(
                     anchor=tk.W, pady=(8, 0))

        ttk.Button(frame, text="Close", width=10,
                   command=dialog.destroy).pack(pady=(12, 0))
        dialog.geometry("+%d+%d" % (
            self.window.win.winfo_x() + 80, self.window.win.winfo_y() + 40))

    def open_user_guide(self, key=None):
        """Open the HTML user guide in the default browser."""
        import webbrowser

        # Locate docs relative to this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        docs_base = os.path.join(script_dir, "docs", "_build", "html")

        # Try Japanese first, then English
        candidates = [
            os.path.join(docs_base, "ja", "gui_guide.html"),
            os.path.join(docs_base, "en", "gui_guide.html"),
            os.path.join(docs_base, "ja", "index.html"),
            os.path.join(docs_base, "en", "index.html"),
        ]

        for path in candidates:
            if os.path.isfile(path):
                webbrowser.open("file://" + path)
                return

        # HTML not built — show instructions
        if messagebox:
            messagebox.showinfo(
                "User Guide",
                "HTML documentation has not been built yet.\n\n"
                "To build it, run:\n"
                "  cd docs\n"
                "  make html\n\n"
                f"Docs source: {os.path.join(script_dir, 'docs')}"
            )

    def show_about(self, key=None):
        """Show the About dialog."""
        if messagebox is None:
            return
        messagebox.showinfo(
            "About Molecular Modeling Toolkit",
            "Molecular Modeling Toolkit\n"
            "Version 0.1.0\n\n"
            "ASE-GUI extension for molecular modeling,\n"
            "Voronoi cell optimization, and structure export.\n\n"
            "Built on ASE (Atomic Simulation Environment)"
        )

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
                M(_('Batch Element Substitution...'),
                  self.batch_substitution_dialog),
                M(_('Build Supercell...'),
                  self.supercell_dialog),
                M(_('Cell Parameters...'),
                  self.cell_parameter_dialog),
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
                M(_('Select All'), self._select_all, 'Ctrl+A'),
                M(_('Deselect All'), self._deselect_all),
                M('---'),
                M(_('Select by Element...'), self.select_by_element),
                M(_('Select by Molecule...'), self.select_by_molecule),
                M(_('Select by Sphere...'), self.select_by_sphere),
                M('---'),
                M(_('Delete Selected Atoms'), self.delete_selected_atoms),
                M(_('Delete Molecule...'), self.delete_molecule),
                M('---'),
                M(_('Move Molecule/Atoms...'), self.move_molecule_dialog),
                M(_('Rotate Molecule/Atoms...'), self.rotate_molecule_dialog),
                M('---'),
                M(_('Focus on Selection'), self.focus_selection, 'Ctrl+F'),
                M(_('Exit Focus Mode'), self.unfocus, 'Ctrl+Shift+F'),
                M('---'),
                M(_('Toggle Drag Mode'), self._toggle_drag_mode_menu),
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
            M(_('Show Bonds'), self.toggle_show_bonds, 'Ctrl+B',
              value=self.config.get('show_bonds', False)),
            M(_('Toggle Hydrogen'), self.toggle_hydrogen),
            M(_('Element Visibility...'), self.toggle_element_visibility),
            M(_('Toggle Side Panel'), self.toggle_side_panel),
            M('---'),
            M(_('Save Screenshot...'), self.save_screenshot, 'Ctrl+Shift+S'),
            M('---'),
        ]
        for label, w, h in self.WINDOW_SIZES:
            view_items.append(M(_(label), self._resize_window_preset(w, h)))
        view_items.append(M(_('Custom Size...'), self._resize_window_custom))
        menu.append((_('View Settings'), view_items))

        # Help menu
        menu.append((
            _('Help'),
            [
                M(_('Keyboard Shortcuts...'), self.show_keyboard_shortcuts),
                M(_('Open User Guide...'), self.open_user_guide),
                M('---'),
                M(_('About...'), self.show_about),
            ]
        ))

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
