# source/gui/add_molecule_dialog.py
"""
Dialog for adding molecules to the system.

Allows users to select a molecule from the registry and specify its position.
Includes a 3D preview of the selected molecule.
"""

import tkinter as tk
from tkinter import ttk
from typing import Dict, Optional, Tuple

import numpy as np

from .base_dialog import BaseDialog, CoordinateEntry, LabeledCombobox
from .molecule_registry import get_registry

# Element colors (CPK coloring scheme)
ELEMENT_COLORS: Dict[str, str] = {
    "H": "#FFFFFF",   # White
    "C": "#909090",   # Gray
    "N": "#3050F8",   # Blue
    "O": "#FF0D0D",   # Red
    "F": "#90E050",   # Green
    "S": "#FFFF30",   # Yellow
    "P": "#FF8000",   # Orange
    "Cl": "#1FF01F",  # Green
    "Br": "#A62929",  # Brown
    "I": "#940094",   # Purple
    "Si": "#F0C8A0",  # Beige
    "Fe": "#E06633",  # Orange-brown
    "default": "#FF1493",  # Pink for unknown
}

# Element radii for display (in Angstroms, scaled)
ELEMENT_RADII: Dict[str, float] = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "S": 1.05,
    "P": 1.07,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
    "Si": 1.11,
    "Fe": 1.26,
    "default": 0.80,
}


class AddMoleculeDialog(BaseDialog):
    """
    Dialog for adding a molecule to the system.

    Allows selection of molecule type and position (Cartesian or fractional).
    Shows a 3D preview of the selected molecule.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    cell : numpy.ndarray, optional
        Unit cell matrix (3x3) for fractional coordinate conversion.
    """

    def __init__(
        self,
        parent: Optional[tk.Widget],
        cell: Optional[np.ndarray] = None,
    ):
        self.cell = cell
        self.registry = get_registry()
        self.canvas = None
        self.preview_enabled = True

        # Try to import matplotlib
        try:
            import matplotlib
            matplotlib.use('TkAgg')
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            from matplotlib.figure import Figure
            self._matplotlib_available = True
        except ImportError:
            self._matplotlib_available = False

        super().__init__(parent, title="Add Molecule")
        # Make dialog larger for preview
        self.geometry("500x700")

    def build_body(self, parent: ttk.Frame) -> None:
        """Build the dialog body."""
        # Molecule selection
        molecule_names = self.registry.get_molecule_names()
        if not molecule_names:
            molecule_names = ["(no molecules found)"]

        self.molecule_combo = LabeledCombobox(
            parent, "Molecule:", molecule_names
        )
        self.molecule_combo.pack(fill=tk.X, pady=2)

        # Molecule info display
        info_frame = ttk.LabelFrame(parent, text="Molecule Info", padding=5)
        info_frame.pack(fill=tk.X, pady=5)

        self.info_label = ttk.Label(info_frame, text="")
        self.info_label.pack(anchor=tk.W)

        # Preview frame
        if self._matplotlib_available:
            preview_frame = ttk.LabelFrame(parent, text="Preview", padding=5)
            preview_frame.pack(fill=tk.BOTH, expand=True, pady=5)
            self._setup_preview(preview_frame)

        # Bind selection change
        self.molecule_combo.combobox.bind(
            "<<ComboboxSelected>>", self._on_molecule_selected
        )
        self._update_molecule_info()
        self._update_preview()

        # Add mode selection (single molecule or full unit cell)
        mode_frame = ttk.LabelFrame(parent, text="Add Mode", padding=5)
        mode_frame.pack(fill=tk.X, pady=5)

        self.add_mode_var = tk.StringVar(value="unitcell")  # Default to unit cell

        self.unitcell_radio = ttk.Radiobutton(
            mode_frame,
            text="All molecules in unit cell",
            variable=self.add_mode_var,
            value="unitcell",
            command=self._on_mode_changed,
        )
        self.unitcell_radio.pack(anchor=tk.W)

        self.single_radio = ttk.Radiobutton(
            mode_frame,
            text="Single molecule at position",
            variable=self.add_mode_var,
            value="single",
            command=self._on_mode_changed,
        )
        self.single_radio.pack(anchor=tk.W)

        # Position/Offset entry frame
        self.position_frame = ttk.LabelFrame(parent, text="Offset", padding=5)
        self.position_frame.pack(fill=tk.X, pady=5)

        self.position_entry = CoordinateEntry(
            self.position_frame, "XYZ:", defaults=(0.0, 0.0, 0.0)
        )
        self.position_entry.pack(fill=tk.X)

        # Fractional coordinates checkbox (default True)
        self.use_fractional_var = tk.BooleanVar(value=True)
        self.fractional_check = ttk.Checkbutton(
            self.position_frame,
            text="Use fractional coordinates",
            variable=self.use_fractional_var,
        )
        self.fractional_check.pack(anchor=tk.W, pady=2)

        # Disable fractional if no cell is provided
        if self.cell is None:
            self.use_fractional_var.set(False)
            self.fractional_check.config(state=tk.DISABLED)

        # Apply initial mode state
        self._on_mode_changed()

    def _on_mode_changed(self) -> None:
        """Handle add mode change."""
        if self.add_mode_var.get() == "single":
            self.position_frame.config(text="Position")
        else:
            self.position_frame.config(text="Offset")

    def _setup_preview(self, parent: ttk.Frame) -> None:
        """Set up the matplotlib preview canvas."""
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        self.fig = Figure(figsize=(4, 3), dpi=100)
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _on_molecule_selected(self, event=None) -> None:
        """Handle molecule selection change."""
        self._update_molecule_info()
        self._update_preview()

    def _update_molecule_info(self) -> None:
        """Update the molecule info display."""
        name = self.molecule_combo.get()
        if name and name != "(no molecules found)":
            num_atoms = self.registry.get_num_atoms(name)
            elements = self.registry.get_unique_elements(name)
            elements_str = ", ".join(elements)
            self.info_label.config(
                text=f"Atoms: {num_atoms}  |  Elements: {elements_str}"
            )
        else:
            self.info_label.config(text="No molecule selected")

    def _update_preview(self) -> None:
        """Update the 3D molecule preview with unit cell and all molecules."""
        if not self._matplotlib_available or self.canvas is None:
            return

        name = self.molecule_combo.get()
        if not name or name == "(no molecules found)":
            self.ax.clear()
            self.ax.set_title("No molecule selected")
            self.canvas.draw()
            return

        # Get molecule config
        config = self.registry.get_molecule(name)
        if config is None:
            self.ax.clear()
            self.ax.set_title("No data available")
            self.canvas.draw()
            return

        # Get molecule data
        symbols = self.registry.get_symbols(name)
        base_positions = self.registry.get_positions(name)

        if len(symbols) == 0 or len(base_positions) == 0:
            self.ax.clear()
            self.ax.set_title("No data available")
            self.canvas.draw()
            return

        # Get lattice vectors
        lattice = np.array(config.get("lattice_vectors", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]), dtype=float)

        # Get center of mass positions for multiple molecules
        center_of_mass = np.array(config.get("center_of_mass", [[0, 0, 0]]), dtype=float)
        num_molecules = int(config.get("num_molecules", 1))

        # Handle fractional coordinates
        use_fractional = config.get("use_fractional", False)
        if use_fractional:
            base_positions = base_positions @ lattice
            center_of_mass = center_of_mass @ lattice

        # Clear and redraw
        self.ax.clear()

        # Center base molecule at origin for placement
        base_com = base_positions.mean(axis=0)
        base_positions_centered = base_positions - base_com

        all_positions = []

        # Plot molecules at each center of mass
        for mol_idx in range(min(num_molecules, len(center_of_mass))):
            com = center_of_mass[mol_idx]
            mol_positions = base_positions_centered + com
            all_positions.extend(mol_positions.tolist())

            # Plot atoms
            for sym, pos in zip(symbols, mol_positions):
                color = ELEMENT_COLORS.get(sym, ELEMENT_COLORS["default"])
                radius = ELEMENT_RADII.get(sym, ELEMENT_RADII["default"])
                size = radius * 80

                self.ax.scatter(
                    pos[0], pos[1], pos[2],
                    c=color,
                    s=size,
                    edgecolors='black',
                    linewidths=0.5,
                    alpha=0.9,
                )

            # Draw bonds within this molecule
            bond_threshold = 1.8
            for i in range(len(mol_positions)):
                for j in range(i + 1, len(mol_positions)):
                    dist = np.linalg.norm(mol_positions[i] - mol_positions[j])
                    if dist < bond_threshold:
                        self.ax.plot(
                            [mol_positions[i][0], mol_positions[j][0]],
                            [mol_positions[i][1], mol_positions[j][1]],
                            [mol_positions[i][2], mol_positions[j][2]],
                            'k-', linewidth=1, alpha=0.5
                        )

        # Draw unit cell box
        self._draw_unit_cell(lattice)

        # Set axis limits based on unit cell
        all_positions = np.array(all_positions) if all_positions else np.zeros((1, 3))

        # Include unit cell corners in range calculation
        corners = np.array([
            [0, 0, 0],
            lattice[0],
            lattice[1],
            lattice[2],
            lattice[0] + lattice[1],
            lattice[0] + lattice[2],
            lattice[1] + lattice[2],
            lattice[0] + lattice[1] + lattice[2],
        ])

        all_points = np.vstack([all_positions, corners])
        margin = 1.0

        self.ax.set_xlim(all_points[:, 0].min() - margin, all_points[:, 0].max() + margin)
        self.ax.set_ylim(all_points[:, 1].min() - margin, all_points[:, 1].max() + margin)
        self.ax.set_zlim(all_points[:, 2].min() - margin, all_points[:, 2].max() + margin)

        self.ax.set_xlabel('X (Å)')
        self.ax.set_ylabel('Y (Å)')
        self.ax.set_zlabel('Z (Å)')
        self.ax.set_title(f'{name} ({num_molecules} mol/cell)')

        try:
            self.ax.set_box_aspect([1, 1, 1])
        except AttributeError:
            pass

        self.canvas.draw()

    def _draw_unit_cell(self, lattice: np.ndarray) -> None:
        """Draw the unit cell box."""
        a, b, c = lattice[0], lattice[1], lattice[2]
        origin = np.array([0, 0, 0])

        # Define the 12 edges of the parallelepiped
        edges = [
            # Bottom face
            (origin, a),
            (origin, b),
            (a, a + b),
            (b, a + b),
            # Top face
            (c, c + a),
            (c, c + b),
            (c + a, c + a + b),
            (c + b, c + a + b),
            # Vertical edges
            (origin, c),
            (a, a + c),
            (b, b + c),
            (a + b, a + b + c),
        ]

        for start, end in edges:
            self.ax.plot(
                [start[0], end[0]],
                [start[1], end[1]],
                [start[2], end[2]],
                'b-', linewidth=1.5, alpha=0.6
            )

    def validate(self) -> bool:
        """Validate dialog input."""
        name = self.molecule_combo.get()
        if not name or name == "(no molecules found)":
            self.show_error("Error", "Please select a molecule.")
            return False

        # Validate position/offset for both modes
        try:
            self.position_entry.get()
        except ValueError:
            error_msg = "Invalid position coordinates." if self.add_mode_var.get() == "single" else "Invalid offset coordinates."
            self.show_error("Error", error_msg)
            return False

        return True

    def apply(self) -> None:
        """Apply dialog result."""
        name = self.molecule_combo.get()
        add_mode = self.add_mode_var.get()

        if add_mode == "single":
            position = self.position_entry.get()
            use_fractional = self.use_fractional_var.get()

            # Convert to Cartesian if using fractional
            if use_fractional and self.cell is not None:
                position = tuple(np.array(position) @ self.cell)

            self.result = {
                "molecule_name": name,
                "position": position,
                "use_fractional": use_fractional,
                "add_mode": "single",
            }
        else:
            # Unit cell mode - return all molecule positions from config with offset
            config = self.registry.get_molecule(name)
            if config:
                lattice = np.array(config.get("lattice_vectors", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]), dtype=float)
                center_of_mass = np.array(config.get("center_of_mass", [[0, 0, 0]]), dtype=float)
                num_molecules = int(config.get("num_molecules", 1))
                use_frac = config.get("use_fractional", False)

                if use_frac:
                    center_of_mass = center_of_mass @ lattice

                # Get offset and convert to Cartesian if needed
                offset = np.array(self.position_entry.get())
                use_fractional_offset = self.use_fractional_var.get()
                if use_fractional_offset and self.cell is not None:
                    offset = offset @ self.cell

                # Apply offset to all positions
                positions = [tuple(center_of_mass[i] + offset) for i in range(min(num_molecules, len(center_of_mass)))]

                self.result = {
                    "molecule_name": name,
                    "positions": positions,  # Multiple positions with offset
                    "lattice": lattice.tolist(),
                    "offset": offset.tolist(),
                    "add_mode": "unitcell",
                }
            else:
                # Fallback if config not found
                self.result = {
                    "molecule_name": name,
                    "position": (0.0, 0.0, 0.0),
                    "add_mode": "single",
                }


class AddMoleculeResult:
    """
    Result from AddMoleculeDialog.

    Attributes
    ----------
    molecule_name : str
        Name of the selected molecule.
    position : tuple[float, float, float]
        Position in Cartesian coordinates.
    use_fractional : bool
        Whether fractional coordinates were used.
    """

    def __init__(self, molecule_name: str, position: Tuple[float, float, float], use_fractional: bool):
        self.molecule_name = molecule_name
        self.position = position
        self.use_fractional = use_fractional


def show_add_molecule_dialog(
    parent: Optional[tk.Widget] = None,
    cell: Optional[np.ndarray] = None,
) -> Optional[dict]:
    """
    Show the add molecule dialog and return the result.

    Parameters
    ----------
    parent : tk.Widget, optional
        Parent widget.
    cell : numpy.ndarray, optional
        Unit cell matrix (3x3) for fractional coordinate conversion.

    Returns
    -------
    dict or None
        Dialog result with keys 'molecule_name', 'position', 'use_fractional',
        or None if cancelled.
    """
    dialog = AddMoleculeDialog(parent, cell=cell)
    dialog.wait_window()
    return dialog.result if not dialog.cancelled else None
