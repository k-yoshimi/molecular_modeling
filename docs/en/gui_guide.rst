GUI Guide
=========

The Molecular Modeling GUI is a graphical interface extending ASE-GUI.
It provides interactive operations for adding molecules, substitutions,
and exporting structure files.

Launching
---------

.. code-block:: bash

   # Launch with empty structure
   python molecular_modeling_gui.py

   # Launch with an existing structure file
   python molecular_modeling_gui.py structure.cif

   # Launch as CIF viewer
   python vis_cif.py structure.cif

Window Layout
-------------

The GUI window consists of:

.. image:: ../images/gui/gui_overview.png
   :align: center
   :width: 100%
   :alt: GUI window overview

1. **File panel (top)**: Icon button groups (described below)
2. **3D viewport**: Molecular structure display (rotatable, translatable, zoomable) with XYZ axis indicator in bottom-left
3. **Toolbar**: Translate / Rotate / Zoom buttons
4. **Bottom panel**: Molecule Tools / Structure Tools icon button groups
5. **Side panel (right)**: View settings, Element Controls, Molecule Controls (collapsible)
6. **Status bar**: Information about selected atoms

Changing Window Size
~~~~~~~~~~~~~~~~~~~~

Use the **View Settings** menu to change the window size:

- **Small (450x450)**: Compact display
- **Medium (700x700)**: Standard size (default)
- **Large (1000x800)**: Large display
- **Wide (1200x600)**: Widescreen display
- **Custom Size...**: Specify any size

You can also resize the window by dragging its borders.

File Panel (top)
~~~~~~~~~~~~~~~~

Icon buttons at the top of the window provide quick access to common operations.

.. image:: ../images/gui/gui_file_panel.png
   :align: center
   :width: 80%
   :alt: File panel (Input / Output / Edit / Pipeline)

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Group
     - Buttons
   * - Input
     - Load Structure / Load Project (TOML project)
   * - Output
     - Save Project / Export Structure / Save Selection as Molecule / Screenshot
   * - Edit
     - Undo / Redo
   * - Pipeline
     - Run Optimize / Run All / Show Results / View PNGs

Bottom Panel
~~~~~~~~~~~~

Icon button groups below the 3D viewport.

.. image:: ../images/gui/gui_bottom_panel.png
   :align: center
   :width: 80%
   :alt: Bottom panel (Molecule Tools / Structure Tools / View Operations)

- **Molecule Tools**: Select / Move / Rotate / Delete / Drag (drag mode)
- **Structure Tools**: Add Molecule / Substitution / Supercell / Cell Parameters

Menu Reference
--------------

MoleculeTools Menu
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Menu Item
     - Description
   * - Undo (Ctrl+Z)
     - Undo the last operation (up to 50 steps)
   * - Redo (Ctrl+Y)
     - Redo the last undone operation
   * - Load Structure... (Ctrl+L)
     - Load a structure file (CIF, XYZ, POSCAR, QE .in, PDB, GEN, etc.)
   * - Load CIF...
     - Load a CIF file and display the structure
   * - Add Molecule... (Ctrl+J)
     - Add a registered molecule to the system
   * - List Molecules/Groups
     - List all molecules and groups in the system
   * - Save Selection as Molecule...
     - Save selected atoms as a new molecule TOML file
   * - Global Element Substitution...
     - Replace all atoms of a specified element
   * - Molecule Substitution...
     - Replace a molecule at a specific site
   * - Molecule Atom Substitution...
     - Replace the atom list within a molecule
   * - Individual Atom Substitution...
     - Replace a single atom
   * - Batch Element Substitution...
     - Define and apply multiple element substitution rules at once (e.g., C→Si, H→F)
   * - Build Supercell...
     - Build a supercell with specified repeat factors
   * - Cell Parameters...
     - Display and edit lattice parameters (a, b, c, α, β, γ). Also shows lattice vectors and volume. Apply rescales atomic positions
   * - Save Project... (Ctrl+K)
     - Save project files (input.toml + data/\*.toml) to a chosen directory
   * - Load Project... (Ctrl+T)
     - Load a project from a directory and build the structure
   * - Measurement Log...
     - Show measurement log (distance/angle/dihedral history)

Recent Files Menu
~~~~~~~~~~~~~~~~~

Displays up to 10 recently opened CIF files.
Click a filename to open it directly.
History is saved in ``~/.molecular_modeling_recent.json``.

Export Menu
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Menu Item
     - Description
   * - Export All (QE, VASP, CIF)... (Ctrl+D)
     - Export to all formats at once
   * - Quantum ESPRESSO...
     - Export as QE input file (.in)
   * - VASP...
     - Export as POSCAR file
   * - CIF...
     - Export as CIF file

Selection Menu
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Menu Item
     - Description
   * - Select All (Ctrl+A)
     - Select all visible atoms
   * - Deselect All (Escape)
     - Clear all selections
   * - Select by Element...
     - Select atoms by element with real-time preview
   * - Select by Molecule...
     - Select atoms by molecule group with color swatches and preview
   * - Select by Sphere...
     - Select atoms within a radius with slider-based real-time preview
   * - Delete Selected Atoms (Delete)
     - Delete selected atoms (with undo support)
   * - Delete Molecule...
     - Delete the entire molecule group of the selected atoms
   * - Move Molecule/Atoms...
     - Translate selected molecule/atoms by a specified vector
   * - Rotate Molecule/Atoms...
     - Rotate selected molecule/atoms with trackball widget
   * - Focus on Selection (Ctrl+F)
     - Open the Focus dialog with 3D preview
   * - Exit Focus Mode (Ctrl+Shift+F)
     - Exit focus mode and restore all atom visibility
   * - Toggle Drag Mode
     - Toggle drag mode for direct atom dragging

Pipeline Menu
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Menu Item
     - Description
   * - Run Optimize (Voronoi) (Ctrl+P)
     - Run Voronoi IoU optimization on the current structure. Results (IoU values, Euler angles) are displayed in a dialog
   * - Run Render (Voronoi PNGs)...
     - Generate 3D Voronoi cell overlay PNG visualizations
   * - Run Export (QE, VASP, CIF)...
     - Export structure files to a selected directory
   * - Run All (Optimize → Render → Export)...
     - Run the full pipeline: optimize → render → export
   * - Show Optimization Results... (Ctrl+H)
     - Display optimization results (IoU values, Euler angles, scales) in a table. Supports CSV export and Voronoi PNG browsing

**Note**: Pipeline operations require a structure with ``mol_id`` information (i.e., built from TOML).
Processing runs in a background thread so the GUI remains responsive.

View Settings Menu
~~~~~~~~~~~~~~~~~~

.. image:: ../images/gui/view_settings_menu.png
   :align: center
   :width: 40%
   :alt: View Settings menu dropdown

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Menu Item
     - Description
   * - Atom Appearance...
     - Open the atom appearance panel to customize colors, color modes (by element/molecule group/tag), and atom sizes
   * - Color by Molecule Group
     - Color atoms by ``mol_id`` using distinct colors per molecule group
   * - Reset Colors
     - Reset colors to ASE default (jmol)
   * - Show Bonds (Ctrl+B)
     - Toggle display of bonds (lines between atoms). Bonds are auto-detected based on covalent radii
   * - Toggle Hydrogen
     - Toggle visibility of hydrogen atoms
   * - Element Visibility...
     - Open dialog to show/hide atoms by element
   * - Toggle Side Panel
     - Show or hide the side panel (quick-access view controls)
   * - Save Screenshot... (Ctrl+Shift+S)
     - Save the current 3D view as a PNG image
   * - Small / Medium / Large / Wide
     - Window size presets
   * - Custom Size...
     - Specify any window size

Help Menu
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Menu Item
     - Description
   * - Keyboard Shortcuts...
     - Show a dialog listing keyboard shortcuts and mouse controls
   * - Open User Guide...
     - Open the built HTML documentation in the browser (requires ``docs/_build/html/``; shows build instructions if not available)
   * - About...
     - Show version information

View Toolbar
------------

The toolbar below the 3D viewport provides view controls.
Each button supports press-and-hold for continuous operation.

.. image:: ../images/gui/gui_toolbar.png
   :align: center
   :width: 80%
   :alt: View Operations toolbar (Pan / Orbit / Zoom / BG)

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Section
     - Buttons
     - Description
   * - Translate
     - ◀ ▶ ▲ ▼
     - Move the view up/down/left/right
   * - Rotate
     - X- X+ Y- Y+ Z- Z+
     - Rotate ±5° around X/Y/Z axis
   * - Zoom
     - \- +
     - Zoom in/out
   * - BG Color
     - White / Black / Dark / Light / Custom
     - Change background color preset or pick a custom color
   * - Actions
     - ↩ ↪ ✖ |move| ↻
     - Undo / Redo / Delete selected / Move / Rotate (each button has a tooltip)
   * - Drag
     - |hand|
     - Toggle drag mode: click to activate, then drag selected atoms to reposition them

Tooltips
~~~~~~~~

All toolbar buttons display tooltips (floating help on mouse hover).
Hover the mouse cursor over any button to see a yellow popup with the operation description and keyboard shortcuts.

Side Panel
~~~~~~~~~~

A side panel is displayed on the right side of the window for quick access to common view settings.
The panel is scrollable — use the **mouse wheel** to scroll up/down when content overflows.

.. image:: ../images/gui/gui_side_panel.png
   :align: center
   :width: 35%
   :alt: Side panel (View Controls)

- **Focus / Unfocus**: Buttons to enter/exit focus mode
- **Color Mode**: Switch color mode via radio buttons (By Element / By Molecule / By Tag)
- **Show Bonds**: Checkbox to toggle bond display (Ctrl+B)
- **3D Spheres**: Enable/disable 3D shaded sphere rendering
- **Element Controls** *(collapsible)*: Click the header to expand/collapse. Per-element row with [visibility checkbox][color swatch][size slider]. Global atom size scale slider (0.3–3.0) included

.. image:: ../images/gui/gui_element_controls.png
   :align: center
   :width: 70%
   :alt: Side panel with Element Controls expanded

- **Molecule Controls** *(collapsible)*: Click the header to expand/collapse. Per-molecule-group row with [visibility checkbox][color swatch][select button]. **In focus mode**, checkboxes act as focus add/remove controls
- **Atom List** *(collapsible, collapsed by default)*: Click the header to expand. Shows a table of all atoms with index (#), element, X/Y/Z coordinates, and molecule group ID. Click a row to select that atom in the 3D view
- **Status info**: Real-time display of atom count, selection count, element breakdown, group count

The side panel can be shown/hidden via **View Settings** → **Toggle Side Panel**.

Status Bar
----------

The status bar at the bottom shows context-dependent information:

.. image:: ../images/gui/gui_atom_selected.png
   :align: center
   :width: 100%
   :alt: Status bar when atom is selected

- **No atoms selected**: System overview — atom count, elements, molecule groups, cell size
  (e.g., ``48 atoms | (C, H) | 4 groups | Cell: 5.5x5.4x18.7 Å``)
- **1 atom selected**: Atom details — index, element name, coordinates, molecule name
  (e.g., ``#3 Carbon (C): 1.114, 0.098, 0.812 Å | mol_id=0 (benzene)``)
- **2 atoms selected**: Inter-atomic distance (MIC corrected) with molecule name
  (e.g., ``#0(C)-#6(H): 1.089 Å | benzene``)
- **3 atoms selected**: Bond angle (MIC corrected) with molecule name
  (e.g., ``#0(C)-#1(C)-#2(C): 120.0° | benzene``)
- **4 atoms selected**: Dihedral angle (MIC corrected)

Measurements are automatically logged and can be reviewed via **MoleculeTools** → **Measurement Log...**.

Keyboard Shortcuts
------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Shortcut
     - Action
   * - Ctrl+Z
     - Undo
   * - Ctrl+Y
     - Redo
   * - Ctrl+A
     - Select all visible atoms
   * - Escape
     - Deselect all
   * - Delete / Backspace
     - Delete selected atoms
   * - Ctrl+L
     - Load structure file (CIF, XYZ, POSCAR, etc.)
   * - Ctrl+J
     - Add molecule
   * - Ctrl+D
     - Export to all formats
   * - Ctrl+P
     - Run Voronoi IoU optimization
   * - Ctrl+H
     - Show optimization results
   * - Ctrl+K
     - Save project
   * - Ctrl+T
     - Load project
   * - Ctrl+B
     - Toggle bond display
   * - Ctrl+F
     - Focus on selection (show only selected atoms)
   * - Ctrl+Shift+F
     - Exit focus mode (restore all atoms)
   * - Ctrl+Shift+S
     - Save screenshot

Standard ASE-GUI keyboard shortcuts (Ctrl+O to open, Ctrl+S to save, etc.) are also available.

Mouse Controls
--------------

Standard ASE-GUI mouse controls are available, plus a right-click context menu:

- **Left click**: Select an atom
- **Shift + left click**: Add to selection
- **Ctrl + left click**: Toggle selection
- **Right-click**: Context menu (Select Molecule Group / Move / Rotate / Delete)
- **Left drag**: Rotate view
- **Right drag** (or Shift + right drag): Translate view
- **Mouse wheel**: Zoom

Axis Indicator
~~~~~~~~~~~~~~

An XYZ axis indicator is displayed in the bottom-left corner of the 3D viewport.
It updates in real-time with the current view rotation so you can always see which direction you're looking.
Axes are color-coded: **X** = red, **Y** = green, **Z** = blue.

Right-Click Context Menu
~~~~~~~~~~~~~~~~~~~~~~~~

Right-clicking in the 3D viewport shows a context menu with quick access to:

.. image:: ../images/gui/gui_context_menu_closeup.png
   :align: center
   :width: 50%
   :alt: Right-click context menu

- **Select Molecule Group**: Select all atoms in the molecule group of the clicked atom
- **Move...**: Open move dialog
- **Rotate...**: Open rotate dialog
- **Delete Selected**: Delete currently selected atoms
- **Delete Molecule**: Delete the molecule of the selected atom

Screenshot
~~~~~~~~~~

Use **Ctrl+Shift+S** or the camera icon in the file panel to save the current 3D view as a PNG file.
If Pillow is installed, the screenshot is saved as PNG; otherwise as EPS.

MoleculeTools Menu Details
--------------------------

Focus Selection Dialog
~~~~~~~~~~~~~~~~~~~~~~

When you run **Focus on Selection** (Ctrl+F) with atoms selected, a preview window appears.

The window is **non-modal**, so clicking atoms in the main view is possible.

- **3D Preview** (top): Focused atoms are shown in color; non-focused atoms are shown as faded gray. The plot can be rotated and zoomed
- **Status**: Shows the number of focused atoms
- **Unit**:

  - **Atom**: Add/delete individual atoms
  - **Molecule**: Add/delete the whole molecule group of the clicked atom

- **Molecules**: Each molecule group is listed with a checkbox and color swatch. Toggling checkboxes adds/removes groups from the focus set

**Auto-toggle behavior**:

The main view highlights currently focused atoms as selected (highlighted).
No Add/Delete mode switching is needed — clicking automatically toggles:

- Click an **unfocused** atom → **adds** it to the focus set
- Click a **focused** atom → **removes** it from the focus set

Regular click, Ctrl+click, and Shift+click all work.

**Example workflow**:

1. Select atoms and click Focus → preview shows molecules of selected atoms, main view highlights the focus set
2. Set Unit to Molecule, click another atom → its molecule is added to the preview, main view highlight updates
3. Set Unit to Atom, click a specific focused atom → only that atom is removed
4. Use the Molecules checkboxes in the dialog to add/remove whole molecule groups
5. Click OK to apply focus mode, or Cancel to restore the original state

**Side panel controls during focus mode**:

While in focus mode, the Molecule Controls checkboxes in the side panel act as focus add/remove controls:

- Check a molecule → adds it to the focused view
- Uncheck a molecule → removes it from the focused view

Add Molecule Dialog
~~~~~~~~~~~~~~~~~~~


- **Molecule**: Select from molecules defined in ``data/*.toml``
- **Position (X, Y, Z)**: Enter placement position in Cartesian coordinates
- **Add Mode**:

  - **Single**: Add one molecule at the specified position
  - **Unit Cell**: Add all molecules from the molecule definition at once (also sets lattice vectors)

- **3D Preview**: Preview of the selected molecule

Save Selection as Molecule Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Select atoms in the 3D viewport (Shift + click for multiple selection)
2. Choose **MoleculeTools** → **Save Selection as Molecule...**
3. Enter the **Molecule name** (one molecule is saved, so ``num_molecules = 1``)
4. If needed, adjust the atoms in the **Selected Atoms** list with **Add Index** /
   **Remove Selected**
5. Click **OK** to save as ``data/<name>.toml``

.. note::

   If the molecule is split across the cell boundary, its atoms are rejoined into one
   whole molecule before saving, so you can select the pieces on either side of the cell.
   The selection should be exactly one whole molecule. If the selected atoms are not
   connected, or are bonded to atoms you did not select, a **Not One Molecule** dialog
   explains the problem and asks whether to save anyway.

Global Element Substitution Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Old Symbol**: Element to replace (selected from existing elements)
- **New Symbol**: Replacement element symbol
- All matching atoms are replaced at once

Molecule Substitution Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Molecule Index**: Index of the molecule to replace
  (auto-filled if atoms are selected)
- **New Molecule**: Select replacement molecule from dropdown

Molecule Atom Substitution Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Replace all atoms of a specific molecule with a new element list.

- **Molecule Index**: Target molecule index (auto-filled when atoms are selected)
- **New Symbols**: New element symbol list. Compact notation (e.g., ``Si6O6``) supported
- The atom count stays the same; only element symbols are changed

Individual Atom Substitution Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Replace a single atom within a molecule with a different element.

- **Molecule Index**: Target molecule index (auto-filled when atoms are selected)
- **Atom Index**: Atom index within the molecule (0-based)
- **New Symbol**: Replacement element symbol

Batch Element Substitution Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Define and apply multiple element substitution rules at once.

- Each row specifies "old element → new element"
- **+** button to add rules, **x** button to remove
- Example: C→Si and H→F can be applied simultaneously
- **OK** applies all rules at once

Build Supercell Dialog
~~~~~~~~~~~~~~~~~~~~~~

Build a supercell by repeating the current structure.

- **a, b, c**: Repetition count along each axis (integer)
- Lattice vectors and atomic positions are expanded accordingly

Move Molecule/Atoms Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~~

Translate selected atoms or molecules.

- **Scope**: "Whole molecule group" or "Selected atoms only"
- **Translation (X, Y, Z)**: Distance in Angstrom for each axis
- **+/- buttons**: Fine adjustment (press-and-hold for continuous operation)
- **Step**: Step size presets (0.01 / 0.1 / 0.5 / 1.0 Angstrom)
- **Apply**: Apply current translation, **Reset**: Return to initial position

Rotate Molecule/Atoms Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rotate selected atoms or molecules around their center of mass.

- **Scope**: "Whole molecule group" or "Selected atoms only"
- **Trackball**: Drag to rotate interactively
- **alpha(X), beta(Y), gamma(Z)**: Euler angles in degrees
- **+/- buttons**: Fine adjustment (press-and-hold for continuous operation)
- **Step**: Rotation step presets (1 / 5 / 15 / 45 / 90 degrees)
- **Apply**: Apply current rotation, **Reset**: Return to initial orientation

Select by Element Dialog
~~~~~~~~~~~~~~~~~~~~~~~~

Select atoms by element using checkboxes.

- Each element shown with a checkbox and atom count
- Selection updates **in real-time** as checkboxes change
- **Mode**: Replace (replace current selection) / Add (add to current selection)

Select by Molecule Dialog
~~~~~~~~~~~~~~~~~~~~~~~~~

Select atoms by molecule group (mol_id).

- Each group shown with a checkbox, color swatch, and atom count
- **Real-time preview** reflects selection changes immediately
- **Mode**: Replace / Add

Select by Sphere Dialog
~~~~~~~~~~~~~~~~~~~~~~~

Select all atoms within a radius from the currently selected atom.

- **Radius slider**: Drag to adjust radius (real-time preview)
- Shows atom count and element breakdown within the selection

Atom Appearance Panel
~~~~~~~~~~~~~~~~~~~~~

Open from **View Settings** → **Atom Appearance...**. The panel has three tabs:

**Color Mode tab**:

- **Default (by element)**: Standard jmol element colors
- **By molecule group**: Color by ``mol_id`` using 20 distinct colors (cycled)
- **By tag**: Color by tag values
- Custom arrays are automatically listed as additional color modes

**Element Colors tab**:

- Change the color of each element individually
- Click **Change...** to open the color picker
- **Reset to Default** restores jmol default colors

**Size tab**:

- **Atom Size Scale**: Global atom size multiplier (0.1–3.0, slider)
- **Per-Element Radius Override**: Adjust covalent radius per element

Cell Parameters Dialog
~~~~~~~~~~~~~~~~~~~~~~

Open from **MoleculeTools** → **Cell Parameters...**.

- **Lattice lengths**: Display and edit a, b, c (in Å)
- **Lattice angles**: Display and edit α, β, γ (in degrees)
- **Lattice vectors**: Current lattice vectors (read-only)
- **Volume**: Cell volume (Å³)
- **Apply**: Update cell with the entered values (atomic positions are scaled)

Export Dialog
~~~~~~~~~~~~~

- Select output directory
- Select output format(s) (QE / VASP / CIF / All)
- Click **OK** to export

.. |move| unicode:: U+2726 .. FOUR DIAMOND ASTERISK
   :trim:
.. |hand| replace:: [Drag]
