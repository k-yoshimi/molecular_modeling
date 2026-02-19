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

   # Launch GUI from pipeline (edit_and_run.py)
   python edit_and_run.py

Window Layout
-------------

.. image:: /_static/placeholder.png
   :width: 600
   :alt: GUI window layout

*Replace this with an actual screenshot*

The GUI window consists of:

1. **Menu bar**: Various operation menus
2. **3D viewport**: Molecular structure display (rotatable, translatable, zoomable)
3. **Toolbar**: Translate / Rotate / Zoom buttons
4. **Status bar**: Information about selected atoms

Changing Window Size
~~~~~~~~~~~~~~~~~~~~

Use the **Window Size** menu to change the window size:

- **Small (450x450)**: Compact display
- **Medium (700x700)**: Standard size (default)
- **Large (1000x800)**: Large display
- **Wide (1200x600)**: Widescreen display
- **Custom Size...**: Specify any size

You can also resize the window by dragging its borders.

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
   * - Select by Element...
     - Select atoms by element symbol. Multiple elements can be chosen. Supports replace/add selection modes
   * - Select by Molecule...
     - Select atoms by molecule group (mol_id). Shows group name and atom count
   * - Select by Sphere...
     - Select atoms within a given radius of the currently selected atom
   * - Delete Selected Atoms
     - Delete selected atoms (with undo support and molecule metadata update)
   * - Delete Molecule...
     - Delete the entire molecule group of the selected atoms
   * - Move Molecule/Atoms...
     - Translate selected molecule/atoms by a specified vector (whole molecule or selected only)
   * - Rotate Molecule/Atoms...
     - Rotate selected molecule/atoms around their center of mass by Euler angles

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
     - Display optimization results (IoU values, Euler angles, scales) in a table

**Note**: Pipeline operations require a structure with ``mol_id`` information (i.e., built from TOML).
Processing runs in a background thread so the GUI remains responsive.

View Settings Menu
~~~~~~~~~~~~~~~~~~

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
   * - Toggle Hydrogen
     - Toggle visibility of hydrogen atoms
   * - Element Visibility...
     - Open dialog to show/hide atoms by element
   * - Toggle Side Panel
     - Show or hide the side panel (quick-access view controls)
   * - Small / Medium / Large / Wide
     - Window size presets
   * - Custom Size...
     - Specify any window size

View Toolbar
------------

The toolbar below the 3D viewport provides view controls.
Each button supports press-and-hold for continuous operation.

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
   * - Actions
     - ↩ ↪ ✖ ✥ ↻
     - Undo / Redo / Delete selected / Move / Rotate (each button has a tooltip)

Side Panel
~~~~~~~~~~

A side panel is displayed on the right side of the window for quick access to common view settings.

- **Show Hydrogen**: Toggle hydrogen atom visibility via checkbox
- **Color Mode**: Switch color mode via radio buttons (By Element / By Molecule / By Tag)
- **Atom Size**: Adjust atom size scale with a slider (0.3–3.0)
- **Element Visibility**: Toggle individual element visibility with checkboxes

The side panel can be shown/hidden via **View Settings** → **Toggle Side Panel**.

Status Bar
----------

The status bar at the bottom shows context-dependent information:

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

Standard ASE-GUI keyboard shortcuts (Ctrl+O to open, Ctrl+S to save, etc.) are also available.

Mouse Controls
--------------

Standard ASE-GUI mouse controls are available:

- **Left click**: Select an atom
- **Shift + left click**: Add to selection
- **Ctrl + left click**: Toggle selection
- **Left drag**: Rotate view
- **Right drag** (or Shift + right drag): Translate view
- **Mouse wheel**: Zoom

Dialog Details
--------------

Add Molecule Dialog
~~~~~~~~~~~~~~~~~~~

.. image:: /_static/placeholder.png
   :width: 400
   :alt: Add Molecule dialog

*Replace this with an actual screenshot*

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
3. Enter molecule name and number of molecules per unit cell
4. Click **OK** to save as ``data/<name>.toml``

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
