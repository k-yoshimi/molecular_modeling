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
   * - Load CIF...
     - Load a CIF file and display the structure
   * - Add Molecule...
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

Export Menu
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Menu Item
     - Description
   * - Export All (QE, VASP, CIF)...
     - Export to all formats at once
   * - Quantum ESPRESSO...
     - Export as QE input file (.in)
   * - VASP...
     - Export as POSCAR file
   * - CIF...
     - Export as CIF file

Pipeline Menu
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Menu Item
     - Description
   * - Run Optimize (Voronoi)
     - Run Voronoi IoU optimization on the current structure. Results (IoU values, Euler angles) are displayed in a dialog
   * - Run Render (Voronoi PNGs)...
     - Generate 3D Voronoi cell overlay PNG visualizations
   * - Run Export (QE, VASP, CIF)...
     - Export structure files to a selected directory
   * - Run All (Optimize → Render → Export)...
     - Run the full pipeline: optimize → render → export

**Note**: Pipeline operations require a structure with ``mol_id`` information (i.e., built from TOML).
Processing runs in a background thread so the GUI remains responsive.

Window Size Menu
~~~~~~~~~~~~~~~~

Change the window size using presets or a custom value.

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

Status Bar
----------

The status bar at the bottom shows context-dependent information:

- **No atoms selected**: System overview — atom count, elements, molecule groups, cell size
  (e.g., ``48 atoms | (C, H) | 4 groups | Cell: 5.5x5.4x18.7 Å``)
- **1 atom selected**: Atom details — index, element name, coordinates, molecule name
  (e.g., ``#3 Carbon (C): 1.114, 0.098, 0.812 Å | mol_id=0 (benzene)``)
- **2 atoms selected**: Inter-atomic distance
- **3 atoms selected**: Bond angles
- **4 atoms selected**: Dihedral angle

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

Export Dialog
~~~~~~~~~~~~~

- Select output directory
- Select output format(s) (QE / VASP / CIF / All)
- Click **OK** to export
