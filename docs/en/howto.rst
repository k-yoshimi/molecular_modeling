How-To Guide
============

A reverse-lookup reference: find what you want to do, then see how to do it.
For details, follow the links to the relevant sections.

.. contents:: Categories
   :local:
   :depth: 1

View / Open Crystal Structures
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Open a CIF file and view the structure
     - GUI: **MoleculeTools > Load CIF...**
       or from the command line: ``python vis_cif.py file.cif``
       (see :doc:`cli_reference` for details)
   * - Open QE / VASP / XYZ or other formats
     - GUI: **MoleculeTools > Load Structure...** (Ctrl+L).
       Any format supported by ASE can be read
   * - Load a TOML project
     - GUI: **MoleculeTools > Load Project...** (Ctrl+T).
       Select ``input.toml`` to restore molecule group information
   * - Quickly reopen a recent file
     - GUI: **Recent Files** menu shows the last 10 files

Build Supercells
----------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Repeat a unit cell to create a supercell
     - GUI: **MoleculeTools > Build Supercell...**
       and specify repetition counts along [a, b, c].
       See :doc:`tutorial` Step 2 for details
   * - View or modify cell parameters (a, b, c, alpha, beta, gamma)
     - GUI: **MoleculeTools > Cell Parameters...**
   * - Configure via TOML
     - Set ``supercell_repeat = [2, 2, 1]`` in ``input.toml``.
       See :doc:`toml_reference` for details

Add and Place Molecules
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Add a molecule to the supercell
     - GUI: **MoleculeTools > Add Molecule...** (Ctrl+J).
       Select from TOML molecule definitions in ``data/`` folder, specify position and orientation.
       See :doc:`gui_guide` Add Molecule dialog for details
   * - Move a molecule to a different position
     - Select atom(s) in the molecule > **Selection > Move Molecule/Atoms...**.
       Adjust displacement with +/- buttons for each X/Y/Z direction
   * - Rotate a molecule
     - Select atom(s) in the molecule > **Selection > Rotate Molecule/Atoms...**.
       Use the trackball or enter Euler angles directly
   * - Drag-move molecules with the mouse
     - Enable **Selection > Toggle Drag Mode**.
       Click and drag to reposition directly
   * - Delete a molecule
     - Select atom(s) in the molecule > **Selection > Delete Molecule...**
   * - Save selected atoms as a new molecule
     - Select atoms > **MoleculeTools > Save Selection as Molecule...**.
       A TOML file is created in ``data/``

Substitute Atoms / Elements
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Replace all atoms of one element globally (e.g., all H to D)
     - GUI: **MoleculeTools > Global Element Substitution...**.
       Or set ``[global_atom_substitutions]`` in ``input.toml``
   * - Replace an entire molecule with a different one
     - GUI: **MoleculeTools > Molecule Substitution...**.
       Orientation is preserved via PCA alignment
   * - Replace a single atom with a different element
     - GUI: **MoleculeTools > Individual Atom Substitution...**
   * - Apply multiple element substitutions at once (e.g., C to Si and H to F)
     - GUI: **MoleculeTools > Batch Element Substitution...**.
       Specify multiple rules in a single dialog
   * - Redefine all atoms in a specific molecule with new elements
     - GUI: **MoleculeTools > Molecule Atom Substitution...**

Optimize Molecular Orientations
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Automatically determine molecular orientations via Voronoi optimization
     - GUI: **Pipeline > Run Optimize (Voronoi)** (Ctrl+P).
       Maximizes Voronoi cell overlap with the reference crystal.
       See :doc:`appendix_voronoi` for algorithm details
   * - Check optimization results (IoU, rotation angles, scale)
     - GUI: **Pipeline > Show Optimization Results...** (Ctrl+H).
       CSV export is also available
   * - Visualize Voronoi cell overlaps
     - GUI: **Pipeline > Run Render (Voronoi PNGs)...**.
       Blue = system, Red = dictionary, Orange = intersection
   * - Run the full pipeline at once
     - GUI: **Pipeline > Run All (Optimize > Render > Export)...**
   * - Run from the command line
     - ``python run.py``.
       See :doc:`cli_reference` for details

Predict Crystal Structures
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Predict candidate crystal structures from one molecule
     - ``molecular-modeling-csp input/csp.toml`` (needs ``pip install -e ".[csp]"``).
       See :doc:`csp_guide`
   * - Make the molecule file for a prediction from a known crystal
     - GUI: select one molecule, then **MoleculeTools > Save Selection as Molecule...**.
       See :doc:`csp_guide`
   * - Compare predicted candidates with an experimental CIF
     - ``python -m molecular_modeling.csp.validate csp_out exp.cif``.
       See :doc:`csp_guide` and :doc:`cli_reference`

Export Structure Files
----------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Create a Quantum ESPRESSO input file
     - GUI: **Export > Quantum ESPRESSO...**
       or **Export > Export All (QE, VASP, CIF)...** (Ctrl+D)
   * - Create VASP files (POSCAR/INCAR/KPOINTS)
     - GUI: **Export > VASP...**
       or **Export > Export All...**
   * - Save as a CIF file
     - GUI: **Export > CIF...**
       or **Export > Export All...**
   * - Export all three formats at once
     - GUI: **Export > Export All (QE, VASP, CIF)...** (Ctrl+D)
   * - Save as a TOML project
     - GUI: **MoleculeTools > Save Project...** (Ctrl+K).
       Saves ``input.toml`` and ``data/*.toml``
   * - Save a screenshot of the 3D view
     - GUI: **View Settings > Save Screenshot...** (Ctrl+Shift+S)

Change Display Settings
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Show bonds
     - **View Settings > Show Bonds** (Ctrl+B).
       Also available from the side panel checkbox
   * - Hide hydrogen atoms
     - **View Settings > Toggle Hydrogen**
   * - Show/hide specific elements
     - **View Settings > Element Visibility...**
   * - Color molecules by group
     - **View Settings > Color by Molecule Group**
   * - Change atom colors and sizes
     - **View Settings > Atom Appearance...**
   * - Change the background color
     - Toolbar **BG Color** button (White/Black/Dark/Light/Custom)
   * - Resize the window
     - **View Settings > Small / Medium / Large / Wide / Custom Size...**
   * - Show/hide the side panel
     - **View Settings > Toggle Side Panel**

Select Atoms
------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Select all atoms
     - **Selection > Select All** (Ctrl+A)
   * - Clear selection
     - **Selection > Deselect All** (Escape)
   * - Select atoms of a specific element (e.g., all C)
     - **Selection > Select by Element...**
   * - Select a specific molecule group
     - **Selection > Select by Molecule...**
   * - Select atoms near a point (sphere selection)
     - Select an atom > **Selection > Select by Sphere...**
   * - Focus on selected molecules only
     - **Selection > Focus on Selection** (Ctrl+F).
       Press Ctrl+Shift+F to return to full view
   * - Measure distance between atoms
     - Select 2 atoms; the distance is shown in the status bar
   * - Measure a bond angle
     - Select 3 atoms; the angle is shown in the status bar
   * - Measure a dihedral angle
     - Select 4 atoms; the dihedral angle is shown in the status bar
   * - View measurement history
     - **MoleculeTools > Measurement Log...**

Undo / Redo
------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - Undo the last operation
     - **MoleculeTools > Undo** (Ctrl+Z). Up to 50 steps
   * - Redo an undone operation
     - **MoleculeTools > Redo** (Ctrl+Y)

TOML Configuration Files
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - I want to...
     - How to do it
   * - See all input.toml parameters
     - :doc:`toml_reference` has full field descriptions
   * - Learn how to write molecule definition files (data/\*.toml)
     - See the "Molecule Definition File" section in :doc:`toml_reference`
   * - Learn how to write substitution rules
     - See the "Substitution Settings" section in :doc:`toml_reference`
   * - See a working example
     - :doc:`tutorial` has a benzene + naphthalene mixed crystal example
