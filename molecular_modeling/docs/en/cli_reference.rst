CLI Reference
=============

run.py — Pipeline Execution
----------------------------

Runs the main pipeline.

.. code-block:: bash

   python run.py

**Processing flow**:

1. Build system from ``input/input.toml`` and ``data/*.toml``
2. Index atoms into molecule groups
3. Optimize molecular orientations via Voronoi IoU maximization
4. Check Voronoi volume tiling consistency
5. Render Voronoi cell overlay visualizations → ``source/voronoi_pngs/``
6. Write structure files → ``source/structure_out/``

**Input files**:

- ``input/input.toml`` — Main configuration
- ``data/<name>.toml`` — Molecule definitions

**Output files**:

- ``source/structure_out/qe.in`` — Quantum ESPRESSO input
- ``source/structure_out/POSCAR`` — VASP structure file
- ``source/structure_out/INCAR`` — VASP settings
- ``source/structure_out/KPOINTS`` — VASP k-point settings
- ``source/structure_out/structure.cif`` — CIF file
- ``source/voronoi_pngs/*.png`` — Voronoi visualizations

molecular_modeling_gui.py — GUI Launch
--------------------------------------

Launches the ASE-GUI extension for molecular modeling.

.. code-block:: bash

   # Launch with empty structure
   python molecular_modeling_gui.py

   # Launch with a structure file
   python molecular_modeling_gui.py [structure_file]

**Arguments**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Argument
     - Description
   * - ``structure_file`` (optional)
     - Path to any structure file readable by ASE

If no file is specified, a 20x20x20 Angstrom box is created.

vis_cif.py — CIF Viewer
------------------------

Displays a CIF file in the GUI.

.. code-block:: bash

   # Open a specific file
   python vis_cif.py <path_to_file.cif>

   # Open file dialog
   python vis_cif.py

**Arguments**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Argument
     - Description
   * - ``path_to_file.cif`` (optional)
     - Path to CIF file. If omitted, a file chooser dialog opens.

**Note**: If ASE's CIF parser fails, the tool retries after stripping
``_citation*`` tags as a fallback.

edit_and_run.py — GUI Edit + Optimization
------------------------------------------

Pipeline that allows manual GUI editing before Voronoi optimization.

.. code-block:: bash

   python edit_and_run.py

**Processing flow**:

1. Build system from TOML
2. **Launch GUI** (user edits the structure)
3. After closing the GUI, continue processing with the edited system
4. Build molecule groups
5. Voronoi IoU optimization
6. Visualization and export

Running Tests
-------------

.. code-block:: bash

   # Run all tests
   pytest

   # Specific test file
   pytest tests/test_voronoi_geom.py

   # Specific test case
   pytest tests/test_voronoi_geom.py::TestPolyhedron::test_volume

   # Skip slow tests
   pytest -m "not slow"
