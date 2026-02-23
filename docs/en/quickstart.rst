Quick Start
===========

This guide explains how to run the pipeline and check results in the shortest steps.

Running the Pipeline
--------------------

With the virtual environment activated, run:

.. code-block:: bash

   cd molecular_modeling
   source .venv/bin/activate
   python run.py

Upon completion, you will see output like:

.. code-block:: text

   === Build system from TOML ===
   ...
   === Optimize Voronoi overlap ===
   ...
   === Write structure files ===
   Written: ./structure_out/qe.in
   Written: ./structure_out/POSCAR
   Written: ./structure_out/INCAR
   Written: ./structure_out/KPOINTS
   Written: ./structure_out/structure.cif

Input Files
-----------

The pipeline reads the following files:

- ``input/input.toml`` — Main configuration file (cell vectors, molecule list, repetitions, etc.)
- ``data/<name>.toml`` — Molecule definition files (atomic positions, lattice vectors, etc.)

Output Files
------------

The following files are generated in ``source/structure_out/``:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Description
   * - ``qe.in``
     - Quantum ESPRESSO input file
   * - ``POSCAR``
     - VASP structure file
   * - ``INCAR``
     - VASP calculation settings
   * - ``KPOINTS``
     - VASP k-point mesh
   * - ``structure.cif``
     - CIF structure file

Voronoi cell visualizations are generated in ``source/voronoi_pngs/``.

Launching the GUI
-----------------

To launch the GUI:

.. code-block:: bash

   python molecular_modeling_gui.py

To open an existing CIF file:

.. code-block:: bash

   python vis_cif.py structure.cif

The GUI provides:

- Adding and removing molecules
- Element substitutions
- Structure file export
- 3D view manipulation (rotation, translation, zoom)

See :doc:`gui_guide` for details.

Next Steps
----------

- :doc:`tutorial` — Detailed tutorial with benzene + naphthalene example
- :doc:`gui_guide` — Full GUI operations guide
- :doc:`cli_reference` — Command-line tools usage
- :doc:`toml_reference` — Complete TOML configuration reference
- :doc:`architecture` — Internal architecture explanation
