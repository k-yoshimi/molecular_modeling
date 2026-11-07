Command-Line Tools
==================

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
5. Render Voronoi cell overlay visualizations → ``voronoi_pngs/``
6. Write structure files → ``structure_out/``

**Input files** (see :doc:`toml_reference` for details):

- ``input/input.toml`` — Main configuration
- ``data/<name>.toml`` — Molecule definitions

**Output files**:

- ``structure_out/qe.in`` — Quantum ESPRESSO input
- ``structure_out/POSCAR`` — VASP structure file
- ``structure_out/INCAR`` — VASP settings
- ``structure_out/KPOINTS`` — VASP k-point settings
- ``structure_out/structure.cif`` — CIF file
- ``voronoi_pngs/*.png`` — Voronoi visualizations

molecular_modeling_gui.py — GUI Launch
--------------------------------------

Launches the ASE-GUI extension for molecular modeling.
See :doc:`gui_guide` for details.

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

molecular-modeling-csp — Crystal Structure Prediction
-----------------------------------------------------

Generates, relaxes and ranks candidate molecular crystals. Needs the ``csp`` extra
(``pip install -e ".[csp]"``; for OrbMol-v2 the ``csp-ml`` extra); see :doc:`csp_guide` for the
settings and the output.

.. code-block:: bash

   molecular-modeling-csp [config] [--relax-method METHOD]

**Arguments**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Argument
     - Description
   * - ``config`` (optional)
     - CSP settings file (TOML). Default: ``input/csp.toml``
   * - ``--relax-method`` (optional)
     - Relaxation method: ``GFN2-xTB``, ``GFN1-xTB`` or ``OrbMol-v2``. Overrides ``relax_method`` in
       the settings file; ``run.log`` records that it came from the command line

**Output** (in ``out_dir``, default ``csp_out/``): ``candidates.csv`` (ranking), ``NNN_<id>.cif`` for
the top 50 candidates, and ``run.log``. If a package of the ``csp`` extra is missing, the command stops
before generating anything and says so. Input errors (bad settings, missing or unreadable files, an
``out_dir`` in use by another run, an OrbMol-v2 model that cannot be loaded, e.g. without network
access on first use) are printed as ``error: <message>`` with exit code 2; Ctrl-C exits
with code 130 and leaves no results. See :doc:`csp_guide`, "Reusing out_dir".

python -m molecular_modeling.csp.validate — Compare with Experiment
-------------------------------------------------------------------

Compares the candidate CIFs of a CSP run with experimental CIFs (cluster RMSD of 15 molecules).
Also needs the ``csp`` extra (pymatgen, spglib); without it the command stops with the install hint.
See :doc:`csp_guide`, "Validation against experiment", for how to read the report.

.. code-block:: bash

   python -m molecular_modeling.csp.validate <out_dir> <exp.cif> [<exp.cif> ...] [options]

**Arguments**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Argument
     - Description
   * - ``out_dir``
     - Output directory of ``molecular-modeling-csp`` (``candidates.csv`` and CIFs). The command
       stops with an error if it does not exist or holds no candidates
   * - ``exp.cif``
     - One or more experimental CIF files (disorder-free; hydrogens may be missing)
   * - ``--n-top N``
     - A match within rank ``N`` counts as success (default 20; at least 1)
   * - ``--rmsd-cut X``
     - RMSD\ :sub:`15` threshold in Å for a match (default 0.8; a positive number)
   * - ``--flip A,B``
     - Two different elements treated as one for the separate flip-match check (default ``O,S``; the
       check is on by default; case does not matter)
   * - ``--no-flip``
     - Skip the flip-match check (useful for molecules without O or S); cannot be combined with ``--flip``

Input errors (bad options, unreadable CIFs, a missing or busy output directory, results of a run that
did not finish, a lock file that cannot be created, more than one MPI process) end with a usage
message and exit code 2. Validation holds the output directory's lock file while it runs, so a CSP run
cannot start in that directory meanwhile; it therefore needs write access to the output directory and
never reads it without the lock.

The report is printed and written to ``<out_dir>/validation.md``.

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
