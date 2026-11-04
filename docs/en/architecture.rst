Architecture
============

Pipeline Flow
-------------

The main pipeline (``run.py`` → ``source/main.py``) processes in the following order:

.. code-block:: text

   input/input.toml + data/*.toml
          |
          v
   +-------------------------+
   |  build_system_from_toml |  TOML -> ASE Atoms construction
   |                         |  Supercell generation, substitutions
   +----------+--------------+
              |  system (Atoms + system.info)
              v
   +-------------------------+
   |  build_molecule_groups  |  Atom grouping by mol_id
   +----------+--------------+
              |  groups: list[list[int]]
              v
   +-------------------------+
   | optimize_voronoi_overlap|  Euler angle grid search +
   |                         |  Nelder-Mead optimization
   +----------+--------------+
              |  optimized_angles_deg, scales, iou
              v
   +-------------------------------+
   | check_system_voronoi_volume   |  Voronoi volume
   | _tiling (voronoi_geom.py)     |  tiling verification
   +----------+--------------------+
              |
              v
   +-------------------------+
   | render_voronoi_overlays |  3D matplotlib visualization
   +----------+--------------+
              |  voronoi_pngs/*.png
              v
   +-------------------------+
   |     write_structure     |  QE / VASP / CIF export
   +----------+--------------+
              |
              v
   structure_out/  (qe.in, POSCAR, INCAR, KPOINTS, structure.cif)

Module Descriptions
-------------------

build_system_from_toml.py
~~~~~~~~~~~~~~~~~~~~~~~~~

**Input**: ``input/input.toml``, ``data/*.toml``

**Output**: ASE ``Atoms`` object

Main operations:

1. Read and parse TOML files
2. Convert fractional → Cartesian coordinates (when ``use_fractional`` is true)
3. Build supercell (repetition tiling)
4. Place each molecule (translation to center-of-mass position)
5. Apply 4 types of substitution tables
6. Generate ``mol_id`` array (maps each atom to its site)
7. Set metadata in ``system.info``

build_molecule_groups.py
~~~~~~~~~~~~~~~~~~~~~~~~

**Input**: ``Atoms`` object (with ``mol_id`` array)

**Output**: ``list[list[int]]`` — atom index list per site

Scans the ``mol_id`` array and groups atom indices with the same value.

voronoi_geom.py
~~~~~~~~~~~~~~~

The core geometry module, providing:

- **Polyhedron class**: Convex polyhedron representation (vertices + faces)
- **Voronoi cell computation**: Periodic Voronoi tessellation via SciPy halfspace intersection
- **IoU computation**: Volume-based Intersection over Union using trimesh and manifold3d
- **Rotation and scaling**: Polyhedron rotation by Euler angles and uniform scaling
- **Periodic boundary handling**: Wrapping at cell boundaries

optimize_voronoi_overlap.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Process**: Search for Euler angles that maximize Voronoi IoU per site

Algorithm:

1. Compute system Voronoi cell volume ``V_sys``
2. Compute dictionary Voronoi cell volume ``V_dict``
3. Set uniform scale ``s = (V_sys / V_dict)^(1/3)``
4. Coarse grid search over Euler angles (alpha, beta, gamma) (default 13x13x13)
5. Refine from best point using Nelder-Mead
6. Apply optimal rotation and scaling to atom positions

Results are stored in ``system.info``:

- ``optimized_angles_deg``: Optimal Euler angles (Nx3)
- ``optimized_scales``: Optimal scales (N,)
- ``optimized_iou``: IoU values (N,)

render_voronoi_overlays.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

Renders 3D matplotlib plots overlaying system and dictionary Voronoi cells.

- Blue translucent: System Voronoi cell
- Red translucent: Dictionary Voronoi cell
- One PNG per site

write_structure.py
~~~~~~~~~~~~~~~~~~

Uses ASE I/O to write structure files in the following formats:

- **Quantum ESPRESSO**: ``qe.in`` (with SCF calculation settings)
- **VASP**: ``POSCAR``, ``INCAR``, ``KPOINTS``
- **CIF**: ``structure.cif``

molecule_unwrap.py
~~~~~~~~~~~~~~~~~~

Rejoins molecules that are split across periodic boundaries (breadth-first walk over the bond graph,
1.2 × natural covalent radii), checks that a set of atoms is one connected molecule, and lists bonded
atom pairs. Used by the GUI (Save Selection as Molecule) and by the CSP package.

csp/ — Crystal Structure Prediction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A separate pipeline (``molecular-modeling-csp`` → ``source/csp/cli.py`` → ``pipeline.run_csp``) that
builds candidate crystals from one isolated molecule. It needs the ``csp`` extra (PyXtal, tblite,
pymatgen); see :doc:`csp_guide`.

.. code-block:: text

   input/csp.toml + molecule file
          |
          v
   generation (gen_random, gen_voronoi, gen_seed) -> prescreen (screen)
          -> relaxation (relax, GFN2-xTB) -> deduplication (dedupe)
          -> csp_out/ (candidates.csv, NNN_<id>.cif, run.log)

- ``config.py``: ``CSPConfig`` and ``load_csp_config`` (TOML settings, blind/informed modes, range checks)
- ``pipeline.py``: ``run_csp`` (stages, time budget, dependency check, outputs); ``cli.py``: command-line entry point
- ``molecule.py``: molecule loading, volume, RMSD and rigid alignment (``match_and_align``)
- ``gen_random.py``, ``gen_voronoi.py``, ``gen_seed.py``: random (PyXtal), Voronoi-oriented and seed-based generators
- ``screen.py``: overlap and packing checks, rigid Lennard-Jones prescreen
- ``relax.py``: parallel relaxation and the deformation / bond checks
- ``dedupe.py``: merging identical structures
- ``candidate.py``: the ``Candidate`` record passed between stages
- ``deps.py``: check that the ``csp`` extra can be imported
- ``validate.py``: comparison with experimental CIFs (``python -m molecular_modeling.csp.validate``)

system.info Metadata
--------------------

The ``info`` dictionary on ASE ``Atoms`` objects is used to pass data
between pipeline stages.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Type
     - Description
   * - ``groups``
     - list[list[int]]
     - Atom indices per site
   * - ``site_mol_name``
     - list[str]
     - Molecule name per site
   * - ``dict_meta``
     - dict
     - Per-molecule metadata (lattice vectors, COMs, voronoi_ref_id)
   * - ``rng``
     - dict
     - Random seeds (global, optimize, render)
   * - ``optimized_angles_deg``
     - ndarray(N,3)
     - Optimized Euler angles (degrees)
   * - ``optimized_scales``
     - ndarray(N,)
     - Optimized scale factors
   * - ``optimized_iou``
     - ndarray(N,)
     - Optimized IoU values

Atom Arrays
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Array
     - Type
     - Description
   * - ``mol_id``
     - int[N_atoms]
     - Site index for each atom

GUI Architecture
----------------

The GUI follows an MVC (Model-View-Controller) pattern.

.. code-block:: text

   +------------------------------------------------------------+
   |              MolecularModelingGUI (View)                   |
   |                                                            |
   |  Menus: MoleculeTools / Recent Files / Export /            |
   |         Selection / Pipeline / View Settings / Help        |
   |                                                            |
   |  +-----------------+  +------------------------------+    |
   |  |  File Panel     |  |  Side Panel (right)          |    |
   |  |  (Input/Output/ |  |  - Color mode                |    |
   |  |   Edit/Pipeline)|  |  - Show Bonds / 3D Spheres   |    |
   |  +-----------------+  |  - Element Controls           |    |
   |  +-----------------+  |  - Molecule Controls          |    |
   |  |  3D Viewport    |  |  - Atom List                  |    |
   |  |  (ASE Canvas)   |  +------------------------------+    |
   |  +-----------------+                                       |
   |  +-----------------+  +------------------------------+    |
   |  |  Toolbar        |  |  Bottom Panel                |    |
   |  |  (Pan/Rotate/   |  |  (Molecule Tools /           |    |
   |  |   Zoom/BG color)|  |   Structure Tools)           |    |
   |  +-----------------+  +------------------------------+    |
   |  Status Bar (selection info / distance / angle / dihedral)|
   +------------------------------------+-----------------------+
                                        | Dialog calls
                                        v
   +------------------------------------------------------------+
   |                     Dialogs (View)                         |
   |  AddMolecule / Substitution (4 types) / BatchSubstitution |
   |  Export / SaveMolecule / MoleculeList / Supercell          |
   |  CellParameters / Focus / Move / Rotate / SelectBy...     |
   +----------------------------+-------------------------------+
                                | Return results
                                v
   +------------------------------------------------------------+
   |           MolecularController (Controller)                 |
   |  add_molecule / substitute_* / build_supercell /           |
   |  run_optimize / run_render / run_export / get_atoms        |
   +----------------------------+-------------------------------+
                                | Modify Atoms
                                v
   +------------------------------------------------------------+
   |               ASE Atoms (Model)                            |
   |  positions / symbols / cell / pbc / info / arrays          |
   +------------------------------------------------------------+

GUI module structure:

- ``molecular_modeling_gui.py``: Main GUI class (extends ASE-GUI). Manages window layout, menus, toolbar, side panel, and status bar
- ``source/gui/controller.py``: Bridge between GUI dialogs and backend operations
- ``source/gui/theme.py``: Theme configuration (colors, fonts, spacing, ttk theme setup)
- ``source/gui/base_dialog.py``: Base class for modal dialogs (OK/Cancel, validation, coordinate input)
- ``source/gui/molecule_registry.py``: Dynamic molecule definition loader from ``data/*.toml``
- ``source/gui/add_molecule_dialog.py``: Add molecule dialog
- ``source/gui/substitution_dialogs.py``: Four substitution dialog types + batch element substitution
- ``source/gui/export_dialog.py``: Export dialog
- ``source/gui/save_molecule_dialog.py``: Save selected atoms as molecule dialog
- ``source/gui/molecule_list_dialog.py``: Molecule list dialog
- ``source/gui/supercell_dialog.py``: Supercell settings dialog
