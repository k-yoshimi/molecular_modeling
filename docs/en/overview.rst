Overview
========

Molecular Modeling Toolkit is a molecular modeling toolkit built on
ASE (Atomic Simulation Environment). It builds supercells from TOML-based
molecular definitions, optimizes molecular orientations via Voronoi cell
IoU maximization, and exports to Quantum ESPRESSO / VASP / CIF formats.

What This Toolkit Can Do
------------------------

- **Crystal structure construction**: Load CIF files, create supercells, add and place molecules
- **Automatic molecular orientation optimization**: Determine optimal orientations by maximizing Voronoi cell overlap with a reference crystal
- **Atom and molecule substitution**: Bulk element replacement, molecule swapping, individual atom changes (4 substitution methods)
- **Structure file export**: Batch export to Quantum ESPRESSO / VASP / CIF formats
- **Interactive GUI**: 3D visualization, mouse operations, intuitive editing via dialogs
- **Crystal structure prediction (CSP)**: Generate candidate crystals from one isolated molecule, relax them with
  GFN2-xTB, rank them by energy and compare them with experimental CIFs (see :doc:`csp_guide`)

Processing Pipeline
-------------------

The core processing consists of a six-stage pipeline:

1. **Build** -- Assemble a supercell from TOML configuration and molecule definitions, applying atomic substitutions
2. **Group** -- Index atoms by molecular site for subsequent per-molecule operations
3. **Optimize** -- Maximize Voronoi cell IoU via Euler angle grid search + Nelder-Mead refinement
4. **Check** -- Verify that Voronoi cell volumes sum to the supercell volume
5. **Render** -- Generate 3D Voronoi cell overlay visualizations (comparison PNGs before and after optimization)
6. **Export** -- Batch export to Quantum ESPRESSO, VASP, and CIF formats

.. code-block:: text

   input.toml + data/*.toml (molecule definitions)
        |
        v
   [1. Build] Supercell generation + atomic substitutions
        |
        v
   [2. Group] Per-site atom indexing
        |
        v
   [3. Optimize] Voronoi IoU maximization (grid search + Nelder-Mead)
        |
        v
   [4. Check] Voronoi volume tiling consistency check
        |
        v
   [5. Render] Voronoi overlay PNG generation
        |
        v
   [6. Export] qe.in / POSCAR / INCAR / KPOINTS / structure.cif

Voronoi Optimization
--------------------

The Voronoi optimization, a unique feature of this toolkit, proceeds as follows:

- **Voronoi cell construction**: Cells are computed from molecular center-of-mass coordinates with periodic boundary conditions
- **Volume scaling**: A uniform scale corrects for volume differences between the system and dictionary Voronoi cells
- **Rotation optimization**: A 13x13x13 Euler angle grid search provides a coarse solution, refined by the Nelder-Mead method
- **Evaluation metric**: IoU (Intersection over Union) = intersection volume / union volume (0--1, higher is better)

See :doc:`appendix_voronoi` for details.

Substitution System
-------------------

Four types of substitution enable flexible crystal structure customization.
Substitutions are applied in the following order:

1. **Molecule substitution**: Replace an entire molecule at a specific cell position with a different molecule
2. **Global element substitution**: Replace all atoms of one element across the entire system (e.g., all H to He)
3. **Molecule-atom substitution**: Redefine the atom list of a specific molecule instance (overrides global substitution results)
4. **Individual atom substitution**: Swap a single atom specified by cell, molecule, and atom index

Input and Output
----------------

**Input files:**

- ``input/input.toml`` -- Lattice vectors, supercell repetition, molecule list, center-of-mass coordinates, substitution rules
- ``data/<name>.toml`` -- Per-molecule atomic positions, lattice vectors, center of mass, element symbols (compact notation ``"C6H6"`` supported)

**Output files:**

- ``structure_out/qe.in`` -- Quantum ESPRESSO input file
- ``structure_out/POSCAR``, ``INCAR``, ``KPOINTS`` -- VASP input file set
- ``structure_out/structure.cif`` -- CIF crystal structure file
- ``voronoi_pngs/*.png`` -- Per-site Voronoi cell overlay visualizations

GUI
---

A Tkinter-based GUI extending ASE-GUI supports both command-line and interactive workflows. Key features include:

- CIF file loading and 3D visualization
- Interactive molecule addition, deletion, and substitution via dialogs
- Cell parameter editing and supercell construction
- One-click optimization, rendering, and export
- Atom selection, focus display, and measurements (distance, angle, dihedral)
- Per-element and per-molecule-group display, color, and size control
