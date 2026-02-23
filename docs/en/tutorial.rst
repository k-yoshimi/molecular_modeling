Tutorial
========

This tutorial walks through all pipeline steps using a benzene-naphthalene
mixed crystal as an example.

Step 1: Understanding Molecule Definitions
------------------------------------------

Let's examine the contents of ``data/benzene.toml``:

.. code-block:: toml

   molecule_name = "benzene"

   symbols = [
     "C", "C", "C", "C", "C", "C",
     "H", "H", "H", "H", "H", "H",
   ]

   use_fractional = false

   positions = [
     [ -1.113710, -0.097913, -0.811833 ],
     [  1.113710,  0.097913,  0.811833 ],
     # ... (coordinates for all 12 atoms)
   ]

   lattice_vectors = [
     [ 5.522000,  0.000000,  0.000000],
     [ 0.000000,  5.439600,  0.000000],
     [-2.693272,  0.000000,  7.184363],
   ]

   num_molecules = 2

   center_of_mass = [
     [ 0.000000,  0.000000,  0.000000],
     [-1.346636,  2.719800,  3.592182],
   ]

Field descriptions:

- ``molecule_name``: Molecule identifier (referenced by ``molecules`` in ``input.toml``)
- ``symbols``: List of element symbols. Compact notation like ``"C6H6"`` is also supported
- ``use_fractional``: If ``true``, coordinates are interpreted as fractional
- ``positions``: Atomic Cartesian coordinates (Angstrom). Fractional if ``use_fractional = true``
- ``lattice_vectors``: Lattice vectors (rows define a, b, c)
- ``num_molecules``: Number of molecules per unit cell
- ``center_of_mass``: Center of mass position for each molecule

See :doc:`toml_reference` for full key descriptions.

Naphthalene (``data/naphthalene.toml``) follows the same structure,
using compact notation ``"C10H8"`` for element symbols.

Step 2: Configuring the Input File
-----------------------------------

``input/input.toml`` specifies how to build the supercell:

.. code-block:: toml

   # Random seed
   rng = 1

   # Lattice vectors (Angstrom)
   cell = [
     [ 5.522000,  0.000000,  0.000000 ],
     [ 0.000000,  5.439600,  0.000000 ],
     [-3.501254,  0.000000,  9.339672 ],
   ]

   # Supercell repetition counts [a, b, c directions]
   supercell_repeat = [1, 1, 2]

   # Molecule names to use
   molecules = ["benzene", "naphthalene"]

   # Whether to use fractional coordinates
   use_fractional = false

   # Center of mass for each molecule
   center_of_mass = [
     [  0.000000,   0.000000,   0.000000 ],
     [ -1.750627,   2.719800,   4.669836 ],
   ]

``supercell_repeat = [1, 1, 2]`` creates a supercell doubled along the c-axis.
Two molecule types (benzene and naphthalene) are placed at their respective
center-of-mass positions.

Step 3: Running the Pipeline
----------------------------

.. code-block:: bash

   python run.py

The pipeline processes in the following order:

1. **Build**: Read TOML files and construct ASE ``Atoms`` object
2. **Group**: Group atoms by ``mol_id``
3. **Optimize**: Maximize Voronoi IoU via Euler angle grid search
4. **Check**: Verify Voronoi volume tiling consistency
5. **Render**: Generate 3D Voronoi cell overlay visualizations
6. **Export**: Write QE, VASP, CIF files

Step 4: Checking Output Files
------------------------------

Structure Files
~~~~~~~~~~~~~~~

The following files are generated in ``structure_out/``:

**qe.in** (Quantum ESPRESSO):

.. code-block:: text

   &CONTROL
     calculation = 'scf'
     prefix = 'molecule'
     outdir = './tmp'
   /
   &SYSTEM
     ecutwfc = 40
     ecutrho = 320
   /
   ...

**POSCAR** (VASP): Atomic positions in fractional coordinates.

**structure.cif**: Output in universal CIF format.

Step 5: Reading Voronoi Visualizations
--------------------------------------

PNG images in ``voronoi_pngs/`` show the Voronoi cell overlay for each site:

- **Blue translucent polyhedron**: System Voronoi cell (actual cell partition)
- **Red translucent polyhedron**: Dictionary Voronoi cell (optimization target)
- Higher IoU (Intersection over Union) values indicate better molecular orientation reproduction

Step 6: Interactive GUI Operations
-----------------------------------

See :doc:`gui_guide` for the full GUI reference.

Loading a CIF file
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python vis_cif.py structure_out/structure.cif

Or from the GUI menu:

**MoleculeTools** → **Load CIF...**

Adding Molecules
~~~~~~~~~~~~~~~~

1. Select **MoleculeTools** → **Add Molecule...**
2. Choose a molecule from the dropdown (e.g., benzene)
3. Enter placement position coordinates
4. Click **OK**

Substitution Operations
~~~~~~~~~~~~~~~~~~~~~~~

The GUI supports 4 types of substitutions:

1. **Global Element Substitution**: Replace all atoms of a specified element
2. **Molecule Substitution**: Replace a molecule at a specific site
3. **Molecule Atom Substitution**: Change the atom list within a molecule
4. **Individual Atom Substitution**: Replace a single atom

Export
~~~~~~

**Export** → **Export All (QE, VASP, CIF)...** exports to all formats at once.

Step 7: Substitution Details
-----------------------------

``input.toml`` supports 4 types of substitution tables.

Molecule Substitutions
~~~~~~~~~~~~~~~~~~~~~~

Replace a molecule at a specific cell position with a different molecule:

.. code-block:: toml

   [[molecule_substitutions]]
   cell_index       = [0, 0, 1]    # Cell position in supercell
   molecule_index   = 1            # Molecule index (0-based)
   new_molecule     = "benzene"    # Replacement molecule name

Global Atom Substitutions
~~~~~~~~~~~~~~~~~~~~~~~~~

Replace all atoms of a specified element throughout the system:

.. code-block:: toml

   [[global_atom_substitutions]]
   old_symbol = "H"     # Original element
   new_symbol = "He"    # Replacement element

Molecule Atom Substitutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Replace the entire atom list of a specific molecule instance:

.. code-block:: toml

   [[molecule_atom_substitutions]]
   cell_index     = [0, 0, 1]
   molecule_index = 1
   new_symbols    = ["Si6O6"]    # Compact notation supported

Individual Atom Substitutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Replace a single atom:

.. code-block:: toml

   [[atom_substitutions]]
   cell_index     = [0, 0, 0]
   molecule_index = 1
   atom_index     = 2           # Atom index within molecule (0-based)
   new_symbol     = "Fe"
