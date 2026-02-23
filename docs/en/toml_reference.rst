Input / Output File Formats
===========================

File Structure Overview
-----------------------

The toolkit uses two types of TOML files:

1. **Input configuration** (``input/input.toml``): Defines supercell lattice, repetition, molecule placement, and substitution rules
2. **Molecule definitions** (``data/<name>.toml``): Defines atomic structure, lattice vectors, and center of mass for each molecule

The relationship between these files is shown below:

.. code-block:: text

   input/input.toml
   ├── molecules = ["benzene", "naphthalene"]
   │       ↓                    ↓
   │   data/benzene.toml    data/naphthalene.toml
   │   (atom coords/lattice) (atom coords/lattice)
   │
   ├── cell = [...]           ← Supercell lattice vectors
   ├── supercell_repeat = [1,1,2]  ← Repetition counts
   ├── center_of_mass = [...]      ← Placement positions
   └── [[substitutions]]           ← Substitution rules (optional)

**Processing flow**:

1. Load ``input.toml`` and read lattice vectors and repetition counts
2. Load ``data/<name>.toml`` for each entry in the ``molecules`` list
3. Generate cells for the number of repetitions specified by ``supercell_repeat``
4. Place molecules at the ``center_of_mass`` positions in each cell
5. Apply substitutions (molecule → global element → molecule atom → individual atom)
6. Build ASE ``Atoms`` object

Input Configuration: input/input.toml
--------------------------------------

Main configuration file controlling the pipeline.

Basic Settings
~~~~~~~~~~~~~~

.. code-block:: toml

   # Random seed (integer)
   # Used as fallback when optimize/render specific seeds are not set
   rng = 1

   # Lattice vectors (3x3 array, unit: Angstrom)
   # Rows define a, b, c vectors
   cell = [
     [ 5.522000,  0.000000,  0.000000 ],   # a-vector
     [ 0.000000,  5.439600,  0.000000 ],   # b-vector
     [-3.501254,  0.000000,  9.339672 ],   # c-vector
   ]

   # Supercell repetition counts [a, b, c directions]
   supercell_repeat = [1, 1, 2]

   # Molecule names (corresponding to data/<name>.toml files)
   molecules = ["benzene", "naphthalene"]

   # Whether to use fractional coordinates (true: fractional, false: Cartesian)
   use_fractional = false

   # Center of mass for each molecule (same order as molecules)
   center_of_mass = [
     [  0.000000,   0.000000,   0.000000 ],
     [ -1.750627,   2.719800,   4.669836 ],
   ]

All Keys
~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Key
     - Type
     - Description
   * - ``rng``
     - int
     - Global random seed. Fallback when ``optimize_rng`` / ``render_rng`` are not set.
   * - ``optimize_rng``
     - int (optional)
     - Random seed for the Voronoi optimization stage. Defaults to ``rng`` if omitted.
   * - ``render_rng``
     - int (optional)
     - Random seed for the visualization stage. Defaults to ``rng`` if omitted.
   * - ``cell``
     - array[3][3]
     - Lattice vectors (rows: a, b, c, unit: Angstrom)
   * - ``supercell_repeat``
     - array[3]
     - Supercell repetition counts [nx, ny, nz]
   * - ``molecules``
     - array[str]
     - Molecule names (corresponding to ``data/<name>.toml``)
   * - ``use_fractional``
     - bool
     - If true, ``center_of_mass`` is interpreted as fractional coordinates
   * - ``center_of_mass``
     - array[N][3]
     - Center of mass coordinates for each molecule (same count as ``molecules``)

.. _cell-index-explanation-en:

cell_index and Supercell Relationship
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When building a supercell with ``supercell_repeat = [nx, ny, nz]``,
each cell is assigned an index ``[ix, iy, iz]``:

- ``ix``: 0 to nx-1
- ``iy``: 0 to ny-1
- ``iz``: 0 to nz-1

For example, with ``supercell_repeat = [1, 1, 2]``:

.. code-block:: text

   Cell [0, 0, 0] — Original unit cell
   Cell [0, 0, 1] — Copy shifted by one c-vector

Each cell contains the molecules defined in ``molecules``, so
``molecule_index = 0`` refers to benzene and ``molecule_index = 1`` refers to naphthalene.

The ``cell_index`` in substitution tables specifies a cell position within this supercell.

Substitution Tables
~~~~~~~~~~~~~~~~~~~

Four types of substitution tables can be defined using ``[[table_name]]`` syntax.
Substitutions are applied in the following order:

1. **Molecule substitutions** (``molecule_substitutions``): Replace entire molecules
2. **Global atom substitutions** (``global_atom_substitutions``): Replace all atoms of an element
3. **Molecule atom substitutions** (``molecule_atom_substitutions``): Replace atom list of a specific molecule
4. **Individual atom substitutions** (``atom_substitutions``): Replace a single atom

molecule_substitutions
^^^^^^^^^^^^^^^^^^^^^^

Replace a molecule at a specific cell position with a different molecule.

.. code-block:: toml

   [[molecule_substitutions]]
   cell_index       = [0, 0, 1]       # Cell position in supercell
   molecule_index   = 1               # Molecule index (0-based)
   new_molecule     = "benzene"       # Replacement molecule name

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Key
     - Type
     - Description
   * - ``cell_index``
     - array[3]
     - Cell position in supercell [ix, iy, iz]
   * - ``molecule_index``
     - int
     - Molecule index within cell (0-based)
   * - ``new_molecule``
     - str
     - Replacement molecule name

global_atom_substitutions
^^^^^^^^^^^^^^^^^^^^^^^^^

Replace all atoms of a specified element throughout the system.

.. code-block:: toml

   [[global_atom_substitutions]]
   old_symbol = "H"      # Original element symbol
   new_symbol = "He"     # Replacement element symbol

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Key
     - Type
     - Description
   * - ``old_symbol``
     - str
     - Original element symbol
   * - ``new_symbol``
     - str
     - Replacement element symbol

molecule_atom_substitutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Replace the atom list of a specific molecule instance with a new list.

.. note::

   This substitution overwrites the entire atom list, so any
   ``global_atom_substitutions`` applied to the same molecule will be overridden.

.. code-block:: toml

   [[molecule_atom_substitutions]]
   cell_index     = [0, 0, 1]
   molecule_index = 1
   new_symbols    = ["Si6O6"]    # Compact notation supported

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Key
     - Type
     - Description
   * - ``cell_index``
     - array[3]
     - Cell position in supercell
   * - ``molecule_index``
     - int
     - Molecule index within cell
   * - ``new_symbols``
     - array[str]
     - New element symbol list (compact notation supported)

atom_substitutions
^^^^^^^^^^^^^^^^^^

Replace a single specific atom.

.. code-block:: toml

   [[atom_substitutions]]
   cell_index     = [0, 0, 0]
   molecule_index = 1
   atom_index     = 2           # Atom index within molecule (0-based)
   new_symbol     = "Fe"

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Key
     - Type
     - Description
   * - ``cell_index``
     - array[3]
     - Cell position in supercell
   * - ``molecule_index``
     - int
     - Molecule index within cell
   * - ``atom_index``
     - int
     - Atom index within molecule (0-based)
   * - ``new_symbol``
     - str
     - Replacement element symbol

Molecule Definition: data/<name>.toml
--------------------------------------

Files defining the structure information for each molecule.
The file name (without ``.toml`` extension) must match the names listed
in ``molecules`` in ``input.toml``.

.. note::

   The ``center_of_mass`` in ``data/*.toml`` is used as **molecule dictionary reference data**.
   The actual placement position in the supercell is determined by ``center_of_mass`` in ``input.toml``.

   The ``lattice_vectors`` in ``data/*.toml`` is the **molecule's own crystal lattice**.
   The overall supercell lattice is defined by ``cell`` in ``input.toml``.
   For mixed crystals (co-crystals of different molecules), these will have different values.

Basic Structure
~~~~~~~~~~~~~~~

.. code-block:: toml

   molecule_name = "benzene"

   # Element symbols: individual list or compact notation
   symbols = ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]
   # or: symbols = ["C6H6"]

   # Coordinate system
   use_fractional = false

   # Atomic coordinates (Nx3 array)
   positions = [
     [-1.113710, -0.097913, -0.811833],
     [ 1.113710,  0.097913,  0.811833],
     # ...
   ]

   # Lattice vectors (3x3 array, unit: Angstrom)
   lattice_vectors = [
     [ 5.522000,  0.000000,  0.000000],
     [ 0.000000,  5.439600,  0.000000],
     [-2.693272,  0.000000,  7.184363],
   ]

   # Number of molecules per unit cell
   num_molecules = 2

   # Center of mass for each molecule
   center_of_mass = [
     [ 0.000000,  0.000000,  0.000000],
     [-1.346636,  2.719800,  3.592182],
   ]

   # Voronoi reference molecule index (default: 0)
   # voronoi_ref_id = 0

All Keys
~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 10 10 55

   * - Key
     - Type
     - Required
     - Description
   * - ``molecule_name``
     - str
     - Yes
     - Molecule identifier
   * - ``symbols``
     - array[str]
     - Yes
     - Element symbols. Compact notation supported (e.g., ``"C6H6"``)
   * - ``use_fractional``
     - bool
     - No
     - If true, positions and center_of_mass are fractional coordinates
   * - ``positions``
     - array[N][3]
     - Yes
     - Atomic coordinates (Cartesian Angstrom or fractional)
   * - ``lattice_vectors``
     - array[3][3]
     - Yes
     - Lattice vectors (rows: a, b, c)
   * - ``num_molecules``
     - int
     - Yes
     - Number of molecules per unit cell
   * - ``center_of_mass``
     - array[K][3]
     - Yes
     - Center of mass for K molecules
   * - ``voronoi_ref_id``
     - int
     - No
     - Reference molecule index for Voronoi cell computation (default: 0)

Compact Notation
~~~~~~~~~~~~~~~~

The ``symbols`` field supports concatenated element-count notation:

.. code-block:: text

   "C6H6"    → ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]
   "C10H8"   → ["C"] x 10 + ["H"] x 8
   "Si6O6"   → ["Si"] x 6 + ["O"] x 6

Coordinate Conversion
~~~~~~~~~~~~~~~~~~~~~

When ``use_fractional = true``, the following conversion is applied internally:

.. code-block:: text

   pos_cart = pos_frac @ lattice_vectors
   pos_frac = pos_cart @ inv(lattice_vectors)

.. note::

   The ``use_fractional`` flags in ``input.toml`` and ``data/*.toml`` are **independent**.
   One can use fractional coordinates while the other uses Cartesian coordinates.

Complete Configuration Examples
-------------------------------

Below are the complete configuration examples for the included benzene-naphthalene
mixed crystal.

input/input.toml Complete Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   # Random seed
   rng = 1

   # Supercell lattice vectors (Angstrom)
   # Lattice for benzene-naphthalene mixed crystal
   cell = [
     [ 5.522000,  0.000000,  0.000000 ],
     [ 0.000000,  5.439600,  0.000000 ],
     [-3.501254,  0.000000,  9.339672 ],
   ]

   # Build a supercell doubled along the c-axis
   supercell_repeat = [1, 1, 2]

   # Use two molecule types
   # References data/benzene.toml and data/naphthalene.toml
   molecules = ["benzene", "naphthalene"]

   # Use Cartesian coordinates
   use_fractional = false

   # Placement position for each molecule (same order as molecules)
   # benzene → origin, naphthalene → offset position
   center_of_mass = [
     [  0.000000,   0.000000,   0.000000 ],
     [ -1.750627,   2.719800,   4.669836 ],
   ]

   # --- Optional substitution rules below ---

   # Replace the 2nd molecule (naphthalene) in cell [0,0,1] with benzene
   [[molecule_substitutions]]
   cell_index       = [0, 0, 1]
   molecule_index   = 1
   new_molecule     = "benzene"

   # Replace all H with He
   [[global_atom_substitutions]]
   old_symbol = "H"
   new_symbol = "He"

   # Replace atom list of the 2nd molecule in cell [0,0,1] with Si6O6
   [[molecule_atom_substitutions]]
   cell_index     = [0, 0, 1]
   molecule_index = 1
   new_symbols    = ["Si6O6"]

   # Replace the 3rd atom of the 2nd molecule in cell [0,0,0] with Fe
   [[atom_substitutions]]
   cell_index     = [0, 0, 0]
   molecule_index = 1
   atom_index     = 2
   new_symbol     = "Fe"

data/benzene.toml Complete Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses the **individual list** format for element symbols.

.. code-block:: toml

   molecule_name = "benzene"

   # 12 atoms listed individually
   symbols = [
     "C", "C", "C", "C", "C", "C",
     "H", "H", "H", "H", "H", "H",
   ]

   use_fractional = false

   # Atomic coordinates of the 1st molecule (Angstrom)
   # Corresponds to symbols in the same order
   positions = [
     # Carbon atoms (6)
     [ -1.113710, -0.097913, -0.811833 ],
     [  1.113710,  0.097913,  0.811833 ],
     [ -0.115564, -1.087920, -0.865716 ],
     [  0.115564,  1.087920,  0.865716 ],
     [ -0.979384,  0.951930,  0.052446 ],
     [  0.979384, -0.951930, -0.052446 ],
     # Hydrogen atoms (6)
     [ -1.861055, -0.158292, -1.362155 ],
     [  1.861055,  0.158292,  1.362155 ],
     [ -0.197159, -1.816826, -1.437591 ],
     [  0.197159,  1.816826,  1.437591 ],
     [ -1.651612,  1.593259,  0.089805 ],
     [  1.651612, -1.593259, -0.089805 ],
   ]

   # Benzene crystal lattice vectors
   lattice_vectors = [
     [ 5.522000,  0.000000,  0.000000],
     [ 0.000000,  5.439600,  0.000000],
     [-2.693272,  0.000000,  7.184363],
   ]

   # 2 molecules per unit cell
   num_molecules = 2

   # Center of mass for each molecule
   center_of_mass = [
     [ 0.000000,  0.000000,  0.000000],
     [-1.346636,  2.719800,  3.592182],
   ]

data/naphthalene.toml Complete Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses **compact notation** for element symbols.

.. code-block:: toml

   molecule_name = "naphthalene"

   # Compact notation: C10H8 → C x 10 + H x 8 = 18 atoms
   symbols = ["C10H8"]

   use_fractional = false

   # Atomic coordinates of the 1st molecule (Angstrom)
   # Corresponds to expanded symbols (C first, then H)
   positions = [
     # Carbon atoms (10)
     [ -0.873883,  0.110730,  2.365899 ],
     [  0.873883, -0.110730, -2.365899 ],
     [ -0.125794,  0.966533,  1.590346 ],
     [  0.125794, -0.966533, -1.590346 ],
     [  0.213579,  0.617262,  0.251858 ],
     [ -0.213579, -0.617262, -0.251858 ],
     [  0.981797,  1.483077, -0.575778 ],
     [ -0.981797, -1.483077,  0.575778 ],
     [ -1.309714, -1.131450,  1.848626 ],
     [  1.309714,  1.131450, -1.848626 ],
     # Hydrogen atoms (8)
     [ -1.094109,  0.351038,  3.237058 ],
     [  1.094109, -0.351038, -3.237058 ],
     [  0.159973,  1.778750,  1.942805 ],
     [ -0.159973, -1.778750, -1.942805 ],
     [  1.263671,  2.303540, -0.240443 ],
     [ -1.263671, -2.303540,  0.240443 ],
     [ -1.818179, -1.705126,  2.375888 ],
     [  1.818179,  1.705126, -2.375888 ],
   ]

   # Naphthalene crystal lattice vectors (different from benzene)
   lattice_vectors = [
     [ 8.034800,  0.000000,  0.000000],
     [ 0.000000,  5.889900,  0.000000],
     [-4.738553,  0.000000,  7.134798],
   ]

   # 2 molecules per unit cell
   num_molecules = 2

   # Center of mass for each molecule
   center_of_mass = [
     [ 0.000000,  0.000000,  0.000000],
     [-1.346636,  2.719800,  3.592182],
   ]

Creating New Molecule Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To add a custom molecule, create a ``data/<name>.toml`` file with the following steps:

1. Set ``molecule_name`` to the molecule name
2. Define element symbols in ``symbols`` (individual list or compact notation)
3. Write the **1st molecule's** atomic coordinates in ``positions`` (center of mass near the origin)
4. Set the molecule's crystal lattice in ``lattice_vectors``
5. Set the number of molecules per unit cell in ``num_molecules``
6. Write the center of mass for each molecule in ``center_of_mass``

.. tip::

   You can use **MoleculeTools** → **Save Selection as Molecule...** in the GUI
   to save selected atoms from the 3D view as a ``data/*.toml`` file.
