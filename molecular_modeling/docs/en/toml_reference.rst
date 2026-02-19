TOML Configuration Reference
============================

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
     - Global random seed. Fallback when optimize/render seeds are not set.
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

Substitution Tables
~~~~~~~~~~~~~~~~~~~

Four types of substitution tables can be defined using ``[[table_name]]`` syntax.

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
