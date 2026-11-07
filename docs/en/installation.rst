Installation
============

Requirements
------------

- Python 3.11 or later (tested with 3.13)
- pip or uv package manager

Dependencies
~~~~~~~~~~~~

Key dependencies (managed in ``pyproject.toml``):

- **ASE** (3.25.0) — Atomic Simulation Environment
- **NumPy** (2.3.1) — Numerical computing
- **SciPy** (1.16.0) — Voronoi computation, optimization
- **matplotlib** (3.10.3) — Visualization
- **trimesh** (4.7.1) — Mesh operations
- **manifold3d** (3.2.1) — Manifold reconstruction
- **Pillow** (11.3.0) — Image processing

Optional dependencies for crystal structure prediction (the ``csp`` extra, see :doc:`csp_guide`):

- **PyXtal** (>= 1.1.5) — Random symmetric molecular crystals
- **tblite** (>= 0.8.0) — GFN2-xTB energies and forces for relaxation
- **pymatgen** (>= 2025.1.9) — CIF reading and writing; brings in **spglib** (space-group detection)

For relaxation with the machine-learned potential OrbMol-v2, the ``csp-ml`` extra adds (Python 3.12 or
later, about 1 GB):

- **orb-models** (>= 0.7, < 0.8) — OrbMol-v2 potential; brings in **PyTorch** and **NVIDIA Warp**

Installation Steps
------------------

1. Clone the repository
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone <repository-url>
   cd molecular_modeling   # リポジトリルート

2. Create a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python -m venv .venv
   source .venv/bin/activate    # Linux / macOS
   .venv\Scripts\activate       # Windows

3. Install the package
~~~~~~~~~~~~~~~~~~~~~~

Using pip:

.. code-block:: bash

   pip install -e .

Using uv:

.. code-block:: bash

   uv sync

For crystal structure prediction, add the ``csp`` extra:

.. code-block:: bash

   pip install -e ".[csp]"     # or: uv sync --extra csp
   pip install -e ".[csp-ml]"  # csp plus OrbMol-v2 (or: uv sync --extra csp-ml)

4. Verify installation
~~~~~~~~~~~~~~~~~~~~~~

Run the test suite to confirm everything is installed correctly:

.. code-block:: bash

   pytest

All tests should pass.

