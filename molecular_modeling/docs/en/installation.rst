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

Installation Steps
------------------

1. Clone the repository
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone <repository-url>
   cd molecular_modeling

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

4. Verify installation
~~~~~~~~~~~~~~~~~~~~~~

Run the test suite to confirm everything is installed correctly:

.. code-block:: bash

   pytest

All 48 tests should pass.

Conda environment
~~~~~~~~~~~~~~~~~

You can also set up a Conda environment using ``environment.yml``:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate molecular_modeling

See ``CONDA_SETUP.md`` for details.
