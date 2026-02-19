Molecular Modeling Toolkit Documentation
=========================================

Molecular Modeling Toolkit is a molecular modeling toolkit built on
ASE (Atomic Simulation Environment). It builds supercells from TOML-based
molecular definitions, optimizes molecular orientations via Voronoi cell
IoU maximization, and exports to Quantum ESPRESSO / VASP / CIF formats.

Key Features
------------

- **TOML-based molecular definitions**: Define molecular structures with intuitive configuration files
- **Supercell construction**: Generate supercells with arbitrary repetition counts
- **Voronoi optimization**: Grid search over Euler angles and Nelder-Mead refinement for molecular orientation optimization
- **4 types of substitutions**: Global element, molecule, molecule-atom, and individual atom substitutions
- **Multi-format export**: Batch export to QE, VASP, CIF formats
- **GUI**: Interactive environment extending ASE-GUI

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   tutorial
   gui_guide

.. toctree::
   :maxdepth: 2
   :caption: Reference

   cli_reference
   toml_reference
   architecture
   api/index

Other Languages
---------------

- `日本語 <../ja/index.html>`_

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
