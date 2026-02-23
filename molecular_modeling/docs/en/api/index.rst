API Reference
=============

This section provides auto-generated API documentation for the main modules.

Module Organization
-------------------

The toolkit consists of the following module groups:

- **Pipeline**: ``main`` is the entry point. Processing flows through ``build_system_from_toml`` → ``build_molecule_groups`` → ``optimize_voronoi_overlap`` → ``render_voronoi_overlays`` → ``write_structure``
- **Geometry**: ``voronoi_geom`` handles Voronoi cell construction, volume calculation, and IoU evaluation. The ``Polyhedron`` class is the core abstraction
- **GUI**: ``gui.controller`` bridges views (dialogs) and the model (ASE Atoms). ``gui.molecule_registry`` manages dynamic loading of ``data/*.toml``

Basic programmatic usage:

.. code-block:: python

   from molecular_modeling.main import main

   # Read input/input.toml and run the full pipeline
   main()

Main Pipeline
-------------

main
~~~~

.. automodule:: molecular_modeling.main
   :members:
   :undoc-members:
   :no-index:

System Construction
-------------------

build_system_from_toml
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.build_system_from_toml
   :members:
   :undoc-members:
   :no-index:

build_molecule_groups
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.build_molecule_groups
   :members:
   :undoc-members:
   :no-index:

Voronoi Geometry
----------------

voronoi_geom
~~~~~~~~~~~~

.. automodule:: molecular_modeling.voronoi_geom
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

Optimization
------------

optimize_voronoi_overlap
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.optimize_voronoi_overlap
   :members:
   :undoc-members:
   :no-index:

Visualization
-------------

render_voronoi_overlays
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.render_voronoi_overlays
   :members:
   :undoc-members:
   :no-index:

Output
------

write_structure
~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.write_structure
   :members:
   :undoc-members:
   :no-index:

GUI
---

controller
~~~~~~~~~~

.. automodule:: molecular_modeling.gui.controller
   :members:
   :undoc-members:
   :no-index:

molecule_registry
~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.gui.molecule_registry
   :members:
   :undoc-members:
   :no-index:
