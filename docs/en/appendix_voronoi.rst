Voronoi Optimization Details
============================

This appendix explains the Voronoi cell IoU optimization algorithm,
a core feature of the toolkit, in detail.

Problem Setting
---------------

When modeling organic molecular crystals and molecular solids, it is necessary
to appropriately determine the **orientation** (direction) of molecules placed
in a unit cell.
This toolkit uses the centers of mass of molecules as generator points for
Voronoi tessellation and automatically determines optimal orientations by
maximizing the **geometric overlap** between the Voronoi cells of the system
and those of the dictionary (reference crystal).

Here, "system" refers to the supercell actually constructed, and "dictionary"
refers to the reference single-crystal structure.

What Are Voronoi Cells?
-----------------------

When Voronoi tessellation is performed using the center of mass of each
molecular site in the crystal as a generator point, space is divided into
regions that are "closest" to each site.
These regions are Voronoi cells (Voronoi polyhedra).

.. code-block:: text

   Properties of Voronoi cells:
   - One convex polyhedron corresponds to each generator point
   - Sum of all cell volumes = unit cell volume (space is filled without gaps)
   - Cell shapes are determined by relative positions between molecules

This toolkit performs Voronoi tessellation with periodic boundary conditions.
Specifically, it generates images for ±1 cell (27 periodic images) and then
computes Voronoi cell vertices using SciPy's halfspace intersection method
(HalfspaceIntersection).

IoU (Intersection over Union)
-----------------------------

IoU is a metric that represents the degree of overlap between two solids,
taking values from 0 to 1.

.. code-block:: text

                Intersection Volume (Overlap)
   IoU = ──────────────────────────────────────────
          V_sys + V_dict - Intersection Volume

   IoU = 1: Two Voronoi cells match perfectly
   IoU = 0: Two Voronoi cells do not overlap at all

Where:

- **V_sys**: Voronoi cell volume of the system (supercell)
- **V_dict**: Voronoi cell volume of the dictionary (reference crystal, after scaling)
- **Intersection volume**: Volume of the overlapping region of the two cells

The intersection volume is computed by performing Voronoi tessellation on the
combined generator points of both the system and dictionary, and obtaining the
volume of the combined cell for the reference point.
This corresponds to the intersection volume of V_sys and V_dict.

Per-Site Optimization Procedure
-------------------------------

Optimization is performed independently for each molecular site.
Each site undergoes the following 5 steps.

Step 1: System Voronoi Cell Volume Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using all molecular centers of mass in the supercell as generator points,
the Voronoi cell volume V_sys at the target site's center of mass is computed.

.. code-block:: text

   Input: Supercell lattice A_sys, all molecular center-of-mass coordinates
   Output: V_sys (angstrom^3)

The system's Voronoi cell does not change during optimization (because
center-of-mass positions remain unchanged even when orientations change).

Step 2: Dictionary Voronoi Cell Volume Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the lattice vectors and centers of mass from the dictionary
(reference crystal), the dictionary-side Voronoi cell volume V_dict is computed.

.. code-block:: text

   Input: Dictionary lattice A_dict, dictionary center-of-mass coordinates
   Output: V_dict (angstrom^3)

Step 3: Volume Scaling
~~~~~~~~~~~~~~~~~~~~~~

When the system and dictionary Voronoi cell volumes differ (e.g., in mixed
crystals), uniform scaling is applied to the dictionary's lattice and centers
of mass to match the volumes.

.. code-block:: text

   Scale factor: s = (V_sys / V_dict)^(1/3)

This scaling makes the dictionary Voronoi cell volume equal to that of the
system. The cube root is used because in 3D uniform scaling, the volume
scales as s^3.

Step 4: Euler Angle Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dictionary Voronoi cell is rotated to find Euler angles (alpha, beta, gamma)
that maximize IoU with the system Voronoi cell.

**Coarse Grid Search (Phase 1)**:

.. code-block:: text

   alpha, beta, gamma each divided into 13 steps over [-180 deg, 180 deg]
   -> 13 x 13 x 13 = 2,197 combinations evaluated
   -> The angles with the highest IoU are recorded as the initial solution

**Nelder-Mead Refinement (Phase 2)**:

.. code-block:: text

   Refine using the Nelder-Mead method with the best grid search result as initial value
   - Objective function: f(alpha, beta, gamma) = 1 - IoU (minimization)
   - Angle tolerance: 0.001 deg
   - IoU tolerance: 0.001
   - Maximum iterations: 2,000
   - Maximum function evaluations: 10,000

The Nelder-Mead method is a derivative-free direct search method, well-suited
for objective functions like IoU where derivatives are difficult to compute.

Euler angle rotation convention:

.. code-block:: text

   Rotation order: XYZ (extrinsic)
   - alpha: rotation about X axis
   - beta: rotation about Y axis
   - gamma: rotation about Z axis

   Rotation matrix R = R_z(gamma) x R_y(beta) x R_x(alpha)

Step 5: Application to Atomic Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimized rotation angles and scale are applied to all atomic coordinates
of the target site.

.. code-block:: text

   New coordinates for each atom:
     pos_new = p_com + s x R x (pos_old - p_com)

   p_com: Center of mass of the site
   s:     Scale factor (Step 3)
   R:     Optimal rotation matrix (Step 4)

Since rotation and scaling are performed about the center of mass,
the molecular center-of-mass positions remain unchanged.

Dictionary Voronoi Cell Construction
-------------------------------------

The dictionary Voronoi cell is constructed by the following procedure:

1. Translate dictionary centers of mass so that the reference molecule (``voronoi_ref_id``, default 0) is at the origin
2. Rotate by Euler angles (alpha, beta, gamma)
3. Uniformly scale by factor s
4. Translate to the center-of-mass position of the target site in the system
5. Apply the same rotation and scale to the dictionary lattice and generate images for ±1 cell
6. Perform combined Voronoi tessellation from all generator points of both system and dictionary

This procedure places the dictionary Voronoi cell in the system's coordinate
space, enabling IoU computation.

Volume Tiling Verification
--------------------------

After optimization, the toolkit verifies that the total Voronoi cell volumes
match the supercell volume.

.. code-block:: text

   Verification condition:
     |sum(V_i) - V_cell| <= abs_tol  or  |sum(V_i) - V_cell| / V_cell <= rel_tol

   Default tolerances:
     abs_tol = 1e-6 angstrom^3
     rel_tol = 1e-3 (0.1%)

If this verification fails, there may be a problem with the Voronoi tessellation
(e.g., molecules placed too closely, extremely distorted lattice).

Interpreting Optimization Results
---------------------------------

Optimization results are stored in ``system.info`` and can be viewed in the GUI
via **Pipeline** > **Show Optimization Results**.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Output Item
     - Meaning
   * - IoU
     - Degree of overlap between system and dictionary Voronoi cells (0--1). Closer to 1 means the dictionary molecular orientation is better reproduced in the system
   * - alpha, beta, gamma
     - Optimized Euler angles (degrees). Amount of rotation relative to the dictionary molecular orientation
   * - scale
     - Volume scaling factor. 1.0 means the system and dictionary Voronoi volumes are equal. Greater than 1 means the system site is larger than the dictionary

**IoU Guidelines**:

- **0.9 or above**: The dictionary orientation is very well reproduced in the system
- **0.7 to 0.9**: Generally good reproduction
- **Below 0.7**: Significant differences between the system and dictionary environments (e.g., effects of different molecules in mixed crystals)

Visualization (Voronoi Overlay PNG)
------------------------------------

Optimization results are visualized as 3D images:

- **Blue semi-transparent polyhedron**: System Voronoi cell (actual spatial partitioning)
- **Red semi-transparent polyhedron**: Dictionary Voronoi cell (reference cell rotated and scaled by optimization)
- **Orange semi-transparent polyhedron**: Intersection of the two cells

The higher the IoU, the more the blue and red polyhedra overlap, becoming nearly
identical in shape.

Libraries Used
--------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Library
     - Purpose
   * - SciPy (``scipy.spatial``)
     - Voronoi tessellation, halfspace intersection, convex hull computation, Nelder-Mead optimization
   * - SciPy (``scipy.spatial.transform``)
     - Rotation matrix generation from Euler angles
   * - trimesh
     - Polyhedral mesh construction, volume computation, intersection calculations
   * - manifold3d
     - Polyhedral boolean operations (intersection, union)
   * - NumPy
     - Numerical array operations, linear algebra
