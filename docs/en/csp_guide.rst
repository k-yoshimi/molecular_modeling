Crystal Structure Prediction (CSP) Guide
========================================

This guide describes how to build candidate molecular crystals from the structure of an isolated
molecule alone, relax them with GFN2-xTB, rank them by energy, and compare them with experimental CIF
files. It is written with OS-2SMe (molecular crystals with O/S exchange, "flip", between molecules)
in mind, but it works for ordinary organic molecular crystals as well.

.. contents:: Contents
   :local:
   :depth: 2

Installation
------------

.. code-block:: bash

   source .venv/bin/activate
   pip install -e ".[csp]"

The ``csp`` extra installs PyXtal, tblite (GFN2-xTB) and pymatgen. They are not part of the core
dependencies.

Workflow
--------

.. code-block:: bash

   molecular-modeling-csp input/csp.toml     # predict (default config: input/csp.toml)
   python -m molecular_modeling.csp.validate csp_out exp1.cif exp2.cif   # compare with experiment

A run goes through these stages:

1. **Generation**: PyXtal builds crystals for each space group (random). A share of them gets its
   molecular orientations re-chosen from Voronoi cells (Voronoi generator). If seed CIFs are given,
   candidates are also derived from them (seed generator).
2. **Overlap check and prescreening**: candidates with atoms too close together, broken molecules or a
   packing fraction outside the range are dropped. The rest are ranked by a cheap rigid-molecule
   intermolecular energy (UFF Lennard-Jones) and the best ``n_prescreen_keep`` are kept.
3. **GFN2-xTB relaxation**: atoms and cell are relaxed (in parallel). Candidates whose molecules
   deform too much from the input are dropped.
4. **Deduplication and ranking**: identical structures are merged and ranked by energy per molecule.
5. **Output**: a ranking table and CIF files for the top 50.

Atomic coordinates of the molecule are never scaled. The only operations on a molecule are rigid
rotation and translation, apart from the GFN2-xTB relaxation.

Input file (input/csp.toml)
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 20 58

   * - Key
     - Default
     - Meaning
   * - ``molecule``
     - (required)
     - Molecule file: ``data/<name>.toml`` or ``.xyz``. Cartesian coordinates of one isolated molecule
       whose bonds form a single connected piece (a fragmented molecule is an error)
   * - ``mode``
     - ``"blind"``
     - ``"blind"`` or ``"informed"``; see "Two modes" below
   * - ``space_groups``
     - 6 groups (below)
     - Space group numbers to try. Required only in ``informed`` mode; in ``blind`` mode it is
       replaced by the default list
   * - ``z``
     - omitted
     - Z per space group (e.g. ``z = { 14 = 4 }``). Default: general-position multiplicity x ``z_prime``
   * - ``z_prime``
     - ``1``
     - Molecules per asymmetric unit. Use ``informed`` mode for Z' = 2 (the cost grows)
   * - ``n_candidates_per_sg``
     - ``150``
     - Candidates generated per space group
   * - ``generators``
     - ``{random=0.6, voronoi=0.4}``
     - Relative weights of the generators; ``seed`` is added automatically when seeds are given
   * - ``seed_structures``
     - omitted
     - List of seed CIFs (``[[seed_structures]]``), see below
   * - ``packing_range``
     - ``[0.4, 0.8]``
     - Allowed packing fraction of the generated, unrelaxed cells
       (molecules x molecular vdW volume / cell volume)
   * - ``overlap_factor``
     - ``0.7``
     - Reject candidates with an interatomic distance below this factor x the sum of vdW radii
   * - ``n_prescreen_keep``
     - ``32``
     - Candidates kept after prescreening, i.e. relaxed with GFN2-xTB
   * - ``relax_steps``
     - ``60``
     - Maximum relaxation steps per candidate
   * - ``relax_fmax``
     - ``0.05``
     - Force convergence threshold (eV/A)
   * - ``n_procs``
     - number of CPU cores
     - Number of worker processes
   * - ``max_mol_rmsd``
     - ``0.5``
     - Reject candidates in which a molecule deviates more than this heavy-atom RMSD (A) from the
       input molecule after relaxation
   * - ``time_budget_min``
     - ``60``
     - Wall-clock target in minutes; see "Timings and defaults"
   * - ``rng``
     - ``1``
     - Random seed (same value, same candidates)
   * - ``out_dir``
     - ``"csp_out"``
     - Output directory

The default space groups are ``[14, 2, 19, 15, 4, 61]`` (P2\ :sub:`1`/c, P-1,
P2\ :sub:`1`\ 2\ :sub:`1`\ 2\ :sub:`1`, C2/c, P2\ :sub:`1`, Pbca). Unknown keys in the file are an error.

Seed CIFs (seed_structures)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   [[seed_structures]]
   cif = "parent.cif"
   substitutions = [["O", "S", 1]]   # (old element, new element, count per molecule)
   volume_scale = "auto"             # or a number

The parent CIF must be **free of disorder** (no partial occupancies): order it, for example in VESTA,
and save it as P1. A CIF with partial occupancy is rejected.

Two modes
---------

``blind``
   A prediction that uses no experimental information: the six default space groups, Z equal to the
   general-position multiplicity (Z' = 1). Use it to test whether the experimental structure can be
   reproduced from the molecule alone.

``informed``
   You give the space groups, ``z`` and ``z_prime`` known from experiment and, optionally,
   ``seed_structures`` (e.g. a known crystal with O replaced by S) as extra starting points. The search
   concentrates on the region around the answer.

How Voronoi is used
-------------------

The earlier Voronoi pipeline (``run.py``) maximised the IoU against Voronoi cells of a centroid lattice
taken from a parent structure. In CSP the Voronoi cell is used to choose the **initial orientation**:

- Take the molecular centroids of a PyXtal crystal and build the Voronoi cell of each molecule.
- Find the orientation that maximises the IoU between the cell and **the molecule's own shape** (the
  convex hull of its atoms inflated by vdW radii). For the comparison only, the shape is scaled to the
  cell volume. Atomic coordinates are not scaled; the molecule is only rotated as a rigid body.
- For elongated, nearly symmetric molecules head and tail are hard to tell apart, so both the best
  orientation and that orientation turned by 180 degrees about the long axis become candidates.
- The space-group symmetry is kept: the orientation is chosen for the asymmetric-unit molecules, and the
  others follow from the symmetry operations.

The ``origin`` column of the validation table shows whether the closest candidate came from the random,
Voronoi or seed generator.

Reading the output
------------------

``csp_out/`` contains:

``candidates.csv``
   The ranking, lowest energy first.

   .. list-table::
      :header-rows: 1
      :widths: 28 72

      * - Column
        - Meaning
      * - ``rank``
        - Rank (1 = most stable)
      * - ``id``
        - Candidate identifier (``<origin>-<space group>-<serial>``)
      * - ``energy_per_mol_eV``
        - GFN2-xTB energy per molecule (eV)
      * - ``rel_energy_kJmol``
        - Difference from rank 1 (kJ/mol per molecule)
      * - ``density_gcm3``
        - Density after relaxation (g/cm3)
      * - ``space_group_in``
        - Space group requested at generation
      * - ``space_group_out``
        - Space group re-detected by spglib after relaxation (0 = detection failed)
      * - ``origin``
        - Generator (``random`` / ``voronoi`` / seed name)
      * - ``max_mol_rmsd``
        - Deformation of the molecules after relaxation (heavy-atom RMSD, A)
      * - ``status``
        - ``relaxed`` (accepted)

``001_<id>.cif`` ... ``050_<id>.cif``
   Structures of the top 50.

``run.log``
   How many candidates survived each stage and how long it took. Example (benzene, P-1)::

      generated 200 candidates {'random': 120, 'voronoi': 80} in 23.0 s
      prescreen {'generated': 200, 'no_overlap': 85, 'packing': 85, 'kept': 85} in 3.4 s
      relax: 2 relaxed, 0 rejected_deformed, 0 failed in 206.1 s
      dedupe: 2 -> 2 in 0.0 s

   ``rejected_deformed`` counts molecules that deformed beyond ``max_mol_rmsd``; ``failed`` counts
   failed calculations. A line starting with ``time budget:`` means the number of relaxations was reduced
   to fit the time target.

Validation against experiment
-----------------------------

.. code-block:: bash

   python -m molecular_modeling.csp.validate csp_out exp_OS-i.cif exp_OS-ii.cif
   # optional: --n-top 20  --rmsd-cut 0.8  --flip O,S  (--no-flip skips the flip check)

The report is printed and written to ``csp_out/validation.md``. Experimental CIFs must be
**disorder-free** (full occupancy).

The comparison follows the COMPACK approach: a cluster of 15 molecules is taken from the experimental
structure and from each candidate and overlaid, giving RMSD\ :sub:`15`. A candidate **matches** if
RMSD\ :sub:`15` < 0.8 A with all 15 molecules paired.

``success``
   yes if a match exists within the top ``n_top`` (default 20).
``best match rank``
   Rank of the best-ranked matching candidate. It is listed even when it is below 20, which separates
   "generated but ranked too low" from "never generated".
``mirror``
   Whether the match used the mirror image of the candidate.
``flip match rank``
   Best rank of a candidate that matches only when O and S are treated as one element: **the same
   packing with O and S exchanged on the molecules**. It is never counted as success. It reads
   ``skipped`` when the experimental structure lacks O or S.
``closest candidates``
   If nothing matches, the ten closest candidates (matched molecules, RMSD) are shown.
``Not comparable``
   Candidates that could not be compared, with the reason, for example: the cluster molecules have
   different compositions, there are no molecules to compare, or the candidate CIF could not be read.

Timings and defaults
--------------------

Measured on an Apple Silicon Mac (8 cores). Naphthalene (18 atoms, 10 heavy atoms) stands in for
OS-2SMe, in space group P2\ :sub:`1`/c with Z = 4 (72 atoms). OS-2SMe has somewhat more heavy atoms, so
expect it to be slower.

.. list-table::
   :header-rows: 1
   :widths: 60 40

   * - Item
     - Measured
   * - PyXtal generation (per candidate)
     - about 0.04 s
   * - Rigid LJ energy (per candidate)
     - about 0.08 s
   * - Voronoi orientation (per skeleton, gives 2 candidates)
     - about 0.7 s
   * - GFN2-xTB relaxation step (one process running)
     - about 4.5 s
   * - Same, 8 processes at once
     - about 10.9 s (1.4 s effective per candidate)
   * - (reference) benzene P-1 Z = 2 (24 atoms), 150 steps
     - about 17 s

GFN2-xTB cost grows nearly with the cube of the atom count, so 72 atoms cost about 40 times more than 24.
Relaxation takes almost all of the time.

**Why these defaults** (fit within 90 % of ``time_budget_min = 60``, i.e. 54 min):

- Generation and prescreening: 6 space groups x 150 = about 900 candidates (about 40 s), Voronoi
  orientation (30 skeletons per space group x 0.7 s, about 2 min in total), LJ ranking (about 20 s):
  about 4 min.
- One relaxation round (8 candidates at once): 60 steps x 10.9 s = about 11 min.
- The remaining 54 - 4 = 50 min allows 4 rounds, about 44 min, i.e. 8 x 4 = **32 candidates**
  (``n_prescreen_keep = 32``).
- Total about 48 min, inside the target.

With ``relax_steps = 150`` one round takes about 27 min and only one round (8 candidates) fits into an
hour. A relaxation that stops before full convergence is still adequate for narrowing down candidates,
so 60 is used. For production runs, raise ``time_budget_min``, reduce the number of space groups
(``informed`` mode), or adjust ``n_procs`` to the number of cores.

The time target works as follows. The first round (``n_procs`` candidates) is timed, and only as many whole
rounds as fit in the remaining time are run. **The first round always runs**, so a ``time_budget_min``
shorter than one round is exceeded.

About the packing range
~~~~~~~~~~~~~~~~~~~~~~~

The unrelaxed cells PyXtal produces are much looser than real crystals (0.6-0.8): about 0.46 for
benzene and 0.52 for naphthalene, almost the same for every space group and random seed. The default
``packing_range`` is therefore ``[0.4, 0.8]``. Denser starting cells are possible by lowering PyXtal's
volume factor from 1.1 to 0.8, but more candidates then fail the overlap check (for naphthalene, the
fraction without overlaps drops from about a fifth to nearly zero), so the default factor was kept
and the range was widened.

Limitations
-----------

- **The molecule is a rigid input conformer.** Torsions are not varied at generation. Relaxation moves
  atoms a little, but a candidate whose molecule deviates beyond ``max_mol_rmsd`` is dropped.
- **GFN2-xTB ranking may not separate close polymorphs**, especially for dipolar molecules, where the
  energy differences are smaller than the method error. A ``best match rank`` below 20 does not mean the
  structure was not generated.
- **Flip disorder is not modelled.** Mixtures of O/S-exchanged molecules (partially occupied sites) are
  not reproduced; they are only reported separately as a "flip match".
- Seed and experimental CIFs must be **free of disorder**.
- Relaxations stop at ``relax_steps`` and may not be fully converged; energies refer to that state.
- Re-evaluation with DFT and large-scale searches on supercomputers are out of scope.
