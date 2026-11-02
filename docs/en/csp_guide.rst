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

Run every command in this guide from the repository root (the directory that contains
``pyproject.toml``). Python 3.11 or later is required.

.. code-block:: bash

   python -m venv .venv          # once, if the virtual environment does not exist yet
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
     - Space group numbers to try. Required in ``informed`` mode; in ``blind`` mode it is
       replaced by the default list
   * - ``z``
     - omitted
     - Z per space group (e.g. ``z = { 14 = 4 }``). Default: general-position multiplicity × ``z_prime``.
       Ignored in ``blind`` mode
   * - ``z_prime``
     - ``1``
     - Molecules per asymmetric unit. Set to 1 in ``blind`` mode (with a warning if the file says
       otherwise); use ``informed`` mode for Z' = 2 (the cost grows)
   * - ``n_candidates_per_sg``
     - ``150``
     - Candidates generated per space group
   * - ``generators``
     - ``{random=0.6, voronoi=0.4}``
     - Relative weights of the generators (non-negative; names ``random``, ``voronoi``, ``seed``).
       ``random`` and ``voronoi`` share ``n_candidates_per_sg`` in this ratio. ``seed`` is only an
       on/off switch: when seeds are given it defaults to 0.5 (on), and ``seed = 0`` ignores them; each
       seed always gives one candidate per substitution choice
   * - ``seed_structures``
     - omitted
     - List of seed CIFs (``[[seed_structures]]``), see below
   * - ``packing_range``
     - ``[0.4, 0.8]``
     - Allowed packing fraction of the generated, unrelaxed cells
       (molecules × molecular vdW volume / cell volume)
   * - ``overlap_factor``
     - ``0.7``
     - Reject candidates with an interatomic distance below this factor × the sum of vdW radii
   * - ``n_prescreen_keep``
     - ``32``
     - Candidates kept after prescreening, i.e. relaxed with GFN2-xTB
   * - ``relax_steps``
     - ``60``
     - Maximum relaxation steps per candidate
   * - ``relax_fmax``
     - ``0.05``
     - Force convergence threshold (eV/Å)
   * - ``n_procs``
     - number of CPU cores
     - Number of worker processes
   * - ``max_mol_rmsd``
     - ``0.5``
     - Reject candidates in which a molecule deviates more than this heavy-atom RMSD (Å) from the
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

Minimal input files
~~~~~~~~~~~~~~~~~~~

A blind search needs only the molecule; every other key takes its default:

.. code-block:: toml

   molecule = "data/benzene.toml"
   mode = "blind"

An informed search gives what is known from experiment, here with a seed CIF:

.. code-block:: toml

   molecule = "data/my_molecule.toml"
   mode = "informed"
   space_groups = [14]
   z = { 14 = 4 }
   z_prime = 1

   [[seed_structures]]
   cif = "parent.cif"
   substitutions = [["O", "S", 1]]

``input/csp.toml`` in the repository is a commented example that lists every key.

The molecule file ``data/<name>.toml`` has the same format as the molecule definitions used by the rest of
the toolkit (see :doc:`toml_reference`, "Molecule Definition"). CSP reads only ``symbols``, ``positions``
and, when ``use_fractional = true``, ``lattice_vectors`` (to convert the positions); other keys such as
``center_of_mass`` are not used. An ``.xyz`` file holding one molecule works as well. The molecule is
shifted so that its plain (unweighted) centroid is at the origin.

Seed CIFs (seed_structures)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   [[seed_structures]]
   cif = "parent.cif"
   substitutions = [["O", "S", 1]]   # (old element, new element, count per molecule)
   volume_scale = "auto"             # or a number

The parent CIF must be **free of disorder** (no partial occupancies): order it, for example in VESTA,
and save it as P1. A CIF with partial occupancy is rejected.

``volume_scale = "auto"`` scales the parent cell isotropically by the ratio of the van der Waals volumes
(input molecule / parent molecule), i.e. each cell length by the cube root of that ratio. A number gives
the volume ratio directly. Molecules are moved rigidly so that their centroids keep their fractional
coordinates; the molecules themselves are never scaled.

Seed candidates come in addition to ``n_candidates_per_sg``: each seed gives one candidate per
substitution choice, all in P1.

Two modes
---------

``blind``
   Uses no information about the target crystal: the six default space groups, Z equal to the
   general-position multiplicity and Z' = 1. ``space_groups`` and ``z`` in the file are ignored and
   ``z_prime`` is set to 1 (with a warning if the file says otherwise). ``seed_structures`` are still
   used if given, because a seed is a related crystal (e.g. the parent compound before O→S substitution),
   not the experimental structure of the target. For a prediction from the molecule alone, leave
   ``seed_structures`` out (or set ``seed = 0`` in ``generators``).

``informed``
   You give the space groups (required), ``z`` and ``z_prime`` known from experiment. The search
   concentrates on the region around the answer. ``seed_structures`` can be added as in ``blind`` mode.

How Voronoi is used
-------------------

The earlier Voronoi pipeline (``run.py``) maximised the IoU against Voronoi cells of a centroid lattice
taken from a parent structure. In CSP the Voronoi cell is used to choose the **initial orientation**:

- Take the molecular centroids of a PyXtal crystal and build the Voronoi cell of each molecule.
- Find the orientation that maximises the IoU between the cell and **the molecule's own shape** (the
  convex hull of its atoms inflated by vdW radii). For the comparison only, the shape is scaled to the
  cell volume. Atomic coordinates are not scaled; the molecule is only rotated as a rigid body.
- For elongated, nearly symmetric molecules head and tail are hard to tell apart, so both the best
  orientation and that orientation turned by 180° about the long axis become candidates.
- The space-group symmetry is kept: the orientation is chosen for the asymmetric-unit molecules, and the
  others follow from the symmetry operations.

Every candidate records its generator in the ``origin`` column of ``candidates.csv``
(``random``, ``voronoi`` or ``seed:<name>``). The validation report ``validation.md`` repeats it for the
best match and for each of the closest candidates, so it shows which generator produced the structure
closest to experiment.

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
        - Density after relaxation (g/cm³)
      * - ``space_group_in``
        - Space group requested at generation
      * - ``space_group_out``
        - Space group re-detected by spglib after relaxation (0 = detection failed)
      * - ``origin``
        - Generator (``random`` / ``voronoi`` / seed name)
      * - ``max_mol_rmsd``
        - Deformation of the molecules after relaxation (heavy-atom RMSD, Å)
      * - ``status``
        - ``relaxed`` (accepted)

``001_<id>.cif`` ... ``050_<id>.cif``
   Structures of the top 50. When ``out_dir`` is reused, ranked CIFs (``NNN_<id>.cif``) of the
   earlier run are deleted first; ``candidates.csv`` and ``run.log`` are overwritten. Validation reads
   only the CIFs listed in ``candidates.csv``.

``run.log``
   How many candidates survived each stage and how long it took. Example from a real run with this
   input (all other keys at their defaults: ``generators = {random = 0.6, voronoi = 0.4}``,
   ``packing_range = [0.4, 0.8]``, ``overlap_factor = 0.7``):

   .. code-block:: toml

      molecule = "data/benzene.toml"
      mode = "informed"
      space_groups = [2, 14]          # P-1 (Z = 2) and P2_1/c (Z = 4)
      n_candidates_per_sg = 100
      n_prescreen_keep = 150
      relax_steps = 60
      n_procs = 2
      time_budget_min = 3
      rng = 1

   The log (per-space-group lines omitted)::

      generated 200 candidates {'random': 120, 'voronoi': 80} in 22.9 s
      prescreen {'generated': 200, 'no_overlap': 85, 'packing': 85, 'kept': 85} in 3.3 s
      time budget: 37.0 s per candidate, 62 s left; relaxing 0 of 83 remaining candidates
      relax: 2 relaxed, 0 rejected_deformed, 0 failed in 73.9 s
      dedupe: 2 -> 2 in 0.0 s

   The counts follow from the input: 2 space groups × 100 = 200 candidates, 60 random and 40 Voronoi
   per space group. 85 have no overlap and all of them are inside the packing range; because
   ``n_prescreen_keep = 150`` is larger, all 85 are kept. The first round (``n_procs`` = 2 candidates,
   74 s) already used most of the 3-minute budget, so the other 83 were not relaxed. Generation takes
   23 s mostly because of the Voronoi orientation search. Timings vary by a few seconds between runs.

   ``rejected_deformed`` counts candidates in which a molecule deformed beyond ``max_mol_rmsd`` or
   became covalently bonded to another molecule or to its own periodic image during relaxation;
   ``failed`` counts
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
structure and from each candidate and overlaid, giving RMSD\ :sub:`15`. RMSD\ :sub:`15` is over heavy
atoms (hydrogens excluded; for a molecule without heavy atoms all atoms are used). A candidate
**matches** if RMSD\ :sub:`15` < 0.8 Å with all 15 molecules paired. The comparison is made with each
structure as the reference and the worse result counts, so extra or missing molecules in either
structure prevent a match. Structures whose molecular densities (molecules per volume) differ by more
than a factor of 1.2 are listed as not comparable; the same polymorph after relaxation differs by only a
few percent.

``success``
   yes if a match exists within the top ``n_top`` (default 20).
``best match rank``
   Rank of the best-ranked matching candidate among the compared CIFs. It is listed even when it is
   below ``n_top``, which shows that the structure was generated but ranked too low. Only the CIFs that
   were written (the top 50 after deduplication) are compared, so "no match" does **not** prove that the
   structure was never generated: it may have been ranked below 50, or dropped by prescreening, the time
   budget or the deformation check.
``mirror``
   Whether the match used the mirror image of the candidate.
``flip match rank``
   Best rank of a candidate that matches only when O and S are treated as one element: **the same
   packing with O and S exchanged on the molecules**. It is never counted as success. It reads
   ``skipped`` when the experimental structure lacks O or S.
``closest candidates``
   Always shown for each experimental structure: up to 10 comparable candidates, ordered by the number of
   matched molecules (most first), then RMSD, then rank.
``Not comparable``
   Candidates that could not be compared, with the reason, for example: the cluster molecules have
   different compositions, the molecular densities differ too much, there are no molecules to compare,
   or the candidate CIF could not be read or is missing.

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

GFN2-xTB cost grows faster than the cube of the atom count (a cubic law would give 27 times); 72 atoms
cost roughly 40 times more per step than 24 atoms in these measurements.
Relaxation takes almost all of the time.

**Why these defaults** (fit within 90 % of ``time_budget_min = 60``, i.e. 54 min):

- Generation and prescreening: 6 space groups × 150 = 900 candidates. Per space group these are
  90 random crystals and 60 Voronoi candidates; the Voronoi candidates come from 30 PyXtal skeletons,
  each giving 2 orientations (best and turned 180°, head/tail). PyXtal builds 6 × (90 + 30) = 720
  crystals (about 30 s), the Voronoi orientation takes 6 × 30 skeletons × 0.7 s (about 2 min), LJ
  ranking about 20 s: about 3-4 min. Seed candidates, if any, come on top.
- One relaxation round (8 candidates at once): 60 steps × 10.9 s = about 11 min.
- The remaining 54 - 4 = 50 min allows 4 rounds, about 44 min, i.e. 8 × 4 = **32 candidates**
  (``n_prescreen_keep = 32``).
- Total about 48 min, inside the target. **This total is calculated from the measured step times; a full
  default run was not timed end to end.** Convergence at 60 steps was not checked: ``relax_fmax`` may not be
  reached.
- A larger molecule than naphthalene, such as OS-2SMe, takes longer per step and therefore gets fewer rounds
  within the same time.

With ``relax_steps = 150`` one round takes about 27 min and only one round (8 candidates) fits into an
hour. A relaxation that stops before full convergence is still adequate for narrowing down candidates,
so 60 is used. For production runs, raise ``time_budget_min``, reduce the number of space groups
(``informed`` mode), or adjust ``n_procs`` to the number of cores.

The time target works as follows. Candidates are relaxed in rounds of ``n_procs``, and every round is
timed. Before each further round the time left (90 % of ``time_budget_min`` minus the time used so far)
is compared with the slowest round so far; if it is shorter, the remaining candidates are skipped.
Because the check is repeated before every round, later rounds that are slower than the first one stop
the run early; a single round that is unexpectedly slow can still overrun the target by up to its own
length.
**The first round always runs** (about 11 min at the defaults for naphthalene, plus about 4 min of
generation), so a ``time_budget_min`` under about 15 min is exceeded.

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
  energy differences are smaller than the method error. The experimental structure can be generated
  yet ranked below ``n_top`` (or below 50, where it is no longer compared).
- **Flip disorder is not modelled.** Mixtures of O/S-exchanged molecules (partially occupied sites) are
  not reproduced; they are only reported separately as a "flip match".
- Seed and experimental CIFs must be **free of disorder**.
- Relaxations stop at ``relax_steps`` and may not be fully converged; energies refer to that state.
- Re-evaluation with DFT and large-scale searches on supercomputers are out of scope.
