Crystal Structure Prediction (CSP) Guide
========================================

This guide describes how to build candidate molecular crystals from the structure of an isolated
molecule alone, relax them with GFN2-xTB (the default) or the machine-learned potential OrbMol-v2,
rank them by energy, and compare them with experimental CIF files. It is written with OS-2SMe (molecular crystals with O/S exchange, "flip", between molecules)
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

The ``csp`` extra installs PyXtal, tblite (GFN2-xTB) and pymatgen; spglib, which CSP also uses,
comes in with them. They are not part of the core dependencies. If one of them cannot be imported,
``molecular-modeling-csp`` stops before generating anything and names the missing package, with the
hint ``pip install -e ".[csp]"``.

To relax with the machine-learned potential OrbMol-v2 (see `Choosing the relaxation method`_), install
the ``csp-ml`` extra instead; it includes everything of ``csp``:

.. code-block:: bash

   pip install -e ".[csp-ml]"    # or: uv sync --extra csp-ml

It adds orb-models with PyTorch and NVIDIA Warp (about 1 GB, all under the Apache 2.0 licence) and
needs Python 3.12 or later. Users of GFN2-xTB alone do not need it.

Workflow
--------

.. code-block:: bash

   molecular-modeling-csp input/csp.toml     # predict (default config: input/csp.toml)
   python -m molecular_modeling.csp.validate csp_out exp1.cif exp2.cif   # compare with experiment

``exp1.cif`` and ``exp2.cif`` stand for experimental CIF files that you supply; none ship with the
repository. Relative paths (the config file, and ``molecule``, ``cif`` and ``out_dir`` inside it) are
resolved against the current directory, not against the directory of the config file. Run both
commands as a single process: under ``mpiexec`` with more than one process they stop with an error
(use ``n_procs`` for parallel relaxation). If the input is
wrong (bad settings, a missing or unreadable file, a missing package, a relaxation method that cannot
run), the command prints
``error: <message>`` and exits with code 2.

A run goes through these stages:

1. **Generation**: PyXtal builds crystals for each space group (random). A share of them gets its
   molecular orientations re-chosen from Voronoi cells (Voronoi generator). If seed CIFs are given,
   candidates are also derived from them (seed generator).
2. **Overlap check and prescreening**: candidates with atoms too close together, broken molecules or a
   packing fraction outside the range are dropped. The rest are ranked by a cheap rigid-molecule
   intermolecular energy (UFF Lennard-Jones) and the best ``n_prescreen_keep`` are kept.
3. **Relaxation**: atoms and cell are relaxed (in parallel) with the method set by ``relax_method``
   (GFN2-xTB by default). Candidates whose molecules
   deform too much from the input are dropped.
4. **Deduplication and ranking**: identical structures are merged and ranked by energy per molecule.
5. **Output**: a ranking table and CIF files for the top 50.

Atomic coordinates of the molecule are never scaled. The only operations on a molecule are rigid
rotation and translation, apart from the relaxation.

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
       replaced by the default list (with a warning if the file gives a different list)
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
     - Candidates requested per space group (random + Voronoi). This is a maximum: PyXtal gives up
       on a space group after 5 × the requested number of attempts, and a Voronoi variant that
       breaks a molecule apart is skipped, so fewer can result. Shortfalls are logged in ``run.log``
   * - ``generators``
     - ``{random=0.6, voronoi=0.4}``
     - Relative weights of the generators (finite, non-negative; names ``random``, ``voronoi``,
       ``seed``). The table is merged over the default: a name you leave out keeps its default weight,
       so ``generators = {seed = 0}`` only switches seeds off. Write 0 to switch a generator off.
       ``random`` and ``voronoi`` share ``n_candidates_per_sg`` in this ratio. ``seed`` is only an
       on/off switch: when seeds are given it defaults to 0.5 (on), any positive value turns them on and
       ``seed = 0`` ignores them; each
       seed always gives one candidate per substitution choice. At least one generator must be on
   * - ``seed_structures``
     - omitted
     - List of seed CIFs (``[[seed_structures]]``), see below
   * - ``packing_range``
     - ``[0.4, 0.8]``
     - Allowed packing fraction of the generated, unrelaxed cells
       (molecules × molecular vdW volume / cell volume)
   * - ``overlap_factor``
     - ``0.7``
     - Reject candidates in which atoms of different molecules (including a molecule and its own
       periodic images) are closer than this factor × the sum of their vdW radii
   * - ``n_prescreen_keep``
     - ``32``
     - Candidates kept after prescreening, i.e. relaxed
   * - ``relax_method``
     - ``"GFN2-xTB"``
     - Relaxation method: ``"GFN2-xTB"``, ``"GFN1-xTB"`` or ``"OrbMol-v2"`` (needs the ``csp-ml``
       extra); see `Choosing the relaxation method`_. The command-line option ``--relax-method``
       overrides it
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
       input molecule after relaxation. Candidates in which a molecule gained or lost a bond within
       itself (a ring closed, a hydrogen moved to another atom) are rejected whatever the RMSD
   * - ``time_budget_min``
     - ``60``
     - Wall-clock target in minutes; see "Timings and defaults"
   * - ``rng``
     - ``1``
     - Random seed. The same value gives the same candidates only with the same input files, settings
       and software versions (PyXtal, NumPy, ...)
   * - ``out_dir``
     - ``"csp_out"``
     - Output directory

The default space groups are ``[14, 2, 19, 15, 4, 61]`` (P2\ :sub:`1`/c, P-1,
P2\ :sub:`1`\ 2\ :sub:`1`\ 2\ :sub:`1`, C2/c, P2\ :sub:`1`, Pbca). Unknown keys in the file are an error,
and so are values outside their range: ``packing_range`` needs 0 < lower < upper ≤ 1, space groups
(and the keys of ``z``) 1-230, ``n_candidates_per_sg``, ``n_prescreen_keep``, ``n_procs``, ``z_prime``
and the values of ``z`` at least 1, ``relax_steps`` at least 0, and ``relax_fmax``, ``max_mol_rmsd``,
``time_budget_min`` and ``overlap_factor`` positive, and ``rng`` a non-negative integer. A value of the
wrong type (for example ``space_groups = 14``) is an error naming the key, and so is an unknown key
inside ``[[seed_structures]]``. Repeated space groups are removed, and ``z`` given for a space group
that is not searched is ignored; both with a warning. Warnings issued while reading the file (such as
keys that ``blind`` mode ignores) are printed and repeated at the top of ``run.log``.

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

``input/csp.toml`` in the repository is a commented example that lists every key. It is a light
benzene example, not the built-in defaults: ``n_candidates_per_sg = 50`` (default 150) and
``n_procs = 4`` (default: all CPU cores), so that it finishes quickly. Its ``space_groups`` line is
commented out because ``blind`` mode ignores it; uncomment it together with ``mode = "informed"``.

The molecule file ``data/<name>.toml`` has the same format as the molecule definitions used by the rest of
the toolkit (see :doc:`toml_reference`, "Molecule Definition"); ``data/benzene.toml`` and
``data/naphthalene.toml`` ship with the repository. CSP needs only the element symbols and the Cartesian
coordinates (Å), for example (an O=C=O molecule):

.. code-block:: toml

   symbols = ["O", "C", "O"]
   positions = [[-1.16, 0.0, 0.0],
                [ 0.00, 0.0, 0.0],
                [ 1.16, 0.0, 0.0]]

CSP reads only ``symbols``, ``positions`` and, when ``use_fractional = true``, ``lattice_vectors`` (to
convert the positions); other keys such as ``center_of_mass`` are not used. An ``.xyz`` file holding one
molecule (exactly one frame; extended XYZ with a ``Properties=`` header is read by its columns) works as well. Symbols must be element symbols with the usual
capitalisation (``C``, ``Cl``), and no two atoms may be closer than 0.5 Å; otherwise loading stops with a
message naming the file. The molecule is shifted so that its plain (unweighted) centroid is at the origin.

.. tip::

   A molecule file can be made from a crystal structure in the GUI: select the atoms of one molecule and
   choose **MoleculeTools** → **Save Selection as Molecule...**, which writes ``data/<name>.toml``
   (see :doc:`gui_guide`). A molecule split across the cell boundary is rejoined before it is saved.
   The selection must be one whole, connected molecule; if it is not (several molecules, or atoms
   missing), the dialog warns and asks before saving.

**Supported elements**: H, C, N, O, F, S, Cl, Se, Br and I. The rigid Lennard-Jones prescreen has UFF
parameters only for these; a molecule with any other element is rejected with an error before
generation starts.

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
   ``z_prime`` is set to 1 (each with a warning if the file says otherwise). ``seed_structures`` are still
   used if given, because a seed is a related crystal (e.g. the parent compound before O→S substitution),
   not the experimental structure of the target. For a prediction from the molecule alone, leave
   ``seed_structures`` out, or switch seeds off with ``generators = {seed = 0}`` (the other weights keep
   their defaults).

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
- The shape cannot tell head from tail, so for every skeleton both the best orientation and that
  orientation turned by 180° about the axis normal to the molecular plane (through the centroid) become
  candidates. That turn exchanges the two ends (for OS-2SMe, the O end and the S end) and keeps the
  molecular plane. If the turn maps every molecule onto itself (atoms of the same element within 0.1 Å,
  as for benzene or naphthalene), the turned copy would be a duplicate and is left out; more skeletons
  are then used, so the Voronoi share of ``n_candidates_per_sg`` is still filled. In ``candidates.csv``
  both variants have origin ``voronoi``.
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
        - Candidate identifier ``<origin>-<space group>-<serial>``, with ``:`` and ``/`` in the origin
          replaced by ``_``; e.g. ``random-14-3``, ``voronoi-2-17``, ``seed_parent-1-0``
      * - ``energy_per_mol_eV``
        - Energy per molecule (eV) from the relaxation method. Each method has its own zero of
          energy: compare energies only between rows of the same method
      * - ``rel_energy_kJmol``
        - Difference from rank 1 (kJ/mol per molecule)
      * - ``density_gcm3``
        - Density after relaxation (g/cm³)
      * - ``space_group_in``
        - Space group requested at generation
      * - ``space_group_out``
        - Space group re-detected by spglib after relaxation (0 = detection failed)
      * - ``origin``
        - Generator: ``random``, ``voronoi`` or ``seed:<name>`` (``<name>`` = seed CIF file name without
          extension, e.g. ``seed:parent``)
      * - ``max_mol_rmsd``
        - Deformation of the molecules after relaxation (heavy-atom RMSD, Å)
      * - ``status``
        - ``relaxed``: the candidate passed the checks after relaxation (deformation, bonds). It does
          not mean that the relaxation converged to ``relax_fmax``; it may have stopped at
          ``relax_steps``
      * - ``relax_method``
        - Relaxation method used for the run (``GFN2-xTB``, ``GFN1-xTB`` or ``OrbMol-v2``)

``001_<id>.cif`` ... ``050_<id>.cif``
   Structures of the top 50. Validation reads only the CIFs listed in ``candidates.csv``.

Reusing ``out_dir``
   The settings, the molecule and the seed CIFs are checked first; if one of them is wrong the run stops
   without touching ``out_dir``, so the results of an earlier run stay as they were. Once the run
   starts, it deletes the earlier ``candidates.csv``, ranked CIFs (``NNN_<id>.cif``) and
   ``validation.md`` and overwrites ``run.log``. The new results are written at the end, with
   ``candidates.csv`` last. If the run stops with an error or is interrupted (Ctrl-C), the reason is
   written to ``run.log`` and no results are left, so an earlier run can never be mistaken for the new
   one. While the results are being replaced the directory holds the marker ``.csp_incomplete``; it
   is written before anything old is deleted and removed only after the new ``candidates.csv`` is in
   place. A run that fails, is interrupted or is killed leaves the marker behind, and validation refuses
   the directory until a run finishes. Do not start two runs with the same ``out_dir``: a running job
   holds the lock file ``.csp_run.lock``, and validation holds the same lock while it reads the results
   and writes ``validation.md``, so a run and a validation of the same directory cannot overlap; the
   second one is refused. Validation therefore needs write access to ``out_dir``: if it cannot create
   the lock file (for example a read-only copy of the results), it stops with an error instead of
   reading without the lock; copy the directory somewhere writable first. If a run was killed so hard
   that the lock file stayed behind, delete the lock file **and rerun**: the marker keeps the partial
   results from being validated.

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

      generated 200 candidates {'random': 120, 'voronoi': 80} in 49.5 s
      prescreen {'generated': 200, 'no_overlap': 86, 'packing': 86, 'kept': 86} in 3.8 s
      time budget: 43.6 s per candidate, 21 s left; relaxing 0 of 84 remaining candidates
      relax: 2 relaxed, 0 rejected_deformed, 0 failed in 87.3 s
      dedupe: 2 -> 2 in 0.0 s

   The counts follow from the input: 2 space groups × 100 = 200 candidates, 60 random and 40 Voronoi
   per space group. Benzene is unchanged by the 180° turn, so each of the 40 Voronoi candidates comes
   from its own PyXtal skeleton. 86 have no overlap and all of them are inside the packing range; because
   ``n_prescreen_keep = 150`` is larger, all 86 are kept. The first round (``n_procs`` = 2 candidates,
   87 s) already used most of the 3-minute budget, so the other 84 were not relaxed. Generation takes
   50 s mostly because of the Voronoi orientation search. The relaxation counts depend on how long the
   first round took, so they are observations of this timed run rather than fixed results; timings vary by
   a few seconds between runs. The log was reproduced with identical counts using Python 3.13.12,
   pyxtal 1.1.5, tblite 0.8.0, ase 3.29.0, numpy 2.5.4, scipy 1.18.1, pymatgen 2026.10.24 and
   spglib 2.8.0.

   ``rejected_deformed`` counts candidates in which, during relaxation, a molecule deformed beyond
   ``max_mol_rmsd``, gained or lost a bond within itself, or became covalently bonded to another
   molecule or to its own periodic image; ``failed`` counts calculations that stopped with an error.
   A missing tblite does not show up here as every candidate failing: the run stops before generation
   (see `Installation`_). After the ``relax:`` line, one line per reason gives the number of failed or
   rejected candidates and the first message, for example
   ``2 failed: InputError: Too close interatomic distances``. A line starting with ``time budget:``
   means the number of relaxations was reduced to fit the time target.

Validation against experiment
-----------------------------

.. code-block:: bash

   python -m molecular_modeling.csp.validate csp_out exp_OS-i.cif exp_OS-ii.cif
   # optional: --n-top 20  --rmsd-cut 0.8  --flip O,S  (--no-flip skips the flip check)

The report is printed and written to ``csp_out/validation.md``. Experimental CIFs must be
**disorder-free**: every site fully occupied (occupancy exactly 1) by one element; for a disordered
structure, keep only the major component and save it as an ordered P1 CIF. They may lack hydrogen atoms, as many experimental CIFs do: only
heavy atoms are compared, so a CIF with or without hydrogens matches the same candidates. If the output
directory does not exist, holds neither ``candidates.csv`` nor candidate CIFs, is in use by a run, or
holds the results of a run that did not finish (see "Reusing ``out_dir``"), the command stops with an
error instead of writing a report; so does an experimental
CIF that cannot be read. Experimental structures are labelled by file name in the report, or by the path
as given when two files have the same name.

The flip check (``--flip O,S``) is **on by default**. For a molecule without O or S it has nothing to
do and is reported as skipped; give ``--no-flip`` to leave it out, or ``--flip`` with two other, different
elements (for example ``--flip Se,S``; case does not matter). ``--flip`` and ``--no-flip`` cannot be
combined, ``--n-top`` must be at least 1 and ``--rmsd-cut`` a positive number.

The comparison follows the COMPACK approach: a cluster of 15 molecules is taken from the experimental
structure and from each candidate and overlaid, giving RMSD\ :sub:`15`. RMSD\ :sub:`15` is over heavy
atoms (hydrogens excluded; for a molecule without heavy atoms all atoms are used). A candidate
**matches** if RMSD\ :sub:`15` < 0.8 Å with all 15 molecules paired. The comparison is made with each
structure as the reference and the worse result counts, so extra or missing molecules in either
structure prevent a match. Structures whose molecular densities (molecules per volume) differ by more
than a factor of 1.2 are listed as not comparable; the same polymorph after relaxation differs by only a
few percent.

The report uses these terms (written below exactly as they appear in ``validation.md``):

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
   packing with O and S exchanged on the molecules**. It is never counted as success. When the
   experimental structure lacks one of the two elements it reads, for example,
   ``skipped (no O, S in the experimental structure)`` (both missing) or
   ``skipped (no S in the experimental structure)``.
``closest candidates``
   Always shown for each experimental structure: up to 10 comparable candidates, ordered by the number of
   matched molecules (most first), then RMSD, then rank.
``Not comparable``
   Candidates that could not be compared, with the reason, for example: the cluster molecules have
   different compositions, the molecular densities differ too much, there are no molecules to compare,
   or the candidate CIF could not be read or is missing.

Choosing the relaxation method
------------------------------

``relax_method`` (or ``--relax-method`` on the command line, which takes precedence) selects how
candidates are relaxed and ranked:

``GFN2-xTB`` (default), ``GFN1-xTB``
   Tight-binding methods from tblite, part of the ``csp`` extra.
``OrbMol-v2``
   A machine-learned interatomic potential (orb-models, trained on OMol25, which covers every element
   CSP supports, Se included). Needs the ``csp-ml`` extra. Molecules are treated as neutral singlets.

In a small comparison (2026-11-06; naphthalene and two Se compounds, the experimental structure plus
6-12 random PyXtal candidates each, up to 300 relaxation steps), OrbMol-v2 ranked the experimental
structure lowest in all three systems and took about half the time of GFN2-xTB per candidate, while
GFN2-xTB ranked one random candidate 20 kJ/mol below the experimental naphthalene structure and failed
(no SCF convergence) on the experimental structure of diphenyl diselenide. **For molecules with Se, and
when GFN2-xTB ranks poorly, OrbMol-v2 is worth trying.** The comparison was small, used looser
candidates than this pipeline produces, and was not checked against DFT; the default stays GFN2-xTB.

Energies of different methods have different zeros: compare ``energy_per_mol_eV`` and rankings only
within one method. The ``relax_method`` column of ``candidates.csv`` and the first lines of
``run.log`` record which method was used (and, in ``run.log``, where the setting came from:
``default``, ``config``, ``cli`` or ``argument``, the package version and, for OrbMol-v2, the weights
file).

**Weights and network access.** On first use OrbMol-v2 downloads its weights (about 100 MB) into
``~/.cache/cached_path`` (or the directory in ``CACHED_PATH_CACHE_ROOT``). Before an OrbMol-v2 run
touches ``out_dir``, a separate process builds the model and evaluates one water molecule (at most
10 min). If that fails, for example on a computer without network access, the run stops with a
one-line ``error:`` that names the cause and the cache directory. To prepare such a computer, run this
on one with network access and copy the whole cache directory over:

.. code-block:: bash

   python -c "from orb_models.forcefield import pretrained; pretrained.orbmol_v2(device='cpu')"

Each worker process loads the model once, before the first round of relaxations (including its slow
first evaluation); that time counts against ``time_budget_min`` but not as part of a round. **Memory:**
each OrbMol-v2 worker used about 1 GB in our tests, and ``n_procs`` defaults to the number of CPU
cores; on a machine with many cores and little memory, set ``n_procs`` lower. The relaxation settings (``relax_steps``,
``relax_fmax``) and the prescreening are the same for every method. The timings below were measured
with GFN2-xTB; OrbMol-v2 has not yet been timed under the conditions of this pipeline (one thread per
worker).

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

**Why these defaults** (fit within 90 % of ``time_budget_min = 60``, i.e. 54 min). The arithmetic
assumes that every requested candidate is generated. A shortfall (see ``n_candidates_per_sg``) leaves
fewer candidates for screening, but it need not shorten the run: PyXtal keeps retrying up to 5 × the
requested number of attempts, which can make generation take longer, and as long as at least
``n_prescreen_keep`` candidates pass prescreening the relaxation time is unchanged:

- Generation and prescreening: 6 space groups × 150 = 900 candidates. Per space group these are
  90 random crystals and 60 Voronoi candidates. For naphthalene the 180° turn gives the same
  orientation, so each PyXtal skeleton gives one Voronoi candidate and 60 skeletons are needed (a
  molecule with different ends, such as OS-2SMe, gives two per skeleton and needs 30). PyXtal builds
  6 × (90 + 60) = 900 crystals (about 36 s), the Voronoi orientation takes 6 × 60 skeletons × 0.7 s
  (about 4 min), LJ ranking about 20 s: about 5-6 min. Seed candidates, if any, come on top.
- One relaxation round (8 candidates at once): 60 steps × 10.9 s = about 11 min.
- The remaining 54 - 6 = 48 min allows 4 rounds, about 44 min, i.e. 8 × 4 = **32 candidates**
  (``n_prescreen_keep = 32``).
- Total about 50 min, inside the target. **This total is calculated from the measured step times; a full
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
**The first round always runs** (about 11 min at the defaults for naphthalene, plus about 6 min of
generation), so a ``time_budget_min`` under about 17 min is exceeded.

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
- **The energy ranking (GFN2-xTB or OrbMol-v2) may not separate close polymorphs**, especially for dipolar molecules, where the
  energy differences are smaller than the method error. The experimental structure can be generated
  yet ranked below ``n_top`` (or below 50, where it is no longer compared).
- **Flip disorder is not modelled.** Mixtures of O/S-exchanged molecules (partially occupied sites) are
  not reproduced; they are only reported separately as a "flip match".
- Seed and experimental CIFs must be **free of disorder**.
- Relaxations stop at ``relax_steps`` and may not be fully converged; energies refer to that state.
- Re-evaluation with DFT and large-scale searches on supercomputers are out of scope.
