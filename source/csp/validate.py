"""Compare CSP candidates with experimental crystal structures (cluster RMSD, COMPACK style).

Usage::

    python -m molecular_modeling.csp.validate <out_dir> <exp.cif> [<exp.cif> ...]

The Markdown report is printed and written to ``<out_dir>/validation.md``.
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation

from .candidate import MAX_CIF
from .gen_seed import read_seed_molecules, split_molecules
from .molecule import AlignResult, _assign, _initial_rotations, _kabsch, match_and_align

logger = logging.getLogger(__name__)

#: Number of molecules in the comparison cluster.
N_CLUSTER = 15
#: A cluster molecule counts as matched if its centroid is within this distance (Å) after overlay.
MATCH_COM_DIST = 2.0
#: The cluster of ``b`` is drawn from this many times ``n_mol`` nearest molecules, so that ties and
#: slightly different distances at the cluster edge do not leave a molecule of ``a`` without a partner.
POOL_FACTOR = 2
#: Overlay refinement rounds (Kabsch on all matched cluster atoms, then re-matching).
N_REFINE = 3
#: Centre-overlay hypotheses (best first after one matching pass) that get the full refinement.
N_KEEP_REFINE = 4
#: A rotation is a (pseudo-)symmetry of the centre molecule if it maps the fitted atoms onto
#: themselves within this RMSD (Å).
SELF_SYM_TOL = 0.3
#: Random starting rotations used to find those symmetries, and the cap on how many are kept.
N_SELF_STARTS = 48
MAX_SELF_ROTATIONS = 120
#: Two rotation matrices closer than this (Frobenius norm) are the same.
SAME_ROT_TOL = 0.1
#: Two molecules of one structure are the same conformer if their heavy-atom RMSD after superposition is
#: below this (Å); overlays are composed only through a molecule of the same conformer.
CONFORMER_TOL = 0.3
#: Tolerance (Å) for finding symmetry-equivalent molecules; only one of each class is used as a cluster
#: centre. Kept small: a structure that is only nearly symmetric simply gets more centres.
CENTRE_SYMPREC = 0.01
#: Structures whose molecular number densities (molecules per Å^3) differ by more than this factor are
#: not compared. A 15-molecule cluster cannot detect a missing or extra molecule outside it, so a
#: structure with far fewer or more molecules per volume could otherwise look identical. Relaxed and
#: experimental cells of the same polymorph differ by a few percent (thermal expansion, method error),
#: well inside this factor.
MAX_DENSITY_RATIO = 1.2


@dataclass
class _Crystal:
    """Whole molecules of a periodic structure, prepared for cluster building."""

    cell: np.ndarray
    mols: list[Atoms]  # fitted atoms only (heavy atoms, or all atoms of a molecule that has none)
    coms: np.ndarray  # (n_mol, 3) centroids of the fitted atoms, wrapped into the cell
    shifts: np.ndarray  # (n_mol, 3) Cartesian shift from the molecule's own coordinates to the wrapped centroid
    heavy_pos: list[np.ndarray]  # positions of fitted atoms (heavy atoms, or all if the molecule has none)
    heavy_sym: list[np.ndarray]
    # (pseudo-)symmetry rotations of each centre molecule, computed lazily
    self_rotations: dict[int, list[np.ndarray]] = field(default_factory=dict)
    centres: list[int] | None = None  # one molecule per symmetry-equivalence class, computed lazily
    references: dict[int, int] = field(default_factory=dict)  # centre -> reference molecule of its conformer


def _prepare(cell: np.ndarray, mols: list[Atoms]) -> _Crystal:
    """Prepare whole molecules for comparison, keeping only the atoms that are fitted.

    Hydrogens are dropped from every molecule that has heavy atoms (a molecule of hydrogens only keeps
    them all), so the comparison, centroids and symmetry search do not depend on whether a structure
    has hydrogens: experimental CIFs often have none.
    """
    cell = np.asarray(cell, dtype=float)
    heavy_pos, heavy_sym, fitted = [], [], []
    for m in mols:
        sym = np.array(m.get_chemical_symbols())
        keep = sym != "H" if np.any(sym != "H") else np.ones(len(sym), dtype=bool)
        heavy_pos.append(m.get_positions()[keep])
        heavy_sym.append(sym[keep])
        fitted.append(Atoms(list(sym[keep]), positions=heavy_pos[-1]))
    raw = np.array([p.mean(axis=0) for p in heavy_pos]).reshape(-1, 3)
    frac = raw @ np.linalg.inv(cell)
    wrapped = (frac - np.floor(frac)) @ cell
    return _Crystal(cell, fitted, wrapped, wrapped - raw, heavy_pos, heavy_sym)


def _inverted(c: _Crystal) -> _Crystal:
    """The mirror-image (inverted through the origin) crystal."""
    mols = [Atoms(m.get_chemical_symbols(), positions=-m.get_positions()) for m in c.mols]
    return _prepare(-c.cell, mols)


def _with_symbols_mapped(c: _Crystal, element_map: dict[str, str]) -> _Crystal:
    """The same crystal with elements renamed (e.g. ``{"S": "O"}`` treats S and O as one element)."""
    mols = [
        Atoms([element_map.get(s, s) for s in m.get_chemical_symbols()], positions=m.get_positions()) for m in c.mols
    ]
    return _prepare(c.cell, mols)


def _self_fit(x: np.ndarray, sym: np.ndarray, rot: np.ndarray) -> tuple[np.ndarray, float]:
    """Refine ``rot`` mapping centred positions ``x`` onto themselves (per-element assignment + Kabsch)."""
    rmsd = math.inf
    for _ in range(20):
        perm = _assign(x @ rot.T, x, sym, sym)
        rot, _t = _kabsch(x[perm], x)
        new = float(np.sqrt(np.mean(np.sum((x[perm] @ rot.T - x) ** 2, axis=1))))
        if abs(rmsd - new) < 1e-10:
            break
        rmsd = new
    return rot, rmsd


def _self_rotations(pos: np.ndarray, sym: np.ndarray, tol: float = SELF_SYM_TOL) -> list[np.ndarray]:
    """Proper rotations about the centroid of ``pos`` that map the molecule onto itself within ``tol`` Å RMSD.

    Includes the identity. Found from many starting rotations, then closed under composition, so exact
    symmetries (benzene: 12, naphthalene: 4) and near-symmetries (e.g. a molecule turned over with O and S
    nearly superposed) are all returned.
    """
    x = pos - pos.mean(axis=0)
    starts = _initial_rotations(x, x) + list(Rotation.random(N_SELF_STARTS, random_state=1).as_matrix())
    found: list[np.ndarray] = [np.eye(3)]

    def known(r: np.ndarray) -> bool:
        return any(np.linalg.norm(r - q) < SAME_ROT_TOL for q in found)

    for r0 in starts:
        r, rmsd = _self_fit(x, sym, r0)
        if rmsd < tol and not known(r):
            found.append(r)
    grew = True
    while grew and len(found) < MAX_SELF_ROTATIONS:
        grew = False
        for r1 in list(found):
            for r2 in list(found):
                r, rmsd = _self_fit(x, sym, r1 @ r2)
                if rmsd < tol and not known(r):
                    found.append(r)
                    grew = True
    return found


def _cluster(c: _Crystal, centre: int, n: int) -> tuple[np.ndarray, np.ndarray]:
    """The ``n`` molecules (periodic images included) whose centroids are nearest to molecule ``centre``.

    Returns ``(indices, shifts)``: molecule indices and the Cartesian shift to add to each molecule's
    own coordinates. Lattice translations are enlarged until no image outside the searched block can
    be closer than the ``n``-th nearest one found.
    """
    vol = abs(np.linalg.det(c.cell))
    a, b, cc = c.cell
    h_min = vol / max(np.linalg.norm(np.cross(b, cc)), np.linalg.norm(np.cross(cc, a)), np.linalg.norm(np.cross(a, b)))
    origin = c.coms[centre]
    big = 1
    while True:
        r = np.arange(-big, big + 1)
        lat = np.stack(np.meshgrid(r, r, r, indexing="ij"), axis=-1).reshape(-1, 3) @ c.cell
        pts = c.coms[:, None, :] + lat[None, :, :]  # (n_mol, n_lat, 3)
        dist = np.linalg.norm(pts - origin, axis=2).ravel()
        if len(dist) >= n:
            order = np.argsort(dist, kind="stable")[:n]
            if dist[order[-1]] <= big * h_min:
                break
        big += 1
    mol_idx, lat_idx = np.divmod(order, len(lat))
    return mol_idx, c.shifts[mol_idx] + lat[lat_idx]


def _match_clusters(
    a_pos: list[np.ndarray],
    a_sym: list[np.ndarray],
    a_com: np.ndarray,
    b_pos: list[np.ndarray],
    b_sym: list[np.ndarray],
    b_com: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int] | None:
    """Pair cluster molecules by centroid, then atoms within each pair by position (per element).

    Returns ``(target, mobile, com_dist, n_extra)``: stacked matched atoms of a and b, centroid
    distances of the molecule pairs, and the number of *unpaired* b molecules whose centroid lies
    inside a's cluster (closer to a's centre molecule than the farthest cluster molecule minus
    :data:`MATCH_COM_DIST`). Such a molecule has no partner in a although a matching structure would
    have one there, so a denser b is not mistaken for a. ``None`` if a pair has different compositions.
    """
    cost = np.linalg.norm(a_com[:, None, :] - b_com[None, :, :], axis=2)
    rows, cols = linear_sum_assignment(cost)
    tgt, mob = [], []
    for i, j in zip(rows, cols):
        if sorted(a_sym[i]) != sorted(b_sym[j]):
            return None
        perm = _assign(b_pos[j], a_pos[i], b_sym[j], a_sym[i])
        tgt.append(a_pos[i])
        mob.append(b_pos[j][perm])
    radius = float(np.max(np.linalg.norm(a_com - a_com[0], axis=1))) - MATCH_COM_DIST
    unpaired = np.setdiff1d(np.arange(len(b_com)), cols)
    n_extra = int(np.sum(np.linalg.norm(b_com[unpaired] - a_com[0], axis=1) < radius))
    return np.vstack(tgt), np.vstack(mob), cost[rows, cols], n_extra


def _n_matched(com_dist: np.ndarray, n_extra: int) -> int:
    """Matched cluster molecules (partner centroid within MATCH_COM_DIST), less unpaired intruders of b."""
    return max(int(np.sum(com_dist < MATCH_COM_DIST)) - n_extra, 0)


@dataclass
class _Result:
    rmsd: float
    n_matched: int
    mirror: bool = False
    note: str | None = None  # why the structures could not be compared (rmsd is then inf)


def _formula(sym: np.ndarray) -> str:
    els, counts = np.unique(sym, return_counts=True)
    return "".join(f"{e}{n}" for e, n in zip(els, counts))


def _centres(c: _Crystal) -> list[int]:
    """Indices of molecules that represent every symmetry-equivalence class of ``c`` (lowest index each).

    Two molecules are equivalent if a space-group operation of the whole structure (found by spglib with
    :data:`CENTRE_SYMPREC`) maps an atom of one onto an atom of the other; operations map whole molecules
    onto whole molecules, so that suffices. If no symmetry is found every molecule is a centre.
    """
    if c.centres is not None:
        return c.centres
    import spglib

    owner = np.concatenate([np.full(len(m), k) for k, m in enumerate(c.mols)])
    pos = np.vstack([m.get_positions() for m in c.mols])
    numbers = np.concatenate([m.get_atomic_numbers() for m in c.mols])
    parent = list(range(len(c.mols)))

    def find(k: int) -> int:
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    ds = None
    if np.linalg.det(c.cell) > 0:
        ds = spglib.get_symmetry_dataset(
            (c.cell, pos @ np.linalg.inv(c.cell), numbers),  # type: ignore[arg-type]
            symprec=CENTRE_SYMPREC,
        )
    if ds is not None:
        equiv = np.asarray(ds["equivalent_atoms"] if isinstance(ds, dict) else ds.equivalent_atoms)
        for atom, rep in enumerate(equiv):
            ra, rb = find(int(owner[atom])), find(int(owner[rep]))
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    c.centres = sorted({find(k) for k in range(len(c.mols))})
    return c.centres


def _reference(a: _Crystal, ci: int) -> int:
    """The first centre of ``a`` whose molecule is the same conformer as molecule ``ci``.

    Same conformer = same heavy-atom composition and heavy-atom RMSD after superposition below
    :data:`CONFORMER_TOL`. Overlays are composed through this reference (pool centre -> reference ->
    centre), which is only valid if the reference and the centre have the same shape; a crystal with
    two different conformers (Z' = 2) gets one reference per conformer.
    """
    if ci in a.references:
        return a.references[ci]
    key = sorted(a.heavy_sym[ci])
    ref = ci
    for k in _centres(a):
        if k == ci:
            break
        if sorted(a.heavy_sym[k]) != key:
            continue
        try:
            if match_and_align(_atoms_of(a, k), _atoms_of(a, ci), heavy_only=True).rmsd < CONFORMER_TOL:
                ref = k
                break
        except ValueError:
            continue
    a.references[ci] = ref
    return ref


def _atoms_of(c: _Crystal, k: int) -> Atoms:
    return Atoms(c.mols[k].get_chemical_symbols(), positions=c.mols[k].get_positions())


_Pool = tuple[bool, int, np.ndarray, list[np.ndarray], list[np.ndarray], np.ndarray, Atoms]


def _pools(b: _Crystal, n_mol: int) -> list[_Pool]:
    """Clusters of ``b`` and of its mirror image around each symmetry-distinct molecule.

    A molecule related to another by a symmetry operation of ``b`` has a congruent cluster (an improper
    operation is covered by the mirror image), so one molecule per class is enough.
    """
    out: list[_Pool] = []
    centres = _centres(b)
    for mirror, bb in ((False, b), (True, _inverted(b))):
        for j in centres:
            b_idx, b_shift = _cluster(bb, j, POOL_FACTOR * n_mol)
            b_raw = [bb.heavy_pos[k] + s for k, s in zip(b_idx, b_shift)]
            b_com_raw = np.array([bb.mols[k].get_positions().mean(axis=0) + s for k, s in zip(b_idx, b_shift)])
            b_sym = [bb.heavy_sym[k] for k in b_idx]
            b_centre = Atoms(bb.mols[j].get_chemical_symbols(), positions=bb.mols[j].get_positions() + b_shift[0])
            out.append((mirror, j, bb.heavy_sym[j], b_raw, b_sym, b_com_raw, b_centre))
    return out


def _compose(first: AlignResult, second: AlignResult) -> tuple[np.ndarray, np.ndarray]:
    """Rotation and translation of ``x -> second(first(x))`` (column-vector convention ``R x + t``)."""
    return second.rotation @ first.rotation, second.rotation @ first.translation + second.translation


def _rmsd_from_centre(
    a: _Crystal, ci: int, pools: list[_Pool], pool_al: list[AlignResult | str], n_mol: int
) -> _Result:
    """Best overlay of the cluster of ``a`` around molecule ``ci`` onto any cluster in ``pools``.

    ``pool_al[k]`` aligns the centre molecule of pool ``k`` onto the reference molecule
    ``_reference(a, ci)`` (or is the reason it cannot be aligned). Combined with one alignment of that
    reference onto molecule ``ci``, it gives the centre overlay; the self-rotations of ``ci`` and the
    cluster refinement then correct the small error the composition can carry.
    """
    a_idx, a_shift = _cluster(a, ci, n_mol)
    a_pos = [a.heavy_pos[k] + s for k, s in zip(a_idx, a_shift)]
    a_sym = [a.heavy_sym[k] for k in a_idx]
    a_com = np.array([a.mols[k].get_positions().mean(axis=0) + s for k, s in zip(a_idx, a_shift)])
    a_centre = Atoms(a.mols[ci].get_chemical_symbols(), positions=a.mols[ci].get_positions() + a_shift[0])
    if ci not in a.self_rotations:
        a.self_rotations[ci] = _self_rotations(a.heavy_pos[ci], a.heavy_sym[ci])
    self_rots = a.self_rotations[ci]
    pivot = a_pos[0].mean(axis=0)  # the self-rotations act about the centroid of the fitted centre atoms

    def key(res: tuple[np.ndarray, np.ndarray, np.ndarray, int]) -> tuple[bool, float]:
        tgt, mob, com_dist, n_extra = res
        rmsd = float(np.sqrt(np.mean(np.sum((tgt - mob) ** 2, axis=1))))
        return _n_matched(com_dist, n_extra) < n_mol, rmsd

    hyps: list[tuple[tuple[bool, float], int, np.ndarray, np.ndarray]] = []
    notes: set[str] = set()
    ref = _reference(a, ci)
    try:
        ref_to_centre = match_and_align(_atoms_of(a, ref), a_centre, heavy_only=True)
    except ValueError as exc:
        return _Result(math.inf, 0, note=f"centre molecules cannot be aligned ({exc})")
    for pi, (_mirror, _j, centre_sym, b_raw, b_sym, b_com_raw, _b_centre) in enumerate(pools):
        if sorted(centre_sym) != sorted(a.heavy_sym[ci]):
            notes.add(
                f"molecule compositions differ ({_formula(a.heavy_sym[ci])} vs {_formula(centre_sym)}, heavy atoms)"
            )
            continue
        al = pool_al[pi]
        if isinstance(al, str):
            notes.add(al)
            continue
        rot0, trans0 = _compose(al, ref_to_centre)
        # every (pseudo-)symmetry-equivalent overlay of the centre molecule is a separate hypothesis
        for s in self_rots:
            rot = s @ rot0
            trans = s @ (trans0 - pivot) + pivot
            cur = _match_clusters(
                a_pos, a_sym, a_com, [p @ rot.T + trans for p in b_raw], b_sym, b_com_raw @ rot.T + trans
            )
            if cur is None:
                notes.add("cluster molecules have different compositions")
                continue
            hyps.append((key(cur), pi, rot, trans))
    if not hyps:
        return _Result(math.inf, 0, note="; ".join(sorted(notes)) or "no molecules to compare")

    hyps.sort(key=lambda h: h[0])
    best: tuple[tuple[bool, float], _Result] | None = None
    for _k, pi, rot, trans in hyps[:N_KEEP_REFINE]:
        mirror, _j, _cs, b_raw, b_sym, b_com_raw, _bc = pools[pi]
        res = None
        for it in range(N_REFINE + 1):
            cur = _match_clusters(
                a_pos, a_sym, a_com, [p @ rot.T + trans for p in b_raw], b_sym, b_com_raw @ rot.T + trans
            )
            if cur is None:
                break
            res = cur
            if it < N_REFINE:
                tgt, mob_moved, _, _ = cur
                # refit on the untransformed b coordinates of the matched atoms
                rot, trans = _kabsch((mob_moved - trans) @ rot, tgt)
        if res is None:
            continue
        k = key(res)
        if best is None or k < best[0]:
            best = (k, _Result(k[1], _n_matched(res[2], res[3]), mirror))
    assert best is not None
    return best[1]


def _rmsd_crystals(a: _Crystal, b: _Crystal, n_mol: int) -> _Result:
    """Cluster comparison with ``a`` as reference, centred on every symmetry-distinct molecule of ``a``.

    For each centre the best overlay onto ``b`` is found; the worst centre counts (fewest matched
    molecules, then largest RMSD), so an extra or missing molecule anywhere in ``a``'s cell is seen
    whatever the molecule order. A centre that cannot be compared makes the whole result not comparable.
    """
    if not a.mols:
        return _Result(math.inf, 0, note="no molecules to compare")
    pools = _pools(b, n_mol)
    # align each pool centre once onto the reference molecule of each conformer
    pool_al: dict[int, list[AlignResult | str]] = {}
    worst: tuple[tuple[bool, float], _Result] | None = None
    for ci in _centres(a):
        ref = _reference(a, ci)
        if ref not in pool_al:
            pool_al[ref] = []
            for _m, _j, centre_sym, _r, _s, _c, b_centre in pools:
                if sorted(centre_sym) != sorted(a.heavy_sym[ref]):
                    pool_al[ref].append("molecule compositions differ")
                    continue
                try:
                    pool_al[ref].append(match_and_align(b_centre, _atoms_of(a, ref), heavy_only=True))
                except ValueError as exc:
                    pool_al[ref].append(f"centre molecules cannot be aligned ({exc})")
        res = _rmsd_from_centre(a, ci, pools, pool_al[ref], n_mol)
        if res.note:
            return res
        k = (res.n_matched < n_mol, res.rmsd)
        if worst is None or k > worst[0]:
            worst = (k, res)
    assert worst is not None
    return worst[1]


def _density(c: _Crystal) -> float:
    return len(c.mols) / abs(float(np.linalg.det(c.cell)))


def _compare(a: _Crystal, b: _Crystal, n_mol: int) -> _Result:
    """Cluster comparison of ``a`` and ``b`` that must hold both ways and at matching density.

    The clusters are compared with ``a`` as reference and with ``b`` as reference, and the worse of the
    two counts (fewer matched molecules, larger RMSD), so neither structure may contain molecules the
    other lacks near the cluster centre. Structures whose molecular densities differ by more than
    :data:`MAX_DENSITY_RATIO` are reported as not comparable (``rmsd`` inf).
    """
    if not a.mols or not b.mols:
        return _Result(math.inf, 0, note="no molecules to compare")
    fwd = _rmsd_crystals(a, b, n_mol)  # first, so that composition problems are reported as such
    if fwd.note:
        return fwd
    ratio = _density(a) / _density(b)
    if max(ratio, 1.0 / ratio) > MAX_DENSITY_RATIO:
        return _Result(math.inf, 0, note=f"molecular densities differ by a factor {max(ratio, 1.0 / ratio):.2f}")
    rev = _rmsd_crystals(b, a, n_mol)
    if rev.note:
        return rev
    return _Result(max(fwd.rmsd, rev.rmsd), min(fwd.n_matched, rev.n_matched), fwd.mirror)


def rmsd_cluster(a: Atoms, b: Atoms, n_mol: int = N_CLUSTER) -> tuple[float, int]:
    """COMPACK-style cluster RMSD between two periodic molecular crystals.

    A cluster of ``n_mol`` molecules nearest (by centroid, periodic images included) to a centre
    molecule of ``a`` is compared with clusters built around the molecules of ``b`` (and of the mirror
    image of ``b``). The centre molecules are superposed (heavy atoms, proper rotation), and that overlay
    is combined with every (pseudo-)symmetry rotation of ``a``'s centre molecule, because a symmetric
    molecule can be overlaid in several ways of which only one fits the rest of the cluster. For each
    overlay the cluster molecules are paired by centroid distance and their atoms by a per-element
    assignment on the overlaid positions, so atom and molecule order do not matter. The best overlays
    are refined by a Kabsch fit over all matched heavy atoms. The cell choice of either structure does
    not matter.

    Every symmetry-distinct molecule of ``a`` is used as a centre in turn and the worst centre counts;
    clusters of ``b`` are built around its symmetry-distinct molecules (equivalent ones give congruent
    clusters). The comparison is then repeated with the roles of ``a`` and ``b`` exchanged and the worse
    result is returned, so that an extra or missing molecule anywhere in either cell prevents a match,
    whatever the molecule order. Structures whose molecular densities differ by more than
    :data:`MAX_DENSITY_RATIO` are not compared.

    Returns
    -------
    rmsd : float
        Heavy-atom RMSD (Å) over the matched cluster.
    n_matched : int
        The number of molecules of ``a``'s cluster whose centroid lies within 2.0 Å of its partner,
        less the number of molecules of ``b`` that are left without a partner although they lie inside
        the cluster (so a denser ``b`` that contains all of ``a``'s molecules is not a match).

    Both are for the best overlay (fully matched clusters first, then lowest RMSD), taking the worse
    of the two directions. ``(inf, 0)`` if the structures cannot be compared (the reason, e.g.
    different densities, is logged as a warning).
    """
    ca = _prepare(np.array(a.get_cell()), split_molecules(a, "a"))
    cb = _prepare(np.array(b.get_cell()), split_molecules(b, "b"))
    res = _compare(ca, cb, n_mol)
    if res.note:
        logger.warning("structures not comparable: %s", res.note)
    return res.rmsd, res.n_matched


@dataclass
class Comparison:
    """RMSD of one experimental structure against one candidate CIF.

    ``mirror`` is true when the best overlay used the mirror image of the candidate; ``note`` says why
    the two could not be compared (``rmsd`` is then infinite).
    """

    rank: int
    id: str
    origin: str
    rmsd: float
    n_matched: int
    mirror: bool = False
    note: str | None = None


@dataclass
class ValidationEntry:
    """Result for one experimental structure.

    ``best_*`` describe the best-ranked candidate that matches (RMSD below the cut with every cluster
    molecule matched), among all candidate CIFs compared; ``None`` if none matches. ``success`` is true
    when that rank is within the top ``n_top``. ``flip_match_rank`` is the best-ranked candidate that
    matches only when the two ``flip`` elements are treated as one element (a molecule turned over so
    that O and S exchange places); it is never counted as success.
    """

    experimental: str
    best_rank: int | None
    best_rmsd: float | None
    best_origin: str | None
    flip_match_rank: int | None
    success: bool = False
    best_id: str | None = None
    best_mirror: bool | None = None
    flip_rmsd: float | None = None
    flip_skipped: str | None = None
    comparisons: list[Comparison] = field(default_factory=list)
    flip_comparisons: list[Comparison] = field(default_factory=list)


@dataclass
class ValidationReport:
    """Validation of a CSP run against experimental structures."""

    success: bool
    entries: list[ValidationEntry]
    n_top: int = 20
    rmsd_cut: float = 0.8
    n_mol: int = N_CLUSTER
    flip: tuple[str, str] | None = None
    n_compared: int = 0

    def to_markdown(self) -> str:
        """Human-readable summary (Markdown)."""

        def f(x: float | None) -> str:
            return "-" if x is None or not math.isfinite(x) else f"{x:.3f}"

        flip_txt = (
            "off"
            if self.flip is None
            else f"{self.flip[0]} and {self.flip[1]} treated as one element (not counted as success)"
        )
        lines = [
            "# CSP validation",
            "",
            f"- Result: **{'SUCCESS' if self.success else 'NOT FOUND'}** "
            f"(match = RMSD{self.n_mol} < {self.rmsd_cut} Å with {self.n_mol}/{self.n_mol} molecules; "
            f"success = a match within the top {self.n_top})",
            f"- Candidate CIFs compared: {self.n_compared}",
            f"- Flip check: {flip_txt}",
            "",
            "| experimental | success | best match rank | id | origin | RMSD (Å) | mirror | flip match rank "
            "| flip RMSD (Å) |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for e in self.entries:
            mirror = "-" if e.best_mirror is None else ("yes" if e.best_mirror else "no")
            flip_rank = f"skipped ({e.flip_skipped})" if e.flip_skipped else (e.flip_match_rank or "-")
            lines.append(
                f"| {e.experimental} | {'yes' if e.success else 'no'} | {e.best_rank or '-'} | {e.best_id or '-'} "
                f"| {e.best_origin or '-'} | {f(e.best_rmsd)} | {mirror} | {flip_rank} | {f(e.flip_rmsd)} |"
            )
        for e in self.entries:
            lines += [
                "",
                f"## {e.experimental}: closest candidates",
                "",
                "| rank | id | origin | RMSD (Å) | matched molecules | mirror |",
                "|---|---|---|---|---|---|",
            ]
            comparable = [c for c in e.comparisons if c.note is None]
            for c in sorted(comparable, key=lambda c: (-c.n_matched, c.rmsd, c.rank))[:10]:
                lines.append(
                    f"| {c.rank} | {c.id} | {c.origin} | {f(c.rmsd)} | {c.n_matched}/{self.n_mol} "
                    f"| {'yes' if c.mirror else 'no'} |"
                )
            failed = [c for c in e.comparisons if c.note is not None]
            if failed:
                lines += ["", f"Not comparable ({len(failed)}):", ""]
                lines += [f"- rank {c.rank} {c.id}: {c.note}" for c in failed]
        return "\n".join(lines) + "\n"


_Cand = tuple[int, str, str, "_Crystal | str"]


def _load_cif(cif: Path) -> _Crystal | str:
    """A candidate CIF as a prepared crystal, or an error message if it cannot be read."""
    try:
        atoms = read(str(cif), format="cif")
        assert isinstance(atoms, Atoms)
        return _prepare(np.array(atoms.get_cell()), split_molecules(atoms, str(cif)))
    except Exception as exc:  # a broken candidate file must not stop the validation
        logger.warning("%s: cannot be read: %s", cif.name, exc)
        return f"cannot be read: {exc}"


def _read_candidates(out_dir: Path) -> list[_Cand]:
    """Candidate CIFs in ``out_dir`` as (rank, id, origin, crystal or error message), ordered by rank.

    If ``candidates.csv`` exists, only the CIFs of its rows (``<rank:03d>_<id>.cif``) are used, so
    files left from another run are ignored. A row whose rank or id was already seen is skipped with
    a warning; a listed CIF that is missing (within the top :data:`MAX_CIF`) is reported as not
    comparable. Without ``candidates.csv`` every ``NNN_<id>.cif`` file is used (origin ``?``).
    """
    csv_path = out_dir / "candidates.csv"
    out: list[_Cand] = []
    if not csv_path.exists():
        for cif in sorted(out_dir.glob("[0-9][0-9][0-9]_*.cif")):
            out.append((int(cif.stem[:3]), cif.stem[4:], "?", _load_cif(cif)))
        return sorted(out, key=lambda t: t[0])

    seen_ranks: set[int] = set()
    seen_ids: set[str] = set()
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rank, cid = int(row["rank"]), row["id"]
            if rank in seen_ranks or cid in seen_ids:
                logger.warning("candidates.csv: duplicate row (rank %d, id %s) skipped", rank, cid)
                continue
            seen_ranks.add(rank)
            seen_ids.add(cid)
            cif = out_dir / f"{rank:03d}_{cid}.cif"
            if cif.exists():
                out.append((rank, cid, row["origin"], _load_cif(cif)))
            elif rank <= MAX_CIF:
                logger.warning("%s: listed in candidates.csv but missing", cif.name)
                out.append((rank, cid, row["origin"], "CIF file missing"))
    return sorted(out, key=lambda t: t[0])


def _compare_all(
    exp: _Crystal, cands: list[_Cand], n_mol: int, element_map: dict[str, str] | None = None
) -> list[Comparison]:
    out = []
    for rank, cid, origin, crystal in cands:
        if isinstance(crystal, str):
            out.append(Comparison(rank, cid, origin, math.inf, 0, note=crystal))
            continue
        if element_map:
            crystal = _with_symbols_mapped(crystal, element_map)
        res = _compare(exp, crystal, n_mol)
        if res.note:
            logger.warning("%03d_%s.cif: not comparable: %s", rank, cid, res.note)
        out.append(Comparison(rank, cid, origin, res.rmsd, res.n_matched, res.mirror, res.note))
    return out


def _is_match(c: Comparison, rmsd_cut: float, n_mol: int) -> bool:
    return c.rmsd < rmsd_cut and c.n_matched == n_mol


def _best_match(comps: list[Comparison], rmsd_cut: float, n_mol: int) -> Comparison | None:
    hits = [c for c in comps if _is_match(c, rmsd_cut, n_mol)]
    return min(hits, key=lambda c: c.rank) if hits else None


def validate_run(
    out_dir: Path,
    experimental: list[Path],
    n_top: int = 20,
    rmsd_cut: float = 0.8,
    flip: tuple[str, str] | None = ("O", "S"),
) -> ValidationReport:
    """Compare every candidate CIF in ``out_dir`` (the top <=50 written by the pipeline) with experiments.

    Parameters
    ----------
    out_dir : Path
        Output directory of :func:`~molecular_modeling.csp.pipeline.run_csp`.
    experimental : list of Path
        Experimental CIFs (ordered, full occupancy; symmetry operators are expanded).
    n_top : int
        A match within this rank counts as success.
    rmsd_cut : float
        RMSD15 threshold (Å) for a match; all 15 cluster molecules must also be matched.
    flip : (str, str) or None
        Two elements treated as one element (``flip[1]`` renamed to ``flip[0]`` in both structures) for
        the separate "flip match" check: a candidate that matches only this way has molecules turned
        over so that the two elements exchange places. Skipped (with a warning) if either element is
        absent from the experimental structure.

    Returns
    -------
    ValidationReport

    Raises
    ------
    FileNotFoundError
        If ``out_dir`` is not a directory, or holds neither ``candidates.csv`` nor candidate CIFs
        (``NNN_<id>.cif``). A ``candidates.csv`` without rows (a run that kept no candidate) is valid.
    """
    out_dir = Path(out_dir)
    if not out_dir.is_dir():
        raise FileNotFoundError(f"CSP output directory not found: {out_dir}")
    if not (out_dir / "candidates.csv").exists() and not any(out_dir.glob("[0-9][0-9][0-9]_*.cif")):
        raise FileNotFoundError(
            f"no candidates in {out_dir}: expected candidates.csv or NNN_<id>.cif files written by molecular-modeling-csp"
        )
    cands = _read_candidates(out_dir)
    entries = []
    for path in experimental:
        cell, mols = read_seed_molecules(Path(path))
        exp = _prepare(cell, mols)
        comps = _compare_all(exp, cands, N_CLUSTER)
        hit = _best_match(comps, rmsd_cut, N_CLUSTER)
        entry = ValidationEntry(
            experimental=Path(path).name,
            best_rank=hit.rank if hit else None,
            best_rmsd=hit.rmsd if hit else None,
            best_origin=hit.origin if hit else None,
            flip_match_rank=None,
            success=hit is not None and hit.rank <= n_top,
            best_id=hit.id if hit else None,
            best_mirror=hit.mirror if hit else None,
            comparisons=comps,
        )
        if flip is not None:
            present = {s for m in mols for s in m.get_chemical_symbols()}
            missing = [el for el in flip if el not in present]
            if missing:
                entry.flip_skipped = f"no {', '.join(missing)} in the experimental structure"
                logger.warning("%s: flip check skipped: %s", Path(path).name, entry.flip_skipped)
            else:
                emap = {flip[1]: flip[0]}
                entry.flip_comparisons = _compare_all(_with_symbols_mapped(exp, emap), cands, N_CLUSTER, emap)
                plain = {c.rank for c in comps if _is_match(c, rmsd_cut, N_CLUSTER)}
                fhit = _best_match([c for c in entry.flip_comparisons if c.rank not in plain], rmsd_cut, N_CLUSTER)
                if fhit is not None:
                    entry.flip_match_rank, entry.flip_rmsd = fhit.rank, fhit.rmsd
        entries.append(entry)
    return ValidationReport(
        success=any(e.success for e in entries),
        entries=entries,
        n_top=n_top,
        rmsd_cut=rmsd_cut,
        flip=flip,
        n_compared=len(cands),
    )


def main(argv: list[str] | None = None) -> int:
    """Command line entry point; prints the Markdown report and writes ``<out_dir>/validation.md``."""
    parser = argparse.ArgumentParser(
        prog="python -m molecular_modeling.csp.validate",
        description="Compare CSP candidates with experimental crystal structures (RMSD15)",
    )
    parser.add_argument("out_dir", type=Path, help="CSP output directory (candidates.csv and CIFs)")
    parser.add_argument("experimental", type=Path, nargs="+", help="experimental CIF file(s)")
    parser.add_argument("--n-top", type=int, default=20, help="success if a match is within this rank (default 20)")
    parser.add_argument("--rmsd-cut", type=float, default=0.8, help="RMSD15 threshold in Å (default 0.8)")
    parser.add_argument(
        "--flip", default="O,S", help="two elements treated as one for the flip-match check (default O,S)"
    )
    parser.add_argument("--no-flip", action="store_true", help="skip the flipped-match check")
    args = parser.parse_args(argv)
    flip: tuple[str, str] | None = None
    if not args.no_flip:
        parts = [p.strip() for p in args.flip.split(",")]
        if len(parts) != 2:
            parser.error("--flip takes two comma-separated elements, e.g. O,S")
        flip = (parts[0], parts[1])
    try:
        report = validate_run(args.out_dir, args.experimental, n_top=args.n_top, rmsd_cut=args.rmsd_cut, flip=flip)
    except FileNotFoundError as exc:
        parser.error(str(exc))
    md = report.to_markdown()
    print(md, end="")
    (args.out_dir / "validation.md").write_text(md, encoding="utf-8")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
