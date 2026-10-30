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

from .gen_seed import read_seed_molecules, split_molecules
from .molecule import _assign, _initial_rotations, _kabsch, match_and_align

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


@dataclass
class _Crystal:
    """Whole molecules of a periodic structure, prepared for cluster building."""

    cell: np.ndarray
    mols: list[Atoms]
    coms: np.ndarray  # (n_mol, 3) centroids wrapped into the cell
    shifts: np.ndarray  # (n_mol, 3) Cartesian shift from the molecule's own coordinates to the wrapped centroid
    heavy_pos: list[np.ndarray]  # positions of fitted atoms (heavy atoms, or all if the molecule has none)
    heavy_sym: list[np.ndarray]
    self_rotations: list[np.ndarray] | None = None  # (pseudo-)symmetry rotations of molecule 0, computed lazily


def _prepare(cell: np.ndarray, mols: list[Atoms]) -> _Crystal:
    cell = np.asarray(cell, dtype=float)
    raw = np.array([m.get_positions().mean(axis=0) for m in mols])
    frac = raw @ np.linalg.inv(cell)
    wrapped = (frac - np.floor(frac)) @ cell
    heavy_pos, heavy_sym = [], []
    for m in mols:
        sym = np.array(m.get_chemical_symbols())
        keep = sym != "H" if np.any(sym != "H") else np.ones(len(sym), dtype=bool)
        heavy_pos.append(m.get_positions()[keep])
        heavy_sym.append(sym[keep])
    return _Crystal(cell, mols, wrapped, wrapped - raw, heavy_pos, heavy_sym)


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
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Pair cluster molecules by centroid, then atoms within each pair by position (per element).

    Returns ``(target, mobile, com_dist)`` (stacked matched atoms of a and b, centroid distances of the
    molecule pairs), or ``None`` if a pair has different compositions.
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
    return np.vstack(tgt), np.vstack(mob), cost[rows, cols]


@dataclass
class _Result:
    rmsd: float
    n_matched: int
    mirror: bool = False
    note: str | None = None  # why the structures could not be compared (rmsd is then inf)


def _formula(sym: np.ndarray) -> str:
    els, counts = np.unique(sym, return_counts=True)
    return "".join(f"{e}{n}" for e, n in zip(els, counts))


def _rmsd_crystals(a: _Crystal, b: _Crystal, n_mol: int) -> _Result:
    a_idx, a_shift = _cluster(a, 0, n_mol)
    a_pos = [a.heavy_pos[k] + s for k, s in zip(a_idx, a_shift)]
    a_sym = [a.heavy_sym[k] for k in a_idx]
    a_com = np.array([a.mols[k].get_positions().mean(axis=0) + s for k, s in zip(a_idx, a_shift)])
    a_centre = Atoms(a.mols[0].get_chemical_symbols(), positions=a.mols[0].get_positions() + a_shift[0])
    if a.self_rotations is None:
        a.self_rotations = _self_rotations(a.heavy_pos[0], a.heavy_sym[0])
    pivot = a_pos[0].mean(axis=0)  # the self-rotations act about the centroid of the fitted centre atoms

    def key(res: tuple[np.ndarray, np.ndarray, np.ndarray]) -> tuple[bool, float]:
        tgt, mob, com_dist = res
        rmsd = float(np.sqrt(np.mean(np.sum((tgt - mob) ** 2, axis=1))))
        return int(np.sum(com_dist < MATCH_COM_DIST)) < n_mol, rmsd

    pools = []
    hyps: list[tuple[tuple[bool, float], int, np.ndarray, np.ndarray]] = []
    notes: set[str] = set()
    for mirror, bb in ((False, b), (True, _inverted(b))):
        for j in range(len(bb.mols)):
            if sorted(bb.heavy_sym[j]) != sorted(a.heavy_sym[0]):
                notes.add(
                    f"molecule compositions differ ({_formula(a.heavy_sym[0])} vs {_formula(bb.heavy_sym[j])},"
                    " heavy atoms)"
                )
                continue
            b_idx, b_shift = _cluster(bb, j, POOL_FACTOR * n_mol)
            b_raw = [bb.heavy_pos[k] + s for k, s in zip(b_idx, b_shift)]
            b_com_raw = np.array([bb.mols[k].get_positions().mean(axis=0) + s for k, s in zip(b_idx, b_shift)])
            b_sym = [bb.heavy_sym[k] for k in b_idx]
            b_centre = Atoms(bb.mols[j].get_chemical_symbols(), positions=bb.mols[j].get_positions() + b_shift[0])
            try:
                al = match_and_align(b_centre, a_centre, heavy_only=True)
            except ValueError as exc:
                notes.add(f"centre molecules cannot be aligned ({exc})")
                continue
            pools.append((mirror, b_raw, b_sym, b_com_raw))
            # every (pseudo-)symmetry-equivalent overlay of the centre molecule is a separate hypothesis
            for s in a.self_rotations:
                rot = s @ al.rotation
                trans = s @ (al.translation - pivot) + pivot
                cur = _match_clusters(
                    a_pos, a_sym, a_com, [p @ rot.T + trans for p in b_raw], b_sym, b_com_raw @ rot.T + trans
                )
                if cur is None:
                    notes.add("cluster molecules have different compositions")
                    continue
                hyps.append((key(cur), len(pools) - 1, rot, trans))
    if not hyps:
        return _Result(math.inf, 0, note="; ".join(sorted(notes)) or "no molecules to compare")

    hyps.sort(key=lambda h: h[0])
    best: tuple[tuple[bool, float], _Result] | None = None
    for _k, pi, rot, trans in hyps[:N_KEEP_REFINE]:
        mirror, b_raw, b_sym, b_com_raw = pools[pi]
        res = None
        for it in range(N_REFINE + 1):
            cur = _match_clusters(
                a_pos, a_sym, a_com, [p @ rot.T + trans for p in b_raw], b_sym, b_com_raw @ rot.T + trans
            )
            if cur is None:
                break
            res = cur
            if it < N_REFINE:
                tgt, mob_moved, _ = cur
                # refit on the untransformed b coordinates of the matched atoms
                rot, trans = _kabsch((mob_moved - trans) @ rot, tgt)
        if res is None:
            continue
        k = key(res)
        if best is None or k < best[0]:
            best = (k, _Result(k[1], int(np.sum(res[2] < MATCH_COM_DIST)), mirror))
    assert best is not None
    return best[1]


def rmsd_cluster(a: Atoms, b: Atoms, n_mol: int = N_CLUSTER) -> tuple[float, int]:
    """COMPACK-style cluster RMSD between two periodic molecular crystals.

    A cluster of ``n_mol`` molecules nearest (by centroid, periodic images included) to molecule 0 of
    ``a`` is compared with clusters built around every molecule of ``b`` (and of the mirror image of
    ``b``). The centre molecules are superposed with :func:`match_and_align` (heavy atoms, proper
    rotation), and that overlay is combined with every (pseudo-)symmetry rotation of ``a``'s centre
    molecule, because a symmetric molecule can be overlaid in several ways of which only one fits the
    rest of the cluster. For each overlay the cluster molecules are paired by centroid distance and their
    atoms by a per-element assignment on the overlaid positions, so atom and molecule order do not
    matter. The best overlays are refined by a Kabsch fit over all matched heavy atoms. The cell choice
    of either structure does not matter.

    Returns
    -------
    (rmsd, n_matched)
        Heavy-atom RMSD (Å) over the matched cluster and the number of molecules of ``a``'s cluster
        whose centroid lies within 2.0 Å of its partner, for the best overlay (fully matched clusters
        first, then lowest RMSD). ``(inf, 0)`` if the structures cannot be compared (the reason is
        logged as a warning).
    """
    ca = _prepare(np.array(a.get_cell()), split_molecules(a, "a"))
    cb = _prepare(np.array(b.get_cell()), split_molecules(b, "b"))
    res = _rmsd_crystals(ca, cb, n_mol)
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


def _read_candidates(out_dir: Path) -> list[_Cand]:
    """Candidate CIFs in ``out_dir`` as (rank, id, origin, crystal or error message), ordered by rank."""
    origins: dict[str, str] = {}
    csv_path = out_dir / "candidates.csv"
    if csv_path.exists():
        with open(csv_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                origins[row["id"]] = row["origin"]
    out: list[_Cand] = []
    for cif in sorted(out_dir.glob("[0-9][0-9][0-9]_*.cif")):
        rank, cid = int(cif.stem[:3]), cif.stem[4:]
        crystal: _Crystal | str
        try:
            atoms = read(str(cif), format="cif")
            assert isinstance(atoms, Atoms)
            crystal = _prepare(np.array(atoms.get_cell()), split_molecules(atoms, str(cif)))
        except Exception as exc:  # a broken candidate file must not stop the validation
            logger.warning("%s: cannot be read: %s", cif.name, exc)
            crystal = f"cannot be read: {exc}"
        out.append((rank, cid, origins.get(cid, "?"), crystal))
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
        res = _rmsd_crystals(exp, crystal, n_mol)
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
    """
    out_dir = Path(out_dir)
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
    report = validate_run(args.out_dir, args.experimental, n_top=args.n_top, rmsd_cut=args.rmsd_cut, flip=flip)
    md = report.to_markdown()
    print(md, end="")
    (args.out_dir / "validation.md").write_text(md, encoding="utf-8")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
