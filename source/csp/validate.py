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

from .gen_seed import read_seed_molecules, split_molecules
from .molecule import _assign, _kabsch, match_and_align

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


@dataclass
class _Crystal:
    """Whole molecules of a periodic structure, prepared for cluster building."""

    cell: np.ndarray
    mols: list[Atoms]
    coms: np.ndarray  # (n_mol, 3) centroids wrapped into the cell
    shifts: np.ndarray  # (n_mol, 3) Cartesian shift from the molecule's own coordinates to the wrapped centroid
    heavy_pos: list[np.ndarray]  # positions of fitted atoms (heavy atoms, or all if the molecule has none)
    heavy_sym: list[np.ndarray]


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


def _with_symbols_swapped(c: _Crystal, pair: tuple[str, str]) -> _Crystal:
    x, y = pair
    swap = {x: y, y: x}
    mols = [Atoms([swap.get(s, s) for s in m.get_chemical_symbols()], positions=m.get_positions()) for m in c.mols]
    return _prepare(c.cell, mols)


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


def _rmsd_crystals(a: _Crystal, b: _Crystal, n_mol: int) -> tuple[float, int]:
    a_idx, a_shift = _cluster(a, 0, n_mol)
    a_pos = [a.heavy_pos[k] + s for k, s in zip(a_idx, a_shift)]
    a_sym = [a.heavy_sym[k] for k in a_idx]
    a_com = np.array([a.mols[k].get_positions().mean(axis=0) + s for k, s in zip(a_idx, a_shift)])
    a_centre = Atoms(a.mols[0].get_chemical_symbols(), positions=a.mols[0].get_positions() + a_shift[0])

    best = (math.inf, 0)
    for bb in (b, _inverted(b)):
        for j in range(len(bb.mols)):
            b_idx, b_shift = _cluster(bb, j, POOL_FACTOR * n_mol)
            b_raw = [bb.heavy_pos[k] + s for k, s in zip(b_idx, b_shift)]
            b_com_raw = np.array([bb.mols[k].get_positions().mean(axis=0) + s for k, s in zip(b_idx, b_shift)])
            b_sym = [bb.heavy_sym[k] for k in b_idx]
            b_centre = Atoms(bb.mols[j].get_chemical_symbols(), positions=bb.mols[j].get_positions() + b_shift[0])
            try:
                al = match_and_align(b_centre, a_centre, heavy_only=True)
            except ValueError:
                continue
            rot, trans = al.rotation, al.translation
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
            tgt, mob, com_dist = res
            rmsd = float(np.sqrt(np.mean(np.sum((tgt - mob) ** 2, axis=1))))
            if rmsd < best[0]:
                best = (rmsd, int(np.sum(com_dist < MATCH_COM_DIST)))
    return best


def rmsd_cluster(a: Atoms, b: Atoms, n_mol: int = N_CLUSTER) -> tuple[float, int]:
    """COMPACK-style cluster RMSD between two periodic molecular crystals.

    A cluster of ``n_mol`` molecules nearest (by centroid, periodic images included) to molecule 0 of
    ``a`` is compared with clusters built around every molecule of ``b`` (and of the mirror image of
    ``b``). The centre molecules are superposed with :func:`match_and_align` (heavy atoms, proper
    rotation); the cluster molecules are paired by centroid distance and their atoms by a per-element
    assignment on the overlaid positions, so atom and molecule order do not matter. The overlay is then
    refined by a Kabsch fit over all matched heavy atoms. The cell choice of either structure does not
    matter.

    Returns
    -------
    (rmsd, n_matched)
        Heavy-atom RMSD (Å) over the matched cluster, minimised over the centre choices of ``b``, and
        the number of molecules of ``a``'s cluster whose centroid lies within 2.0 Å of its partner.
        ``(inf, 0)`` if no centre of ``b`` has the composition of ``a``'s centre molecule.
    """
    ca = _prepare(np.array(a.get_cell()), split_molecules(a, "a"))
    cb = _prepare(np.array(b.get_cell()), split_molecules(b, "b"))
    return _rmsd_crystals(ca, cb, n_mol)


@dataclass
class Comparison:
    """RMSD of one experimental structure against one candidate CIF."""

    rank: int
    id: str
    origin: str
    rmsd: float
    n_matched: int


@dataclass
class ValidationEntry:
    """Result for one experimental structure.

    ``best_*`` describe the best-ranked candidate that matches (RMSD below the cut with every cluster
    molecule matched), among all candidate CIFs compared; ``None`` if none matches. ``success`` is true
    when that rank is within the top ``n_top``. ``flip_match_rank`` is the best-ranked match after
    swapping the ``flip`` elements in the experimental structure (never counted as success).
    """

    experimental: str
    best_rank: int | None
    best_rmsd: float | None
    best_origin: str | None
    flip_match_rank: int | None
    success: bool = False
    best_id: str | None = None
    flip_rmsd: float | None = None
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

        lines = [
            "# CSP validation",
            "",
            f"- Result: **{'SUCCESS' if self.success else 'NOT FOUND'}** "
            f"(match = RMSD{self.n_mol} < {self.rmsd_cut} Å with {self.n_mol}/{self.n_mol} molecules; "
            f"success = a match within the top {self.n_top})",
            f"- Candidate CIFs compared: {self.n_compared}",
            f"- Flip check: {'off' if self.flip is None else f'{self.flip[0]} <-> {self.flip[1]} (not counted as success)'}",
            "",
            "| experimental | success | best match rank | id | origin | RMSD (Å) | flip match rank | flip RMSD (Å) |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for e in self.entries:
            lines.append(
                f"| {e.experimental} | {'yes' if e.success else 'no'} | {e.best_rank or '-'} | {e.best_id or '-'} "
                f"| {e.best_origin or '-'} | {f(e.best_rmsd)} | {e.flip_match_rank or '-'} | {f(e.flip_rmsd)} |"
            )
        for e in self.entries:
            lines += [
                "",
                f"## {e.experimental}: closest candidates",
                "",
                "| rank | id | origin | RMSD (Å) | matched molecules |",
                "|---|---|---|---|---|",
            ]
            for c in sorted(e.comparisons, key=lambda c: (-c.n_matched, c.rmsd, c.rank))[:10]:
                lines.append(f"| {c.rank} | {c.id} | {c.origin} | {f(c.rmsd)} | {c.n_matched}/{self.n_mol} |")
        return "\n".join(lines) + "\n"


def _read_candidates(out_dir: Path) -> list[tuple[int, str, str, _Crystal]]:
    """Candidate CIFs in ``out_dir`` as (rank, id, origin, crystal), ordered by rank."""
    origins: dict[str, str] = {}
    csv_path = out_dir / "candidates.csv"
    if csv_path.exists():
        with open(csv_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                origins[row["id"]] = row["origin"]
    out = []
    for cif in sorted(out_dir.glob("[0-9][0-9][0-9]_*.cif")):
        rank, cid = int(cif.stem[:3]), cif.stem[4:]
        try:
            atoms = read(str(cif), format="cif")
            assert isinstance(atoms, Atoms)
            crystal = _prepare(np.array(atoms.get_cell()), split_molecules(atoms, str(cif)))
        except Exception as exc:  # a broken candidate file must not stop the validation
            logger.warning("skipping %s: %s", cif, exc)
            continue
        out.append((rank, cid, origins.get(cid, "?"), crystal))
    return sorted(out, key=lambda t: t[0])


def _compare_all(exp: _Crystal, cands: list[tuple[int, str, str, _Crystal]], n_mol: int) -> list[Comparison]:
    out = []
    for rank, cid, origin, crystal in cands:
        rmsd, n = _rmsd_crystals(exp, crystal, n_mol)
        out.append(Comparison(rank, cid, origin, rmsd, n))
    return out


def _best_match(comps: list[Comparison], rmsd_cut: float, n_mol: int) -> Comparison | None:
    hits = [c for c in comps if c.rmsd < rmsd_cut and c.n_matched == n_mol]
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
        Elements swapped in the experimental structure for the separate "flipped match" check.

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
            comparisons=comps,
        )
        present = {s for m in mols for s in m.get_chemical_symbols()}
        if flip is not None and present & set(flip):
            entry.flip_comparisons = _compare_all(_with_symbols_swapped(exp, flip), cands, N_CLUSTER)
            fhit = _best_match(entry.flip_comparisons, rmsd_cut, N_CLUSTER)
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
    parser.add_argument("--flip", default="O,S", help="two elements swapped for the flipped-match check (default O,S)")
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
