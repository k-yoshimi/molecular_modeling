"""Duplicate removal for relaxed candidate structures."""
from __future__ import annotations

import warnings
from dataclasses import replace

from .candidate import Candidate


def deduplicate(
    cands: list[Candidate], ltol: float = 0.2, stol: float = 0.3, angle_tol: float = 5.0
) -> list[Candidate]:
    """Merge candidates that pymatgen's ``StructureMatcher`` considers the same crystal.

    Candidates are visited from the lowest energy up; each one is compared with the
    representatives kept so far (O(n^2) comparisons). A matching candidate is dropped and its
    ``origin`` appended to the representative's ``meta["duplicates"]``.

    Parameters
    ----------
    cands : list of Candidate
        Relaxed candidates. Those without energy are treated as highest in energy.
    ltol, stol, angle_tol : float
        ``StructureMatcher`` tolerances.

    Returns
    -------
    list of Candidate
        Representatives, lowest energy first. Their ``meta`` is copied before ``duplicates`` is added.
    """
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor

    matcher = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol)
    order = sorted(cands, key=lambda c: float("inf") if c.energy_per_mol is None else c.energy_per_mol)
    kept: list[Candidate] = []
    structures: list[Structure] = []
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", module=r"pymatgen\..*")
        for cand in order:
            struct = AseAtomsAdaptor.get_structure(cand.atoms)
            for i, other in enumerate(structures):
                if matcher.fit(struct, other):
                    dups = list(kept[i].meta.get("duplicates", []))
                    dups.append(cand.origin)
                    kept[i] = replace(kept[i], meta={**kept[i].meta, "duplicates": dups})
                    break
            else:
                kept.append(cand)
                structures.append(struct)
    return kept
