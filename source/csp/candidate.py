"""Candidate crystal structure container."""

from __future__ import annotations

from dataclasses import dataclass, field

from ase import Atoms

#: Number of top-ranked candidates written as CIF files (``001_<id>.cif`` ...).
MAX_CIF = 50


@dataclass
class Candidate:
    """A candidate crystal structure.

    Parameters
    ----------
    atoms : ase.Atoms
        Periodic structure. ``atoms.arrays["mol_id"]`` is required.
    origin : str
        ``"random"``, ``"voronoi"`` or ``"seed:<cif name>"``.
    space_group : int
        Space group number the candidate was generated in.
    meta : dict
        Free-form metadata.
    energy_per_mol : float or None
        Energy per molecule (eV) after relaxation.
    density : float or None
        Density in g/cm^3.
    status : str
        Pipeline status, ``"generated"`` initially.
    """

    atoms: Atoms
    origin: str
    space_group: int
    meta: dict = field(default_factory=dict)
    energy_per_mol: float | None = None
    density: float | None = None
    status: str = "generated"
