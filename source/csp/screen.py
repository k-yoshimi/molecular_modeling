"""Cheap prescreening of crystal candidates: overlap, packing fraction and rigid LJ energy."""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.data import atomic_masses, chemical_symbols, vdw_radii
from ase.neighborlist import neighbor_list
from molecular_modeling.csp.candidate import Candidate
from molecular_modeling.csp.config import CSPConfig
from molecular_modeling.molecule_unwrap import is_connected, unwrap_molecule

KCAL_PER_MOL_TO_EV = 0.0433641
AMU_PER_A3_TO_G_PER_CM3 = 1.66054
DEFAULT_VDW_RADIUS = 1.7

#: UFF Lennard-Jones parameters (x_i in Angstrom, D_i in kcal/mol), Rappe et al., JACS 114, 10024 (1992).
UFF_LJ: dict[str, tuple[float, float]] = {
    "H": (2.886, 0.044),
    "C": (3.851, 0.105),
    "N": (3.660, 0.069),
    "O": (3.500, 0.060),
    "F": (3.364, 0.050),
    "S": (4.035, 0.274),
    "Cl": (3.947, 0.227),
    "Se": (4.205, 0.291),
    "Br": (4.189, 0.251),
    "I": (4.500, 0.339),
}


def _mol_ids(atoms: Atoms) -> np.ndarray:
    return atoms.arrays["mol_id"]


def _whole_molecules(atoms: Atoms) -> Atoms:
    """Copy of ``atoms`` in which every ``mol_id`` molecule is unwrapped into one connected piece.

    After this, two atoms of the same molecule at the same lattice shift (``S == 0`` in a neighbour
    list) are the intramolecular pair, whether or not the input had the molecule wrapped across the
    cell boundary.

    Raises
    ------
    ValueError
        If a molecule is not connected even after unwrapping.
    """
    out = atoms.copy()
    mol = _mol_ids(atoms)
    cell = np.asarray(atoms.cell)
    symbols = atoms.get_chemical_symbols()
    pos = atoms.get_positions()
    for m in np.unique(mol):
        idx = np.nonzero(mol == m)[0]
        sym = [symbols[i] for i in idx]
        whole = unwrap_molecule(pos[idx], sym, cell, atoms.pbc)
        if not is_connected(whole, sym):
            raise ValueError(f"molecule {m} is not connected")
        pos[idx] = whole
    out.set_positions(pos)
    return out


def has_intermolecular_overlap(atoms: Atoms, factor: float = 0.7) -> bool:
    """Return True if atoms of different molecules are closer than ``factor`` x (sum of vdW radii).

    A molecule and its own periodic image count as different molecules.
    Periodic boundaries are respected, and the result does not depend on whether molecules are
    wrapped into the cell. Radii come from ``ase.data.vdw_radii`` (1.7 A when undefined). A
    molecule that is not one connected piece counts as an overlap (the candidate is unusable).
    """
    try:
        atoms = _whole_molecules(atoms)
    except ValueError:
        return True
    radii = vdw_radii[atoms.numbers].astype(float)
    radii = np.where(np.isnan(radii), DEFAULT_VDW_RADIUS, radii)
    # per-atom cutoff r_i * factor: the pair cutoff is then factor * (r_i + r_j)
    i, j, s, d = neighbor_list("ijSd", atoms, radii * factor)
    if len(i) == 0:
        return False
    mol = _mol_ids(atoms)
    inter = (mol[i] != mol[j]) | np.any(s != 0, axis=1)
    return bool(np.any(inter & (d < factor * (radii[i] + radii[j]))))


def packing_fraction(atoms: Atoms, mol_volume: float) -> float:
    """Return ``n_mol * mol_volume / cell_volume``."""
    n_mol = len(np.unique(_mol_ids(atoms)))
    return float(n_mol * mol_volume / atoms.get_volume())


def rigid_lj_energy(atoms: Atoms, cutoff: float = 10.0) -> float:
    """Return the intermolecular UFF Lennard-Jones energy per molecule in eV.

    Sums over unique atom pairs of different molecules (a molecule and its own
    periodic image are different) within ``cutoff``, using geometric-mean mixing
    of the UFF parameters.

    Raises
    ------
    KeyError
        If an element has no UFF parameters.
    ValueError
        If a molecule is not one connected piece.
    """
    x = np.empty(len(atoms))
    dep = np.empty(len(atoms))
    for k, z in enumerate(atoms.numbers):
        sym = chemical_symbols[z]
        if sym not in UFF_LJ:
            raise KeyError(f"no UFF Lennard-Jones parameters for element {sym!r}")
        x[k], dep[k] = UFF_LJ[sym]
    atoms = _whole_molecules(atoms)
    i, j, s, d = neighbor_list("ijSd", atoms, cutoff)
    mol = _mol_ids(atoms)
    inter = (mol[i] != mol[j]) | np.any(s != 0, axis=1)
    i, j, d = i[inter], j[inter], d[inter]
    xij = np.sqrt(x[i] * x[j])
    dij = np.sqrt(dep[i] * dep[j])
    t = (xij / d) ** 6
    e = float(np.sum(dij * (t * t - 2.0 * t))) * KCAL_PER_MOL_TO_EV
    # neighbor_list lists each pair in both directions
    return e / 2.0 / len(np.unique(mol))


def crystal_density(atoms: Atoms) -> float:
    return float(atomic_masses[atoms.numbers].sum() * AMU_PER_A3_TO_G_PER_CM3 / atoms.get_volume())


def prescreen(cands: list[Candidate], mol_volume: float, cfg: CSPConfig) -> tuple[list[Candidate], dict[str, int]]:
    """Filter candidates by overlap, then packing fraction, then rank by rigid LJ energy.

    Sets ``meta["lj"]``, ``meta["packing"]`` and ``density`` on candidates that
    pass the first two stages; survivors get ``status = "prescreened"``.

    Returns
    -------
    kept : list of Candidate
        Up to ``cfg.n_prescreen_keep`` candidates, lowest LJ energy first.
    counts : dict
        Number remaining after each stage: ``generated``, ``no_overlap``,
        ``packing``, ``kept``.
    """
    counts = {"generated": len(cands), "no_overlap": 0, "packing": 0, "kept": 0}
    survivors = [c for c in cands if not has_intermolecular_overlap(c.atoms, cfg.overlap_factor)]
    counts["no_overlap"] = len(survivors)

    lo, hi = cfg.packing_range
    packed: list[Candidate] = []
    for c in survivors:
        pf = packing_fraction(c.atoms, mol_volume)
        c.meta["packing"] = pf
        c.density = crystal_density(c.atoms)
        if lo <= pf <= hi:
            packed.append(c)
    counts["packing"] = len(packed)

    for c in packed:
        c.meta["lj"] = rigid_lj_energy(c.atoms)
    kept = sorted(packed, key=lambda c: c.meta["lj"])[: cfg.n_prescreen_keep]
    for c in kept:
        c.status = "prescreened"
    counts["kept"] = len(kept)
    return kept, counts
