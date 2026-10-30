"""Random molecular-crystal candidate generation with PyXtal."""
from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING, Any

import numpy as np
from ase import Atoms
from molecular_modeling.molecule_unwrap import is_connected

from .candidate import Candidate

if TYPE_CHECKING:  # pragma: no cover
    from pyxtal.molecule import pyxtal_molecule

logger = logging.getLogger(__name__)


def to_pyxtal_molecule(mol: Atoms, seed: int | None = None) -> pyxtal_molecule:
    """Convert an ASE molecule to a PyXtal molecule (via ``pymatgen.core.Molecule``).

    The geometry is kept as given (``symmetrize=False``): PyXtal would otherwise snap it to
    its detected point group, changing bond lengths slightly. ``seed`` fixes the molecule's
    internal random state for reproducibility.
    """
    from pymatgen.core import Molecule
    from pyxtal.molecule import pyxtal_molecule as _pm

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", module=r"pymatgen\..*")
        pmg = Molecule(mol.get_chemical_symbols(), mol.get_positions())
        return _pm(pmg, symmetrize=False, seed=seed, random_state=seed)


def default_z(space_group: int, z_prime: int) -> int:
    """General-position multiplicity of ``space_group`` times ``z_prime``."""
    from pyxtal.symmetry import Group

    return int(Group(space_group)[0].multiplicity) * int(z_prime)


def _to_atoms(crystal: Any, n_per_mol: int, z: int) -> Atoms | None:
    """Convert a PyXtal crystal to ASE atoms with ``mol_id``; None if molecules are not contiguous."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", module=r"pymatgen\..*")
        atoms = crystal.to_ase(resort=False)
    if len(atoms) != n_per_mol * z:
        return None
    pos = atoms.get_positions()
    sym = atoms.get_chemical_symbols()
    for k in range(z):
        sl = slice(k * n_per_mol, (k + 1) * n_per_mol)
        if not is_connected(pos[sl], sym[sl]):
            return None
    atoms.set_array("mol_id", np.repeat(np.arange(z), n_per_mol).astype(int))
    return atoms


def generate_random(
    mol: Atoms,
    space_group: int,
    z: int,
    n: int,
    rng: np.random.Generator,
    max_attempts_factor: int = 5,
) -> list[Candidate]:
    """Generate up to ``n`` random rigid-molecule crystals in ``space_group`` with ``z`` molecules per cell.

    Molecules are placed as rigid bodies (no scaling). PyXtal seeds are drawn from ``rng``
    so runs are reproducible. Failed attempts (invalid crystal, exception, or non-contiguous
    molecules) count towards the ``n * max_attempts_factor`` budget; if fewer than ``n``
    candidates result, those found are returned and a warning is logged.
    """
    from pyxtal import pyxtal

    n_per_mol = len(mol)
    out: list[Candidate] = []
    max_attempts = n * max_attempts_factor
    attempts = 0
    while len(out) < n and attempts < max_attempts:
        attempts += 1
        seed = int(rng.integers(0, 2**31 - 1))
        try:
            pmol = to_pyxtal_molecule(mol, seed=seed)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", module=r"pymatgen\..*")
                warnings.filterwarnings("ignore", module=r"pyxtal\..*")
                crystal = pyxtal(molecular=True)
                crystal.from_random(3, space_group, [pmol], [z], random_state=seed)
            if not crystal.valid:
                continue
            atoms = _to_atoms(crystal, n_per_mol, z)
        except Exception as exc:  # PyXtal raises many kinds of errors for impossible setups
            logger.debug("PyXtal failed for space group %d: %s", space_group, exc)
            continue
        if atoms is None:
            continue
        meta = {"pyxtal": crystal, "seed": seed}
        out.append(Candidate(atoms=atoms, origin="random", space_group=space_group, meta=meta))
    if len(out) < n:
        logger.warning(
            "space group %d: generated only %d of %d candidates in %d attempts", space_group, len(out), n, attempts
        )
    return out
