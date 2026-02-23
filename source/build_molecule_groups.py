from typing import List
import numpy as np
from ase.atoms import Atoms

def build_molecule_groups(system: Atoms) -> List[List[int]]:
    """
    Build atom index groups keyed by the `mol_id` array.

    Parameters
    ----------
    system : ase.Atoms
        System holding `system.arrays['mol_id']`, where each value marks the
        placement/site identifier for an atom.

    Returns
    -------
    list[list[int]]
        List of atom index lists, ordered by `mol_id` (0..K-1).

    Raises
    ------
    ValueError
        If `mol_id` is missing from `system.arrays`.
    """
    if "mol_id" not in system.arrays:
        raise ValueError("system has no 'mol_id' array. Please set it in build_system_from_toml.")
    mol_id = system.get_array("mol_id")
    K = int(mol_id.max()) + 1 if mol_id.size else 0
    return [np.nonzero(mol_id == k)[0].tolist() for k in range(K)]
