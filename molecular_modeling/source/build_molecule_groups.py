from typing import List
import numpy as np
from ase.atoms import Atoms

def build_molecule_groups(system: Atoms) -> List[List[int]]:
    """Group atoms by system.arrays['mol_id'] (placement order preserved)."""
    if "mol_id" not in system.arrays:
        raise ValueError("system has no 'mol_id' array. Please set it in build_system_from_toml.")
    mol_id = system.get_array("mol_id")
    K = int(mol_id.max()) + 1 if mol_id.size else 0
    return [np.nonzero(mol_id == k)[0].tolist() for k in range(K)]