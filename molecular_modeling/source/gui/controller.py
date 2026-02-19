# source/gui/controller.py
"""
Controller for GUI-backend integration.

Provides high-level operations that bridge GUI dialogs with the molecular modeling backend.
"""

import re
from typing import List, Dict, Optional, Tuple
import numpy as np
from ase import Atoms
from scipy.spatial.transform import Rotation as R

from .molecule_registry import get_registry, MoleculeRegistry


class MolecularController:
    """
    Controller for molecular modeling operations.

    Bridges GUI dialogs with the molecular modeling backend operations.

    Parameters
    ----------
    atoms : ase.Atoms, optional
        Initial atoms object to operate on.
    """

    def __init__(self, atoms: Optional[Atoms] = None):
        self.atoms = atoms
        self.registry = get_registry()
        self._supercell_repeat = [1, 1, 1]
        self._use_fractional = False

    def set_atoms(self, atoms: Atoms) -> None:
        """Set the current atoms object."""
        self.atoms = atoms

    def get_atoms(self) -> Optional[Atoms]:
        """Get the current atoms object."""
        return self.atoms

    # =========================================================================
    # Molecule Addition
    # =========================================================================

    def add_molecule(
        self,
        molecule_name: str,
        position: Tuple[float, float, float],
        rotation: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> bool:
        """
        Add a molecule to the current system.

        Parameters
        ----------
        molecule_name : str
            Name of the molecule to add.
        position : tuple[float, float, float]
            Center of mass position in Cartesian coordinates.
        rotation : tuple[float, float, float], optional
            Euler angles (XYZ, degrees) for rotation.

        Returns
        -------
        bool
            True if molecule was added successfully.
        """
        config = self.registry.get_molecule(molecule_name)
        if config is None:
            return False

        # Get molecule data
        symbols = self.registry.get_symbols(molecule_name)
        positions = self.registry.get_positions(molecule_name)

        if len(symbols) == 0 or len(positions) == 0:
            return False

        # Recenter to origin
        com = positions.mean(axis=0)
        positions_centered = positions - com

        # Apply rotation
        rot = R.from_euler('xyz', rotation, degrees=True)
        positions_rotated = rot.apply(positions_centered)

        # Translate to target position
        positions_final = positions_rotated + np.array(position)

        # Create new atoms object for the molecule
        new_atoms = Atoms(symbols=symbols, positions=positions_final)

        # Add to system
        if self.atoms is None:
            self.atoms = new_atoms
        else:
            self.atoms += new_atoms

        # Update groups info
        self._update_groups_after_add(molecule_name, len(symbols))

        return True

    def _update_groups_after_add(self, molecule_name: str, n_atoms: int) -> None:
        """Update groups metadata after adding a molecule."""
        if self.atoms is None:
            return

        total_atoms = len(self.atoms)

        # Get or create groups list
        groups = self.atoms.info.get("groups", [])
        site_mol_name = self.atoms.info.get("site_mol_name", [])

        # New group is the last n_atoms
        start_idx = total_atoms - n_atoms
        new_group = list(range(start_idx, total_atoms))
        groups.append(new_group)
        site_mol_name.append(molecule_name)

        self.atoms.info["groups"] = groups
        self.atoms.info["site_mol_name"] = site_mol_name

        # Rebuild mol_id array from groups to ensure consistency
        mol_ids = np.zeros(total_atoms, dtype=int)
        for group_idx, group in enumerate(groups):
            for atom_idx in group:
                if atom_idx < total_atoms:
                    mol_ids[atom_idx] = group_idx

        # Delete existing array if present, then create new one
        if "mol_id" in self.atoms.arrays:
            del self.atoms.arrays["mol_id"]
        self.atoms.new_array("mol_id", mol_ids)

    # =========================================================================
    # Molecule Substitution
    # =========================================================================

    def substitute_molecule(self, molecule_index: int, new_molecule_name: str) -> bool:
        """
        Replace a molecule at a specific index with a different molecule.

        Preserves the orientation of the original molecule using principal component analysis.

        Parameters
        ----------
        molecule_index : int
            Index of the molecule group to replace.
        new_molecule_name : str
            Name of the new molecule to insert.

        Returns
        -------
        bool
            True if substitution was successful.
        """
        if self.atoms is None:
            return False

        groups = self.atoms.info.get("groups", [])
        site_mol_name = self.atoms.info.get("site_mol_name", [])

        if molecule_index < 0 or molecule_index >= len(groups):
            return False

        old_group = groups[molecule_index]
        if not old_group:
            return False

        # Get new molecule data
        new_config = self.registry.get_molecule(new_molecule_name)
        if new_config is None:
            return False

        new_symbols = self.registry.get_symbols(new_molecule_name)
        new_positions = self.registry.get_positions(new_molecule_name)

        if len(new_symbols) == 0 or len(new_positions) == 0:
            return False

        # Calculate center of mass of old molecule
        old_positions = self.atoms.get_positions()[old_group]
        old_com = old_positions.mean(axis=0)

        # Center both molecules at origin
        old_centered = old_positions - old_com
        new_com = new_positions.mean(axis=0)
        new_centered = new_positions - new_com

        # Calculate principal axes using SVD
        old_axes = self._get_principal_axes(old_centered)
        new_axes = self._get_principal_axes(new_centered)

        # Calculate rotation matrix to align new molecule with old orientation
        rotation_matrix = old_axes @ new_axes.T

        # Apply rotation and translate to old COM position
        new_positions_rotated = (new_centered @ rotation_matrix.T) + old_com

        # Get current data
        all_symbols = self.atoms.get_chemical_symbols()
        all_positions = self.atoms.get_positions().tolist()
        cell = self.atoms.get_cell()
        pbc = self.atoms.get_pbc()

        # Remove old atoms (in reverse order to preserve indices)
        old_group_sorted = sorted(old_group, reverse=True)
        for idx in old_group_sorted:
            del all_symbols[idx]
            del all_positions[idx]

        # Calculate insertion point (where old group started)
        insert_idx = min(old_group)

        # Insert new atoms at the same position with preserved orientation
        for i, (sym, pos) in enumerate(zip(new_symbols, new_positions_rotated)):
            all_symbols.insert(insert_idx + i, sym)
            all_positions.insert(insert_idx + i, pos.tolist())

        # Rebuild atoms object
        from ase import Atoms as AseAtoms
        new_atoms = AseAtoms(
            symbols=all_symbols,
            positions=all_positions,
            cell=cell,
            pbc=pbc,
        )

        # Update groups: recalculate all indices
        old_count = len(old_group)
        new_count = len(new_symbols)
        diff = new_count - old_count

        new_groups = []
        new_site_mol_name = []

        for i, group in enumerate(groups):
            if i == molecule_index:
                # This is the replaced molecule
                new_group = list(range(insert_idx, insert_idx + new_count))
                new_groups.append(new_group)
                new_site_mol_name.append(new_molecule_name)
            else:
                # Adjust indices for other groups
                adjusted_group = []
                for idx in group:
                    if idx < insert_idx:
                        adjusted_group.append(idx)
                    elif idx >= insert_idx + old_count:
                        adjusted_group.append(idx + diff)
                    # Indices within old_group are removed, skip them
                new_groups.append(adjusted_group)
                if i < len(site_mol_name):
                    new_site_mol_name.append(site_mol_name[i])
                else:
                    new_site_mol_name.append(f"mol_{i}")

        new_atoms.info["groups"] = new_groups
        new_atoms.info["site_mol_name"] = new_site_mol_name

        # Rebuild mol_id array
        total_atoms = len(new_atoms)
        mol_ids = np.zeros(total_atoms, dtype=int)
        for group_idx, group in enumerate(new_groups):
            for atom_idx in group:
                if 0 <= atom_idx < total_atoms:
                    mol_ids[atom_idx] = group_idx
        new_atoms.new_array("mol_id", mol_ids)

        # Update controller's atoms reference
        self.atoms = new_atoms

        return True

    @staticmethod
    def _get_principal_axes(positions: np.ndarray) -> np.ndarray:
        """
        Calculate principal axes of a molecule using SVD.

        Parameters
        ----------
        positions : np.ndarray
            Atom positions centered at origin (N x 3).

        Returns
        -------
        np.ndarray
            Principal axes as a 3x3 rotation matrix (rows are axes).
        """
        if len(positions) < 2:
            # For single atom or empty, return identity
            return np.eye(3)

        # Use SVD to find principal axes
        _, _, Vt = np.linalg.svd(positions, full_matrices=False)

        # Ensure right-handed coordinate system
        if np.linalg.det(Vt) < 0:
            Vt[2, :] *= -1

        return Vt

    # =========================================================================
    # Atom Substitution
    # =========================================================================

    def global_atom_substitution(self, old_symbol: str, new_symbol: str) -> int:
        """
        Replace all atoms of one element with another.

        Parameters
        ----------
        old_symbol : str
            Element symbol to replace.
        new_symbol : str
            New element symbol.

        Returns
        -------
        int
            Number of atoms replaced.
        """
        if self.atoms is None:
            return 0

        symbols = self.atoms.get_chemical_symbols()
        count = 0

        for i, sym in enumerate(symbols):
            if sym == old_symbol:
                symbols[i] = new_symbol
                count += 1

        if count > 0:
            self.atoms.set_chemical_symbols(symbols)

        return count

    def substitute_atom(self, atom_index: int, new_symbol: str) -> bool:
        """
        Replace a single atom with a new element.

        Parameters
        ----------
        atom_index : int
            Index of the atom to replace.
        new_symbol : str
            New element symbol.

        Returns
        -------
        bool
            True if substitution was successful.
        """
        if self.atoms is None:
            return False

        if atom_index < 0 or atom_index >= len(self.atoms):
            return False

        symbols = self.atoms.get_chemical_symbols()
        symbols[atom_index] = new_symbol
        self.atoms.set_chemical_symbols(symbols)

        return True

    def substitute_molecule_atoms(
        self,
        molecule_index: int,
        new_symbols: List[str],
    ) -> Tuple[bool, str]:
        """
        Replace all atoms in a molecule with new symbols.

        Parameters
        ----------
        molecule_index : int
            Index of the molecule group.
        new_symbols : List[str]
            New element symbols (can be compact notation like "C6H6").

        Returns
        -------
        tuple[bool, str]
            (success, error_message). error_message is empty if success.
        """
        if self.atoms is None:
            return False, "No atoms loaded"

        groups = self.atoms.info.get("groups", [])
        if not groups:
            return False, "No molecule groups found. Add molecules using 'Add Molecule' first."

        if molecule_index < 0 or molecule_index >= len(groups):
            return False, f"Invalid molecule index {molecule_index}. Valid range: 0-{len(groups)-1}"

        group = groups[molecule_index]
        expanded_symbols = self._expand_symbols(new_symbols)

        if len(expanded_symbols) != len(group):
            return False, f"Atom count mismatch: molecule has {len(group)} atoms, but '{' '.join(new_symbols)}' expands to {len(expanded_symbols)} atoms"

        symbols = self.atoms.get_chemical_symbols()
        for i, atom_idx in enumerate(group):
            symbols[atom_idx] = expanded_symbols[i]
        self.atoms.set_chemical_symbols(symbols)

        return True, ""

    @staticmethod
    def _expand_symbols(tokens: List[str]) -> List[str]:
        """Expand compact chemical tokens into elemental symbols."""
        pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
        result: List[str] = []
        for token in tokens:
            for element, count_str in pattern.findall(token):
                count = int(count_str) if count_str else 1
                result.extend([element] * count)
        return result

    # =========================================================================
    # Cell Operations
    # =========================================================================

    def set_supercell_repeat(self, repeats: List[int]) -> None:
        """Set the supercell repeat factors."""
        self._supercell_repeat = list(repeats)

    def get_supercell_repeat(self) -> List[int]:
        """Get the current supercell repeat factors."""
        return self._supercell_repeat.copy()

    def set_use_fractional(self, use_fractional: bool) -> None:
        """Set whether to use fractional coordinates."""
        self._use_fractional = use_fractional

    def get_use_fractional(self) -> bool:
        """Get whether fractional coordinates are used."""
        return self._use_fractional

    # =========================================================================
    # Molecule Information
    # =========================================================================

    def get_molecule_groups(self) -> List[Dict]:
        """
        Get information about all molecule groups in the system.

        Returns
        -------
        list[dict]
            List of dictionaries with molecule group information.
        """
        if self.atoms is None:
            return []

        groups = self.atoms.info.get("groups", [])
        site_mol_name = self.atoms.info.get("site_mol_name", [])
        positions = self.atoms.get_positions()
        symbols = self.atoms.get_chemical_symbols()

        result = []
        for i, group in enumerate(groups):
            if not group:
                continue

            mol_name = site_mol_name[i] if i < len(site_mol_name) else f"mol_{i}"
            mol_symbols = [symbols[idx] for idx in group]
            mol_positions = positions[group]
            com = mol_positions.mean(axis=0)

            result.append({
                "index": i,
                "name": mol_name,
                "atom_count": len(group),
                "atom_indices": group,
                "elements": sorted(set(mol_symbols)),
                "center_of_mass": com.tolist(),
            })

        return result

    def get_selected_molecule_index(self, selected_atoms: List[int]) -> Optional[int]:
        """
        Get the molecule index for a set of selected atoms.

        Parameters
        ----------
        selected_atoms : list[int]
            List of selected atom indices.

        Returns
        -------
        int or None
            Molecule index if all selected atoms belong to the same molecule.
        """
        if self.atoms is None or not selected_atoms:
            return None

        if "mol_id" not in self.atoms.arrays:
            return None

        mol_ids = self.atoms.get_array("mol_id")
        selected_mol_ids = set(mol_ids[i] for i in selected_atoms)

        if len(selected_mol_ids) == 1:
            return selected_mol_ids.pop()

        return None


# Global controller instance
_controller: Optional[MolecularController] = None


def get_controller() -> MolecularController:
    """
    Get the global controller instance.

    Returns
    -------
    MolecularController
        The global controller instance.
    """
    global _controller
    if _controller is None:
        _controller = MolecularController()
    return _controller


def reset_controller() -> MolecularController:
    """
    Reset and return a fresh global controller instance.

    Use this when starting a new GUI session to avoid stale state.

    Returns
    -------
    MolecularController
        A fresh controller instance.
    """
    global _controller
    _controller = MolecularController()
    return _controller
