# source/gui/molecule_registry.py
"""
Molecule registry for dynamic loading of molecule definitions from TOML files.

Scans the data/*.toml directory and provides access to molecule configurations.
"""

import re
import tomllib
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np


class MoleculeRegistry:
    """
    Registry for molecule definitions loaded from TOML files.

    Scans the data/ directory for .toml files and loads molecule configurations.

    Parameters
    ----------
    data_dir : str or Path, optional
        Path to the directory containing molecule TOML files.
        Defaults to "data" relative to the package root.

    Attributes
    ----------
    molecules : dict[str, dict]
        Dictionary mapping molecule names to their configurations.
    """

    def __init__(self, data_dir: Optional[str | Path] = None):
        if data_dir is None:
            # Default to data/ relative to the molecular_modeling root
            package_root = Path(__file__).parent.parent.parent
            data_dir = package_root / "data"
        self.data_dir = Path(data_dir)
        self.molecules: Dict[str, dict] = {}
        self._load_all()

    def _load_all(self) -> None:
        """Load all molecule TOML files from the data directory."""
        if not self.data_dir.exists():
            return

        for toml_path in sorted(self.data_dir.glob("*.toml")):
            try:
                with open(toml_path, "rb") as f:
                    config = tomllib.load(f)
                # Use filename (without extension) as key if molecule_name not present
                name = config.get("molecule_name", toml_path.stem)
                self.molecules[name] = config
            except Exception as e:
                print(f"Warning: Failed to load {toml_path}: {e}")

    def reload(self) -> None:
        """Reload all molecule definitions from disk."""
        self.molecules.clear()
        self._load_all()

    def get_molecule_names(self) -> List[str]:
        """
        Get list of available molecule names.

        Returns
        -------
        list[str]
            Sorted list of molecule names.
        """
        return sorted(self.molecules.keys())

    def get_molecule(self, name: str) -> Optional[dict]:
        """
        Get molecule configuration by name.

        Parameters
        ----------
        name : str
            Molecule name.

        Returns
        -------
        dict or None
            Molecule configuration, or None if not found.
        """
        return self.molecules.get(name)

    def get_symbols(self, name: str) -> List[str]:
        """
        Get expanded atomic symbols for a molecule.

        Parameters
        ----------
        name : str
            Molecule name.

        Returns
        -------
        list[str]
            List of atomic symbols (expanded from compact notation).
        """
        config = self.get_molecule(name)
        if config is None:
            return []
        return self._expand_symbols(config.get("symbols", []))

    def get_positions(self, name: str) -> np.ndarray:
        """
        Get atomic positions for a molecule.

        Parameters
        ----------
        name : str
            Molecule name.

        Returns
        -------
        numpy.ndarray
            Atomic positions as (N, 3) array.
        """
        config = self.get_molecule(name)
        if config is None:
            return np.array([])
        return np.array(config.get("positions", []), dtype=float)

    def get_num_atoms(self, name: str) -> int:
        """
        Get number of atoms in a molecule.

        Parameters
        ----------
        name : str
            Molecule name.

        Returns
        -------
        int
            Number of atoms.
        """
        return len(self.get_symbols(name))

    def get_unique_elements(self, name: str) -> List[str]:
        """
        Get unique elements in a molecule.

        Parameters
        ----------
        name : str
            Molecule name.

        Returns
        -------
        list[str]
            Sorted list of unique element symbols.
        """
        return sorted(set(self.get_symbols(name)))

    @staticmethod
    def _expand_symbols(tokens: List[str]) -> List[str]:
        """
        Expand compact chemical tokens into elemental symbols.

        Parameters
        ----------
        tokens : list[str]
            Strings such as "C6H6" or "Si6O6".

        Returns
        -------
        list[str]
            Flattened elemental symbols expanded by count.
        """
        pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
        result: List[str] = []
        for token in tokens:
            for element, count_str in pattern.findall(token):
                count = int(count_str) if count_str else 1
                result.extend([element] * count)
        return result


# Global registry instance (lazy initialization)
_registry: Optional[MoleculeRegistry] = None


def get_registry() -> MoleculeRegistry:
    """
    Get the global molecule registry instance.

    Returns
    -------
    MoleculeRegistry
        The global registry instance.
    """
    global _registry
    if _registry is None:
        _registry = MoleculeRegistry()
    return _registry
