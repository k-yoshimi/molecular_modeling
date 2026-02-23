"""
Unit tests for build_molecule_groups module.
"""
import numpy as np
import pytest
from ase import Atoms
from molecular_modeling.build_molecule_groups import build_molecule_groups


class TestBuildMoleculeGroups:
    """Test cases for build_molecule_groups function."""

    def test_single_molecule(self):
        """
        Test grouping with a single molecule.
        """
        atoms = Atoms(
            symbols=['C', 'H', 'H'],
            positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
            cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            pbc=True,
        )
        atoms.new_array('mol_id', np.array([0, 0, 0], dtype=int))
        groups = build_molecule_groups(atoms)
        assert len(groups) == 1
        assert groups[0] == [0, 1, 2]

    def test_multiple_molecules(self):
        """
        Test grouping with multiple molecules.
        """
        atoms = Atoms(
            symbols=['C', 'H', 'C', 'H'],
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [6.0, 0.0, 0.0],
            ],
            cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            pbc=True,
        )
        atoms.new_array('mol_id', np.array([0, 0, 1, 1], dtype=int))
        groups = build_molecule_groups(atoms)
        assert len(groups) == 2
        assert groups[0] == [0, 1]
        assert groups[1] == [2, 3]

    def test_non_sequential_mol_ids(self):
        """
        Test grouping with non-sequential molecule IDs.
        """
        atoms = Atoms(
            symbols=['C', 'H', 'C', 'H'],
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [6.0, 0.0, 0.0],
            ],
            cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            pbc=True,
        )
        atoms.new_array('mol_id', np.array([0, 0, 2, 2], dtype=int))
        groups = build_molecule_groups(atoms)
        assert len(groups) == 3
        assert groups[0] == [0, 1]
        assert groups[1] == []
        assert groups[2] == [2, 3]

    def test_empty_system(self):
        """
        Test grouping with empty system.
        """
        atoms = Atoms(
            symbols=[],
            positions=[],
            cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            pbc=True,
        )
        atoms.new_array('mol_id', np.array([], dtype=int))
        groups = build_molecule_groups(atoms)
        assert len(groups) == 0

    def test_missing_mol_id(self):
        """
        Test that ValueError is raised when mol_id array is missing.
        """
        atoms = Atoms(
            symbols=['C', 'H'],
            positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            pbc=True,
        )
        with pytest.raises(ValueError, match="system has no 'mol_id' array"):
            build_molecule_groups(atoms)


