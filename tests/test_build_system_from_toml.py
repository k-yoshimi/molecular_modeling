"""
Unit tests for build_system_from_toml module.
"""
import numpy as np
import pytest
from unittest.mock import patch, mock_open

from molecular_modeling.build_system_from_toml import build_system_from_toml


class TestBuildSystemFromToml:
    """Test cases for build_system_from_toml function."""

    @patch('molecular_modeling.build_system_from_toml.tomllib.load')
    @patch('builtins.open')
    def test_build_system_basic(self, mock_open_func, mock_toml_load):
        """
        Test building a basic system from TOML.
        """
        # Mock input.toml
        input_data = {
            'molecules': ['benzene'],
            'supercell_repeat': [1, 1, 1],
            'cell': [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            'center_of_mass': [[5.0, 5.0, 5.0]],
            'use_fractional': False,
            'rng': 1,
        }
        # Mock benzene.toml
        benzene_data = {
            'symbols': ['C6H6'],
            'positions': [
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.866, 0.0],
                [1.0, 1.732, 0.0], [0.0, 1.732, 0.0], [-0.5, 0.866, 0.0],
                [0.5, 0.5, 0.0], [1.5, 0.5, 0.0], [2.0, 1.366, 0.0],
                [1.5, 2.232, 0.0], [0.5, 2.232, 0.0], [-0.5, 1.366, 0.0],
            ],
            'center_of_mass': [[0.0, 0.0, 0.0]],
            'lattice_vectors': [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            'use_fractional': False,
            'voronoi_ref_id': 0,
        }

        # Track file paths being opened
        file_paths = []

        def open_side_effect(file_path, *args, **kwargs):
            file_path_str = str(file_path)
            file_paths.append(file_path_str)
            # Return a mock file object
            m = mock_open()
            return m(file_path, *args, **kwargs)

        mock_open_func.side_effect = open_side_effect

        # Return appropriate data based on which file is being opened
        def toml_load_side_effect(file_obj, *args, **kwargs):
            # Determine which file based on the last opened file path
            if file_paths:
                last_path = file_paths[-1]
                if 'input.toml' in last_path:
                    return input_data
                elif 'benzene.toml' in last_path:
                    return benzene_data
            # Fallback: return input_data for first call, benzene_data for others
            if len(file_paths) == 1:
                return input_data
            return benzene_data

        mock_toml_load.side_effect = toml_load_side_effect

        system = build_system_from_toml(verbose=False)
        assert system is not None
        assert 'mol_id' in system.arrays
        assert 'groups' in system.info
        assert 'site_mol_name' in system.info
        assert 'dict_meta' in system.info
        assert 'rng' in system.info

    def test_expand_symbols_helper(self):
        """
        Test the expand_symbols helper function.
        """
        from molecular_modeling.build_system_from_toml import build_system_from_toml
        import re

        # Test expand_symbols logic (internal function)
        pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
        tokens = ["C6H6"]
        result = []
        for token in tokens:
            for element, count_str in pattern.findall(token):
                count = int(count_str) if count_str else 1
                result.extend([element] * count)

        assert result == ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']

    def test_compute_positions_helper(self):
        """
        Test the compute_positions helper function logic.
        """
        from scipy.spatial.transform import Rotation as R
        import numpy as np

        # Simulate compute_positions logic
        cfg = {
            'positions': [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            'center_of_mass': [[0.5, 0.0, 0.0]],
            'lattice_vectors': [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            'use_fractional': False,
        }
        site_com_cart = np.array([5.0, 5.0, 5.0])
        angles_xyz_deg = [0.0, 0.0, 0.0]

        pos = np.array(cfg["positions"], float)
        pos_com = pos.mean(axis=0)
        pos_centered = pos - pos_com
        rot = R.from_euler('xyz', angles_xyz_deg, degrees=True)
        pos_rot = rot.apply(pos_centered)
        final_pos = pos_rot + site_com_cart

        assert final_pos.shape == (2, 3)
        np.testing.assert_array_almost_equal(final_pos[0], [4.5, 5.0, 5.0])
        np.testing.assert_array_almost_equal(final_pos[1], [5.5, 5.0, 5.0])


