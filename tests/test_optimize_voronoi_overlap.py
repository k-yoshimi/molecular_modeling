"""
Unit tests for optimize_voronoi_overlap module.
"""
import numpy as np
import pytest
from ase import Atoms
from molecular_modeling.build_molecule_groups import build_molecule_groups
from molecular_modeling.optimize_voronoi_overlap import optimize_voronoi_overlap


class TestOptimizeVoronoiOverlap:
    """Test cases for optimize_voronoi_overlap function."""

    def test_optimize_voronoi_overlap_basic(self, simple_atoms):
        """
        Test basic optimization with a simple system.
        """
        # Add required metadata
        simple_atoms.info['site_mol_name'] = ['test_mol']
        simple_atoms.info['dict_meta'] = {
            'test_mol': {
                'lattice_vectors': np.array([
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [0.0, 0.0, 10.0],
                ], dtype=float),
                'center_of_mass': np.array([[5.0, 5.0, 5.0]], dtype=float),
                'voronoi_ref_id': 0,
            },
        }
        simple_atoms.info['rng'] = {'optimize': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        result = optimize_voronoi_overlap(
            simple_atoms,
            groups=groups,
            grid_num=(3, 3, 3),  # Small grid for fast testing
            rng=1,
        )

        assert result is not None
        assert 'optimized_angles_deg' in result.info
        assert 'optimized_scales' in result.info
        assert 'optimized_iou' in result.info

    def test_optimize_voronoi_overlap_missing_metadata(self, simple_atoms):
        """
        Test that ValueError is raised when required metadata is missing.
        """
        groups = build_molecule_groups(simple_atoms)
        with pytest.raises(ValueError, match="system.info must contain"):
            optimize_voronoi_overlap(simple_atoms, groups=groups)

    def test_optimize_voronoi_overlap_with_initial_angles(self, simple_atoms):
        """
        Test optimization with initial Euler angles.
        """
        simple_atoms.info['site_mol_name'] = ['test_mol']
        simple_atoms.info['dict_meta'] = {
            'test_mol': {
                'lattice_vectors': np.array([
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [0.0, 0.0, 10.0],
                ], dtype=float),
                'center_of_mass': np.array([[5.0, 5.0, 5.0]], dtype=float),
                'voronoi_ref_id': 0,
            },
        }
        simple_atoms.info['rng'] = {'optimize': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        x0 = np.array([[10.0, 20.0, 30.0]], dtype=float)
        result = optimize_voronoi_overlap(
            simple_atoms,
            groups=groups,
            x0_euler_deg=x0,
            grid_num=(3, 3, 3),
            rng=1,
        )

        assert result is not None
        assert 'optimized_angles_deg' in result.info
        angles = result.info['optimized_angles_deg']
        assert angles.shape == (1, 3)

    def test_optimize_voronoi_overlap_wrong_angle_shape(self, simple_atoms):
        """
        Test that ValueError is raised when x0_euler_deg has wrong shape.
        """
        simple_atoms.info['site_mol_name'] = ['test_mol']
        simple_atoms.info['dict_meta'] = {
            'test_mol': {
                'lattice_vectors': np.array([
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [0.0, 0.0, 10.0],
                ], dtype=float),
                'center_of_mass': np.array([[5.0, 5.0, 5.0]], dtype=float),
                'voronoi_ref_id': 0,
            },
        }
        simple_atoms.info['rng'] = {'optimize': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        x0 = np.array([[10.0, 20.0]], dtype=float)  # Wrong shape
        with pytest.raises(ValueError, match="x0_euler_deg must have shape"):
            optimize_voronoi_overlap(
                simple_atoms,
                groups=groups,
                x0_euler_deg=x0,
                grid_num=(3, 3, 3),
                rng=1,
            )

    def test_optimize_groups_none_fallback(self, simple_atoms):
        """Test that groups=None falls back to build_molecule_groups."""
        simple_atoms.info['site_mol_name'] = ['test_mol']
        simple_atoms.info['dict_meta'] = {
            'test_mol': {
                'lattice_vectors': np.array([
                    [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0],
                ], dtype=float),
                'center_of_mass': np.array([[5.0, 5.0, 5.0]], dtype=float),
                'voronoi_ref_id': 0,
            },
        }
        simple_atoms.info['rng'] = {'optimize': 1, 'global': 1}

        result = optimize_voronoi_overlap(
            simple_atoms, groups=None, grid_num=(3, 3, 3), rng=1,
        )
        assert result is not None
        assert 'optimized_angles_deg' in result.info

    def test_optimize_rng_none_fallback(self, simple_atoms):
        """Test that rng=None falls back to system.info['rng']."""
        simple_atoms.info['site_mol_name'] = ['test_mol']
        simple_atoms.info['dict_meta'] = {
            'test_mol': {
                'lattice_vectors': np.array([
                    [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0],
                ], dtype=float),
                'center_of_mass': np.array([[5.0, 5.0, 5.0]], dtype=float),
                'voronoi_ref_id': 0,
            },
        }
        simple_atoms.info['rng'] = {'optimize': 2, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        result = optimize_voronoi_overlap(
            simple_atoms, groups=groups, grid_num=(3, 3, 3), rng=None,
        )
        assert result is not None

    def test_optimize_site_angles_from_info(self, simple_atoms):
        """Test that x0_euler_deg falls back to system.info['site_angles_deg']."""
        simple_atoms.info['site_mol_name'] = ['test_mol']
        simple_atoms.info['dict_meta'] = {
            'test_mol': {
                'lattice_vectors': np.array([
                    [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0],
                ], dtype=float),
                'center_of_mass': np.array([[5.0, 5.0, 5.0]], dtype=float),
                'voronoi_ref_id': 0,
            },
        }
        simple_atoms.info['rng'] = {'optimize': 1, 'global': 1}
        simple_atoms.info['site_angles_deg'] = [[10.0, 20.0, 30.0]]

        groups = build_molecule_groups(simple_atoms)
        result = optimize_voronoi_overlap(
            simple_atoms, groups=groups, grid_num=(3, 3, 3), rng=1,
        )
        assert result is not None
        assert 'optimized_angles_deg' in result.info

    def test_optimize_log_dir(self, simple_atoms):
        """Test that log_dir writes optimization results to files."""
        import tempfile
        from pathlib import Path

        simple_atoms.info['site_mol_name'] = ['test_mol']
        simple_atoms.info['dict_meta'] = {
            'test_mol': {
                'lattice_vectors': np.array([
                    [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0],
                ], dtype=float),
                'center_of_mass': np.array([[5.0, 5.0, 5.0]], dtype=float),
                'voronoi_ref_id': 0,
            },
        }
        simple_atoms.info['rng'] = {'optimize': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            result = optimize_voronoi_overlap(
                simple_atoms, groups=groups, grid_num=(3, 3, 3), rng=1,
                log_dir=tmpdir,
            )
            res_file = Path(tmpdir) / "site_0" / "res.txt"
            assert res_file.exists()
            content = res_file.read_text()
            assert "fx =" in content
            assert "alpha =" in content
            assert "scale =" in content

    def test_optimize_voronoi_overlap_select_sites(self, simple_atoms):
        """
        Test optimization with site selection.
        """
        # Create system with 2 sites
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
        atoms.info['site_mol_name'] = ['test_mol', 'test_mol']
        atoms.info['dict_meta'] = {
            'test_mol': {
                'lattice_vectors': np.array([
                    [10.0, 0.0, 0.0],
                    [0.0, 10.0, 0.0],
                    [0.0, 0.0, 10.0],
                ], dtype=float),
                'center_of_mass': np.array([[5.0, 5.0, 5.0]], dtype=float),
                'voronoi_ref_id': 0,
            },
        }
        atoms.info['rng'] = {'optimize': 1, 'global': 1}

        groups = build_molecule_groups(atoms)
        result = optimize_voronoi_overlap(
            atoms,
            groups=groups,
            select_sites=[0],  # Only optimize site 0
            grid_num=(3, 3, 3),
            rng=1,
        )

        assert result is not None
        assert 'optimized_angles_deg' in result.info
        angles = result.info['optimized_angles_deg']
        assert angles.shape == (2, 3)


