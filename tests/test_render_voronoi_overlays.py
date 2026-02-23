"""
Unit tests for render_voronoi_overlays module.
"""
import tempfile
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms

from molecular_modeling.render_voronoi_overlays import render_voronoi_overlays
from molecular_modeling.build_molecule_groups import build_molecule_groups


class TestRenderVoronoiOverlays:
    """Test cases for render_voronoi_overlays function."""

    def test_render_voronoi_overlays_basic(self, simple_atoms):
        """
        Test basic rendering with a simple system.
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            render_voronoi_overlays(
                simple_atoms,
                groups=groups,
                out_dir=tmpdir,
                rng=1,
            )
            outdir = Path(tmpdir)
            # Check that PNG files were created (at least one)
            png_files = list(outdir.glob("*.png"))
            assert len(png_files) > 0

    def test_render_voronoi_overlays_creates_directory(self, simple_atoms):
        """
        Test that render_voronoi_overlays creates the output directory.
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "new_render_dir"
            render_voronoi_overlays(
                simple_atoms,
                groups=groups,
                out_dir=outdir,
                rng=1,
            )
            assert outdir.exists()

    def test_render_rng_none_fallback(self, simple_atoms):
        """Test rng=None falls back to system.info['rng']."""
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            render_voronoi_overlays(
                simple_atoms, groups=groups, out_dir=tmpdir, rng=None,
            )
            assert len(list(Path(tmpdir).glob("*.png"))) > 0

    def test_render_groups_none_fallback(self, simple_atoms):
        """Test groups=None falls back to build_molecule_groups."""
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}

        with tempfile.TemporaryDirectory() as tmpdir:
            render_voronoi_overlays(
                simple_atoms, groups=None, out_dir=tmpdir, rng=1,
            )
            assert len(list(Path(tmpdir).glob("*.png"))) > 0

    def test_render_missing_metadata_raises(self, simple_atoms):
        """Test that missing site_mol_name/dict_meta raises ValueError."""
        # Remove required metadata
        simple_atoms.info.pop('site_mol_name', None)
        simple_atoms.info.pop('dict_meta', None)
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValueError, match="system.info must contain"):
                render_voronoi_overlays(
                    simple_atoms, groups=groups, out_dir=tmpdir, rng=1,
                )

    def test_render_angles_from_system_info(self, simple_atoms):
        """Test angles_deg falls back to system.info['optimized_angles_deg']."""
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}
        simple_atoms.info['optimized_angles_deg'] = np.array([[10.0, 20.0, 30.0]])

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            render_voronoi_overlays(
                simple_atoms, groups=groups, out_dir=tmpdir, rng=1,
                angles_deg=None,
            )
            assert len(list(Path(tmpdir).glob("*.png"))) > 0

    def test_render_angles_wrong_shape_raises(self, simple_atoms):
        """Test that angles_deg with wrong shape raises ValueError."""
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValueError, match="angles_deg must have shape"):
                render_voronoi_overlays(
                    simple_atoms, groups=groups, out_dir=tmpdir, rng=1,
                    angles_deg=np.array([[1.0, 2.0]]),  # wrong shape
                )

    def test_render_optimized_scales_from_info(self, simple_atoms):
        """Test scales loaded from system.info['optimized_scales']."""
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}
        simple_atoms.info['optimized_scales'] = np.array([1.5])

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            render_voronoi_overlays(
                simple_atoms, groups=groups, out_dir=tmpdir, rng=1,
            )
            assert len(list(Path(tmpdir).glob("*.png"))) > 0

    def test_render_optimized_scales_wrong_shape(self, simple_atoms):
        """Test that wrong-shaped scales fallback to ones."""
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}
        simple_atoms.info['optimized_scales'] = np.array([1.5, 2.0])  # wrong shape (N=1)

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            render_voronoi_overlays(
                simple_atoms, groups=groups, out_dir=tmpdir, rng=1,
            )
            assert len(list(Path(tmpdir).glob("*.png"))) > 0

    def test_render_explicit_sites(self, simple_atoms):
        """Test sites parameter converts to int list."""
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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            render_voronoi_overlays(
                simple_atoms, groups=groups, out_dir=tmpdir, rng=1,
                sites=[0],
            )
            assert len(list(Path(tmpdir).glob("*.png"))) > 0

    def test_render_voronoi_overlays_multiple_sites(self):
        """
        Test rendering with multiple sites.
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
        atoms.info['rng'] = {'render': 1, 'global': 1}

        groups = build_molecule_groups(atoms)
        with tempfile.TemporaryDirectory() as tmpdir:
            render_voronoi_overlays(
                atoms,
                groups=groups,
                out_dir=tmpdir,
                rng=1,
            )
            outdir = Path(tmpdir)
            png_files = list(outdir.glob("*.png"))
            assert len(png_files) >= 2  # At least one per site

    def test_render_set_box_aspect_exception(self, simple_atoms):
        """Test that set_box_aspect exception is silently caught (lines 110-111)."""
        from unittest.mock import patch

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
        simple_atoms.info['rng'] = {'render': 1, 'global': 1}

        groups = build_molecule_groups(simple_atoms)

        # Wrap the real _set_axes_equal_from_points to force the exception path
        from molecular_modeling import render_voronoi_overlays as rvo
        original_fn = rvo._set_axes_equal_from_points

        def patched_fn(ax, pts):
            """Call original but make set_box_aspect raise."""
            original_set = ax.set_box_aspect
            ax.set_box_aspect = lambda *a, **kw: (_ for _ in ()).throw(
                AttributeError("no set_box_aspect"))
            try:
                original_fn(ax, pts)
            finally:
                ax.set_box_aspect = original_set

        with patch.object(rvo, '_set_axes_equal_from_points', side_effect=patched_fn):
            with tempfile.TemporaryDirectory() as tmpdir:
                render_voronoi_overlays(
                    simple_atoms, groups=groups, out_dir=tmpdir, rng=1,
                )
                assert len(list(Path(tmpdir).glob("*.png"))) > 0


