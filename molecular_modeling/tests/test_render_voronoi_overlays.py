"""
Unit tests for render_voronoi_overlays module.
"""
import tempfile
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms

from molecular_modeling.source.render_voronoi_overlays import render_voronoi_overlays
from molecular_modeling.source.build_molecule_groups import build_molecule_groups


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


