"""
Unit tests for write_structure module.
"""
import tempfile
import shutil
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms

from molecular_modeling.source.write_structure import write_structure


class TestWriteStructure:
    """Test cases for write_structure function."""

    def test_write_structure_basic(self, simple_atoms):
        """
        Test writing structure files to a directory.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            write_structure(simple_atoms, outdir=tmpdir, verbose=False)
            outdir = Path(tmpdir)

            # Check that files were created
            assert (outdir / "qe.in").exists()
            assert (outdir / "POSCAR").exists()
            assert (outdir / "INCAR").exists()
            assert (outdir / "KPOINTS").exists()
            assert (outdir / "POSCAR.vasp").exists()
            assert (outdir / "structure.cif").exists()

    def test_write_structure_creates_directory(self, simple_atoms):
        """
        Test that write_structure creates the output directory if it doesn't exist.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "new_dir"
            write_structure(simple_atoms, outdir=outdir, verbose=False)
            assert outdir.exists()
            assert (outdir / "structure.cif").exists()

    def test_write_structure_qe_in_content(self, simple_atoms):
        """
        Test that QE input file has expected content.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            write_structure(simple_atoms, outdir=tmpdir, verbose=False)
            qe_file = Path(tmpdir) / "qe.in"
            assert qe_file.exists()
            content = qe_file.read_text()
            assert "&CONTROL" in content or "calculation" in content.lower()

    def test_write_structure_incar_content(self, simple_atoms):
        """
        Test that INCAR file has expected content.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            write_structure(simple_atoms, outdir=tmpdir, verbose=False)
            incar_file = Path(tmpdir) / "INCAR"
            assert incar_file.exists()
            content = incar_file.read_text()
            assert "PREC" in content
            assert "ENCUT" in content

    def test_write_structure_kpoints_content(self, simple_atoms):
        """
        Test that KPOINTS file has expected content.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            write_structure(simple_atoms, outdir=tmpdir, verbose=False)
            kpoints_file = Path(tmpdir) / "KPOINTS"
            assert kpoints_file.exists()
            content = kpoints_file.read_text()
            assert "Gamma" in content or "gamma" in content.lower()

    def test_write_structure_cif_exists(self, simple_atoms):
        """
        Test that CIF file is created.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            write_structure(simple_atoms, outdir=tmpdir, verbose=False)
            cif_file = Path(tmpdir) / "structure.cif"
            assert cif_file.exists()
            content = cif_file.read_text()
            assert "data_" in content or "_cell" in content

    def test_write_structure_with_multiple_species(self):
        """
        Test writing structure with multiple chemical species.
        """
        atoms = Atoms(
            symbols=['C', 'H', 'O', 'N'],
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            cell=[[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            pbc=True,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            write_structure(atoms, outdir=tmpdir, verbose=False)
            outdir = Path(tmpdir)
            assert (outdir / "structure.cif").exists()
            assert (outdir / "POSCAR").exists()

    def test_write_structure_verbose(self, simple_atoms, capsys):
        """
        Test that verbose mode prints output.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            write_structure(simple_atoms, outdir=tmpdir, verbose=True)
            captured = capsys.readouterr()
            assert "Structure input" in captured.out or "written" in captured.out.lower()


