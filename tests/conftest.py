"""
Pytest configuration and shared fixtures.
"""
import numpy as np
import pytest
from ase import Atoms


@pytest.fixture
def simple_cell():
    """
    Simple cubic cell for testing.

    Returns
    -------
    numpy.ndarray
        Cell vectors (3, 3) in Cartesian coordinates.
    """
    return np.array([
        [10.0, 0.0, 0.0],
        [0.0, 10.0, 0.0],
        [0.0, 0.0, 10.0],
    ], dtype=float)


@pytest.fixture
def simple_atoms(simple_cell):
    """
    Simple Atoms object for testing.

    Parameters
    ----------
    simple_cell : numpy.ndarray
        Cell vectors.

    Returns
    -------
    ase.Atoms
        Atoms object with 2 atoms.
    """
    atoms = Atoms(
        symbols=['C', 'H'],
        positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        cell=simple_cell,
        pbc=True,
    )
    atoms.new_array('mol_id', np.array([0, 0], dtype=int))
    atoms.info['groups'] = [[0, 1]]
    atoms.info['site_mol_name'] = ['test_mol']
    return atoms


@pytest.fixture
def simple_polyhedron():
    """
    Simple cube polyhedron for testing.

    Returns
    -------
    Polyhedron
        Cube polyhedron.
    """
    from molecular_modeling.voronoi_geom import Polyhedron
    # Unit cube vertices
    verts = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float)
    # Cube faces (triangulated)
    faces = np.array([
        [0, 1, 2], [0, 2, 3],  # bottom
        [4, 7, 6], [4, 6, 5],  # top
        [0, 4, 5], [0, 5, 1],  # front
        [2, 6, 7], [2, 7, 3],  # back
        [0, 3, 7], [0, 7, 4],  # left
        [1, 5, 6], [1, 6, 2],  # right
    ], dtype=int)
    return Polyhedron(verts=verts, faces=faces)


