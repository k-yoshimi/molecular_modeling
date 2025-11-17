"""
Unit tests for render_voronoi_overlays helper functions.
"""
import numpy as np
import pytest

from molecular_modeling.source.voronoi_geom import Polyhedron
from molecular_modeling.source.render_voronoi_overlays import (
    _poly_faces_to_verts,
    _convex_hull_to_polyhedron,
)


class TestPolyFacesToVerts:
    """Test cases for _poly_faces_to_verts function."""

    def test_poly_faces_to_verts_basic(self, simple_polyhedron):
        """
        Test converting polyhedron faces to vertex arrays.
        """
        verts_list = _poly_faces_to_verts(simple_polyhedron)
        assert len(verts_list) == len(simple_polyhedron.faces)
        for i, face_verts in enumerate(verts_list):
            assert face_verts.shape == (3, 3)  # Triangle with 3 vertices
            # Check that vertices match the face indices
            face_indices = simple_polyhedron.faces[i]
            np.testing.assert_array_equal(
                face_verts, simple_polyhedron.verts[face_indices]
            )

    def test_poly_faces_to_verts_empty(self):
        """
        Test with empty polyhedron.
        """
        poly = Polyhedron(
            verts=np.zeros((0, 3), dtype=float),
            faces=np.zeros((0, 3), dtype=int),
        )
        verts_list = _poly_faces_to_verts(poly)
        assert len(verts_list) == 0


class TestConvexHullToPolyhedron:
    """Test cases for _convex_hull_to_polyhedron function."""

    def test_convex_hull_to_polyhedron_cube(self):
        """
        Test creating polyhedron from cube points.
        """
        # Unit cube vertices
        points = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
        ], dtype=float)
        poly = _convex_hull_to_polyhedron(points)
        assert poly is not None
        assert poly.verts.shape[0] == 8
        assert poly.faces.shape[1] == 3  # Triangular faces

    def test_convex_hull_to_polyhedron_tetrahedron(self):
        """
        Test creating polyhedron from tetrahedron points.
        """
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, np.sqrt(3)/2, 0.0],
            [0.5, np.sqrt(3)/6, np.sqrt(2/3)],
        ], dtype=float)
        poly = _convex_hull_to_polyhedron(points)
        assert poly is not None
        assert poly.verts.shape[0] == 4
        assert poly.faces.shape[0] >= 4  # At least 4 triangular faces

    def test_convex_hull_to_polyhedron_insufficient_points(self):
        """
        Test with insufficient points for convex hull.
        """
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ], dtype=float)
        # Should handle gracefully or raise appropriate error
        with pytest.raises((ValueError, RuntimeError)):
            _convex_hull_to_polyhedron(points)


