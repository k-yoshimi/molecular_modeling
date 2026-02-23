"""
Unit tests for voronoi_geom module.
"""
import numpy as np
import pytest

from molecular_modeling.voronoi_geom import (
    Polyhedron,
    wrap_angles_deg,
    polyhedron_volume,
    rotate_polyhedron,
    _triangulate_faces,
    _build_unit_cell_box_mesh,
)


class TestPolyhedron:
    """Test cases for Polyhedron class."""

    def test_polyhedron_creation(self):
        """
        Test creating a Polyhedron.
        """
        verts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=int)
        poly = Polyhedron(verts=verts, faces=faces)
        assert poly.verts.shape == (4, 3)
        assert poly.faces.shape == (4, 3)
        np.testing.assert_array_equal(poly.verts, verts)
        np.testing.assert_array_equal(poly.faces, faces)


class TestWrapAnglesDeg:
    """Test cases for wrap_angles_deg function."""

    def test_wrap_positive_angles(self):
        """
        Test wrapping positive angles.
        """
        angles = np.array([200.0, 370.0, 90.0])
        wrapped = wrap_angles_deg(angles)
        expected = np.array([-160.0, 10.0, 90.0])
        np.testing.assert_array_almost_equal(wrapped, expected)

    def test_wrap_negative_angles(self):
        """
        Test wrapping negative angles.
        """
        angles = np.array([-200.0, -370.0, -90.0])
        wrapped = wrap_angles_deg(angles)
        expected = np.array([160.0, -10.0, -90.0])
        np.testing.assert_array_almost_equal(wrapped, expected)

    def test_wrap_180_degree(self):
        """
        Test that 180 degrees wraps to -180.
        """
        angles = np.array([180.0, -180.0])
        wrapped = wrap_angles_deg(angles)
        expected = np.array([-180.0, -180.0])
        np.testing.assert_array_almost_equal(wrapped, expected)

    def test_wrap_in_range(self):
        """
        Test that angles already in range are unchanged.
        """
        angles = np.array([-90.0, 0.0, 90.0])
        wrapped = wrap_angles_deg(angles)
        np.testing.assert_array_almost_equal(wrapped, angles)

    def test_wrap_scalar(self):
        """
        Test wrapping a scalar angle.
        """
        angle = 200.0
        wrapped = wrap_angles_deg(np.array([angle]))
        assert abs(wrapped[0] - (-160.0)) < 1e-10


class TestPolyhedronVolume:
    """Test cases for polyhedron_volume function."""

    def test_unit_cube_volume(self, simple_polyhedron):
        """
        Test volume calculation for a unit cube.
        """
        volume = polyhedron_volume(simple_polyhedron)
        assert abs(volume - 1.0) < 1e-6

    def test_empty_polyhedron(self):
        """
        Test volume calculation for empty polyhedron.
        """
        poly = Polyhedron(
            verts=np.zeros((0, 3), dtype=float),
            faces=np.zeros((0, 3), dtype=int),
        )
        volume = polyhedron_volume(poly)
        assert volume == 0.0

    def test_tetrahedron_volume(self):
        """
        Test volume calculation for a tetrahedron.
        """
        # Regular tetrahedron with edge length sqrt(2)
        verts = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, np.sqrt(3)/2, 0.0],
            [0.5, np.sqrt(3)/6, np.sqrt(2/3)],
        ], dtype=float)
        faces = np.array([
            [0, 1, 2],
            [0, 1, 3],
            [0, 2, 3],
            [1, 2, 3],
        ], dtype=int)
        poly = Polyhedron(verts=verts, faces=faces)
        volume = polyhedron_volume(poly)
        # Expected volume for regular tetrahedron with edge sqrt(2)
        expected = np.sqrt(2) / 12.0
        assert abs(volume - expected) < 1e-6


class TestRotatePolyhedron:
    """Test cases for rotate_polyhedron function."""

    def test_rotate_about_origin(self, simple_polyhedron):
        """
        Test rotating a polyhedron about the origin.
        """
        center = np.array([0.5, 0.5, 0.5], dtype=float)
        angles = np.array([90.0, 0.0, 0.0])  # 90 degrees around x-axis
        rotated = rotate_polyhedron(simple_polyhedron, angles, center)
        assert rotated.verts.shape == simple_polyhedron.verts.shape
        assert rotated.faces.shape == simple_polyhedron.faces.shape

    def test_rotate_about_center(self, simple_polyhedron):
        """
        Test rotating a polyhedron about its center.
        """
        center = np.array([0.5, 0.5, 0.5], dtype=float)
        angles = np.array([0.0, 90.0, 0.0])  # 90 degrees around y-axis
        rotated = rotate_polyhedron(simple_polyhedron, angles, center)
        # Center should remain unchanged
        rotated_center = rotated.verts.mean(axis=0)
        np.testing.assert_array_almost_equal(rotated_center, center, decimal=6)

    def test_no_rotation(self, simple_polyhedron):
        """
        Test that zero rotation preserves the polyhedron.
        """
        center = np.array([0.5, 0.5, 0.5], dtype=float)
        angles = np.array([0.0, 0.0, 0.0])
        rotated = rotate_polyhedron(simple_polyhedron, angles, center)
        np.testing.assert_array_almost_equal(
            rotated.verts, simple_polyhedron.verts, decimal=10
        )


class TestTriangulateFaces:
    """Test cases for _triangulate_faces function."""

    def test_already_triangulated(self):
        """
        Test that already triangulated faces are unchanged.
        """
        faces = np.array([[0, 1, 2], [0, 2, 3], [1, 2, 3]], dtype=int)
        result = _triangulate_faces(faces)
        np.testing.assert_array_equal(result, faces)

    def test_quad_to_triangles(self):
        """
        Test converting a quad to triangles.
        """
        faces = np.array([[0, 1, 2, 3]], dtype=int)
        result = _triangulate_faces(faces)
        assert result.shape == (2, 3)
        assert result[0, 0] == 0
        assert result[1, 0] == 0

    def test_empty_faces(self):
        """
        Test handling empty faces array.
        """
        faces = np.zeros((0, 3), dtype=int)
        result = _triangulate_faces(faces)
        assert result.shape == (0, 3)

    def test_pentagon_to_triangles(self):
        """
        Test converting a pentagon to triangles.
        """
        faces = np.array([[0, 1, 2, 3, 4]], dtype=int)
        result = _triangulate_faces(faces)
        assert result.shape == (3, 3)
        # All triangles should start with vertex 0
        assert np.all(result[:, 0] == 0)


class TestBuildUnitCellBoxMesh:
    """Test cases for _build_unit_cell_box_mesh function."""

    def test_cubic_cell(self, simple_cell):
        """
        Test building mesh for cubic cell.
        """
        mesh = _build_unit_cell_box_mesh(simple_cell)
        assert mesh is not None
        assert mesh.vertices.shape[0] == 8  # Cube has 8 vertices
        assert mesh.faces.shape[1] == 3  # Triangular faces

    def test_orthorhombic_cell(self):
        """
        Test building mesh for orthorhombic cell.
        """
        cell = np.array([
            [5.0, 0.0, 0.0],
            [0.0, 10.0, 0.0],
            [0.0, 0.0, 15.0],
        ], dtype=float)
        mesh = _build_unit_cell_box_mesh(cell)
        assert mesh is not None
        assert mesh.vertices.shape[0] == 8

    def test_triclinic_cell(self):
        """
        Test building mesh for triclinic cell.
        """
        cell = np.array([
            [5.0, 0.0, 0.0],
            [2.0, 4.0, 0.0],
            [1.0, 1.0, 3.0],
        ], dtype=float)
        mesh = _build_unit_cell_box_mesh(cell)
        assert mesh is not None
        assert mesh.vertices.shape[0] == 8


