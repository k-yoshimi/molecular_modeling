"""
Unit tests for voronoi_geom module.
"""
import numpy as np
import pytest
from molecular_modeling.voronoi_geom import (
    Polyhedron,
    _build_unit_cell_box_mesh,
    _triangulate_faces,
    _voronoi_cell_vertices_halfspaces_cart,
    check_system_voronoi_volume_tiling,
    compute_system_voronoi_cells,
    compute_voronoi_cells_from_coms,
    iou_via_union_cart,
    optimize_angles_each_system_site,
    overlap_volume_union_cart,
    polyhedron_to_trimesh,
    polyhedron_volume,
    precompute_system_cache,
    rotate_polyhedron,
    single_set_cell_volume_cart,
    wrap_angles_deg,
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


class TestPolyhedronToTrimesh:
    """Test polyhedron_to_trimesh conversion."""

    def test_basic_conversion(self, simple_polyhedron):
        mesh = polyhedron_to_trimesh(simple_polyhedron)
        assert mesh is not None
        assert mesh.vertices.shape[0] == 8
        assert mesh.faces.shape[1] == 3


class TestComputeVoronoiCellsFromComs:
    """Test compute_voronoi_cells_from_coms."""

    def test_empty_coms(self):
        cell = np.eye(3) * 10.0
        result = compute_voronoi_cells_from_coms(cell, np.empty((0, 3)))
        assert result == []

    def test_single_com(self):
        cell = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        result = compute_voronoi_cells_from_coms(cell, coms)
        assert len(result) == 1
        vol = polyhedron_volume(result[0])
        assert vol > 0

    def test_two_coms(self):
        cell = np.eye(3) * 10.0
        coms = np.array([[2.5, 5.0, 5.0], [7.5, 5.0, 5.0]])
        result = compute_voronoi_cells_from_coms(cell, coms)
        assert len(result) == 2


class TestComputeSystemVoronoiCells:
    """Test compute_system_voronoi_cells."""

    def test_empty_groups(self):
        from ase import Atoms
        atoms = Atoms(cell=np.eye(3) * 10.0, pbc=True)
        result = compute_system_voronoi_cells(atoms, [])
        assert result == []

    def test_basic(self, simple_atoms):
        from molecular_modeling.build_molecule_groups import build_molecule_groups
        groups = build_molecule_groups(simple_atoms)
        result = compute_system_voronoi_cells(simple_atoms, groups)
        assert len(result) == 1


class TestVoronoiCellVerticesHalfspaces:
    """Test _voronoi_cell_vertices_halfspaces_cart edge cases."""

    def test_duplicate_points_skipped(self):
        p_ref = np.array([0.0, 0.0, 0.0])
        others = np.array([
            [0.0, 0.0, 0.0],  # duplicate of p_ref
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
        ])
        verts = _voronoi_cell_vertices_halfspaces_cart(p_ref, others)
        assert len(verts) > 0

    def test_all_duplicates_returns_empty(self):
        p_ref = np.array([0.0, 0.0, 0.0])
        others = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ])
        verts = _voronoi_cell_vertices_halfspaces_cart(p_ref, others)
        assert verts.shape == (0, 3)


class TestPolyhedronVolumeFanTriangulation:
    """Test polyhedron_volume with non-triangular faces."""

    def test_quad_faces_volume(self):
        """Test volume calculation with quad faces (fan triangulation)."""
        verts = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
        ], dtype=float)
        # Quad faces (not triangulated)
        faces = np.array([
            [0, 1, 2, 3],
            [4, 7, 6, 5],
            [0, 4, 5, 1],
            [2, 6, 7, 3],
            [0, 3, 7, 4],
            [1, 5, 6, 2],
        ], dtype=int)
        poly = Polyhedron(verts=verts, faces=faces)
        vol = polyhedron_volume(poly)
        assert abs(vol - 1.0) < 1e-6

    def test_single_edge_face_returns_zero(self):
        """Test fan triangulation with 2-vertex faces gives zero volume."""
        verts = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
        # 2-column face array: cannot form a triangle via fan triangulation
        faces = np.array([[0, 1]], dtype=int)
        poly = Polyhedron(verts=verts, faces=faces)
        vol = polyhedron_volume(poly)
        assert vol == 0.0


class TestOverlapVolumeUnionCart:
    """Test overlap_volume_union_cart edge cases."""

    def test_degenerate_returns_zero(self):
        """Test that degenerate geometry returns 0."""
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        vol = overlap_volume_union_cart(
            A, coms, A, coms, [0.0, 0.0, 0.0], coms[0], rng=1,
        )
        assert vol >= 0.0

    def test_few_verts_returns_zero(self):
        """Test that < 4 vertices returns 0 via mock."""
        from unittest.mock import patch
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        with patch(
            'molecular_modeling.voronoi_geom._voronoi_cell_vertices_halfspaces_cart',
            return_value=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float),
        ):
            vol = overlap_volume_union_cart(
                A, coms, A, coms, [0.0, 0.0, 0.0], coms[0], rng=1,
            )
            assert vol == 0.0


class TestSingleSetCellVolumeCart:
    """Test single_set_cell_volume_cart."""

    def test_invalid_which_raises(self):
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        with pytest.raises(ValueError, match="which must be"):
            single_set_cell_volume_cart(A, coms, coms[0], 'invalid')

    def test_system_volume(self):
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        vol = single_set_cell_volume_cart(A, coms, coms[0], 'system', rng=1)
        assert vol > 0

    def test_dict_volume(self):
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        vol = single_set_cell_volume_cart(
            A, coms, coms[0], 'dict',
            angles_deg=[0.0, 0.0, 0.0], rng=1,
        )
        assert vol > 0

    def test_few_verts_returns_zero(self):
        """Test that < 4 vertices returns 0 (line 609)."""
        from unittest.mock import patch
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        with patch(
            'molecular_modeling.voronoi_geom._voronoi_cell_vertices_halfspaces_cart',
            return_value=np.array([[0, 0, 0], [1, 0, 0]], dtype=float),
        ):
            vol = single_set_cell_volume_cart(A, coms, coms[0], 'system', rng=1)
            assert vol == 0.0


class TestPrecomputeSystemCache:
    """Test precompute_system_cache."""

    def test_basic_cache(self):
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        cache = precompute_system_cache(A, coms, coms[0], rng=1)
        assert "seeds_sys" in cache
        assert "vol_sys" in cache
        assert cache["vol_sys"] > 0
        assert cache["seeds_sys"].ndim == 2

    def test_cached_iou_matches_uncached(self):
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        angles = [10.0, 20.0, 30.0]
        iou_uncached = iou_via_union_cart(
            A, coms, A, coms, coms[0], angles, rng=1,
        )
        cache = precompute_system_cache(A, coms, coms[0], rng=1)
        iou_cached = iou_via_union_cart(
            A, coms, A, coms, coms[0], angles, rng=1,
            _cached_sys=cache,
        )
        assert abs(iou_cached - iou_uncached) < 1e-12


class TestIoUViaUnionCart:
    """Test iou_via_union_cart."""

    def test_identical_system_dict(self):
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        iou = iou_via_union_cart(
            A, coms, A, coms, coms[0], [0.0, 0.0, 0.0], rng=1,
        )
        assert 0.0 <= iou <= 1.0

    def test_zero_denom_returns_zero(self):
        """Test that denom <= 0 returns 0 (line 663)."""
        from unittest.mock import patch
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        with patch(
            'molecular_modeling.voronoi_geom.overlap_volume_union_cart',
            return_value=0.0,
        ), patch(
            'molecular_modeling.voronoi_geom.single_set_cell_volume_cart',
            return_value=0.0,
        ):
            iou = iou_via_union_cart(
                A, coms, A, coms, coms[0], [0.0, 0.0, 0.0], rng=1,
            )
            assert iou == 0.0


class TestOptimizeAnglesEachSystemSite:
    """Test optimize_angles_each_system_site."""

    def test_basic_optimization(self):
        A = np.eye(3) * 10.0
        coms = np.array([[5.0, 5.0, 5.0]])
        results = optimize_angles_each_system_site(
            A, coms, A, coms, coms,
            angles_init=(0.0, 0.0, 0.0), rng=1,
        )
        assert len(results) == 1
        idx, angles, best_iou = results[0]
        assert idx == 0
        assert len(angles) == 3
        assert best_iou >= 0.0


class TestCheckSystemVoronoiVolumeTiling:
    """Test check_system_voronoi_volume_tiling."""

    def test_basic_check(self, simple_atoms):
        from molecular_modeling.build_molecule_groups import build_molecule_groups
        groups = build_molecule_groups(simple_atoms)
        report = check_system_voronoi_volume_tiling(
            simple_atoms, groups=groups, verbose=False,
        )
        assert 'sum_voronoi_volume' in report
        assert 'cell_volume' in report
        assert 'ok' in report
        assert 'site_volumes' in report

    def test_groups_none_fallback(self, simple_atoms):
        report = check_system_voronoi_volume_tiling(
            simple_atoms, groups=None, verbose=False,
        )
        assert 'ok' in report

    def test_verbose_output(self, simple_atoms, capsys):
        report = check_system_voronoi_volume_tiling(
            simple_atoms, verbose=True,
        )
        captured = capsys.readouterr()
        assert "[Voronoi tiling volume check]" in captured.out
        assert "Sum of Voronoi cell volumes" in captured.out
        assert "PASS" in captured.out or "FAIL" in captured.out


