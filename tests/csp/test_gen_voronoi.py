import itertools
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from molecular_modeling.csp.candidate import Candidate
from molecular_modeling.csp.config import SeedSpec
from molecular_modeling.csp.gen_voronoi import (
    best_orientation,
    convex_iou,
    generate_voronoi,
    molecule_shape_points,
    orient_candidate,
    site_voronoi_vertices,
)
from molecular_modeling.csp.molecule import load_molecule, molecule_rmsd
from scipy.spatial import ConvexHull

pytest.importorskip("pyxtal")

# ASE 3.x still assigns to ndarray.shape, which NumPy 2.5 deprecates (third-party, harmless).
pytestmark = pytest.mark.filterwarnings("ignore:Setting the shape on a NumPy array:DeprecationWarning")

DATA = Path(__file__).parent / "data"


def spglib_number(atoms: Atoms) -> int:
    import spglib

    cell = (atoms.get_cell().array, atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    return int(spglib.get_symmetry_dataset(cell, symprec=1e-3).number)


def rigid_rmsd(ref: Atoms, mol: Atoms) -> float:
    """RMSD to ``ref`` or to its mirror image (improper symmetry copies are mirror images)."""
    mirror = ref.copy()
    mirror.positions[:, 0] *= -1
    return min(molecule_rmsd(ref, mol), molecule_rmsd(mirror, mol))


def test_molecule_shape_points_centered_and_sized():
    mol = load_molecule("data/benzene.toml")
    pts = molecule_shape_points(mol, n_sphere=26)
    assert pts.shape == (12 * 26, 3)
    np.testing.assert_allclose(pts.mean(axis=0), 0.0, atol=1e-10)
    # every point lies on a vdW sphere: about 1.2 Å (H) or 1.7 Å (C) from its atom
    d = np.linalg.norm(pts.reshape(12, 26, 3) - mol.get_positions()[:, None, :], axis=2)
    np.testing.assert_allclose(d[:6], 1.7)
    np.testing.assert_allclose(d[6:], 1.2)
    # 26 points per atom already give a hull within 10 % of a finely sampled one
    fine = ConvexHull(molecule_shape_points(mol, n_sphere=400)).volume
    assert 0.9 * fine < ConvexHull(pts).volume <= fine


def test_site_voronoi_vertices_single_site_is_cell_box():
    cell = np.diag([4.0, 5.0, 6.0])
    v = site_voronoi_vertices(cell, np.array([[1.0, 1.0, 1.0]]), 0)
    np.testing.assert_allclose(v.mean(axis=0), 0.0, atol=1e-10)
    assert ConvexHull(v).volume == pytest.approx(120.0)


@pytest.mark.parametrize("n_sites", [1, 2])
def test_site_voronoi_vertices_independent_of_skewed_basis(n_sites):
    # the same lattice in a reduced and a strongly skewed basis must give the same cells
    reduced = np.diag([10.0, 1.0, 10.0])
    skewed = np.array([[10.0, 0.0, 0.0], [20.0, 1.0, 0.0], [0.0, 0.0, 10.0]])
    frac = np.array([[0.1, 0.5, 0.1], [0.6, 0.2, 0.55]])[:n_sites]
    coms = frac @ skewed  # Cartesian centres inside the skewed cell; same points for both bases
    for k in range(n_sites):
        v_red = ConvexHull(site_voronoi_vertices(reduced, coms, k)).volume
        v_skew = ConvexHull(site_voronoi_vertices(skewed, coms, k)).volume
        assert v_skew == pytest.approx(v_red)
    if n_sites == 1:
        assert v_red == pytest.approx(100.0)


def test_convex_iou_identity_and_disjoint():
    cube = np.array(list(itertools.product([-1, 1], repeat=3)), float)
    assert convex_iou(cube, cube) == pytest.approx(1.0)
    assert convex_iou(cube, cube * 0.5) == pytest.approx(1 / 8)


def test_convex_iou_no_common_interior_is_zero():
    cube = np.array(list(itertools.product([-1, 1], repeat=3)), float)
    assert convex_iou(cube, cube + 5.0) == 0.0


def test_best_orientation_aligns_elongated_shape_with_elongated_cell():
    # a rod along x in an orthorhombic cell elongated along z -> rod should end up along z
    rod = np.array([[x, y, z] for x in (-3, 3) for y in (-0.5, 0.5) for z in (-0.5, 0.5)], float)
    cell = np.diag([4.0, 4.0, 12.0])
    coms = np.zeros((1, 3))
    rot, iou = best_orientation(rod, site_voronoi_vertices(cell, coms, 0))
    assert abs(rot.apply([1, 0, 0])[2]) > 0.95
    # the volume-matched rod (19 x 3.2 x 3.2 Å) overlaps the 4 x 4 x 12 Å cell by 12 x 3.2 x 3.2 Å
    assert iou == pytest.approx(0.47, abs=0.03)


def test_best_orientation_is_deterministic():
    rod = np.array([[x, y, z] for x in (-3, 3) for y in (-0.5, 0.5) for z in (-0.5, 0.5)], float)
    cv = site_voronoi_vertices(np.diag([4.0, 6.0, 12.0]), np.zeros((1, 3)), 0)
    r1, i1 = best_orientation(rod, cv)
    r2, i2 = best_orientation(rod, cv)
    np.testing.assert_allclose(r1.as_matrix(), r2.as_matrix())
    assert i1 == i2


def test_orient_candidate_keeps_symmetry_and_rigidity():
    from molecular_modeling.csp.gen_random import generate_random

    mol = load_molecule("data/benzene.toml")
    skel = generate_random(mol, 14, 4, n=1, rng=np.random.default_rng(2))[0]
    before = skel.atoms.get_positions().copy()
    out = orient_candidate(skel)
    # benzene is unchanged by the 180° turn, so only the best orientation is returned
    assert len(out) == 1 and all(c.origin == "voronoi" for c in out)
    np.testing.assert_allclose(skel.atoms.get_positions(), before)  # input untouched
    for c in out:
        assert c.space_group == 14
        assert 0.0 < c.meta["voronoi_iou"] <= 1.0
        for k in range(4):
            assert rigid_rmsd(mol, c.atoms[c.atoms.arrays["mol_id"] == k]) < 1e-3
        assert spglib_number(c.atoms) == 14
    assert molecule_rmsd(mol, out[0].atoms[out[0].atoms.arrays["mol_id"] == 0]) < 1e-3


def test_orient_candidate_orients_every_wyckoff_site():
    from molecular_modeling.csp.gen_random import generate_random

    mol = load_molecule("data/benzene.toml")
    skel = generate_random(mol, 2, 4, n=1, rng=np.random.default_rng(3))[0]  # Z'=2 in P-1
    sites = skel.meta["pyxtal"].mol_sites
    assert len(sites) >= 2
    old = [s.orientation.matrix.copy() for s in sites]
    out = orient_candidate(skel)
    assert len(out[0].meta["voronoi_iou_per_site"]) == len(sites)
    new_sites = out[0].meta["pyxtal"].mol_sites
    for s_old, s_new in zip(old, new_sites):
        assert not np.allclose(s_old, s_new.orientation.matrix)
    # original pyxtal object unchanged
    for s, m in zip(sites, old):
        np.testing.assert_allclose(s.orientation.matrix, m)
    for c in out:
        for k in range(4):
            assert rigid_rmsd(mol, c.atoms[c.atoms.arrays["mol_id"] == k]) < 1e-3
        assert spglib_number(c.atoms) == 2


def test_orient_seed_candidate_rotates_each_molecule_about_its_centre():
    from molecular_modeling.csp.gen_seed import generate_from_seed

    mol = load_molecule("data/benzene.toml")
    skel = generate_from_seed(mol, SeedSpec(cif=DATA / "seed_benzene_p21c.cif", substitutions=[], volume_scale=1.0))[0]
    out = orient_candidate(skel)
    assert len(out) == 1 and all(c.origin == "voronoi" for c in out)  # benzene: the turned copy is identical
    ids = skel.atoms.arrays["mol_id"]
    for c in out:
        assert "pyxtal" not in c.meta
        np.testing.assert_array_equal(c.atoms.arrays["mol_id"], ids)
        assert len(c.meta["voronoi_iou_per_site"]) == len(set(ids))
        for k in sorted(set(ids)):
            sel = ids == k
            np.testing.assert_allclose(
                c.atoms.positions[sel].mean(axis=0), skel.atoms.positions[sel].mean(axis=0), atol=1e-8
            )
            assert rigid_rmsd(mol, c.atoms[sel]) < 1e-3


def test_generate_voronoi_concatenates():
    cube_mol = Atoms("C2", positions=[[-0.7, 0, 0], [0.7, 0, 0]])
    atoms = Atoms("C2", positions=cube_mol.positions + 2.0, cell=np.diag([4.0, 4.0, 8.0]), pbc=True)
    atoms.set_array("mol_id", np.array([0, 0]))
    skel = Candidate(atoms=atoms, origin="seed:x", space_group=1)
    out = generate_voronoi([skel, skel])
    # one candidate per skeleton: turning the C2 rod end over end gives the same orientation
    assert len(out) == 2 and all(c.origin == "voronoi" for c in out)
    # the C2 rod ends up along the long cell axis
    d = out[0].atoms.positions[1] - out[0].atoms.positions[0]
    assert abs(d[2]) / np.linalg.norm(d) > 0.95


# Planar zigzag O-C-C-C-C-S: a 180° turn about the plane normal maps every atom onto the position of its
# partner at the other end (O <-> S), so the swapped variant puts O where S was and vice versa.
OS_ROD = Atoms(
    "OCCCCS",
    positions=[[-3.0, 0.5, 0], [-1.7, -0.2, 0], [-0.5, 0.5, 0], [0.5, -0.5, 0], [1.7, 0.2, 0], [3.0, -0.5, 0]],
)


@pytest.mark.parametrize("with_pyxtal", [False, True])
def test_second_variant_swaps_head_and_tail(with_pyxtal):
    if with_pyxtal:
        from molecular_modeling.csp.gen_random import generate_random

        skel = generate_random(OS_ROD, 2, 2, n=1, rng=np.random.default_rng(0))[0]
    else:
        atoms = Atoms(
            OS_ROD.symbols, positions=OS_ROD.positions + [4.0, 4.5, 5.0], cell=np.diag([8.0, 9.0, 10.0]), pbc=True
        )
        atoms.set_array("mol_id", np.zeros(6, dtype=int))
        skel = Candidate(atoms=atoms, origin="seed:x", space_group=1)
    best, swapped = orient_candidate(skel)
    assert swapped.meta["voronoi_variant"] == "flip"
    ids = best.atoms.arrays["mol_id"]
    for k in sorted(set(ids)):
        sel = np.nonzero(ids == k)[0]
        cell = np.asarray(best.atoms.cell)
        p_best = best.atoms.positions[sel]
        p_swap = swapped.atoms.positions[sel]
        # rejoin across the boundary along the bonds (the rod is longer than half the cell)
        from ase.geometry import find_mic
        from molecular_modeling.molecule_unwrap import unwrap_molecule

        syms = list(OS_ROD.symbols)
        p_best = unwrap_molecule(p_best, syms, cell)
        p_swap = unwrap_molecule(p_swap, syms, cell)
        c_best, c_swap = p_best.mean(axis=0), p_swap.mean(axis=0)
        np.testing.assert_allclose(find_mic(c_swap - c_best, cell)[1], 0.0, atol=1e-6)  # same centroid
        # O of the swapped variant sits where S was, and S where O was; the molecular plane is kept
        np.testing.assert_allclose(p_swap[0] - c_swap, p_best[5] - c_best, atol=1e-3)
        np.testing.assert_allclose(p_swap[5] - c_swap, p_best[0] - c_best, atol=1e-3)
        # rigid: all pair distances preserved
        d = lambda p: np.linalg.norm(p[:, None] - p[None], axis=2)  # noqa: E731
        np.testing.assert_allclose(d(p_swap), d(OS_ROD.positions), atol=1e-6)


def _single_site(mol: Atoms) -> Candidate:
    atoms = Atoms(mol.symbols, positions=mol.positions + [4.0, 4.5, 5.0], cell=np.diag([8.0, 9.0, 10.0]), pbc=True)
    atoms.set_array("mol_id", np.zeros(len(mol), dtype=int))
    return Candidate(atoms=atoms, origin="seed:x", space_group=1)


def test_turned_copy_skipped_when_it_is_the_same_orientation():
    # benzene is unchanged by a 180° turn about its plane normal: the second variant would be a duplicate
    out = orient_candidate(_single_site(load_molecule("data/benzene.toml")))
    assert [c.meta["voronoi_variant"] for c in out] == ["best"]


def test_turned_copy_kept_when_head_and_tail_differ():
    out = orient_candidate(_single_site(OS_ROD))
    assert [c.meta["voronoi_variant"] for c in out] == ["best", "flip"]
