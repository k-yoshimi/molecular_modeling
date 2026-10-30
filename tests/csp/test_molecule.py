import numpy as np
import pytest
from ase import Atoms
from molecular_modeling.csp.molecule import load_molecule, match_and_align, molecule_rmsd
from scipy.spatial.transform import Rotation


def test_load_molecule_centers_benzene():
    mol = load_molecule("data/benzene.toml")
    assert len(mol) == 12
    np.testing.assert_allclose(mol.get_positions().mean(axis=0), 0.0, atol=1e-10)


def test_load_molecule_rejects_fragmented(tmp_path):
    (tmp_path / "frag.xyz").write_text("2\n\nC 0 0 0\nC 6 0 0\n")
    with pytest.raises(ValueError, match="fragmented"):
        load_molecule(tmp_path / "frag.xyz")


def test_load_molecule_fractional(tmp_path):
    p = tmp_path / "m.toml"
    p.write_text(
        'symbols = ["C", "C"]\nuse_fractional = true\npositions = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]\n'
        "lattice_vectors = [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]]\n"
    )
    mol = load_molecule(p)
    assert abs(np.linalg.norm(mol.positions[0] - mol.positions[1]) - 1.0) < 1e-9


def _moved(mol, seed=0):
    rot = Rotation.from_euler("xyz", [30, 60, 10], degrees=True)
    moved = mol.copy()
    moved.positions = rot.apply(mol.positions) + [1, 2, 3]
    return moved[np.random.default_rng(seed).permutation(len(moved))]


def test_match_and_align_recovers_rotation_and_permutation():
    mol = load_molecule("data/benzene.toml")
    moved = _moved(mol)
    res = match_and_align(moved, mol)
    assert res.rmsd < 1e-6
    assert molecule_rmsd(mol, moved) < 1e-6
    # contract: R @ mobile[perm[i]] + t ~ target[i]
    out = moved.positions[res.perm] @ res.rotation.T + res.translation
    np.testing.assert_allclose(out, mol.positions, atol=1e-5)
    assert np.isclose(np.linalg.det(res.rotation), 1.0)


def test_match_and_align_asymmetric_molecule():
    mol = Atoms("CCOH", positions=[[0, 0, 0], [1.5, 0, 0], [2.2, 1.1, 0.3], [-0.5, 0.9, 0.4]])
    res = match_and_align(_moved(mol, 3), mol)
    assert res.rmsd < 1e-6


def test_match_and_align_element_map():
    mol = load_molecule("data/benzene.toml")
    sub = _moved(mol)
    sub.symbols[[i for i, s in enumerate(sub.get_chemical_symbols()) if s == "C"][0]] = "S"
    with pytest.raises(ValueError):
        match_and_align(sub, mol)
    res = match_and_align(sub, mol, element_map={"S": "C"})
    assert res.rmsd < 1e-6


def _cubane():
    a = 1.57 / 2
    corners = np.array([[x, y, z] for x in (-a, a) for y in (-a, a) for z in (-a, a)])
    h = corners * (1 + 1.09 / (a * np.sqrt(3)))
    return Atoms("C8H8", positions=np.vstack([corners, h]))


def test_match_and_align_degenerate_inertia_cube():
    mol = _cubane()
    bad = 0
    for s in range(50):
        moved = mol.copy()
        moved.positions = Rotation.random(random_state=s).apply(mol.positions) + [1, -2, 3]
        moved = moved[np.random.default_rng(s).permutation(len(moved))]
        bad += match_and_align(moved, mol).rmsd >= 1e-6
    assert bad == 0


def test_perm_maps_hydrogens_with_heavy_only():
    mol = load_molecule("data/benzene.toml")
    moved = _moved(mol, 5)
    res = match_and_align(moved, mol, heavy_only=True)
    out = moved.positions[res.perm] @ res.rotation.T + res.translation
    np.testing.assert_allclose(out, mol.positions, atol=1e-5)
