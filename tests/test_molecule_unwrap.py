import numpy as np
from ase import Atoms
from molecular_modeling.molecule_unwrap import bond_pairs, extract_molecule, is_connected, unwrap_molecule


def test_unwrap_rejoins_molecule_split_across_boundary():
    cell = np.diag([10.0, 10.0, 10.0])
    pos = np.array([[9.4, 5.0, 5.0], [0.9, 5.0, 5.0]])
    out = unwrap_molecule(pos, ["C", "C"], cell)
    assert np.isclose(np.linalg.norm(out[1] - out[0]), 1.5)
    assert is_connected(out, ["C", "C"])
    assert not is_connected(pos, ["C", "C"])


def test_extract_molecule_uses_plain_mean_and_unwraps():
    atoms = Atoms("CO", positions=[[9.4, 5, 5], [0.7, 5, 5]], cell=np.diag([10.0] * 3), pbc=True)
    syms, centered, com = extract_molecule(atoms, [0, 1])
    assert syms == ["C", "O"]
    np.testing.assert_allclose(centered.mean(axis=0), 0.0, atol=1e-12)
    np.testing.assert_allclose(com, [(9.4 + 10.7) / 2, 5, 5])


def test_unwrap_chain_spanning_multiple_boundaries():
    cell = np.diag([10.0, 10.0, 10.0])
    # L-shaped chain crossing both the x and y boundaries
    true = np.array([[8.0 + 1.4 * k, 9.0, 5.0] for k in range(3)] + [[10.8, 9.0 + 1.4 * k, 5.0] for k in range(1, 3)])
    out = unwrap_molecule(true % 10.0, ["C"] * 5, cell)
    expected = true - true[0] + (true[0] % 10.0)
    np.testing.assert_allclose(out, expected, atol=1e-9)


def test_single_atom_and_non_periodic_axis():
    assert is_connected(np.zeros((1, 3)), ["C"])
    cell = np.diag([10.0, 10.0, 10.0])
    pos = np.array([[9.4, 5.0, 1.0], [0.9, 5.0, 1.0]])
    out = unwrap_molecule(pos, ["C", "C"], cell, pbc=[False, True, True])
    assert not np.isclose(np.linalg.norm(out[1] - out[0]), 1.5)


def test_bond_pairs_lists_each_bond_once():
    pos = np.array([[-1.4, 0, 0], [0, 0, 0], [0, 1.4, 0]])
    assert bond_pairs(pos, ["C", "C", "C"]) == {(0, 1), (1, 2)}
    assert bond_pairs(pos[:1], ["C"]) == set()
