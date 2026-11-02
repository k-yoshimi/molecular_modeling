"""Tests for seed-crystal candidate generation."""

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.io import write

pytest.importorskip("pyxtal")

from molecular_modeling.csp.config import SeedSpec  # noqa: E402
from molecular_modeling.csp.gen_random import generate_random  # noqa: E402
from molecular_modeling.csp.gen_seed import generate_from_seed, read_seed_molecules  # noqa: E402
from molecular_modeling.csp.molecule import load_molecule, molecule_rmsd, molecule_volume  # noqa: E402

DATA = Path(__file__).parent / "data"


def rigid_rmsd(ref, mol) -> float:
    """RMSD to ``ref`` or to its mirror image (symmetry ops of the parent can be improper)."""
    mirror = ref.copy()
    mirror.positions = -mirror.positions
    return min(molecule_rmsd(ref, mol), molecule_rmsd(mirror, mol))


def write_seed_cif(path: Path, occupancy_of_first_site: float = 1.0) -> Path:
    """Write a wrapped benzene P2_1/c crystal (from the random generator) as a P1 CIF with occupancies."""
    mol = load_molecule("data/benzene.toml")
    cand = generate_random(mol, 14, 4, n=1, rng=np.random.default_rng(0))[0]
    atoms = cand.atoms.copy()
    atoms.wrap()
    a, b, c, al, be, ga = atoms.cell.cellpar()
    frac = atoms.get_scaled_positions()
    lines = [
        "data_seed",
        f"_cell_length_a {a:.6f}",
        f"_cell_length_b {b:.6f}",
        f"_cell_length_c {c:.6f}",
        f"_cell_angle_alpha {al:.6f}",
        f"_cell_angle_beta {be:.6f}",
        f"_cell_angle_gamma {ga:.6f}",
        "_symmetry_space_group_name_H-M 'P 1'",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    for i, (s, f) in enumerate(zip(atoms.get_chemical_symbols(), frac)):
        occ = occupancy_of_first_site if i == 0 else 1.0
        lines.append(f"{s}{i + 1} {s} {f[0]:.8f} {f[1]:.8f} {f[2]:.8f} {occ}")
    path.write_text("\n".join(lines) + "\n")
    return path


def test_read_seed_molecules_splits_and_rejoins():
    cell, mols = read_seed_molecules(DATA / "seed_benzene_p21c.cif")
    assert cell.shape == (3, 3)
    assert len(mols) == 4
    ref = load_molecule("data/benzene.toml")
    assert all(len(m) == 12 and molecule_rmsd(ref, m) < 1e-3 for m in mols)


def test_seed_rejoins_and_keeps_rigid_molecules():
    mol = load_molecule("data/benzene.toml")
    cands = generate_from_seed(mol, SeedSpec(cif=DATA / "seed_benzene_p21c.cif", substitutions=[], volume_scale=1.0))
    assert len(cands) == 1
    assert cands[0].origin == "seed:seed_benzene_p21c"
    assert cands[0].space_group == 1
    atoms = cands[0].atoms
    for k in set(atoms.arrays["mol_id"]):
        assert rigid_rmsd(mol, atoms[atoms.arrays["mol_id"] == k]) < 1e-3


def test_seed_substitution_enumerates_uniform_choices():
    # substitute 1 of benzene's 6 H by F: 6 choices, symmetric copies uniform across molecules
    mol = load_molecule("data/benzene.toml")
    mol_f = mol.copy()
    mol_f.symbols[6] = "F"
    spec = SeedSpec(cif=DATA / "seed_benzene_p21c.cif", substitutions=[("H", "F", 1)], volume_scale=1.0)
    cands = generate_from_seed(mol_f, spec)
    assert len(cands) == 6
    assert all(c.atoms.get_chemical_symbols().count("F") == 4 for c in cands)
    assert sorted(c.meta["substitution_choice"][0] for c in cands) == list(range(6, 12))
    for c in cands:
        at = c.atoms
        for k in range(4):
            sub = at[at.arrays["mol_id"] == k]
            assert rigid_rmsd(mol_f, sub) < 1e-3


def test_seed_substituted_site_is_same_in_every_molecule():
    mol = load_molecule("data/benzene.toml")
    mol_f = mol.copy()
    mol_f.symbols[6] = "F"
    spec = SeedSpec(cif=DATA / "seed_benzene_p21c.cif", substitutions=[("H", "F", 1)], volume_scale=1.0)
    cands = generate_from_seed(mol_f, spec)
    for c in cands:
        at = c.atoms
        sites = []
        for k in range(4):
            sym = at[at.arrays["mol_id"] == k].get_chemical_symbols()
            sites.append(sym.index("F"))
        assert len(set(sites)) == 1, (c.meta["substitution_choice"], sites)
        assert sites[0] == c.meta["substitution_choice"][0]


def test_seed_auto_volume_scales_cell_not_molecules():
    mol = load_molecule("data/benzene.toml")
    spec = SeedSpec(cif=DATA / "seed_benzene_p21c.cif", substitutions=[], volume_scale="auto")
    cell0, _ = read_seed_molecules(spec.cif)
    big = mol.copy()
    big.symbols[6] = "Cl"
    spec_cl = SeedSpec(cif=spec.cif, substitutions=[("H", "Cl", 1)], volume_scale="auto")
    c_auto = generate_from_seed(mol, spec)[0]
    assert abs(c_auto.atoms.get_volume() / abs(np.linalg.det(cell0)) - 1.0) < 0.02
    c_cl = generate_from_seed(big, spec_cl)[0]
    ratio = c_cl.atoms.get_volume() / abs(np.linalg.det(cell0))
    assert ratio > 1.05
    for k in range(4):
        sub = c_cl.atoms[c_cl.atoms.arrays["mol_id"] == k]
        assert rigid_rmsd(big, sub) < 1e-3


def test_seed_rejects_partial_occupancy(tmp_path):
    cif = write_seed_cif(tmp_path / "dis.cif", occupancy_of_first_site=0.7)
    with pytest.raises(ValueError, match="occupancy") as ei:
        read_seed_molecules(cif)
    assert "C1" in str(ei.value)


def test_seed_mismatched_composition_raises():
    mol = load_molecule("data/benzene.toml")
    mol.symbols[6] = "F"  # F in the input but no substitution requested
    with pytest.raises(ValueError):
        generate_from_seed(mol, SeedSpec(cif=DATA / "seed_benzene_p21c.cif", substitutions=[], volume_scale=1.0))


def test_molecule_volume_benzene():
    v = molecule_volume(load_molecule("data/benzene.toml"))
    assert 75.0 < v < 90.0  # ~78.8 with ASE's Bondi radii (converged grid)
    assert v == molecule_volume(load_molecule("data/benzene.toml"))


def test_load_molecule_accepts_compact_symbols():
    from molecular_modeling.molecule_unwrap import is_connected

    nap = load_molecule("data/naphthalene.toml")
    assert len(nap) == 18
    assert is_connected(nap.get_positions(), nap.get_chemical_symbols())


# Asymmetric molecule (no internal symmetry) for the approximate-symmetry tests.
ASYM_MOL = Atoms("CNOS", positions=[[0, 0, 0], [1.3, 0, 0], [-0.4, 1.2, 0.3], [0.2, -0.6, 1.6]])


def _pair_distances(pos: np.ndarray) -> np.ndarray:
    return np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2)


@pytest.mark.parametrize(
    "w_mat",
    [np.diag([1, -1, 1]), np.diag([-1, 1, -1])],  # mirror (improper) and 2-fold axis (proper)
    ids=["mirror", "twofold"],
)
def test_seed_approximate_symmetry_keeps_molecules_rigid(tmp_path, w_mat):
    # gamma is 0.23 degrees off 90, so the operation is a symmetry only within the detection tolerance
    # and its Cartesian matrix is not exactly orthogonal.
    cell = np.array([[10.0, 0, 0], [0.04, 10.0, 0], [0, 0, 10.0]])
    frac0 = (ASYM_MOL.positions + [2.0, 2.5, 3.0]) @ np.linalg.inv(cell)
    frac1 = frac0 @ w_mat.T + [0.5, 0.5, 0.5]
    seed = Atoms("CNOS" * 2, scaled_positions=np.vstack([frac0, frac1]) % 1.0, cell=cell, pbc=True)
    cif = tmp_path / "seed.cif"
    write(cif, seed)
    mol = ASYM_MOL.copy()
    mol.positions -= mol.positions.mean(axis=0)
    cands = generate_from_seed(mol, SeedSpec(cif=cif, volume_scale=1.0))
    ref = _pair_distances(mol.positions)
    atoms = cands[0].atoms
    for m in np.unique(atoms.arrays["mol_id"]):
        idx = np.nonzero(atoms.arrays["mol_id"] == m)[0]
        assert atoms.get_chemical_symbols()[idx[0] : idx[-1] + 1] == list("CNOS")
        np.testing.assert_allclose(_pair_distances(atoms.positions[idx]), ref, atol=1e-9)
    # each molecule still sits where (and with the handedness) the parent has it, up to the tolerance
    _, parents = read_seed_molecules(cif)
    for m, parent in enumerate(parents):
        pos = atoms.positions[atoms.arrays["mol_id"] == m]
        np.testing.assert_allclose(pos - pos.mean(axis=0), parent.positions - parent.positions.mean(axis=0), atol=0.02)
        np.testing.assert_allclose(pos.mean(axis=0), parent.positions.mean(axis=0), atol=1e-9)
