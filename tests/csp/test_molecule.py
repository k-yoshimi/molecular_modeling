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


BENZ_SYM = '["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]'


def _benz_positions(n: int = 12) -> str:
    from molecular_modeling.csp.molecule import load_molecule as _lm

    pos = _lm("data/benzene.toml").positions[:n]
    return "[" + ", ".join(f"[{x:.6f}, {y:.6f}, {z:.6f}]" for x, y, z in pos) + "]"


@pytest.mark.parametrize(
    ("content", "needle"),
    [
        ("this is not toml at all ===\n", "invalid TOML"),
        (f"symbols = {BENZ_SYM}\n", "positions"),
        (f"positions = {_benz_positions()}\n", "symbols"),
        (f"symbols = {BENZ_SYM}\npositions = {_benz_positions(11)}\n", "12 symbols but 11 positions"),
        (f"symbols = {BENZ_SYM.replace('C', 'c')}\npositions = {_benz_positions()}\n", "'c'"),
        (
            f"symbols = {BENZ_SYM.replace(chr(34) + 'C' + chr(34), chr(34) + 'Xx' + chr(34), 1)}\npositions = {_benz_positions()}\n",
            "Xx",
        ),
        (f"symbols = {BENZ_SYM}\nuse_fractional = true\npositions = {_benz_positions()}\n", "lattice_vectors"),
        ("symbols = []\npositions = []\n", "no atoms"),
        ('symbols = ["C", "C"]\npositions = [[0, 0, 0], [0.1, 0, 0]]\n', "overlapping"),
        ('symbols = ["C", "C"]\npositions = [[0, 0], [1.4, 0]]\n', "positions"),
    ],
)
def test_load_molecule_toml_errors_name_the_problem(tmp_path, content, needle):
    p = tmp_path / "mol.toml"
    p.write_text(content)
    with pytest.raises(ValueError) as info:
        load_molecule(p)
    assert needle in str(info.value) and "mol.toml" in str(info.value)


def test_load_molecule_missing_file_names_it(tmp_path):
    with pytest.raises(FileNotFoundError, match="nope.toml"):
        load_molecule(tmp_path / "nope.toml")


@pytest.mark.parametrize(
    ("content", "needle"),
    [
        ("12\ncomment\nC 0 0 0\nC 1.4 0\nthis is garbage\n", "XYZ"),
        ("2\nc\nC 0 0 0\nO 1.2 0 0\n2\nc\nC 0 0 0\nO 1.2 0 0\n", "2 frames"),
        ("2\nc\nC 0 0 0\nC 0.2 0 0\n", "overlapping"),
    ],
)
def test_load_molecule_xyz_errors(tmp_path, content, needle):
    p = tmp_path / "mol.xyz"
    p.write_text(content)
    with pytest.raises(ValueError) as info:
        load_molecule(p)
    assert needle in str(info.value) and "mol.xyz" in str(info.value)


def test_load_molecule_single_frame_xyz_still_loads(tmp_path):
    p = tmp_path / "co.xyz"
    p.write_text("2\nc\nC 0 0 0\nO 1.2 0 0\n")
    assert load_molecule(p).get_chemical_symbols() == ["C", "O"]


def test_load_molecule_extended_xyz_single_frame(tmp_path):
    from ase.io import write as ase_write

    p = tmp_path / "co_ext.xyz"
    ase_write(p, Atoms("CO", positions=[[0, 0, 0], [1.2, 0, 0]]), format="extxyz")
    assert "Properties=" in p.read_text()
    mol = load_molecule(p)
    assert mol.get_chemical_symbols() == ["C", "O"]
    np.testing.assert_allclose(np.linalg.norm(mol.positions[1] - mol.positions[0]), 1.2)


def test_load_molecule_extended_xyz_reads_positions_not_other_columns(tmp_path):
    # forces come before positions in this file; a plain-XYZ reader would take them as coordinates
    p = tmp_path / "co_forces.xyz"
    p.write_text(
        "2\n"
        'Properties=species:S:1:forces:R:3:pos:R:3 pbc="F F F"\n'
        "C 9.0 9.0 9.0 0.0 0.0 0.0\n"
        "O 7.0 7.0 7.0 1.2 0.0 0.0\n"
    )
    mol = load_molecule(p)
    assert mol.get_chemical_symbols() == ["C", "O"]
    np.testing.assert_allclose(mol.positions, [[-0.6, 0, 0], [0.6, 0, 0]], atol=1e-12)


@pytest.mark.parametrize(
    "comment", ["", "benzene from a lab notebook, conformer 1", "Lattice optimization", "stress=0"]
)
def test_load_molecule_plain_xyz_with_any_comment_line(tmp_path, comment):
    p = tmp_path / "co.xyz"
    p.write_text(f"2\n{comment}\nC 0 0 0\nO 1.2 0 0\n")
    mol = load_molecule(p)
    assert mol.get_chemical_symbols() == ["C", "O"]
    np.testing.assert_allclose(mol.positions, [[-0.6, 0, 0], [0.6, 0, 0]], atol=1e-12)


def test_load_molecule_broken_extended_xyz_is_not_read_as_plain(tmp_path):
    p = tmp_path / "bad_ext.xyz"
    p.write_text("2\nProperties=species:S:1:pos:R:3:forces:R:3\nC 0 0 0\nO 1.2 0 0\n")  # forces missing
    with pytest.raises(ValueError, match="bad_ext.xyz"):
        load_molecule(p)
