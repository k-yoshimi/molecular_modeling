"""Tests for comparing CSP candidates with experimental structures (RMSD_15)."""

import shutil
from dataclasses import replace

import numpy as np
import pytest
from ase import Atoms
from ase.build import make_supercell
from ase.io import write
from scipy.spatial.transform import Rotation

pytest.importorskip("pyxtal")
pytest.importorskip("tblite")
pytest.importorskip("pymatgen")

from molecular_modeling.csp import validate as validate_mod  # noqa: E402
from molecular_modeling.csp.gen_random import generate_random  # noqa: E402
from molecular_modeling.csp.gen_seed import split_molecules  # noqa: E402
from molecular_modeling.csp.molecule import load_molecule  # noqa: E402
from molecular_modeling.csp.pipeline import CSV_HEADER, run_csp  # noqa: E402
from molecular_modeling.csp.validate import ValidationReport, rmsd_cluster, validate_run  # noqa: E402

from .conftest import default_cfg  # noqa: E402


def _benzene(sg: int, z: int, seed: int) -> Atoms:
    return generate_random(load_molecule("data/benzene.toml"), sg, z, n=1, rng=np.random.default_rng(seed))[0].atoms


def test_rmsd_cluster_zero_for_same_structure_and_rigid_motion():
    c = generate_random(load_molecule("data/benzene.toml"), 14, 4, n=1, rng=np.random.default_rng(5))[0].atoms
    moved = c.copy()
    moved.translate([0.3, -0.2, 0.1])
    moved.wrap()
    rmsd, n = rmsd_cluster(c, moved)
    assert rmsd < 1e-3 and n == 15


def test_rmsd_cluster_large_for_different_packing():
    mol = load_molecule("data/benzene.toml")
    a = generate_random(mol, 14, 4, n=1, rng=np.random.default_rng(5))[0].atoms
    b = generate_random(mol, 2, 2, n=1, rng=np.random.default_rng(6))[0].atoms
    rmsd, n = rmsd_cluster(a, b)
    assert rmsd > 0.8 or n < 15


def test_rmsd_cluster_ignores_atom_and_molecule_order():
    a = _benzene(14, 4, 5)
    rng = np.random.default_rng(1)
    mol_ids = a.get_array("mol_id")
    order = []
    for m in rng.permutation(np.unique(mol_ids)):
        order.extend(rng.permutation(np.nonzero(mol_ids == m)[0]))
    b = a[order]
    rmsd, n = rmsd_cluster(a, b)
    assert rmsd < 1e-3 and n == 15


@pytest.mark.parametrize(
    "matrix",
    [[[1, 1, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 0], [0, 0, 1], [1, 0, 0]], [[2, 0, 0], [0, 1, 0], [0, 0, 1]]],
)
def test_rmsd_cluster_independent_of_cell_choice(matrix):
    a = _benzene(14, 4, 5)
    b = make_supercell(a, np.array(matrix))
    b.wrap()
    rmsd, n = rmsd_cluster(a, b)
    assert rmsd < 1e-3 and n == 15


def test_rmsd_cluster_rotated_crystal():
    a = _benzene(2, 2, 3)
    b = a.copy()
    b.rotate(37, "x", rotate_cell=True)
    b.rotate(-71, "z", rotate_cell=True)
    rmsd, n = rmsd_cluster(a, b)
    assert rmsd < 1e-3 and n == 15


def test_rmsd_cluster_mirror_image_crystal():
    a = _benzene(19, 4, 2)
    b = a.copy()
    b.positions *= -1
    b.set_cell(-np.array(b.get_cell()))
    b.wrap()
    rmsd, n = rmsd_cluster(a, b)
    assert rmsd < 1e-3 and n == 15


def test_rmsd_cluster_small_for_slightly_perturbed():
    a = _benzene(14, 4, 5)
    b = a.copy()
    b.positions += np.random.default_rng(0).normal(scale=0.05, size=b.positions.shape)
    rmsd, n = rmsd_cluster(a, b)
    assert 0.01 < rmsd < 0.2 and n == 15


def test_validate_run_reports_rank_and_origin(tmp_path):
    cfg = replace(
        default_cfg(),
        space_groups=[2],
        n_candidates_per_sg=6,
        n_prescreen_keep=2,
        relax_steps=3,
        n_procs=1,
        out_dir=tmp_path,
    )
    out = run_csp(cfg)
    exp = tmp_path / "exp.cif"
    shutil.copy(next(tmp_path.glob("001_*.cif")), exp)
    report = validate_run(tmp_path, [exp], flip=None)
    assert isinstance(report, ValidationReport)
    entry = report.entries[0]
    assert report.success and entry.best_rank == 1 and entry.best_origin == out[0].origin
    assert entry.best_rmsd is not None and entry.best_rmsd < 1e-3
    assert entry.flip_match_rank is None
    md = report.to_markdown()
    assert "exp.cif" in md


def _fake_out_dir(tmp_path, structures):
    """Write a minimal out_dir (candidates.csv + CIFs) for the given (id, origin, atoms) list."""
    lines = [CSV_HEADER]
    for rank, (cid, origin, atoms) in enumerate(structures, start=1):
        lines.append(f"{rank},{cid},-1.0,0.0,1.0,2,2,{origin},0.0,relaxed")
        write(tmp_path / f"{rank:03d}_{cid}.cif", atoms)
    (tmp_path / "candidates.csv").write_text("\n".join(lines) + "\n")


def test_validate_run_symmetric_experimental_cif_and_rank_outside_top(tmp_path):
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.io.cif import CifWriter

    target = _benzene(14, 4, 5)
    others = [_benzene(2, 2, s) for s in (6, 7)]
    _fake_out_dir(
        tmp_path,
        [("random-2-0", "random", others[0]), ("random-2-1", "random", others[1]), ("voronoi-14-0", "voronoi", target)],
    )
    exp = tmp_path / "exp_p21c.cif"
    CifWriter(AseAtomsAdaptor.get_structure(target), symprec=0.01).write_file(str(exp))
    assert "P 1" not in exp.read_text().split("_symmetry_space_group_name_H-M")[1].splitlines()[0]

    report = validate_run(tmp_path, [exp], n_top=2, flip=None)
    entry = report.entries[0]
    assert entry.best_rank == 3 and entry.best_origin == "voronoi"
    assert entry.best_rmsd is not None and entry.best_rmsd < 0.1
    assert not entry.success and not report.success  # rank 3 is outside the top 2
    assert len(entry.comparisons) == 3


def _swap_first(atoms: Atoms, x: str, y: str) -> Atoms:
    """Copy of ``atoms`` in which, in every molecule, the first ``x`` and the first ``y`` exchange elements."""
    out = atoms.copy()
    sym = np.array(out.get_chemical_symbols())
    mol_ids = out.get_array("mol_id")
    for m in np.unique(mol_ids):
        idx = np.nonzero(mol_ids == m)[0]
        ix = idx[sym[idx] == x][0]
        iy = idx[sym[idx] == y][0]
        sym[ix], sym[iy] = y, x
    out.set_chemical_symbols(list(sym))
    return out


# Non-planar molecule with unequal O/S counts (1 O, 2 S), like OS-2SMe (1 O, 3 S).
# O and the first S are both 1.62 Å from C so that the molecule stays bonded when they exchange elements.
OSS_MOL = Atoms("CSOSN", positions=[[0, 0, 0], [1.62, 0, 0], [0, 1.62, 0], [-1.2, -1.2, 0.3], [-0.4, 0.2, 1.3]])


def test_validate_run_flip_match_is_recorded_separately(tmp_path):
    exp_atoms = generate_random(OSS_MOL, 2, 2, n=1, rng=np.random.default_rng(4))[0].atoms
    flipped = _swap_first(exp_atoms, "O", "S")  # same composition, O and one S exchanged in every molecule
    _fake_out_dir(tmp_path, [("random-2-0", "random", flipped), ("random-2-1", "random", exp_atoms)])
    exp = tmp_path / "exp.cif"
    write(exp, exp_atoms)
    report = validate_run(tmp_path, [exp])
    entry = report.entries[0]
    assert entry.best_rank == 2  # plain match only with the identical structure
    assert entry.flip_match_rank == 1  # the flipped candidate is a flip match, not a plain match
    assert entry.flip_rmsd is not None and entry.flip_rmsd < 1e-3
    assert not entry.success or entry.best_rank <= report.n_top
    no_flip = validate_run(tmp_path, [exp], flip=None)
    assert no_flip.entries[0].flip_match_rank is None


def test_validate_run_flip_not_recorded_for_plain_match(tmp_path):
    exp_atoms = generate_random(OSS_MOL, 2, 2, n=1, rng=np.random.default_rng(4))[0].atoms
    _fake_out_dir(tmp_path, [("random-2-0", "random", exp_atoms)])
    exp = tmp_path / "exp.cif"
    write(exp, exp_atoms)
    entry = validate_run(tmp_path, [exp]).entries[0]
    assert entry.best_rank == 1 and entry.flip_match_rank is None


def test_validate_run_flip_skipped_with_warning_when_element_absent(tmp_path, caplog):
    a = _benzene(14, 4, 5)
    _fake_out_dir(tmp_path, [("random-14-0", "random", a)])
    exp = tmp_path / "exp.cif"
    write(exp, a)
    with caplog.at_level("WARNING", logger="molecular_modeling.csp.validate"):
        entry = validate_run(tmp_path, [exp]).entries[0]
    assert entry.flip_match_rank is None
    assert "flip" in caplog.text


def _reoriented(atoms: Atoms, angle: float, axis_index: int = 0) -> Atoms:
    """Turn every molecule about one of its principal axes (0 = plane normal for a planar molecule)."""
    out = Atoms(cell=atoms.get_cell(), pbc=True)
    for m in split_molecules(atoms):
        c = m.positions.mean(axis=0)
        x = m.positions - c
        _, v = np.linalg.eigh(x.T @ x)
        m.positions = x @ Rotation.from_rotvec(angle * v[:, axis_index]).as_matrix().T + c
        out += m
    return out


@pytest.mark.parametrize(
    ("name", "sg", "z", "seed", "angle"),
    [
        ("benzene", 14, 4, 5, np.pi / 3),
        ("benzene", 2, 2, 3, np.pi / 3),
        ("benzene", 19, 4, 2, np.pi / 3),
        ("naphthalene", 14, 4, 1, np.pi),
        ("naphthalene", 2, 2, 3, np.pi),
    ],
)
def test_rmsd_cluster_symmetry_equivalent_molecule_orientation(name, sg, z, seed, angle):
    mol = load_molecule(f"data/{name}.toml")
    a = generate_random(mol, sg, z, n=1, rng=np.random.default_rng(seed))[0].atoms
    b = _reoriented(a, angle)
    for x, y in ((a, b), (b, a)):
        rmsd, n = rmsd_cluster(x, y)
        assert rmsd < 0.1 and n == 15


@pytest.mark.parametrize("name", ["benzene", "naphthalene"])
def test_rmsd_cluster_noise_on_both_sides(name):
    mol = load_molecule(f"data/{name}.toml")
    # (14, 4, 1) is avoided: its naphthalene cell has H...H contacts that noise turns into bonds.
    for sg, z, seed in [(14, 4, 2), (2, 2, 3), (19, 4, 2), (14, 4, 5)]:
        a0 = generate_random(mol, sg, z, n=1, rng=np.random.default_rng(seed))[0].atoms
        for t in range(3):
            r = np.random.default_rng(100 + t)
            a = a0.copy()
            a.positions += r.normal(scale=0.05, size=a.positions.shape)
            b = a0.copy()
            b.positions += r.normal(scale=0.05, size=b.positions.shape)
            rmsd, n = rmsd_cluster(a, b)
            assert rmsd < 0.8 and n == 15, (sg, seed, t, rmsd, n)


def test_validate_run_reports_reason_when_not_comparable(tmp_path, caplog):
    nap = generate_random(load_molecule("data/naphthalene.toml"), 2, 2, n=1, rng=np.random.default_rng(3))[0].atoms
    _fake_out_dir(tmp_path, [("random-2-0", "random", nap)])
    exp = tmp_path / "exp.cif"
    write(exp, _benzene(14, 4, 5))
    with caplog.at_level("WARNING", logger="molecular_modeling.csp.validate"):
        report = validate_run(tmp_path, [exp], flip=None)
    comp = report.entries[0].comparisons[0]
    assert comp.note is not None and "composition" in comp.note
    assert "composition" in caplog.text and "001_random-2-0.cif" in caplog.text
    assert "composition" in report.to_markdown()


def test_validate_run_records_mirror_match(tmp_path):
    # OSS_MOL is chiral, so in P2_12_12_1 the inverted crystal is not superposable by a proper motion.
    a = generate_random(OSS_MOL, 19, 4, n=1, rng=np.random.default_rng(2))[0].atoms
    b = a.copy()
    b.positions *= -1  # inverted through the origin; same (right-handed) cell so the CIF keeps the handedness
    b.wrap()
    _fake_out_dir(tmp_path, [("random-19-0", "random", b), ("random-19-1", "random", a)])
    exp = tmp_path / "exp.cif"
    write(exp, a)
    entry = validate_run(tmp_path, [exp], flip=None).entries[0]
    assert entry.best_rank == 1 and entry.best_mirror is True
    assert [c.mirror for c in entry.comparisons] == [True, False]


def test_cli_writes_markdown(tmp_path, capsys):
    a = _benzene(14, 4, 5)
    _fake_out_dir(tmp_path, [("random-14-0", "random", a)])
    exp = tmp_path / "exp.cif"
    write(exp, a)
    assert validate_mod.main([str(tmp_path), str(exp)]) == 0
    printed = capsys.readouterr().out
    assert "exp.cif" in printed
    assert "exp.cif" in (tmp_path / "validation.md").read_text()


# Asymmetric three-atom molecule used for the extra-molecule (density) tests.
CNO_MOL = Atoms("CNO", positions=[[0, 0, 0], [1.3, 0, 0], [-0.4, 1.2, 0.3]])


def _cno_crystal(body_centred: bool) -> Atoms:
    """One CNO molecule in a 20 Å cubic cell, optionally with a second copy at the body centre."""
    pos = [CNO_MOL.positions + 5.0]
    if body_centred:
        pos.append(CNO_MOL.positions + 15.0)
    n = len(pos)
    a = Atoms("CNO" * n, positions=np.vstack(pos), cell=np.eye(3) * 20.0, pbc=True)
    a.set_array("mol_id", np.repeat(np.arange(n), 3))
    return a


@pytest.mark.parametrize("reverse", [False, True])
def test_rmsd_cluster_extra_candidate_molecules_are_not_a_match(reverse):
    # The body-centred crystal contains every molecule of the simple cubic one plus as many again;
    # it must not count as the same structure in either direction.
    simple, centred = _cno_crystal(False), _cno_crystal(True)
    a, b = (centred, simple) if reverse else (simple, centred)
    rmsd, n = rmsd_cluster(a, b)
    assert n < 15


def test_validate_run_extra_candidate_molecules_not_matched(tmp_path):
    _fake_out_dir(tmp_path, [("random-1-0", "random", _cno_crystal(True))])
    exp = tmp_path / "exp.cif"
    write(exp, _cno_crystal(False))
    entry = validate_run(tmp_path, [exp], flip=None).entries[0]
    assert entry.best_rank is None and not entry.success


def test_validate_run_flipped_only_candidate(tmp_path):
    exp_atoms = generate_random(OSS_MOL, 2, 2, n=1, rng=np.random.default_rng(4))[0].atoms
    _fake_out_dir(tmp_path, [("random-2-0", "random", _swap_first(exp_atoms, "O", "S"))])
    exp = tmp_path / "exp.cif"
    write(exp, exp_atoms)
    report = validate_run(tmp_path, [exp])
    entry = report.entries[0]
    assert not entry.success and not report.success
    assert entry.best_rank is None
    assert entry.flip_match_rank == 1


def test_validate_run_flip_check_rejects_extra_molecules(tmp_path):
    # Flip (element-mapped) comparisons get the same extra-molecule check as plain ones.
    simple = _cno_crystal(False)
    simple.set_chemical_symbols(["C", "S", "O"])
    centred = _cno_crystal(True)
    centred.set_chemical_symbols(["C", "O", "S"] * 2)  # flipped and doubled
    _fake_out_dir(tmp_path, [("random-1-0", "random", centred)])
    exp = tmp_path / "exp.cif"
    write(exp, simple)
    entry = validate_run(tmp_path, [exp]).entries[0]
    assert entry.best_rank is None and entry.flip_match_rank is None
