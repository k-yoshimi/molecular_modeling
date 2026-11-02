import time
from dataclasses import replace

import numpy as np
import pytest
from ase import Atoms
from molecular_modeling.csp.molecule import load_molecule, molecule_volume
from molecular_modeling.csp.screen import (
    has_intermolecular_overlap,
    packing_fraction,
    prescreen,
    rigid_lj_energy,
)

from tests.csp.conftest import default_cfg

pytest.importorskip("pyxtal")
from molecular_modeling.csp.gen_random import generate_random  # noqa: E402

pytestmark = pytest.mark.filterwarnings("ignore:Setting the shape on a NumPy array:DeprecationWarning")


def expand_coms(atoms, factor):
    """Scale the cell by ``factor`` and move each molecule rigidly with its centroid (fractional kept)."""
    new = atoms.copy()
    new.set_cell(atoms.cell * factor, scale_atoms=False)
    mol_id = atoms.arrays["mol_id"]
    pos = atoms.get_positions()
    out = pos.copy()
    for m in np.unique(mol_id):
        sel = mol_id == m
        com = pos[sel].mean(axis=0)
        out[sel] = pos[sel] - com + com * factor
    new.set_positions(out)
    return new


def _benzene_cands(n, sg=2, z=2, seed=0):
    mol = load_molecule("data/benzene.toml")
    return mol, generate_random(mol, sg, z, n=n, rng=np.random.default_rng(seed))


def test_overlap_detected_for_squeezed_cell():
    _, cands = _benzene_cands(1)
    squeezed = cands[0].atoms.copy()
    squeezed.set_cell(squeezed.cell * 0.7, scale_atoms=False)
    assert has_intermolecular_overlap(squeezed)


def test_no_overlap_for_sparse_cell():
    _, cands = _benzene_cands(1)
    assert not has_intermolecular_overlap(expand_coms(cands[0].atoms, 1.5))


def test_overlap_same_molecule_same_image_ignored():
    a = Atoms("HH", positions=[[0, 0, 0], [0.5, 0, 0]], cell=[20, 20, 20], pbc=True)
    a.set_array("mol_id", np.array([0, 0]))
    assert not has_intermolecular_overlap(a)
    a.set_array("mol_id", np.array([0, 1]), dtype=int)
    assert has_intermolecular_overlap(a)


def test_overlap_own_periodic_image_counts():
    a = Atoms("C", positions=[[0, 0, 0]], cell=[1.5, 20, 20], pbc=True)
    a.set_array("mol_id", np.array([0]))
    assert has_intermolecular_overlap(a)


def test_packing_fraction():
    mol, cands = _benzene_cands(1)
    a = cands[0].atoms
    expected = 2 * molecule_volume(mol) / a.get_volume()
    assert packing_fraction(a, molecule_volume(mol)) == pytest.approx(expected)


def test_lj_prefers_reasonable_density():
    _, cands = _benzene_cands(1)
    a = cands[0].atoms
    assert rigid_lj_energy(expand_coms(a, 1.5)) > rigid_lj_energy(a)


def test_lj_two_atom_analytic():
    # two isolated H atoms of different molecules at r = x_HH: E = -D per pair, 2 molecules
    x, d = 2.886, 0.044 * 0.0433641
    a = Atoms("HH", positions=[[0, 0, 0], [x, 0, 0]], cell=[50, 50, 50], pbc=True)
    a.set_array("mol_id", np.array([0, 1]))
    assert rigid_lj_energy(a, cutoff=10.0) == pytest.approx(-d / 2)


def test_lj_unknown_element():
    a = Atoms("Xe", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
    a.set_array("mol_id", np.array([0]))
    with pytest.raises(KeyError, match="Xe"):
        rigid_lj_energy(a)


def test_prescreen_counts_and_keep():
    mol, cands = _benzene_cands(10)
    order = [id(c) for c in cands]
    kept, counts = prescreen(cands, molecule_volume(mol), replace(default_cfg(), n_prescreen_keep=3))
    assert [id(c) for c in cands] == order
    assert len(kept) <= 3
    assert set(counts) == {"generated", "no_overlap", "packing", "kept"} and counts["generated"] == 10
    assert counts["generated"] >= counts["no_overlap"] >= counts["packing"] >= counts["kept"] == len(kept)
    ljs = [k.meta["lj"] for k in kept]
    assert ljs == sorted(ljs)
    for k in kept:
        assert k.status == "prescreened" and k.density and "packing" in k.meta


def test_prescreen_empty():
    kept, counts = prescreen([], 78.8, default_cfg())
    assert kept == [] and counts == {"generated": 0, "no_overlap": 0, "packing": 0, "kept": 0}


@pytest.mark.slow
def test_timing_benzene_z4(capsys):
    _, cands = _benzene_cands(1, sg=14, z=4)
    a = cands[0].atoms
    t = time.perf_counter()
    rigid_lj_energy(a)
    t1 = time.perf_counter() - t
    t = time.perf_counter()
    has_intermolecular_overlap(a)
    t2 = time.perf_counter() - t
    print(f"lj {t1:.4f}s overlap {t2:.4f}s")
    assert t1 < 1.0 and t2 < 1.0


def _c2_box(x0: float) -> Atoms:
    """One C2 molecule (bond 1.4 Å along x) starting at ``x0`` in a 10 Å cubic cell, then wrapped."""
    a = Atoms("C2", positions=[[x0, 5, 5], [x0 + 1.4, 5, 5]], cell=[10, 10, 10], pbc=True)
    a.set_array("mol_id", np.array([0, 0]))
    a.wrap()
    return a


def test_screen_same_result_for_wrapped_molecule():
    whole, split = _c2_box(4.0), _c2_box(9.4)
    assert split.positions[1, 0] == pytest.approx(0.8)  # the bond crosses the cell boundary
    assert not has_intermolecular_overlap(whole)
    assert not has_intermolecular_overlap(split)
    assert rigid_lj_energy(split) == pytest.approx(rigid_lj_energy(whole), abs=1e-9)


def test_screen_crystal_invariant_under_wrapping():
    c = generate_random(load_molecule("data/benzene.toml"), 14, 4, n=1, rng=np.random.default_rng(5))[0].atoms
    wrapped = c.copy()
    wrapped.translate([3.1, -2.2, 1.7])
    wrapped.wrap()
    unwrapped_pos = c.positions + [3.1, -2.2, 1.7]
    assert np.abs(wrapped.positions - unwrapped_pos).max() > 1.0  # some molecules really are split
    assert has_intermolecular_overlap(wrapped) == has_intermolecular_overlap(c)
    assert rigid_lj_energy(wrapped) == pytest.approx(rigid_lj_energy(c), rel=1e-9)


def test_screen_broken_molecule_rejected():
    a = Atoms("C2", positions=[[2, 5, 5], [6, 5, 5]], cell=[10, 10, 10], pbc=True)  # 4 Å apart: not bonded
    a.set_array("mol_id", np.array([0, 0]))
    assert has_intermolecular_overlap(a)
    with pytest.raises(ValueError, match="not connected"):
        rigid_lj_energy(a)
