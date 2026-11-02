"""Tests for cell relaxation of candidates."""
import importlib.util
from dataclasses import replace

import numpy as np
import pytest
from ase import Atoms
from molecular_modeling.csp.candidate import Candidate
from molecular_modeling.csp.gen_random import generate_random
from molecular_modeling.csp.molecule import load_molecule, molecule_rmsd
from molecular_modeling.csp.relax import (
    make_calculator,
    max_deformation,
    relax_candidate,
    relax_many,
)

from tests.csp.conftest import default_cfg


def test_make_calculator_rejects_unknown():
    with pytest.raises(ValueError):
        make_calculator("GFN-FF")


def test_mirror_image_molecules_not_flagged():
    mol = load_molecule("data/benzene.toml")
    # a chiral-ish reference: displace one atom out of plane so the mirror differs from any rotation
    mol = mol.copy()
    mol.positions[0] += [0.0, 0.0, 0.3]
    mirror = mol.copy()
    mirror.positions = -mirror.positions
    cell = np.diag([30.0, 30.0, 30.0])
    n = len(mol)
    crystal = Atoms(
        symbols=list(mol.symbols) + list(mirror.symbols),
        positions=np.vstack([mol.positions + [7, 7, 7], mirror.positions + [20, 20, 20]]),
        cell=cell,
        pbc=True,
    )
    crystal.set_array("mol_id", np.array([0] * n + [1] * n))
    assert max_deformation(crystal, mol) < 1e-3
    # the mirror branch is needed: the plain rmsd of the mirrored molecule is clearly nonzero
    assert molecule_rmsd(mol, mirror) > 0.05


def test_max_deformation_across_cell_boundary():
    mol = load_molecule("data/benzene.toml")
    n = len(mol)
    atoms = Atoms(symbols=list(mol.symbols), positions=mol.positions + [0.5, 0.5, 5.0], cell=np.diag([12.0, 12.0, 12.0]), pbc=True)
    atoms.set_array("mol_id", np.zeros(n, dtype=int))
    atoms.wrap()
    assert np.any(np.abs(atoms.positions - (mol.positions + [0.5, 0.5, 5.0])) > 1.0)  # really wrapped
    assert max_deformation(atoms, mol) < 1e-3


def _emt_cands(n, gold=()):
    out = []
    for k in range(n):
        el = "Au" if k in gold else "Cu"
        atoms = Atoms(el + "4", positions=[[0, 0, 0], [0, 0, 2.5], [5, 5, 0], [5, 5, 2.5]], cell=[10, 10, 10], pbc=True)
        atoms.set_array("mol_id", np.array([0, 0, 1, 1]))
        out.append(Candidate(atoms=atoms, origin="random", space_group=1, meta={"pyxtal": object(), "tag": k}))
    return out


def test_relax_many_parallel_order_and_meta():
    mol = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, 2.5]])
    cfg = replace(default_cfg(), n_procs=2, relax_steps=2, max_mol_rmsd=5.0)
    out = relax_many(_emt_cands(4), mol, cfg, method="EMT")
    assert [c.meta["tag"] for c in out] == [0, 1, 2, 3]
    assert all("pyxtal" not in c.meta for c in out)
    assert all(c.status in {"relaxed", "rejected_deformed"} for c in out)


def test_relax_many_survives_worker_crash():
    mol = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, 2.5]])
    cfg = replace(default_cfg(), n_procs=2, relax_steps=2, max_mol_rmsd=5.0)
    out = relax_many(_emt_cands(4, gold={1}), mol, cfg, method="EMT-CRASH-ON-AU")
    assert [c.meta["tag"] for c in out] == [0, 1, 2, 3]
    assert out[1].status == "failed" and "error" in out[1].meta
    assert [out[i].status for i in (0, 2, 3)] == ["relaxed"] * 3


def test_relax_candidate_emt_fills_fields():
    mol = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, 2.5]])
    atoms = Atoms("Cu4", positions=[[0, 0, 0], [0, 0, 2.5], [5, 5, 0], [5, 5, 2.5]],
                  cell=[10, 10, 10], pbc=True)
    atoms.set_array("mol_id", np.array([0, 0, 1, 1]))
    c = Candidate(atoms=atoms, origin="random", space_group=1, meta={"pyxtal": object(), "lj": 1.0})
    out = relax_candidate(c, mol, "EMT", fmax=0.5, steps=3, max_mol_rmsd=5.0)
    assert out.status in {"relaxed", "rejected_deformed"}
    assert out.energy_per_mol is not None and out.density > 0
    assert "max_mol_rmsd" in out.meta and isinstance(out.meta["converged"], bool)


needs_tblite = pytest.mark.skipif(importlib.util.find_spec("tblite") is None, reason="tblite not installed")


@needs_tblite
def test_relax_candidate_gfn2_benzene_p1_short():
    mol = load_molecule("data/benzene.toml")
    c = generate_random(mol, 2, 2, n=1, rng=np.random.default_rng(3))[0]
    out = relax_candidate(c, mol, "GFN2-xTB", fmax=0.1, steps=5, max_mol_rmsd=0.5)
    assert out.status in {"relaxed", "rejected_deformed"}
    assert out.energy_per_mol is not None and out.density > 0


@needs_tblite
def test_relax_many_survives_failure():
    mol = load_molecule("data/benzene.toml")
    cands = generate_random(mol, 2, 2, n=2, rng=np.random.default_rng(3))
    cands[0].atoms.positions[0] = cands[0].atoms.positions[1]  # coincident atoms
    cfg = replace(default_cfg(), n_procs=1, relax_steps=2)
    out = relax_many(cands, mol, cfg, method="GFN2-xTB")
    assert len(out) == 2
    assert [c.status for c in out].count("failed") >= 1
    assert out[1].status != "failed"
    failed = next(c for c in out if c.status == "failed")
    assert "error" in failed.meta


class _ControlledPool:
    """In-process stand-in for ``ProcessPoolExecutor`` in which chosen jobs "kill their worker".

    Jobs run synchronously at submission. A job whose ``meta["tag"]`` is in ``crash`` breaks the
    pool: its future fails with ``BrokenProcessPool`` and every later ``submit`` raises it, as a real
    pool does once a worker has died.
    """

    crash: set[int] = set()
    submitted: list[list[int]] = []

    def __init__(self, max_workers=None, initializer=None):
        self.broken = False
        self.log: list[int] = []
        _ControlledPool.submitted.append(self.log)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def shutdown(self, wait=True):
        pass

    def submit(self, fn, job):
        from concurrent.futures import Future
        from concurrent.futures.process import BrokenProcessPool

        if self.broken:
            raise BrokenProcessPool("A child process terminated abruptly")
        cand = job[0]
        self.log.append(cand.meta["tag"])
        fut: Future = Future()
        if cand.meta["tag"] in self.crash:
            self.broken = True
            fut.set_exception(BrokenProcessPool("A child process terminated abruptly"))
        else:
            cand.status = "relaxed"
            fut.set_result(cand)
        return fut


@pytest.mark.parametrize(("n_procs", "crash"), [(4, {1}), (4, {1, 2}), (1, {0, 1, 3}), (3, {2, 3, 4})])
def test_relax_many_recovers_from_broken_pool_at_submit(monkeypatch, n_procs, crash):
    from molecular_modeling.csp import relax as relax_mod

    _ControlledPool.crash = set(crash)
    _ControlledPool.submitted = []
    monkeypatch.setattr(relax_mod, "ProcessPoolExecutor", _ControlledPool)
    mol = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, 2.5]])
    cfg = replace(default_cfg(), n_procs=n_procs)
    out = relax_many(_emt_cands(6), mol, cfg, method="EMT")
    assert [c.meta["tag"] for c in out] == list(range(6))
    for c in out:
        if c.meta["tag"] in crash:
            assert c.status == "failed" and "died" in c.meta["error"]
        else:
            assert c.status == "relaxed", c.meta.get("error")


def test_relax_many_survives_consecutive_worker_crashes():
    mol = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, 2.5]])
    cfg = replace(default_cfg(), n_procs=2, relax_steps=2, max_mol_rmsd=5.0)
    out = relax_many(_emt_cands(5, gold={1, 2}), mol, cfg, method="EMT-CRASH-ON-AU")
    assert [c.meta["tag"] for c in out] == [0, 1, 2, 3, 4]
    assert [c.status for c in out] == ["relaxed", "failed", "failed", "relaxed", "relaxed"]


def _c2_pair(dy: float, cell: float = 10.0, shift=(3.0, 3.0, 3.0)) -> Candidate:
    """Two C2 molecules (bond 1.4 Å) side by side, ``dy`` Å apart."""
    pos = np.array([[0, 0, 0], [1.4, 0, 0], [0, dy, 0], [1.4, dy, 0]]) + shift
    atoms = Atoms("C4", positions=pos, cell=[cell] * 3, pbc=True)
    atoms.set_array("mol_id", np.array([0, 0, 1, 1]))
    atoms.wrap()
    return Candidate(atoms=atoms, origin="random", space_group=1)


C2 = Atoms("C2", positions=[[0, 0, 0], [1.4, 0, 0]])


@pytest.mark.parametrize("shift", [(3.0, 3.0, 3.0), (9.0, 9.5, 3.0)])  # second: bonds across the boundary
def test_relax_candidate_rejects_bonded_molecules(shift):
    out = relax_candidate(_c2_pair(1.4, shift=shift), C2, "EMT", fmax=0.05, steps=0, max_mol_rmsd=0.5)
    assert out.status == "rejected_deformed"
    assert "bond" in out.meta["reject_reason"]


def test_relax_candidate_accepts_separate_molecules():
    out = relax_candidate(_c2_pair(4.0), C2, "EMT", fmax=0.05, steps=0, max_mol_rmsd=0.5)
    assert out.status == "relaxed" and "reject_reason" not in out.meta
