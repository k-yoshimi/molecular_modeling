"""Tests for cell relaxation of candidates."""
import importlib.util
from dataclasses import replace

import numpy as np
import pytest
from ase import Atoms
from molecular_modeling.csp.candidate import Candidate
from molecular_modeling.csp.gen_random import generate_random
from molecular_modeling.csp.molecule import load_molecule
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
    failed = next(c for c in out if c.status == "failed")
    assert "error" in failed.meta
