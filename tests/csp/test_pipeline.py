"""Tests for deduplication, the CSP pipeline and the CLI."""
import time
from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("pyxtal")
pytest.importorskip("tblite")
pytest.importorskip("pymatgen")

from molecular_modeling.csp import cli  # noqa: E402
from molecular_modeling.csp.dedupe import deduplicate  # noqa: E402
from molecular_modeling.csp.gen_random import generate_random  # noqa: E402
from molecular_modeling.csp.molecule import load_molecule  # noqa: E402
from molecular_modeling.csp.pipeline import CSV_HEADER, run_csp  # noqa: E402

from .conftest import default_cfg  # noqa: E402


def test_deduplicate_merges_identical():
    c = generate_random(load_molecule("data/benzene.toml"), 2, 2, n=1, rng=np.random.default_rng(0))[0]
    a, b = replace(c, energy_per_mol=-1.0), replace(c, atoms=c.atoms.copy(), energy_per_mol=-0.5)
    out = deduplicate([a, b])
    assert len(out) == 1
    assert out[0].energy_per_mol == -1.0
    assert out[0].meta["duplicates"] == [b.origin]


def test_run_csp_benzene_end_to_end(tmp_path):
    cfg = replace(
        default_cfg(), space_groups=[2], n_candidates_per_sg=6, n_prescreen_keep=2,
        relax_steps=3, n_procs=1, out_dir=tmp_path,
    )
    out = run_csp(cfg)
    rows = (tmp_path / "candidates.csv").read_text().splitlines()
    assert rows[0].startswith("rank,id,energy_per_mol_eV")
    assert len(rows) - 1 == len(out) and len(out) >= 1
    assert {c.origin for c in out} <= {"random", "voronoi"}
    assert (tmp_path / "run.log").exists()
    assert len(list(tmp_path.glob("001_*.cif"))) == 1


def test_run_csp_handles_zero_survivors(tmp_path):
    cfg = replace(default_cfg(), space_groups=[2], n_candidates_per_sg=3, packing_range=(0.99, 1.0), out_dir=tmp_path)
    assert run_csp(cfg) == []
    assert (tmp_path / "candidates.csv").read_text().splitlines() == [CSV_HEADER]


def test_cli_returns_zero(tmp_path):
    cfg_path = tmp_path / "csp.toml"
    cfg_path.write_text(
        'molecule = "data/benzene.toml"\nmode = "informed"\nspace_groups = [2]\n'
        f'n_candidates_per_sg = 4\nn_prescreen_keep = 1\nrelax_steps = 2\nn_procs = 1\nout_dir = "{tmp_path / "out"}"\n'
    )
    assert cli.main([str(cfg_path)]) == 0
    assert (tmp_path / "out" / "candidates.csv").exists()


@pytest.mark.slow
def test_run_csp_time_budget_respected(tmp_path):
    """With far more candidates than a 3-minute budget allows, the run stops near the budget.

    The budget is a target for the whole run; the first round of ``n_procs`` relaxations always runs,
    so ``relax_steps`` is chosen such that one round fits well inside the budget.
    """
    cfg = replace(
        default_cfg(), space_groups=[2, 14], n_candidates_per_sg=100, n_prescreen_keep=150,
        relax_steps=60, n_procs=2, time_budget_min=3, out_dir=tmp_path,
    )
    t0 = time.monotonic()
    run_csp(cfg)
    assert time.monotonic() - t0 < 3 * 60 * 1.2
    assert "time budget" in (tmp_path / "run.log").read_text()
