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
    # ~85 candidates survive prescreening but at most a few rounds fit, so truncation is certain.
    assert "time budget" in (tmp_path / "run.log").read_text()


class _FakeClock:
    """Stand-in for the ``time`` module whose clock only moves when the fake relaxation runs."""

    def __init__(self) -> None:
        self.now = 0.0

    def perf_counter(self) -> float:
        return self.now


def _run_budget(monkeypatch, budget_min, n_cands=40, n_procs=4, round_s=100.0):
    from molecular_modeling.csp import pipeline

    clock = _FakeClock()
    sizes: list[int] = []

    def fake_relax_many(cands, mol, cfg, method):
        sizes.append(len(cands))
        clock.now += round_s
        return list(cands)

    monkeypatch.setattr(pipeline, "time", clock)
    monkeypatch.setattr(pipeline, "relax_many", fake_relax_many)
    cfg = replace(default_cfg(), n_procs=n_procs, time_budget_min=budget_min)
    done = pipeline._relax_with_budget(list(range(n_cands)), None, cfg, "EMT", 0.0)  # type: ignore[arg-type]
    return done, sizes


def test_relax_budget_floors_to_whole_rounds(monkeypatch):
    # 10 min * 0.9 = 540 s; first round takes 100 s -> 440 s left -> 4 whole rounds of 4 candidates.
    done, sizes = _run_budget(monkeypatch, budget_min=10)
    assert sizes == [4, 16]
    assert len(done) == 20


def test_relax_budget_first_round_always_runs(monkeypatch):
    # 1 min budget is shorter than one round (100 s): the first round still runs, nothing more does.
    done, sizes = _run_budget(monkeypatch, budget_min=1)
    assert sizes == [4]
    assert len(done) == 4
