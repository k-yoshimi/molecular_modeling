"""Tests for deduplication, the CSP pipeline and the CLI."""

import math
import time
from dataclasses import replace

import numpy as np
import pytest
from ase.io import write

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
        default_cfg(),
        space_groups=[2],
        n_candidates_per_sg=6,
        n_prescreen_keep=2,
        relax_steps=3,
        n_procs=1,
        out_dir=tmp_path,
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
        default_cfg(),
        space_groups=[2, 14],
        n_candidates_per_sg=100,
        n_prescreen_keep=150,
        relax_steps=60,
        n_procs=2,
        time_budget_min=3,
        out_dir=tmp_path,
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
    """Run the budgeted relaxation with a fake clock.

    ``round_s`` is the wall time of one round of ``n_procs`` parallel relaxations (a number, or a
    function of the 0-based round index); a call relaxing ``n`` candidates costs ``ceil(n / n_procs)``
    rounds.
    """
    from molecular_modeling.csp import pipeline

    clock = _FakeClock()
    sizes: list[int] = []
    cost = round_s if callable(round_s) else (lambda k: round_s)

    def fake_relax_many(cands, mol, cfg, method):
        for _ in range(math.ceil(len(cands) / cfg.n_procs)):
            clock.now += cost(len(sizes))
        sizes.append(len(cands))
        return list(cands)

    monkeypatch.setattr(pipeline, "time", clock)
    monkeypatch.setattr(pipeline, "relax_many", fake_relax_many)
    cfg = replace(default_cfg(), n_procs=n_procs, time_budget_min=budget_min)
    done = pipeline._relax_with_budget(list(range(n_cands)), None, cfg, "EMT", 0.0)  # type: ignore[arg-type]
    return done, sizes, clock.now


def test_relax_budget_runs_whole_rounds_while_time_is_left(monkeypatch):
    # 10 min * 0.9 = 540 s; each round of 4 takes 100 s. After 5 rounds (500 s) only 40 s are left,
    # less than a round, so 20 candidates are relaxed.
    done, sizes, t = _run_budget(monkeypatch, budget_min=10)
    assert sizes == [4] * 5
    assert len(done) == 20 and t == 500


def test_relax_budget_first_round_always_runs(monkeypatch):
    # 1 min budget is shorter than one round (100 s): the first round still runs, nothing more does.
    done, sizes, _ = _run_budget(monkeypatch, budget_min=1)
    assert sizes == [4]
    assert len(done) == 4


def test_relax_budget_stops_when_rounds_get_slower(monkeypatch):
    # The first round takes 1 s, later ones 100 s. A 60 s budget (54 s usable) must not commit all
    # remaining candidates on the strength of the first fast round.
    done, sizes, t = _run_budget(monkeypatch, budget_min=1, round_s=lambda k: 1.0 if k == 0 else 100.0)
    assert sizes == [4, 4]
    assert len(done) == 8 and t == 101


def test_relax_budget_relaxes_everything_when_time_allows(monkeypatch):
    done, sizes, _ = _run_budget(monkeypatch, budget_min=60, n_cands=10)
    assert sizes == [4, 4, 2]
    assert done == list(range(10))


@pytest.mark.parametrize(("weight", "expect_seed"), [(0.0, False), (0.5, True), (None, True)])
def test_seed_generator_weight_switches_seed_candidates(weight, expect_seed):
    from pathlib import Path

    from molecular_modeling.csp import pipeline
    from molecular_modeling.csp.config import SeedSpec

    gens = {"random": 1.0} if weight is None else {"random": 1.0, "seed": weight}
    cfg = replace(
        default_cfg(),
        space_groups=[],
        generators=gens,
        seed_structures=[SeedSpec(cif=Path("tests/csp/data/seed_benzene_p21c.cif"), volume_scale=1.0)],
    )
    cands = pipeline._generate(cfg, load_molecule("data/benzene.toml"), np.random.default_rng(0))
    assert bool(cands) == expect_seed
    assert all(c.origin.startswith("seed:") for c in cands)


def test_write_outputs_quotes_awkward_seed_names(tmp_path):
    import csv

    from molecular_modeling.csp import pipeline

    c = generate_random(load_molecule("data/benzene.toml"), 2, 2, n=1, rng=np.random.default_rng(0))[0]
    c.origin = 'seed:parent,form "II"'
    c.energy_per_mol, c.density, c.status = -1.25, 1.1, "relaxed"
    c.meta.update(id=f"{pipeline._tag(c.origin)}-1-0", max_mol_rmsd=0.01)
    pipeline._write_outputs(replace(default_cfg(), out_dir=tmp_path), [c])
    text = (tmp_path / "candidates.csv").read_text()
    assert text.splitlines()[0] == CSV_HEADER
    rows = list(csv.DictReader(text.splitlines()))
    assert len(rows) == 1 and None not in rows[0]  # no extra columns
    assert rows[0]["origin"] == 'seed:parent,form "II"'
    assert rows[0]["id"] == c.meta["id"]
    assert rows[0]["energy_per_mol_eV"] == "-1.250000" and rows[0]["density_gcm3"] == "1.1000"
    assert rows[0]["max_mol_rmsd"] == "0.0100" and rows[0]["status"] == "relaxed"
    assert (tmp_path / f"001_{c.meta['id']}.cif").exists()


def test_reused_output_dir_with_fewer_survivors(tmp_path):
    """A second run into the same directory leaves CSV and ranked CIFs in agreement."""
    import csv

    from molecular_modeling.csp import pipeline
    from molecular_modeling.csp.validate import validate_run

    cands = generate_random(load_molecule("data/benzene.toml"), 2, 2, n=3, rng=np.random.default_rng(0))
    for k, c in enumerate(cands):
        c.energy_per_mol, c.status = -1.0 + k, "relaxed"
        c.meta["id"] = f"random-2-{k}"
    cfg = replace(default_cfg(), out_dir=tmp_path)
    pipeline._write_outputs(cfg, cands)  # first run: 3 survivors
    assert len(list(tmp_path.glob("[0-9][0-9][0-9]_*.cif"))) == 3
    (tmp_path / "notes.cif").write_text("not a ranked file")

    second = [replace(cands[2], energy_per_mol=-2.0)]  # second run: 1 survivor, a different structure
    pipeline._write_outputs(cfg, second)
    with open(tmp_path / "candidates.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    listed = {f"{int(r['rank']):03d}_{r['id']}.cif" for r in rows}
    ranked = {p.name for p in tmp_path.glob("[0-9][0-9][0-9]_*.cif")}
    assert listed == ranked == {"001_random-2-2.cif"}
    assert (tmp_path / "notes.cif").exists()  # files that are not ranked CIFs are left alone

    exp = tmp_path / "exp.cif"
    write(exp, cands[0].atoms)  # matches only the stale first-run rank 1
    report = validate_run(tmp_path, [exp], flip=None)
    assert report.n_compared == 1 and report.entries[0].best_rank is None


def test_run_csp_rejects_elements_without_lj_parameters(tmp_path, monkeypatch):
    from molecular_modeling.csp import pipeline

    xyz = tmp_path / "silane.xyz"
    xyz.write_text("5\n\nSi 0 0 0\nH 0.85 0.85 0.85\nH -0.85 -0.85 0.85\nH -0.85 0.85 -0.85\nH 0.85 -0.85 -0.85\n")

    def no_generation(*args, **kwargs):
        raise AssertionError("generation must not start")

    monkeypatch.setattr(pipeline, "_generate", no_generation)
    cfg = replace(default_cfg(), molecule=xyz, out_dir=tmp_path / "out")
    with pytest.raises(ValueError, match="Si") as exc:
        run_csp(cfg)
    assert "H, C, N, O, F, S, Cl, Se, Br, I" in str(exc.value)


def test_config_warnings_reach_run_log(tmp_path, monkeypatch):
    from molecular_modeling.csp import pipeline
    from molecular_modeling.csp.config import load_csp_config

    p = tmp_path / "csp.toml"
    p.write_text(
        f'molecule = "data/benzene.toml"\nmode = "blind"\nspace_groups = [2]\nz_prime = 2\nout_dir = "{tmp_path / "out"}"\n'
    )
    cfg = load_csp_config(p)  # warnings are issued here, before run.log exists
    monkeypatch.setattr(pipeline, "_generate", lambda *a, **k: [])
    run_csp(cfg, method="EMT")
    log = (tmp_path / "out" / "run.log").read_text()
    assert "WARNING" in log and "space_groups" in log and "z_prime" in log


def test_run_csp_reports_missing_dependency_before_generation(tmp_path, monkeypatch):
    import importlib

    from molecular_modeling.csp import pipeline

    real = importlib.import_module

    def fake_import(name, *args, **kwargs):
        if name == "tblite":
            raise ImportError("No module named 'tblite'")
        return real(name, *args, **kwargs)

    def no_generation(*args, **kwargs):
        raise AssertionError("generation must not start")

    monkeypatch.setattr(pipeline.importlib, "import_module", fake_import)
    monkeypatch.setattr(pipeline, "_generate", no_generation)
    cfg = replace(default_cfg(), out_dir=tmp_path / "out")
    with pytest.raises(ImportError, match=r"tblite.*pip install -e \"\.\[csp\]\""):
        run_csp(cfg)
    # the test calculator does not need tblite
    monkeypatch.setattr(pipeline, "_generate", lambda *a, **k: [])
    assert run_csp(cfg, method="EMT") == []
