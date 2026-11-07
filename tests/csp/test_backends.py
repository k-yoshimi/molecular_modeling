"""Tests for the relaxation calculators and their dependency checks."""

import importlib
import importlib.util
import subprocess
import sys
import time

import pytest
from molecular_modeling.csp.backends import (
    CONFIG_METHODS,
    METHODS,
    ML_METHODS,
    TEST_METHODS,
    RelaxBackendError,
    backend_info,
    make_calculator,
    preflight,
)
from molecular_modeling.csp.pipeline import check_dependencies


def test_methods_lists():
    assert CONFIG_METHODS == ("GFN2-xTB", "GFN1-xTB", "OrbMol-v2")
    assert set(TEST_METHODS) <= set(METHODS) and set(CONFIG_METHODS) <= set(METHODS)
    assert "OrbMol-v2" in ML_METHODS


def test_make_calculator_test_methods():
    assert make_calculator("EMT").__class__.__name__ == "EMT"
    with pytest.raises(RuntimeError, match="test backend is broken"):
        make_calculator("EMT-BROKEN")
    with pytest.raises(ValueError, match="GFN-FF"):
        make_calculator("GFN-FF")


def test_importing_backends_does_not_import_torch():
    # a fresh interpreter, so that torch loaded by other tests does not count
    code = (
        "import sys, molecular_modeling.csp.backends, molecular_modeling.csp.pipeline, molecular_modeling.csp.cli; "
        "print('torch' in sys.modules)"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"


def _fake_find_spec(monkeypatch, present: bool):
    real = importlib.util.find_spec

    def fake(name, *args, **kwargs):
        if name in ("orb_models", "torch"):
            return object() if present else None
        return real(name, *args, **kwargs)

    monkeypatch.setattr(importlib.util, "find_spec", fake)


def test_check_dependencies_orbmol_names_ml_extra(monkeypatch):
    _fake_find_spec(monkeypatch, present=False)
    with pytest.raises(ImportError, match=r"orb_models, torch.*pip install -e \"\.\[csp-ml\]\""):
        check_dependencies("OrbMol-v2")


def test_check_dependencies_orbmol_does_not_need_tblite(monkeypatch):
    real = importlib.import_module

    def fake_import(name, *args, **kwargs):
        if name.startswith("tblite"):
            raise ImportError(f"No module named {name!r}")
        return real(name, *args, **kwargs)

    monkeypatch.setattr(importlib, "import_module", fake_import)
    _fake_find_spec(monkeypatch, present=True)
    check_dependencies("OrbMol-v2")
    with pytest.raises(ImportError, match="tblite"):
        check_dependencies("GFN2-xTB")


def test_preflight_ok_for_working_backend():
    assert preflight("EMT") == {}


def test_preflight_failure_is_one_line_with_advice():
    with pytest.raises(RelaxBackendError) as ei:
        preflight("EMT-BROKEN")
    msg = str(ei.value)
    assert "\n" not in msg and "EMT-BROKEN" in msg and "RuntimeError: test backend is broken" in msg
    assert "cached_path" in msg and "pretrained.orbmol_v2" in msg


def test_preflight_cache_location_follows_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("CACHED_PATH_CACHE_ROOT", str(tmp_path / "weights"))
    with pytest.raises(RelaxBackendError, match=str(tmp_path / "weights")):
        preflight("EMT-BROKEN")


def test_preflight_times_out():
    t0 = time.monotonic()
    with pytest.raises(RelaxBackendError, match="TimeoutError: no answer within 2 s"):
        preflight("EMT-SLOW", timeout_s=2)
    assert time.monotonic() - t0 < 15


def test_backend_info_names_package_and_version():
    from importlib.metadata import version

    assert backend_info("GFN2-xTB") == f"tblite {version('tblite')}"
    assert backend_info("EMT") == "test calculator"


needs_orb = pytest.mark.skipif(importlib.util.find_spec("orb_models") is None, reason="orb-models not installed")


@pytest.mark.slow
@needs_orb
def test_preflight_orbmol_real():
    info = preflight("OrbMol-v2")
    assert info["weights"].startswith("orbmol-v2") and info["weights"].endswith(".ckpt") and info["cache"]


@pytest.mark.slow
@needs_orb
def test_run_csp_orbmol_small(tmp_path):
    from dataclasses import replace

    from molecular_modeling.csp.pipeline import run_csp

    from tests.csp.conftest import default_cfg

    cfg = replace(default_cfg(), n_candidates_per_sg=4, n_prescreen_keep=2, relax_steps=3, n_procs=2, out_dir=tmp_path)
    out = run_csp(cfg, method="OrbMol-v2")
    log = (tmp_path / "run.log").read_text()
    assert "relax method OrbMol-v2 (from argument); orb-models 0.7" in log and "worker pool" in log
    rows = (tmp_path / "candidates.csv").read_text().splitlines()[1:]
    assert len(out) >= 1 and len(rows) == len(out) and all(r.endswith(",OrbMol-v2") for r in rows)
