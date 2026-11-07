"""Tests for the relaxation calculators and their dependency checks."""

import importlib
import importlib.util
import subprocess
import sys

import pytest
from molecular_modeling.csp.backends import (
    CONFIG_METHODS,
    METHODS,
    ML_METHODS,
    TEST_METHODS,
    make_calculator,
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
