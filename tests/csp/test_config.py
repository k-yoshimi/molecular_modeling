from pathlib import Path

import pytest
from molecular_modeling.csp.candidate import Candidate
from molecular_modeling.csp.config import DEFAULT_SPACE_GROUPS, load_csp_config

from tests.csp.conftest import default_cfg


def test_blind_mode_uses_default_space_groups(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nmode = "blind"\nspace_groups = [1]\n')
    cfg = load_csp_config(p)
    assert cfg.space_groups == [14, 2, 19, 15, 4, 61]
    assert cfg.space_groups == DEFAULT_SPACE_GROUPS
    assert cfg.z is None


def test_informed_mode_requires_space_groups(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nmode = "informed"\n')
    with pytest.raises(ValueError):
        load_csp_config(p)


def test_informed_mode_with_z_and_seed(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text(
        'molecule = "data/benzene.toml"\nmode = "informed"\nspace_groups = [14]\n'
        "[z]\n14 = 4\n"
        '[[seed_structures]]\ncif = "a.cif"\nsubstitutions = [["O", "S", 2]]\n'
    )
    cfg = load_csp_config(p)
    assert cfg.z == {14: 4}
    assert cfg.seed_structures[0].cif == Path("a.cif")
    assert cfg.seed_structures[0].substitutions == [("O", "S", 2)]
    assert cfg.seed_structures[0].volume_scale == "auto"
    assert "seed" in cfg.generators and cfg.generators["random"] == 0.6


def test_defaults_and_no_seed(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\n')
    cfg = load_csp_config(p)
    assert cfg.mode == "blind"
    assert cfg.generators == {"random": 0.6, "voronoi": 0.4}
    assert cfg.n_candidates_per_sg == 150 and cfg.rng == 1


def test_unknown_key_rejected(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nbogus = 1\n')
    with pytest.raises(ValueError, match="bogus"):
        load_csp_config(p)


def test_example_input_loads():
    cfg = load_csp_config("input/csp.toml")
    assert cfg.n_candidates_per_sg == 50


def test_default_cfg_and_candidate():
    from ase import Atoms

    cfg = default_cfg()
    assert cfg.space_groups == [2] and cfg.n_procs == 1
    c = Candidate(atoms=Atoms("H"), origin="random", space_group=2)
    assert c.status == "generated" and c.meta == {}
