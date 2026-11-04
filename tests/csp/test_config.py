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


def _seed_toml(tmp_path, generators: str = "") -> Path:
    p = tmp_path / "csp.toml"
    p.write_text(
        'molecule = "data/benzene.toml"\n'
        + generators
        + '[[seed_structures]]\ncif = "tests/csp/data/seed_benzene_p21c.cif"\n'
    )
    return p


def test_seed_weight_added_when_seeds_given(tmp_path):
    assert load_csp_config(_seed_toml(tmp_path)).generators["seed"] > 0


def test_seed_weight_zero_kept(tmp_path):
    cfg = load_csp_config(_seed_toml(tmp_path, "generators = {random = 1.0, seed = 0}\n"))
    assert cfg.generators["seed"] == 0


@pytest.mark.parametrize(
    "generators", ["generators = {random = 1.0, seed = -1}\n", "generators = {random = 1.0, sede = 1}\n"]
)
def test_bad_generator_weights_rejected(tmp_path, generators):
    with pytest.raises(ValueError, match="generators"):
        load_csp_config(_seed_toml(tmp_path, generators))


@pytest.mark.parametrize("value", ["nan", "inf", "-inf"])
@pytest.mark.parametrize("name", ["seed", "random"])
def test_non_finite_generator_weights_rejected(tmp_path, name, value):
    other = "voronoi = 1.0" if name == "random" else "random = 1.0"
    with pytest.raises(ValueError, match="generators"):
        load_csp_config(_seed_toml(tmp_path, f"generators = {{{other}, {name} = {value}}}\n"))


@pytest.mark.parametrize("value", [0, 0.0, 2, 0.25])
def test_finite_generator_weights_accepted(tmp_path, value):
    cfg = load_csp_config(_seed_toml(tmp_path, f"generators = {{random = 1.0, seed = {value}}}\n"))
    assert cfg.generators["seed"] == value


def test_blind_mode_forces_z_prime_one_and_keeps_seeds(tmp_path, caplog):
    p = tmp_path / "csp.toml"
    p.write_text(
        'molecule = "data/benzene.toml"\nmode = "blind"\nz_prime = 2\n'
        '[[seed_structures]]\ncif = "tests/csp/data/seed_benzene_p21c.cif"\n'
    )
    with caplog.at_level("WARNING", logger="molecular_modeling.csp.config"):
        cfg = load_csp_config(p)
    assert cfg.z_prime == 1 and cfg.z is None
    assert "z_prime" in caplog.text
    assert len(cfg.seed_structures) == 1  # seeds are a parent compound, not experimental data: kept


def test_informed_mode_keeps_z_prime(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nmode = "informed"\nspace_groups = [2]\nz_prime = 2\n')
    assert load_csp_config(p).z_prime == 2


def test_partial_generators_table_keeps_other_defaults(tmp_path):
    # switching seeds off must not silently switch the random and Voronoi generators off too
    cfg = load_csp_config(_seed_toml(tmp_path, "generators = {seed = 0}\n"))
    assert cfg.generators == {"random": 0.6, "voronoi": 0.4, "seed": 0}
    p = tmp_path / "plain.toml"
    p.write_text('molecule = "data/benzene.toml"\ngenerators = {random = 1.0}\n')
    assert load_csp_config(p).generators == {"random": 1.0, "voronoi": 0.4}


def test_all_generators_off_rejected(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\ngenerators = {random = 0, voronoi = 0}\n')
    with pytest.raises(ValueError, match="no generator"):
        load_csp_config(p)
    # seeds alone are enough
    cfg = load_csp_config(_seed_toml(tmp_path, "generators = {random = 0, voronoi = 0}\n"))
    assert cfg.generators["seed"] > 0


def test_blind_mode_warns_when_space_groups_are_overridden(tmp_path, caplog):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nmode = "blind"\nspace_groups = [14, 2]\n')
    with caplog.at_level("WARNING", logger="molecular_modeling.csp.config"):
        cfg = load_csp_config(p)
    assert cfg.space_groups == DEFAULT_SPACE_GROUPS
    assert "space_groups" in caplog.text
    assert any("space_groups" in w for w in cfg.config_warnings)


def test_blind_mode_no_warning_without_space_groups(tmp_path, caplog):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nmode = "blind"\n')
    with caplog.at_level("WARNING", logger="molecular_modeling.csp.config"):
        cfg = load_csp_config(p)
    assert caplog.text == "" and cfg.config_warnings == []


def test_example_input_loads_without_warnings(caplog):
    with caplog.at_level("WARNING", logger="molecular_modeling.csp.config"):
        cfg = load_csp_config("input/csp.toml")
    assert cfg.config_warnings == [] and caplog.text == ""


def test_config_warnings_is_not_an_input_key(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nconfig_warnings = ["x"]\n')
    with pytest.raises(ValueError, match="config_warnings"):
        load_csp_config(p)


@pytest.mark.parametrize(
    ("line", "key"),
    [
        ("packing_range = [0.8, 0.4]", "packing_range"),
        ("packing_range = [0.0, 0.5]", "packing_range"),
        ("packing_range = [0.5, 1.2]", "packing_range"),
        ("packing_range = [0.5, 0.5]", "packing_range"),
        ('mode = "informed"\nspace_groups = [0]', "space_groups"),
        ('mode = "informed"\nspace_groups = [231]', "space_groups"),
        ('mode = "informed"\nspace_groups = []', "space_groups"),
        ("n_candidates_per_sg = 0", "n_candidates_per_sg"),
        ("n_candidates_per_sg = 2.5", "n_candidates_per_sg"),
        ("n_prescreen_keep = -1", "n_prescreen_keep"),
        ("relax_steps = -5", "relax_steps"),
        ("max_mol_rmsd = 0", "max_mol_rmsd"),
        ("time_budget_min = -1", "time_budget_min"),
        ("time_budget_min = nan", "time_budget_min"),
        ("n_procs = 0", "n_procs"),
        ("relax_fmax = 0", "relax_fmax"),
        ("overlap_factor = -0.1", "overlap_factor"),
        ('mode = "informed"\nspace_groups = [2]\nz_prime = 0', "z_prime"),
        ('mode = "informed"\nspace_groups = [2]\n[z]\n2 = 0', "z"),
    ],
)
def test_out_of_range_values_rejected(tmp_path, line, key):
    p = tmp_path / "csp.toml"
    p.write_text(f'molecule = "data/benzene.toml"\n{line}\n')
    with pytest.raises(ValueError, match=key):
        load_csp_config(p)


def test_boundary_values_accepted(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text(
        'molecule = "data/benzene.toml"\nmode = "informed"\nspace_groups = [1, 230]\n'
        "packing_range = [0.01, 1.0]\nn_candidates_per_sg = 1\nn_prescreen_keep = 1\nrelax_steps = 0\n"
        "n_procs = 1\nmax_mol_rmsd = 0.01\ntime_budget_min = 0.5\n"
    )
    cfg = load_csp_config(p)
    assert cfg.space_groups == [1, 230] and cfg.relax_steps == 0
