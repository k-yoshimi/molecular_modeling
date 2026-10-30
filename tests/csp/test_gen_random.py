import numpy as np
import pytest
from molecular_modeling.csp.gen_random import default_z, generate_random
from molecular_modeling.csp.molecule import load_molecule, molecule_rmsd

pytest.importorskip("pyxtal")

# ASE 3.x still assigns to ndarray.shape, which NumPy 2.5 deprecates (third-party, harmless).
pytestmark = pytest.mark.filterwarnings("ignore:Setting the shape on a NumPy array:DeprecationWarning")


def test_generate_random_benzene_p21c():
    mol = load_molecule("data/benzene.toml")
    cands = generate_random(mol, 14, 4, n=3, rng=np.random.default_rng(1))
    assert len(cands) == 3
    for c in cands:
        assert c.origin == "random" and c.space_group == 14
        assert len(c.atoms) == 48
        assert sorted(set(c.atoms.arrays["mol_id"])) == [0, 1, 2, 3]
        assert molecule_rmsd(mol, c.atoms[c.atoms.arrays["mol_id"] == 0]) < 1e-3  # rigid
        assert "pyxtal" in c.meta


def test_generate_random_reproducible():
    mol = load_molecule("data/benzene.toml")
    a = generate_random(mol, 14, 4, n=1, rng=np.random.default_rng(7))
    b = generate_random(mol, 14, 4, n=1, rng=np.random.default_rng(7))
    np.testing.assert_allclose(a[0].atoms.get_positions(), b[0].atoms.get_positions())


def test_generate_random_skips_impossible_group(caplog):
    mol = load_molecule("data/benzene.toml")
    cands = generate_random(mol, 230, 1, n=2, rng=np.random.default_rng(1), max_attempts_factor=1)
    assert cands == []
    assert "230" in caplog.text


def test_default_z():
    assert default_z(14, 1) == 4 and default_z(2, 2) == 4
