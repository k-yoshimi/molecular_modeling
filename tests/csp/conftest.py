"""Shared helpers for the CSP tests."""
from pathlib import Path

from molecular_modeling.csp.config import CSPConfig


def default_cfg() -> CSPConfig:
    """Return a small benzene configuration for tests (P-1, 10 candidates, single process)."""
    return CSPConfig(
        molecule=Path("data/benzene.toml"),
        mode="blind",
        space_groups=[2],
        z=None,
        n_candidates_per_sg=10,
        n_procs=1,
        relax_steps=5,
        out_dir=Path("csp_out_test"),
    )
