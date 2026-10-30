"""Configuration for crystal structure prediction."""
from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Literal

DEFAULT_SPACE_GROUPS: list[int] = [14, 2, 19, 15, 4, 61]
DEFAULT_GENERATORS: dict[str, float] = {"random": 0.6, "voronoi": 0.4}
#: Relative weight given to the seed generator when seed structures are supplied.
SEED_GENERATOR_WEIGHT = 0.5


@dataclass
class SeedSpec:
    """A seed crystal (CIF) with element substitutions.

    Parameters
    ----------
    cif : Path
        CIF file of the parent crystal.
    substitutions : list of (old, new, count)
        Element substitutions: original element, new element, number of
        substituted atoms per molecule.
    volume_scale : float or "auto"
        Volume scaling applied to the seed cell.
    """

    cif: Path
    substitutions: list[tuple[str, str, int]] = field(default_factory=list)
    volume_scale: float | Literal["auto"] = "auto"


@dataclass
class CSPConfig:
    """Settings for a CSP run (see ``input/csp.toml``)."""

    molecule: Path
    mode: Literal["blind", "informed"] = "blind"
    space_groups: list[int] = field(default_factory=lambda: list(DEFAULT_SPACE_GROUPS))
    z: dict[int, int] | None = None
    z_prime: int = 1
    n_candidates_per_sg: int = 500
    generators: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_GENERATORS))
    seed_structures: list[SeedSpec] = field(default_factory=list)
    packing_range: tuple[float, float] = (0.6, 0.8)
    overlap_factor: float = 0.7
    n_prescreen_keep: int = 150
    relax_steps: int = 150
    relax_fmax: float = 0.05
    n_procs: int = field(default_factory=lambda: os.cpu_count() or 1)
    max_mol_rmsd: float = 0.5
    time_budget_min: float = 60
    rng: int = 1
    out_dir: Path = Path("csp_out")


def load_csp_config(path: str | Path) -> CSPConfig:
    """Load a :class:`CSPConfig` from a TOML file.

    In ``blind`` mode ``space_groups`` is replaced by the default list and ``z``
    is ignored. In ``informed`` mode ``space_groups`` is required.

    Raises
    ------
    ValueError
        On unknown keys, bad mode, or ``informed`` mode without ``space_groups``.
    """
    with open(path, "rb") as f:
        raw: dict[str, Any] = tomllib.load(f)

    known = {f.name for f in fields(CSPConfig)}
    unknown = sorted(set(raw) - known)
    if unknown:
        raise ValueError(f"unknown key(s) in {path}: {', '.join(unknown)}")
    if "molecule" not in raw:
        raise ValueError(f"'molecule' is required in {path}")

    mode = raw.get("mode", "blind")
    if mode not in ("blind", "informed"):
        raise ValueError(f"mode must be 'blind' or 'informed', got {mode!r}")

    kwargs: dict[str, Any] = {k: v for k, v in raw.items() if k not in ("seed_structures", "z")}
    kwargs["molecule"] = Path(raw["molecule"])
    if "out_dir" in kwargs:
        kwargs["out_dir"] = Path(kwargs["out_dir"])
    if "packing_range" in kwargs:
        lo, hi = kwargs["packing_range"]
        kwargs["packing_range"] = (float(lo), float(hi))

    if mode == "blind":
        kwargs["space_groups"] = list(DEFAULT_SPACE_GROUPS)
        kwargs["z"] = None
    else:
        if "space_groups" not in raw:
            raise ValueError("informed mode requires 'space_groups'")
        kwargs["z"] = {int(k): int(v) for k, v in raw["z"].items()} if "z" in raw else None

    seeds = [
        SeedSpec(
            cif=Path(s["cif"]),
            substitutions=[(str(a), str(b), int(n)) for a, b, n in s.get("substitutions", [])],
            volume_scale=s.get("volume_scale", "auto"),
        )
        for s in raw.get("seed_structures", [])
    ]
    kwargs["seed_structures"] = seeds

    generators = dict(kwargs.get("generators", DEFAULT_GENERATORS))
    if seeds and "seed" not in generators:
        generators["seed"] = SEED_GENERATOR_WEIGHT
    kwargs["generators"] = generators

    return CSPConfig(**kwargs)
