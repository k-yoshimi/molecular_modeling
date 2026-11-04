"""Configuration for crystal structure prediction."""

from __future__ import annotations

import logging
import math
import os
import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Literal

logger = logging.getLogger(__name__)

DEFAULT_SPACE_GROUPS: list[int] = [14, 2, 19, 15, 4, 61]
DEFAULT_GENERATORS: dict[str, float] = {"random": 0.6, "voronoi": 0.4}
#: Relative weight given to the seed generator when seed structures are supplied. For the seed
#: generator the weight is an on/off switch: > 0 uses the seeds, 0 ignores them.
SEED_GENERATOR_WEIGHT = 0.5
#: Generator names accepted in ``generators``.
GENERATOR_NAMES = ("random", "voronoi", "seed")


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
    n_candidates_per_sg: int = 150
    generators: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_GENERATORS))
    seed_structures: list[SeedSpec] = field(default_factory=list)
    packing_range: tuple[float, float] = (0.4, 0.8)
    overlap_factor: float = 0.7
    n_prescreen_keep: int = 32
    relax_steps: int = 60
    relax_fmax: float = 0.05
    n_procs: int = field(default_factory=lambda: os.cpu_count() or 1)
    max_mol_rmsd: float = 0.5
    time_budget_min: float = 60
    rng: int = 1
    out_dir: Path = Path("csp_out")
    #: Warnings issued while loading the file (not an input key); ``run_csp`` repeats them in ``run.log``.
    config_warnings: list[str] = field(default_factory=list, repr=False, compare=False)


#: Fields of :class:`CSPConfig` that are set by the loader, not read from the TOML file.
_NOT_INPUT_KEYS = {"config_warnings"}


def _is_int(x: Any) -> bool:
    return isinstance(x, int) and not isinstance(x, bool)


def _is_real(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _check_ranges(cfg: CSPConfig) -> None:
    """Raise ``ValueError`` naming the first setting outside its allowed range."""
    lo, hi = cfg.packing_range
    if not (_is_real(lo) and _is_real(hi) and 0 < lo < hi <= 1):
        raise ValueError(f"packing_range must be [lo, hi] with 0 < lo < hi <= 1, got {list(cfg.packing_range)}")
    if not cfg.space_groups or not all(_is_int(g) and 1 <= g <= 230 for g in cfg.space_groups):
        raise ValueError(f"space_groups must be a non-empty list of integers 1..230, got {cfg.space_groups}")
    if cfg.z is not None and not all(1 <= g <= 230 and v >= 1 for g, v in cfg.z.items()):
        raise ValueError(f"z must map space groups 1..230 to positive integers, got {cfg.z}")
    for name in ("z_prime", "n_candidates_per_sg", "n_prescreen_keep", "n_procs"):
        v = getattr(cfg, name)
        if not (_is_int(v) and v >= 1):
            raise ValueError(f"{name} must be a positive integer, got {v!r}")
    if not (_is_int(cfg.relax_steps) and cfg.relax_steps >= 0):
        raise ValueError(f"relax_steps must be a non-negative integer, got {cfg.relax_steps!r}")
    for name in ("relax_fmax", "max_mol_rmsd", "time_budget_min", "overlap_factor"):
        v = getattr(cfg, name)
        if not (_is_real(v) and v > 0):
            raise ValueError(f"{name} must be a positive finite number, got {v!r}")


def load_csp_config(path: str | Path) -> CSPConfig:
    """Load a :class:`CSPConfig` from a TOML file.

    In ``blind`` mode ``space_groups`` is replaced by the default list, ``z`` is ignored and
    ``z_prime`` is set to 1 (with a warning if the file says otherwise); ``seed_structures`` are
    still used, since a seed is a related (parent) crystal, not the target's experimental structure.
    In ``informed`` mode ``space_groups`` is required. Warnings are logged and also kept in
    ``config_warnings``.

    Raises
    ------
    ValueError
        On unknown keys, bad mode, ``informed`` mode without ``space_groups``, unknown generator
        names / negative or non-finite generator weights, no generator enabled, or a setting outside
        its range: ``0 < packing_range[0] < packing_range[1] <= 1``; space groups (and ``z`` keys)
        1..230; positive integers ``z`` values, ``z_prime``, ``n_candidates_per_sg``,
        ``n_prescreen_keep``, ``n_procs``; non-negative integer ``relax_steps``; positive finite
        ``relax_fmax``, ``max_mol_rmsd``, ``time_budget_min``, ``overlap_factor``.
    """
    warnings: list[str] = []

    def warn(msg: str) -> None:
        logger.warning("%s", msg)
        warnings.append(msg)

    with open(path, "rb") as f:
        raw: dict[str, Any] = tomllib.load(f)

    known = {f.name for f in fields(CSPConfig)} - _NOT_INPUT_KEYS
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
        pr = kwargs["packing_range"]
        # checked before conversion: float(True) would silently give 1.0
        if not (isinstance(pr, list) and len(pr) == 2 and all(_is_real(x) for x in pr)):
            raise ValueError(f"packing_range must be two finite numbers [lo, hi], got {pr!r}")
        kwargs["packing_range"] = (float(pr[0]), float(pr[1]))

    if mode == "blind":
        # blind uses no information about the target crystal: default space groups, Z' = 1
        kwargs["space_groups"] = list(DEFAULT_SPACE_GROUPS)
        kwargs["z"] = None
        if "space_groups" in raw and list(raw["space_groups"]) != DEFAULT_SPACE_GROUPS:
            warn(
                f"blind mode: space_groups = {raw['space_groups']} is ignored and the default list "
                f"{DEFAULT_SPACE_GROUPS} is used; use informed mode to choose space groups"
            )
        if "z" in raw:
            warn("blind mode: 'z' is ignored (Z is the general-position multiplicity)")
        if raw.get("z_prime", 1) != 1:
            warn(f"blind mode: z_prime = {raw['z_prime']} is ignored and set to 1; use informed mode for Z' > 1")
        kwargs["z_prime"] = 1
    else:
        if "space_groups" not in raw:
            raise ValueError("informed mode requires 'space_groups'")
        if "z" in raw:
            # values checked before conversion: int(1.9) and int(True) would silently give 1
            bad = {k: v for k, v in raw["z"].items() if not _is_int(v)}
            if bad:
                raise ValueError(f"z values must be integers, got {bad}")
            try:
                kwargs["z"] = {int(k): v for k, v in raw["z"].items()}
            except ValueError:
                raise ValueError(f"z keys must be space group numbers, got {list(raw['z'])}") from None
        else:
            kwargs["z"] = None

    seeds = [
        SeedSpec(
            cif=Path(s["cif"]),
            substitutions=[(str(a), str(b), int(n)) for a, b, n in s.get("substitutions", [])],
            volume_scale=s.get("volume_scale", "auto"),
        )
        for s in raw.get("seed_structures", [])
    ]
    kwargs["seed_structures"] = seeds

    # the file's table is merged over the defaults: a name left out keeps its default weight, so a
    # partial table (e.g. only ``seed = 0``) cannot switch other generators off by omission
    generators = {**DEFAULT_GENERATORS, **kwargs.get("generators", {})}
    bad_names = sorted(set(generators) - set(GENERATOR_NAMES))
    if bad_names:
        raise ValueError(f"unknown generators {bad_names}; use {', '.join(GENERATOR_NAMES)}")
    for name, w in generators.items():
        if not isinstance(w, (int, float)) or isinstance(w, bool) or not math.isfinite(w) or w < 0:
            raise ValueError(f"generators: weight of {name!r} must be a finite non-negative number, got {w!r}")
    if seeds and "seed" not in generators:
        generators["seed"] = SEED_GENERATOR_WEIGHT
    if generators.get("random", 0) + generators.get("voronoi", 0) <= 0 and not (
        seeds and generators.get("seed", 0) > 0
    ):
        raise ValueError("generators: no generator is enabled (random and voronoi are 0 and no seeds are used)")
    kwargs["generators"] = generators

    cfg = CSPConfig(**kwargs)
    _check_ranges(cfg)
    cfg.config_warnings = warnings
    return cfg
