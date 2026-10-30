"""Command line entry point: ``molecular-modeling-csp [config]``."""
from __future__ import annotations

import argparse

from .config import load_csp_config
from .pipeline import run_csp


def main(argv: list[str] | None = None) -> int:
    """Run CSP from a TOML config (default ``input/csp.toml``); return the exit code."""
    parser = argparse.ArgumentParser(prog="molecular-modeling-csp", description="Molecular crystal structure prediction")
    parser.add_argument("config", nargs="?", default="input/csp.toml", help="CSP config TOML (default: input/csp.toml)")
    args = parser.parse_args(argv)
    cfg = load_csp_config(args.config)
    run_csp(cfg)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
