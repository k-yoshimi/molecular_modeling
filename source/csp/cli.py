"""Command line entry point: ``molecular-modeling-csp [config] [--relax-method METHOD]``."""

from __future__ import annotations

import argparse
import sys
import tomllib

from .config import CONFIG_RELAX_METHODS, load_csp_config
from .pipeline import run_csp

#: Errors caused by the input or the environment: reported as one line, without a traceback.
EXPECTED_ERRORS = (
    ValueError,
    ImportError,
    FileNotFoundError,
    FileExistsError,
    PermissionError,
    IsADirectoryError,
    NotADirectoryError,
)


def main(argv: list[str] | None = None) -> int:
    """Run CSP from a TOML config (default ``input/csp.toml``); return the exit code.

    ``--relax-method`` overrides ``relax_method`` of the config file.

    Errors caused by the input (bad settings, missing or unreadable files, missing packages) are
    printed as ``error: <message>`` and give exit code 2; other exceptions keep their traceback.
    An interrupted run (Ctrl-C) exits with code 130.
    """
    parser = argparse.ArgumentParser(
        prog="molecular-modeling-csp", description="Molecular crystal structure prediction"
    )
    parser.add_argument("config", nargs="?", default="input/csp.toml", help="CSP config TOML (default: input/csp.toml)")
    parser.add_argument(
        "--relax-method",
        choices=CONFIG_RELAX_METHODS,
        help="relaxation method; overrides relax_method in the config file",
    )
    args = parser.parse_args(argv)
    try:
        cfg = load_csp_config(args.config)
        if args.relax_method is not None:
            cfg.relax_method = args.relax_method
            cfg.relax_method_source = "cli"
        run_csp(cfg)
    except tomllib.TOMLDecodeError as exc:
        print(f"error: {args.config}: invalid TOML: {exc}", file=sys.stderr)
        return 2
    except EXPECTED_ERRORS as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("error: interrupted", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
