"""Command line entry point: ``molecular-modeling-csp [config]``."""

from __future__ import annotations

import argparse
import sys
import tomllib

from .config import load_csp_config
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

    Errors caused by the input (bad settings, missing or unreadable files, missing packages) are
    printed as ``error: <message>`` and give exit code 2; other exceptions keep their traceback.
    An interrupted run (Ctrl-C) exits with code 130.
    """
    parser = argparse.ArgumentParser(
        prog="molecular-modeling-csp", description="Molecular crystal structure prediction"
    )
    parser.add_argument("config", nargs="?", default="input/csp.toml", help="CSP config TOML (default: input/csp.toml)")
    args = parser.parse_args(argv)
    try:
        cfg = load_csp_config(args.config)
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
