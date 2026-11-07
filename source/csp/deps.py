"""Check that the optional packages of the ``csp`` extra can be used."""

from __future__ import annotations

import importlib

INSTALL_HINT = 'install the CSP dependencies with: pip install -e ".[csp]"'


def require(
    modules: list[str],
    attributes: dict[str, str] | None = None,
    purpose: str = "this command",
    hint: str = INSTALL_HINT,
) -> None:
    """Raise ``ImportError`` naming every module in ``modules`` that cannot be imported.

    ``attributes`` maps a module to a name it must provide (e.g. ``{"tblite.ase": "TBLite"}``); a module
    without it counts as missing. The message ends with ``hint`` (by default the ``pip install -e ".[csp]"`` hint).
    """
    attributes = attributes or {}
    missing = []
    for name in modules:
        try:
            mod = importlib.import_module(name)
        except ImportError:
            missing.append(name)
            continue
        attr = attributes.get(name)
        if attr is not None and not hasattr(mod, attr):
            missing.append(f"{name}.{attr}")
    if missing:
        raise ImportError(f"{purpose} needs {', '.join(missing)}, which cannot be imported; {hint}")
