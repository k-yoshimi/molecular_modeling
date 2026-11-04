"""Markers and the lock that keep readers of a CSP output directory away from incomplete results."""

from __future__ import annotations

import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

#: Lock file held by a running CSP job or a validation of the directory.
LOCK_NAME = ".csp_run.lock"
#: Present while the results in the directory are being replaced. It is written before old results are
#: deleted and removed only after ``candidates.csv`` of a finished run is in place; a run that fails or
#: is interrupted leaves it behind, so the directory stays marked as incomplete.
INCOMPLETE_NAME = ".csp_incomplete"


@contextmanager
def hold_lock(out_dir: Path, holder: str) -> Iterator[None]:
    """Hold ``out_dir/.csp_run.lock`` (created atomically) while the block runs.

    Raises
    ------
    FileExistsError
        If a CSP run or a validation already holds the directory.
    OSError
        If the lock file cannot be created (e.g. the directory is read-only). There is no unlocked
        fallback: without the lock a concurrent run could change the results while they are read.
    """
    lock = Path(out_dir) / LOCK_NAME
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise FileExistsError(
            f"{out_dir} is in use by a CSP run or a validation (lock file {lock}); wait until it finishes. "
            "If none is active (a run was killed), delete the lock file and rerun molecular-modeling-csp"
        ) from None
    except OSError as exc:
        # same subtype (e.g. PermissionError), so callers that report input errors still recognise it
        raise type(exc)(
            f"cannot create the lock file {lock} ({exc.strerror or exc}); the output directory must be writable"
        ) from exc
    try:
        try:
            os.write(fd, f"{holder} {os.getpid()}\n".encode())
        finally:
            os.close(fd)
        yield
    finally:
        lock.unlink(missing_ok=True)


def check_single_process() -> None:
    """Raise ``ValueError`` when running under MPI with more than one process.

    The lock is taken per process while ASE reads and writes files collectively across MPI ranks, so a
    multi-process run would deadlock. CSP parallelises with worker processes (``n_procs``) instead.
    """
    import ase.parallel

    size = int(ase.parallel.world.size)
    if size > 1:
        raise ValueError(
            f"CSP must not be started under MPI (mpiexec with {size} processes); run it as one process and "
            "use n_procs for parallel relaxation"
        )


def mark_incomplete(out_dir: Path) -> None:
    """Mark the results in ``out_dir`` as incomplete (before anything in them is changed)."""
    (Path(out_dir) / INCOMPLETE_NAME).write_text("results are being replaced; rerun if this file stays\n")


def is_incomplete(out_dir: Path) -> bool:
    """True if the last run into ``out_dir`` did not finish publishing its results."""
    return (Path(out_dir) / INCOMPLETE_NAME).exists()
