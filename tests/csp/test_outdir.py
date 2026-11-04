"""Tests for the output-directory lock and markers."""

import os

import pytest
from molecular_modeling.csp import outdir


def test_lock_descriptor_closed_and_lock_removed_when_metadata_write_fails(tmp_path, monkeypatch):
    closed = []
    real_close = os.close

    def failing_write(fd, data):
        raise OSError("disk full")

    def spy_close(fd):
        closed.append(fd)
        return real_close(fd)

    monkeypatch.setattr(os, "write", failing_write)
    monkeypatch.setattr(os, "close", spy_close)
    with pytest.raises(OSError, match="disk full"):
        with outdir.hold_lock(tmp_path, "test"):
            pass
    assert len(closed) == 1
    assert not (tmp_path / outdir.LOCK_NAME).exists()


def test_lock_creation_failure_is_reported(tmp_path, monkeypatch):
    def no_create(path, flags, *args):
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr(os, "open", no_create)
    with pytest.raises(OSError, match="lock file"):
        with outdir.hold_lock(tmp_path, "test"):
            pass
