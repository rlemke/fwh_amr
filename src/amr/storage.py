"""Backend-aware paths for the amr domain.

⚠️ A THIN SHIM over ``facetwork.domains.storage``. This module was one of 21
near-copies across the fwh_* repos; the shared layer owns the behaviour now and
this file keeps the import path, the public names, and anything genuinely
specific to amr.

⚠️ The arguments below pin where this domain's data ALREADY sits in the object
store. They are not style — changing one orphans that data rather than moving
it — and they were verified against the previous module across local, s3:// and
hdfs:// before the switch.
"""
from __future__ import annotations
import contextlib
import os
import tempfile
from collections.abc import Iterator
from typing import IO
from facetwork.config import get_output_base
from facetwork.runtime import storage as _fws

from facetwork.domains.storage import domain_storage, is_remote, join  # noqa: F401

_S = domain_storage("amr")


def cache_root() -> str:
    return _S.cache_root()


def output_root() -> str:
    return _S.output_root()


def exists(path: str) -> bool:
    return _S.exists(path)


def localize(path: str) -> str:
    return _S.localize(path)


def list_dir(path: str) -> list[str]:
    return _S.list_dir(path)


def open_read(path: str, mode: str = "r", **kw):
    return _S.open_read(path, mode, **kw)


def open_write(path: str, mode: str = "w", **kw):
    return _S.open_write(path, mode, **kw)


def _data_root() -> str:
    return os.environ.get("FW_DATA_ROOT") or get_output_base()


def organisms_root() -> str:
    """Shared per-organism aggregate dir (fan-out writes here, fan-in reads)."""
    return join(cache_root(), "organisms")
