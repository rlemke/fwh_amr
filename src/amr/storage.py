"""Backend-aware paths for the AMR cache + outputs.

On the fleet (``FW_STORAGE=s3`` / ``FW_DATA_ROOT=s3://afl-cache``) the cached
per-organism aggregates + rendered artifacts land in the shared MinIO object
store. A thin wrapper over ``facetwork.runtime.storage`` (the same shape
census-us / save-earth / migration use), so terminal use and fleet runs share
one cache rooted at ``$FW_DATA_ROOT/cache/amr/``.

Per-organism aggregates live under a **shared** namespace
(``cache/amr/cache/organisms/<group>.json``) so a fan-out of one
``DownloadOrganism`` per pathogen writes side-by-side artifacts that the fan-in
map/table steps read together — the same shared-national-cache pattern the
county atlas uses.
"""

from __future__ import annotations

import contextlib
import os
import tempfile
from collections.abc import Iterator
from typing import IO

from facetwork.config import get_output_base
from facetwork.runtime import storage as _fws


def is_remote(path: str) -> bool:
    return "://" in (path or "")


def _data_root() -> str:
    return os.environ.get("FW_DATA_ROOT") or get_output_base()


def join(*parts: str) -> str:
    parts = [p for p in parts if p]
    if not parts:
        return ""
    base = parts[0].rstrip("/")
    rest = [p.strip("/") for p in parts[1:]]
    return "/".join([base, *[p for p in rest if p]])


def cache_root() -> str:
    ov = os.environ.get("FW_AMR_CACHE_DIR")
    if ov:
        return ov
    r = _data_root()
    return join(r, "cache", "amr", "cache") if is_remote(r) else join(r, "amr-cache")


def organisms_root() -> str:
    """Shared per-organism aggregate dir (fan-out writes here, fan-in reads)."""
    return join(cache_root(), "organisms")


def output_root() -> str:
    ov = os.environ.get("FW_AMR_OUTPUT_DIR")
    if ov:
        return ov
    r = _data_root()
    return join(r, "cache", "amr", "output") if is_remote(r) else join(r, "amr-output")


def exists(path: str) -> bool:
    return _fws.get_storage_backend(path).exists(path)


def localize(path: str) -> str:
    if not is_remote(path):
        return path
    return _fws.localize(path)


def list_dir(path: str) -> list[str]:
    """Best-effort list of immediate entries under ``path`` (full paths)."""
    backend = _fws.get_storage_backend(path)
    lister = getattr(backend, "list", None) or getattr(backend, "listdir", None)
    if lister is not None:
        try:
            return [join(path, os.path.basename(p.rstrip("/"))) for p in lister(path)]
        except Exception:
            pass
    if not is_remote(path) and os.path.isdir(path):
        return [join(path, n) for n in os.listdir(path)]
    return []


def open_read(path: str, mode: str = "r", **kw) -> IO:
    return open(localize(path), mode, **kw)


@contextlib.contextmanager
def open_write(path: str, mode: str = "w", **kw) -> Iterator[IO]:
    if not is_remote(path):
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        f = open(path, mode, **kw)
        try:
            yield f
        finally:
            f.close()
        return

    fd, tmp = tempfile.mkstemp(suffix="_" + os.path.basename(path))
    os.close(fd)
    f = open(tmp, mode, **kw)
    try:
        yield f
        f.close()
        with open(tmp, "rb") as src, _fws.get_storage_backend(path).open(path, "wb") as dst:
            dst.write(src.read())
    finally:
        if not f.closed:
            f.close()
        try:
            os.unlink(tmp)
        except OSError:
            pass
