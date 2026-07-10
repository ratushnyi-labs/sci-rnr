#!/usr/bin/env python3
"""Uniform block-iterator API over the RNR corpus suite.

All experiment code accesses corpora through this module, never through raw
paths, so that (a) every access is checked against MANIFEST.json and (b) the
block segmentation used by the coders (exp-design sections 2 and 5.3,
N in {256, ..., 16384}; engineering-spec sub-block size K) is uniform across
corpora.

    from loader import iter_blocks, corpus_names, corpus_bytes

    for block in iter_blocks("enwik8", K=4096):
        ...                      # successive 4096-byte blocks (last may be short)

    iter_blocks(name, K, drop_last=True)   # exact-K blocks only
    read_range(name, offset, length)       # random-access byte range
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterator

DATA_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = DATA_DIR / "MANIFEST.json"

_manifest_cache: dict | None = None


def load_manifest(refresh: bool = False) -> dict:
    """Parsed MANIFEST.json (cached)."""
    global _manifest_cache
    if _manifest_cache is None or refresh:
        with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
            _manifest_cache = json.load(f)
    return _manifest_cache


def corpus_names(include_deferred: bool = False) -> list[str]:
    """Names of all manifest entries (fetched-and-deferred if requested)."""
    entries = load_manifest()["entries"]
    return [name for name, e in entries.items()
            if include_deferred or not _is_unfetched_deferred(e)]


def entry(name: str) -> dict:
    entries = load_manifest()["entries"]
    if name not in entries:
        raise KeyError(
            f"unknown corpus '{name}'; known: {sorted(entries)}")
    return entries[name]


def _is_unfetched_deferred(e: dict) -> bool:
    return bool(e.get("deferred")) and not (DATA_DIR / e["file"]).exists()


def payload_path(name: str) -> Path:
    """Absolute path of a corpus payload; raises if not materialized."""
    e = entry(name)
    path = DATA_DIR / e["file"]
    if not path.exists():
        hint = (" (deferred entry; run fetch.py --include-deferred)"
                if e.get("deferred") else " (run fetch.py)")
        raise FileNotFoundError(f"payload for '{name}' not materialized{hint}")
    return path


def corpus_bytes(name: str) -> int:
    """Payload size in bytes as recorded in the manifest."""
    return entry(name)["bytes"]


def iter_blocks(name: str, K: int, drop_last: bool = False) -> Iterator[bytes]:
    """Yield successive K-byte blocks of the named corpus.

    The final block may be shorter than K unless drop_last is set.
    Streaming: memory use is O(K) regardless of corpus size.
    """
    if K <= 0:
        raise ValueError(f"block size K must be positive, got {K}")
    path = payload_path(name)
    with open(path, "rb") as f:
        while True:
            block = f.read(K)
            if not block:
                return
            if len(block) < K and drop_last:
                return
            yield block


def num_blocks(name: str, K: int, drop_last: bool = False) -> int:
    """Number of blocks iter_blocks will yield, from manifest size alone."""
    if K <= 0:
        raise ValueError(f"block size K must be positive, got {K}")
    n = corpus_bytes(name)
    return n // K if drop_last else (n + K - 1) // K


def read_range(name: str, offset: int, length: int) -> bytes:
    """Random-access read of [offset, offset+length) from a corpus payload."""
    if offset < 0 or length < 0:
        raise ValueError("offset and length must be non-negative")
    path = payload_path(name)
    with open(path, "rb") as f:
        f.seek(offset)
        return f.read(length)


if __name__ == "__main__":
    # Smoke demonstration: enumerate corpora with sizes and block counts.
    m = load_manifest()
    print(f"{'name':<16} {'class':<28} {'bytes':>12}  blocks@4KiB")
    for name in corpus_names():
        e = entry(name)
        print(f"{name:<16} {e['class']:<28} {e['bytes']:>12}  "
              f"{num_blocks(name, 4096)}")
