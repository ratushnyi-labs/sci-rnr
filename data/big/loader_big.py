#!/usr/bin/env python3
"""Tier-aware block-iterator API over the scale-benchmark corpus (data/big/).

Mirrors data/loader.py, with one addition: every access names a TIER
(decimal MB: 500 / 1500 / 5000; text: 500 / 1000).  Tiers are byte-prefixes
of a single per-type master payload (prefix nesting), so the loader simply
limits reads to the first tier_mb * 10^6 bytes of the master.

    from loader_big import iter_blocks_big, big_names, tier_bytes

    for block in iter_blocks_big("video", tier_mb=1500, K=4096):
        ...

    read_range_big("white_noise", 500, offset, length)   # random access
    tiers_available("text")                              # -> [500, 1000]

All access is checked against MANIFEST-BIG.json; payloads that are not
materialized (status 'deferred'/'partial') raise FileNotFoundError with the
manifest status, so benchmarks fail loudly instead of measuring air.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterator

BIG_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = BIG_DIR / "MANIFEST-BIG.json"
MB = 10**6

_manifest_cache: dict | None = None


def load_manifest(refresh: bool = False) -> dict:
    global _manifest_cache
    if _manifest_cache is None or refresh:
        with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
            _manifest_cache = json.load(f)
    return _manifest_cache


def big_names() -> list[str]:
    """All content-type names in the big manifest."""
    return list(load_manifest()["entries"])


def entry(name: str) -> dict:
    entries = load_manifest()["entries"]
    if name not in entries:
        raise KeyError(f"unknown big corpus '{name}'; known: {sorted(entries)}")
    return entries[name]


def tiers_available(name: str, materialized_only: bool = True) -> list[int]:
    """Tier sizes (decimal MB) defined -- or materialized -- for a type."""
    return [t["mb"] for t in entry(name)["tiers"]
            if t["materialized"] or not materialized_only]


def tier_bytes(name: str, tier_mb: int) -> int:
    """Exact byte length of a tier; raises if the tier is not defined."""
    for t in entry(name)["tiers"]:
        if t["mb"] == tier_mb:
            return t["bytes"]
    raise KeyError(f"tier {tier_mb} MB not defined for '{name}'; "
                   f"defined: {tiers_available(name, materialized_only=False)}")


def payload_path(name: str, tier_mb: int) -> Path:
    """Path of the master payload backing a tier; validates materialization."""
    e = entry(name)
    n = tier_bytes(name, tier_mb)
    path = BIG_DIR / e["file"]
    if not path.exists():
        raise FileNotFoundError(
            f"big corpus '{name}' not materialized (manifest status: "
            f"{e['status']}); run fetch_big.py")
    if path.stat().st_size < n:
        raise FileNotFoundError(
            f"big corpus '{name}' master is {path.stat().st_size} bytes; "
            f"tier {tier_mb} MB needs {n} (manifest status: {e['status']})")
    return path


def iter_blocks_big(name: str, tier_mb: int, K: int,
                    drop_last: bool = False) -> Iterator[bytes]:
    """Successive K-byte blocks of the first tier_mb*10^6 bytes of a master.

    Streaming; memory is O(K).  The final block may be short unless
    drop_last is set; use num_blocks_big for the exact count."""
    if K <= 0:
        raise ValueError(f"block size K must be positive, got {K}")
    path = payload_path(name, tier_mb)
    remaining = tier_bytes(name, tier_mb)
    with open(path, "rb") as f:
        while remaining > 0:
            block = f.read(min(K, remaining))
            if not block:
                return
            remaining -= len(block)
            if len(block) < K and drop_last:
                return
            yield block


def num_blocks_big(name: str, tier_mb: int, K: int,
                   drop_last: bool = False) -> int:
    if K <= 0:
        raise ValueError(f"block size K must be positive, got {K}")
    n = tier_bytes(name, tier_mb)
    return n // K if drop_last else (n + K - 1) // K


def read_range_big(name: str, tier_mb: int, offset: int, length: int) -> bytes:
    """Random-access read of [offset, offset+length) within a tier prefix."""
    if offset < 0 or length < 0:
        raise ValueError("offset and length must be non-negative")
    n = tier_bytes(name, tier_mb)
    if offset + length > n:
        raise ValueError(f"range [{offset}, {offset+length}) exceeds tier "
                         f"{tier_mb} MB = {n} bytes")
    path = payload_path(name, tier_mb)
    with open(path, "rb") as f:
        f.seek(offset)
        return f.read(length)


if __name__ == "__main__":
    m = load_manifest()
    print(f"{'name':<14} {'status':<34} tiers(MB): materialized/defined")
    for name in big_names():
        e = entry(name)
        mat = tiers_available(name)
        allt = tiers_available(name, materialized_only=False)
        print(f"{name:<14} {e['status']:<34} {mat} / {allt}")
