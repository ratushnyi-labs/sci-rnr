#!/usr/bin/env python3
"""PASS/FAIL gate for the data/big/ scale-benchmark corpus.

Checks (default = full; --quick skips full-payload hashing):
  1. MANIFEST-BIG.json parses and has the required structure.
  2. Size-convention sanity: every tier is exactly tier_mb * 10^6 bytes
     (decimal MB), ladder is 500/1500/5000 (text: 500/1000, capped).
  3. Source downloads that the manifest marks as fetched exist with the
     exact recorded byte size (and, unless --quick, the recorded sha256).
  4. Every materialized master matches its manifest byte size; its sha256
     and ALL tier-prefix sha256s are recomputed from the payload in one
     streaming pass (this verifies the prefix-nesting property directly:
     the tier hash IS the hash of the master's first tier_mb*10^6 bytes).
  5. white_noise determinism: two 64 MiB chunks are regenerated from the
     recorded seed formula and compared byte-for-byte against the payload.
  6. loader_big smoke test on every materialized (type, tier).

Deferred/partial entries are reported, not failed: the manifest honestly
records what exists, and the gate verifies exactly that.  A manifest claim
that is contradicted by the filesystem IS a failure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

BIG_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BIG_DIR))

import fetch_big                                    # noqa: E402
import loader_big                                   # noqa: E402

MB = 10**6
FAILURES: list[str] = []
WARNINGS: list[str] = []


def check(ok: bool, label: str, detail: str = "") -> bool:
    tag = "ok  " if ok else "FAIL"
    print(f"  [{tag}] {label}" + (f" -- {detail}" if detail else ""), flush=True)
    if not ok:
        FAILURES.append(label)
    return ok


def warn(label: str) -> None:
    print(f"  [warn] {label}", flush=True)
    WARNINGS.append(label)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="skip full-payload hashing (size checks only)")
    args = ap.parse_args()

    print("== 1. manifest structure")
    man_path = BIG_DIR / "MANIFEST-BIG.json"
    if not check(man_path.exists(), "MANIFEST-BIG.json exists"):
        return finish()
    man = json.loads(man_path.read_text())
    for key in ("suite", "size_convention", "prefix_nesting", "seed",
                "sources", "entries"):
        check(key in man, f"manifest field '{key}'")
    check(man.get("seed") == fetch_big.SEED, "seed matches fetch_big.SEED")
    check(set(man["entries"]) == set(fetch_big.MASTERS),
          "entry set matches fetch_big.MASTERS",
          ", ".join(sorted(man["entries"])))

    print("== 2. tier ladder / size convention")
    for name, e in man["entries"].items():
        want = [500, 1000] if name == "text" else [500, 1500, 5000]
        check(e["tiers_mb"] == want, f"{name}: ladder {want}")
        for t in e["tiers"]:
            check(t["bytes"] == t["mb"] * MB,
                  f"{name}: tier {t['mb']} MB == {t['mb']}*10^6 bytes")

    print("== 3. sources")
    for s in man["sources"]:
        p = BIG_DIR / "downloads" / s["file"]
        if not s.get("fetched"):
            warn(f"source {s['file']} not fetched (deferred)")
            continue
        if not check(p.exists(), f"source {s['file']} exists"):
            continue
        check(p.stat().st_size == s["bytes"],
              f"source {s['file']} size", f"{p.stat().st_size}")
        if not args.quick and s.get("sha256"):
            got = fetch_big.sha256_file(p)
            check(got == s["sha256"], f"source {s['file']} sha256",
                  got[:16] + "...")

    print("== 4. masters: size + sha256 + prefix nesting")
    materialized: list[tuple[str, int]] = []
    for name, e in man["entries"].items():
        p = BIG_DIR / e["file"]
        if e["bytes"] is None:
            if p.exists():
                check(False, f"{name}: manifest says deferred but payload exists")
            else:
                warn(f"{name}: deferred ({e['status']})")
            continue
        if not check(p.exists(), f"{name}: payload exists ({e['status']})"):
            continue
        check(p.stat().st_size == e["bytes"], f"{name}: size {e['bytes']}",
              str(p.stat().st_size))
        mat_tiers = [t for t in e["tiers"] if t["materialized"]]
        for t in mat_tiers:
            check(t["sha256"] is not None,
                  f"{name}: tier {t['mb']} MB has recorded sha256")
        for t in e["tiers"]:
            if not t["materialized"]:
                warn(f"{name}: tier {t['mb']} MB not materialized "
                     f"({e['status']})")
        if not args.quick:
            cuts = [t["bytes"] for t in mat_tiers] + [e["bytes"]]
            got = fetch_big.prefix_hashes(p, sorted(set(cuts)))
            check(got.get(e["bytes"]) == e["sha256"],
                  f"{name}: master sha256 recomputed")
            for t in mat_tiers:
                check(got.get(t["bytes"]) == t["sha256"],
                      f"{name}: tier {t['mb']} MB prefix sha256 recomputed "
                      "(prefix nesting)")
        materialized.extend((name, t["mb"]) for t in mat_tiers)

    print("== 5. white_noise determinism (seeded regeneration)")
    wn = BIG_DIR / fetch_big.MASTERS["white_noise"]
    if wn.exists():
        with open(wn, "rb") as f:
            for ci in (0, 57):
                f.seek(ci * fetch_big.CHUNK)
                disk = f.read(fetch_big.CHUNK)
                regen = fetch_big.white_noise_chunk(ci, len(disk))
                check(disk == regen, f"white_noise chunk {ci} regenerates")
    else:
        warn("white_noise payload absent; determinism check skipped")

    print("== 6. content spot checks")
    def head(path: Path, n: int, off: int = 0) -> bytes:
        with open(path, "rb") as f:
            f.seek(off)
            return f.read(n)

    m_text = BIG_DIR / fetch_big.MASTERS["text"]
    if m_text.exists():
        check(head(m_text, 10) == b"<mediawiki",
              "text: starts with <mediawiki (enwik9)")
    m_img = BIG_DIR / fetch_big.MASTERS["images"]
    if m_img.exists():
        check(head(m_img, 2) == b"\xff\xd8",
              "images: starts with JPEG SOI marker")
    m_pre = BIG_DIR / fetch_big.MASTERS["precompressed"]
    if m_pre.exists():
        check(head(m_pre, 4) == b"\x28\xb5\x2f\xfd",
              "precompressed: starts with zstd magic")
    m_vid = BIG_DIR / fetch_big.MASTERS["video"]
    if m_vid.exists():
        # concat order: master must equal each source at its recorded offset
        off = 0
        ok = True
        for src_name in fetch_big.VIDEO_CONCAT_ORDER:
            p = BIG_DIR / "downloads" / src_name
            if not p.exists() or off >= fetch_big.MASTER_BYTES:
                break
            n = min(1 << 20, fetch_big.MASTER_BYTES - off)
            ok &= head(m_vid, n, off) == head(p, n)
            off += p.stat().st_size
        check(ok, "video: every source file starts at its concat offset")
    m_tar = BIG_DIR / fetch_big.MASTERS["archive_mix"]
    if m_tar.exists():
        import tarfile
        names = []
        try:
            with tarfile.open(m_tar, mode="r|") as tf:
                for i, ti in enumerate(tf):
                    names.append(ti.name)
                    if i >= 10:
                        break
        except tarfile.TarError:
            pass  # truncated tail is expected; header region must parse
        check(len(names) >= 9 and names[0].startswith("mix/small/"),
              "archive_mix: tar header region lists small payloads",
              f"{len(names)} members seen")

    print("== 7. loader_big smoke")
    loader_big.load_manifest(refresh=True)
    for name, mb in materialized:
        try:
            it = loader_big.iter_blocks_big(name, mb, 4096)
            first = next(it)
            ok = len(first) == 4096
            rr = loader_big.read_range_big(name, mb, 0, 4096)
            ok &= rr == first
            n = loader_big.tier_bytes(name, mb)
            ok &= n == mb * MB
            check(ok, f"loader_big {name}@{mb}MB block/range/size")
        except Exception as exc:
            check(False, f"loader_big {name}@{mb}MB", repr(exc))

    return finish()


def finish() -> int:
    print()
    if WARNINGS:
        print(f"deferred/partial (not failures): {len(WARNINGS)}")
        for w in WARNINGS:
            print(f"  - {w}")
    if FAILURES:
        print(f"RESULT: FAIL ({len(FAILURES)} failed checks)")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
