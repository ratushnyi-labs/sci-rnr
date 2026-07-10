#!/usr/bin/env python3
"""Executable checks for the RNR corpus suite (data/).

Verifies, in the style of scripts/verify/:
  C1  manifest schema: MANIFEST.json parses and every entry carries the
      fields required by exp-design section 3 / section 12 (name, class,
      bytes, sha256, source URL, license note, hypothesis qualification).
  C2  payload integrity: every non-deferred entry exists on disk with the
      byte size and sha256 recorded in the manifest.
  C3  loader round-trip: iter_blocks() reassembles each payload exactly
      (total size, block count, and prefix content), for two block sizes,
      including drop_last accounting.
  C4  class coverage: at least one materialized corpus in each exp-design
      class 3.1 (text/source), 3.2 (structured binary), 3.3 (image/
      numeric), 3.4 (heterogeneous + negative controls).
  C5  deferred entries: enwik9 is manifest-only (URL + deferred flag) and
      the loader refuses it cleanly while unfetched.
  C6  content sanity: each payload looks like what its class claims
      (tar readability, SQLite magic + row counts, telemetry spatial
      correlation, negative controls incompressible, base64 decodes to its
      seeded source, enwik8 header), skipped for stand-ins.

Exit code 0 iff every check passes.
"""

from __future__ import annotations

import base64
import hashlib
import json
import random
import sqlite3
import sys
import tarfile
import zlib
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(DATA_DIR))

import loader  # noqa: E402

MiB = 1 << 20
RESULTS: list[tuple[str, bool, str]] = []


def check(cid: str, ok: bool, detail: str) -> bool:
    RESULTS.append((cid, ok, detail))
    print(f"{'PASS' if ok else 'FAIL'}  {cid}: {detail}")
    return ok


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(MiB)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def zlib_ratio(data: bytes, level: int = 6) -> float:
    return len(zlib.compress(data, level)) / len(data)


REQUIRED_FIELDS = ("class", "file", "source_url", "license_note",
                   "qualifies_for", "provenance")
KNOWN_CLASSES = ("3.1-text-source", "3.2-structured-binary",
                 "3.3-image-numeric", "3.4-heterogeneous-negative")


def c1_manifest_schema(manifest: dict) -> None:
    entries = manifest.get("entries", {})
    check("C1.count", len(entries) >= 10,
          f"manifest has {len(entries)} entries (expect >= 10)")
    for name, e in entries.items():
        missing = [f for f in REQUIRED_FIELDS if f not in e]
        deferred_unfetched = bool(e.get("deferred")) and "sha256" not in e
        if not deferred_unfetched:
            missing += [f for f in ("bytes", "sha256") if f not in e]
        check(f"C1.fields[{name}]", not missing,
              "all required fields present" if not missing
              else f"missing fields: {missing}")
        check(f"C1.class[{name}]", e.get("class") in KNOWN_CLASSES,
              f"class '{e.get('class')}' is a known section-3 class")
        check(f"C1.qualifies[{name}]",
              isinstance(e.get("qualifies_for"), list)
              and len(e["qualifies_for"]) >= 1
              and all(q.startswith("H") for q in e["qualifies_for"]),
              f"qualifies_for = {e.get('qualifies_for')}")


def c2_payload_integrity(manifest: dict) -> None:
    for name, e in manifest["entries"].items():
        path = DATA_DIR / e["file"]
        if e.get("deferred") and not path.exists():
            continue  # handled by C5
        if not check(f"C2.exists[{name}]", path.exists(),
                     f"{e['file']} exists"):
            continue
        size = path.stat().st_size
        check(f"C2.bytes[{name}]", size == e["bytes"],
              f"size {size} == manifest {e['bytes']}")
        digest = sha256_file(path)
        check(f"C2.sha256[{name}]", digest == e["sha256"],
              f"sha256 {digest[:16]}... matches manifest"
              if digest == e["sha256"]
              else f"sha256 {digest[:16]}... != manifest {e['sha256'][:16]}...")


def c3_loader_roundtrip(manifest: dict) -> None:
    for name, e in manifest["entries"].items():
        path = DATA_DIR / e["file"]
        if e.get("deferred") and not path.exists():
            continue
        n = e["bytes"]
        for K in (65536, 1 << 20):
            total, count, first = 0, 0, None
            for block in loader.iter_blocks(name, K):
                if first is None:
                    first = block
                total += len(block)
                count += 1
            ok_size = total == n
            ok_count = count == loader.num_blocks(name, K)
            with open(path, "rb") as f:
                ok_prefix = first == f.read(min(K, n))
            check(f"C3.roundtrip[{name},K={K}]",
                  ok_size and ok_count and ok_prefix,
                  f"{count} blocks, {total} bytes reassembled"
                  + ("" if ok_size and ok_count and ok_prefix else
                     f" (size={ok_size} count={ok_count} prefix={ok_prefix})"))
        # drop_last accounting
        K = 65536
        kept = sum(1 for _ in loader.iter_blocks(name, K, drop_last=True))
        check(f"C3.drop_last[{name}]",
              kept == n // K == loader.num_blocks(name, K, drop_last=True),
              f"drop_last yields {kept} == floor({n}/{K})")


def c4_class_coverage(manifest: dict) -> None:
    materialized: dict[str, list[str]] = {c: [] for c in KNOWN_CLASSES}
    for name, e in manifest["entries"].items():
        if (DATA_DIR / e["file"]).exists():
            materialized.setdefault(e["class"], []).append(name)
    for cls in KNOWN_CLASSES:
        names = materialized.get(cls, [])
        check(f"C4.coverage[{cls}]", len(names) >= 1,
              f"materialized corpora: {names or 'NONE'}")


def c5_deferred(manifest: dict) -> None:
    e = manifest["entries"].get("enwik9")
    ok = (e is not None and e.get("deferred") == "large"
          and e.get("source_url", "").startswith("http"))
    check("C5.enwik9_manifest", ok,
          "enwik9 present as deferred=large with source URL")
    if e is not None and not (DATA_DIR / e["file"]).exists():
        try:
            loader.payload_path("enwik9")
            check("C5.enwik9_loader", False,
                  "loader returned a path for an unfetched deferred entry")
        except FileNotFoundError as exc:
            check("C5.enwik9_loader", "include-deferred" in str(exc),
                  "loader raises FileNotFoundError with fetch hint")


def _entry_live(manifest, name):
    e = manifest["entries"].get(name)
    if e is None or not (DATA_DIR / e["file"]).exists():
        return None
    return e


def c6_content_sanity(manifest: dict) -> None:
    # enwik8: Wikipedia XML dump prefix.
    e = _entry_live(manifest, "enwik8")
    if e and not e.get("standin"):
        head = loader.read_range("enwik8", 0, 32)
        check("C6.enwik8_header", head.startswith(b"<mediawiki"),
              f"payload starts with {head[:12]!r}")

    # Calgary / Canterbury: readable tar with the classic member counts.
    for name, min_members in (("calgary", 14), ("canterbury", 11)):
        e = _entry_live(manifest, name)
        if e and not e.get("standin"):
            with tarfile.open(DATA_DIR / e["file"]) as tar:
                members = [m for m in tar.getmembers() if m.isfile()]
            check(f"C6.{name}_tar", len(members) >= min_members,
                  f"tar readable with {len(members)} files "
                  f"(expect >= {min_members})")

    # scripts_src: deterministic tar of the repo's scripts/ directory.
    e = _entry_live(manifest, "scripts_src")
    if e:
        with tarfile.open(DATA_DIR / e["file"]) as tar:
            names = tar.getnames()
        py = [n for n in names if n.endswith(".py")]
        ok = len(py) >= 50 and all(n.startswith("rnr_scripts/") for n in names)
        check("C6.scripts_tar", ok,
              f"{len(py)} .py members under rnr_scripts/ (expect >= 50)")

    # sqlite_synth: magic + row counts match the seeded parameters.
    e = _entry_live(manifest, "sqlite_synth")
    if e:
        head = loader.read_range("sqlite_synth", 0, 16)
        check("C6.sqlite_magic", head == b"SQLite format 3\x00",
              f"header {head!r}")
        con = sqlite3.connect(f"file:{DATA_DIR / e['file']}?mode=ro", uri=True)
        try:
            counts = {t: con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
                      for t in ("sensors", "readings", "logs")}
        finally:
            con.close()
        p = e["params"]
        ok = (counts["sensors"] == p["sensors"]
              and counts["readings"] == p["readings"]
              and counts["logs"] == p["logs"])
        check("C6.sqlite_rows", ok, f"row counts {counts} match params")

    # telemetry_grid: exact size, finite values, strong spatial correlation.
    e = _entry_live(manifest, "telemetry_grid")
    if e:
        import numpy as np
        n0, n1 = e["params"]["shape"]
        expected = n0 * n1 * 4
        check("C6.telemetry_size", e["bytes"] == expected,
              f"{e['bytes']} bytes == {n0}x{n1} float32")
        grid = np.frombuffer(
            (DATA_DIR / e["file"]).read_bytes(), dtype="<f4").reshape(n0, n1)
        finite = bool(np.isfinite(grid).all())
        check("C6.telemetry_finite", finite, "all values finite")
        g = grid - grid.mean()
        lag1 = float((g[:, :-1] * g[:, 1:]).sum()
                     / np.sqrt((g[:, :-1] ** 2).sum() * (g[:, 1:] ** 2).sum()))
        check("C6.telemetry_correlated", lag1 > 0.9,
              f"lag-1 spatial autocorrelation {lag1:.4f} > 0.9")

    # ctrl_urandom: incompressible under zlib.
    e = _entry_live(manifest, "ctrl_urandom")
    if e:
        ratio = zlib_ratio(loader.read_range("ctrl_urandom", 0, MiB))
        check("C6.urandom_incompressible", ratio > 0.99,
              f"zlib ratio {ratio:.4f} > 0.99 on first MiB")

    # ctrl_zstd: zstd magic, incompressible residual, exact round-trip to
    # the enwik8 prefix it was derived from.
    e = _entry_live(manifest, "ctrl_zstd")
    if e:
        head = loader.read_range("ctrl_zstd", 0, 4)
        check("C6.zstd_magic", head == b"\x28\xb5\x2f\xfd",
              f"magic {head.hex()}")
        ratio = zlib_ratio(loader.read_range("ctrl_zstd", 0, MiB))
        check("C6.zstd_incompressible", ratio > 0.95,
              f"zlib ratio {ratio:.4f} > 0.95 on first MiB")
        src = _entry_live(manifest, e["params"]["input"])
        if src:
            import zstandard
            blob = (DATA_DIR / e["file"]).read_bytes()
            raw = zstandard.ZstdDecompressor().decompress(
                blob, max_output_size=2 * e["params"]["input_prefix_bytes"])
            want = loader.read_range(e["params"]["input"], 0,
                                     e["params"]["input_prefix_bytes"])
            check("C6.zstd_roundtrip", raw == want,
                  f"decompresses to the exact {len(want)}-byte source prefix")

    # ctrl_base64: alphabet, and decodes to the seeded PRNG stream.
    e = _entry_live(manifest, "ctrl_base64")
    if e:
        text = (DATA_DIR / e["file"]).read_bytes()
        alphabet = set(b"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                       b"abcdefghijklmnopqrstuvwxyz0123456789+/=\n")
        check("C6.base64_alphabet", set(text) <= alphabet,
              "payload uses only the base64 alphabet + newlines")
        decoded = base64.b64decode(text.replace(b"\n", b""))
        expect = random.Random(e["params"]["seed"]).randbytes(
            e["params"]["raw_bytes"])
        check("C6.base64_seeded_source", decoded == expect,
              f"decodes to the seeded {len(expect)}-byte PRNG stream")

    # Negative-control separation: real text compresses far better than the
    # controls (sanity for downstream rate measurements).
    e8 = _entry_live(manifest, "enwik8")
    eu = _entry_live(manifest, "ctrl_urandom")
    if e8 and eu and not e8.get("standin"):
        r_text = zlib_ratio(loader.read_range("enwik8", 0, MiB))
        r_rand = zlib_ratio(loader.read_range("ctrl_urandom", 0, MiB))
        check("C6.control_separation", r_rand - r_text > 0.3,
              f"zlib ratio urandom {r_rand:.3f} vs enwik8 {r_text:.3f} "
              f"(gap > 0.3)")


def main() -> int:
    print("RNR corpus suite checks (exp-design section 3)")
    print("=" * 70)
    if not (DATA_DIR / "MANIFEST.json").exists():
        print("FAIL  C0.manifest: MANIFEST.json missing (run fetch.py first)")
        print("\nOVERALL: FAIL")
        return 1
    with open(DATA_DIR / "MANIFEST.json", encoding="utf-8") as f:
        manifest = json.load(f)

    c1_manifest_schema(manifest)
    c2_payload_integrity(manifest)
    c3_loader_roundtrip(manifest)
    c4_class_coverage(manifest)
    c5_deferred(manifest)
    c6_content_sanity(manifest)

    n_pass = sum(1 for _, ok, _ in RESULTS if ok)
    n_fail = len(RESULTS) - n_pass
    standins = [name for name, e in manifest["entries"].items()
                if e.get("standin")]
    print("-" * 70)
    print(f"{n_pass} passed, {n_fail} failed"
          + (f"; STAND-INS in use: {standins}" if standins else ""))
    print(f"OVERALL: {'PASS' if n_fail == 0 else 'FAIL'}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
