#!/usr/bin/env python3
"""Corpus suite fetcher/builder for the RNR empirical program.

Materializes the corpus suite described in tex/rnr_experimental_design.tex
Section 3 (classes 3.1 text/source, 3.2 structured binary, 3.3 image/numeric,
3.4 heterogeneous + negative controls) at "small classics" scale:

  downloads   : Calgary corpus, Canterbury corpus, enwik8
  derived     : pinned tar of this repository's scripts/ directory
                (source-code sample), zstd-compressed blob (negative control)
  generated   : synthetic SQLite database, spatially correlated float32
                telemetry grid, os.urandom block, base64 text
  deferred    : enwik9 (1 GB; manifest entry only, fetched with
                --include-deferred)

Behaviour:
  * Idempotent: an entry whose payload already exists with the sha256
    recorded in MANIFEST.json is skipped.
  * Checksums: sha256 and byte size are computed after every
    download/generation and written to MANIFEST.json.
  * Network failure: a deterministic synthetic stand-in is produced and
    the manifest entry is marked "standin": true.

Usage:
  python fetch.py                 # fetch/build all non-deferred entries
  python fetch.py --only enwik8   # single entry
  python fetch.py --include-deferred   # also fetch enwik9 (large!)
  python fetch.py --list          # print registry and exit
"""

from __future__ import annotations

import argparse
import base64
import datetime
import gzip
import hashlib
import json
import os
import random
import shutil
import sqlite3
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent
PAYLOAD_DIR = DATA_DIR / "payloads"
MANIFEST_PATH = DATA_DIR / "MANIFEST.json"
REPO_ROOT = DATA_DIR.parent

SEED = 20260710          # global suite seed; per-entry seeds derive from it
CHUNK = 1 << 20          # 1 MiB streaming chunk
HTTP_TIMEOUT = 60        # seconds
MiB = 1 << 20

# --------------------------------------------------------------------------
# Registry: one record per corpus entry.  MANIFEST.json is the materialized
# form of this registry plus computed byte sizes and sha256 digests.
# --------------------------------------------------------------------------

REGISTRY = [
    {
        "name": "enwik8",
        "class": "3.1-text-source",
        "file": "payloads/enwik8",
        "source_url": "https://mattmahoney.net/dc/enwik8.zip",
        "license_note": ("English Wikipedia XML dump excerpt (2006-03-03); "
                         "text under GFDL / CC BY-SA per Wikipedia terms; "
                         "standard Hutter Prize benchmark file."),
        "qualifies_for": ["H1", "H3", "H4", "H6", "H8"],
        "provenance": "download",
        "deferred": False,
        "notes": ("First 10^8 bytes of English Wikipedia. Small-scale "
                  "stand-in for the enwik9 entry of exp-design section 3.1."),
    },
    {
        "name": "enwik9",
        "class": "3.1-text-source",
        "file": "payloads/enwik9",
        "source_url": "https://mattmahoney.net/dc/enwik9.zip",
        "license_note": ("English Wikipedia XML dump excerpt (2006-03-03); "
                         "text under GFDL / CC BY-SA per Wikipedia terms; "
                         "standard Hutter Prize benchmark file."),
        "qualifies_for": ["H1", "H3", "H4", "H6", "H15"],
        "provenance": "download",
        "deferred": "large",
        "expected_bytes": 1000000000,
        "notes": ("Deferred: ~0.3 GB download, 1 GB payload. Required for "
                  "H15 (Chinchilla-on-enwik9 synthesis). Fetch with "
                  "--include-deferred; sha256 is recorded on first fetch."),
    },
    {
        "name": "calgary",
        "class": "3.4-heterogeneous-negative",
        "file": "payloads/calgary.tar",
        "source_url": "https://corpus.canterbury.ac.nz/resources/calgary.tar.gz",
        "license_note": ("Calgary corpus (Bell, Cleary, Witten 1987-1989); "
                         "distributed freely for compression research by the "
                         "University of Calgary / University of Canterbury."),
        "qualifies_for": ["H1", "H5", "H8"],
        "provenance": "download",
        "deferred": False,
        "notes": ("Classic 14-file mixed corpus (English text, program "
                  "sources, geophysical binary, bitmap picture, object "
                  "files). Stored as the decompressed .tar so one payload "
                  "file carries the heterogeneous mix; exercises "
                  "Type-III-A scheduling (class 3.4)."),
    },
    {
        "name": "canterbury",
        "class": "3.1-text-source",
        "file": "payloads/canterbury.tar",
        "source_url": "https://corpus.canterbury.ac.nz/resources/cantrbry.tar.gz",
        "license_note": ("Canterbury corpus (Arnold, Bell 1997); distributed "
                         "freely for compression research by the University "
                         "of Canterbury."),
        "qualifies_for": ["H1", "H3", "H4", "H8"],
        "provenance": "download",
        "deferred": False,
        "notes": ("Classic 11-file corpus, predominantly text and source "
                  "code. Stored as the decompressed .tar."),
    },
    {
        "name": "scripts_src",
        "class": "3.1-text-source",
        "file": "payloads/rnr_scripts_src.tar",
        "source_url": "local:scripts/",
        "license_note": ("This repository's own verification scripts; same "
                         "license as the repository."),
        "qualifies_for": ["H1", "H3", "H4"],
        "provenance": "derived",
        "deferred": False,
        "notes": ("Deterministic tar (sorted members, zeroed timestamps and "
                  "ownership, ustar format) of the repo's scripts/ "
                  "directory: a pinned Python source-code sample standing "
                  "in for the Linux-kernel-source entry of section 3.1. "
                  "Pinned at generation time; the manifest sha256 is the "
                  "reference."),
    },
    {
        "name": "sqlite_synth",
        "class": "3.2-structured-binary",
        "file": "payloads/sqlite_synth.db",
        "source_url": "generated:sqlite",
        "license_note": "Synthetic data generated by fetch.py; public domain.",
        "qualifies_for": ["H2", "H3", "H4", "H8"],
        "provenance": "generated",
        "deferred": False,
        "params": {"seed": SEED, "sensors": 500, "readings": 250000,
                   "logs": 20000, "page_size": 4096},
        "notes": ("SQLite database built from seeded synthetic rows "
                  "(sensor registry, time-series readings with an index, "
                  "log lines). B-tree pages give the positional regularity "
                  "targeted by Type-II (class 3.2, database-pages "
                  "stand-in). Byte-exact layout can vary across SQLite "
                  "library versions; the manifest sha256 pins the "
                  "generated instance."),
    },
    {
        "name": "telemetry_grid",
        "class": "3.3-image-numeric",
        "file": "payloads/telemetry_grid.f32",
        "source_url": "generated:telemetry",
        "license_note": "Synthetic data generated by fetch.py; public domain.",
        "qualifies_for": ["H4", "H6"],
        "provenance": "generated",
        "deferred": False,
        "params": {"seed": SEED, "shape": [2048, 2048], "dtype": "<f4",
                   "order": "C", "corr_length_px": 24},
        "notes": ("Raw headerless 2048x2048 float32 (little-endian, "
                  "C-order) sensor field: Gaussian-filtered white noise "
                  "(correlation length 24 px) plus a smooth deterministic "
                  "trend. Stand-in for the HDF5/scientific-array entry of "
                  "section 3.3; the float32 numeric class of H6."),
    },
    {
        "name": "ctrl_urandom",
        "class": "3.4-heterogeneous-negative",
        "file": "payloads/ctrl_urandom.bin",
        "source_url": "generated:os.urandom",
        "license_note": "Locally generated randomness; public domain.",
        "qualifies_for": ["H16"],
        "provenance": "generated",
        "deferred": False,
        "params": {"bytes": 16 * MiB},
        "notes": ("Negative control: incompressible. 16 MiB os.urandom "
                  "block (reduced-scale stand-in for the 64 MB /dev/urandom "
                  "entry of section 3.4). Inherently non-reproducible; the "
                  "manifest sha256 pins the drawn instance."),
    },
    {
        "name": "ctrl_zstd",
        "class": "3.4-heterogeneous-negative",
        "file": "payloads/ctrl_zstd.zst",
        "source_url": "derived:enwik8[:32MiB] | zstd -19",
        "license_note": ("Derived from enwik8 (see enwik8 license note); "
                         "compressed form only."),
        "qualifies_for": ["H16"],
        "provenance": "derived",
        "deferred": False,
        "params": {"input": "enwik8", "input_prefix_bytes": 32 * MiB,
                   "zstd_level": 19},
        "notes": ("Negative control: high-entropy residual. First 32 MiB "
                  "of the enwik8 payload compressed with zstandard level "
                  "19 (reduced-scale stand-in for the 128 MB "
                  "zstd-precompressed-text entry of section 3.4)."),
    },
    {
        "name": "ctrl_base64",
        "class": "3.4-heterogeneous-negative",
        "file": "payloads/ctrl_base64.txt",
        "source_url": "generated:base64(seeded PRNG bytes)",
        "license_note": "Synthetic data generated by fetch.py; public domain.",
        "qualifies_for": ["H1"],
        "provenance": "generated",
        "deferred": False,
        "params": {"seed": SEED ^ 0xB64B64, "raw_bytes": 12 * MiB,
                   "line_length": 76},
        "notes": ("Negative/calibration control: base64 text over seeded "
                  "PRNG bytes. Looks like ASCII text but carries exactly 6 "
                  "bits per payload byte, so the ideal coding rate is "
                  "known in closed form (~6.06 bpb including newlines); "
                  "calibrates text-rate measurements (H1)."),
    },
]


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(CHUNK)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def utcnow() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")


def load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    return {
        "manifest_version": 1,
        "suite": "RNR corpus suite (small classics scale)",
        "spec_reference": "tex/rnr_experimental_design.tex section 3",
        "seed": SEED,
        "entries": {},
    }


def save_manifest(manifest: dict) -> None:
    manifest["updated_utc"] = utcnow()
    tmp = MANIFEST_PATH.with_suffix(".json.tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=False)
        f.write("\n")
    tmp.replace(MANIFEST_PATH)


def download(url: str, dest: Path) -> None:
    """Stream url to dest via a temp file; raise on any network error."""
    req = urllib.request.Request(url, headers={"User-Agent": "rnr-corpus-fetch/1"})
    with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as tmp:
        tmp_path = Path(tmp.name)
        try:
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp:
                shutil.copyfileobj(resp, tmp, length=CHUNK)
        except BaseException:
            tmp_path.unlink(missing_ok=True)
            raise
    tmp_path.replace(dest)


def gunzip_file(src: Path, dest: Path) -> None:
    with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=CHUNK)


def extract_zip_member(src: Path, member: str, dest: Path) -> None:
    with zipfile.ZipFile(src) as z, z.open(member) as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=CHUNK)


# --------------------------------------------------------------------------
# Builders (one per provenance kind)
# --------------------------------------------------------------------------

def fetch_download(entry: dict, dest: Path) -> None:
    """Download-and-unpack for calgary / canterbury / enwik8 / enwik9."""
    url = entry["source_url"]
    staging = dest.parent / (dest.name + ".download")
    download(url, staging)
    try:
        if url.endswith(".tar.gz"):
            gunzip_file(staging, dest)
        elif url.endswith(".zip"):
            member = dest.name  # enwik8.zip contains 'enwik8', etc.
            extract_zip_member(staging, member, dest)
        else:
            staging.replace(dest)
            return
    finally:
        staging.unlink(missing_ok=True)


_STANDIN_WORDS = (
    "the of and to in a is that for it as was with be by on not he this are "
    "or his from at which but have an they you were her all she there would "
    "their we him been has when who will no more if out so said what up its "
    "about into than them can only other new some could time these two may "
    "then do first any my now such like our over man me even most made after "
    "also did many before must through years where much your way well down "
    "should because each just those people how too little state good very "
    "make world still own see men work long get here between both life being "
    "under never day same another know while last might us great old year "
    "off come since against go came right used take three").split()


def synthesize_standin(entry: dict, dest: Path, size: int = 8 * MiB) -> None:
    """Deterministic English-like filler used only when a download fails."""
    name_tag = int.from_bytes(
        hashlib.sha256(entry["name"].encode()).digest()[:4], "big")
    rnd = random.Random(SEED ^ name_tag)
    with open(dest, "wb") as f:
        written = 0
        while written < size:
            n = rnd.randint(8, 15)
            line = " ".join(rnd.choice(_STANDIN_WORDS) for _ in range(n))
            data = (line + "\n").encode("ascii")
            f.write(data)
            written += len(data)


def build_scripts_tar(entry: dict, dest: Path) -> None:
    """Deterministic tar of the repository's scripts/ directory."""
    src_root = REPO_ROOT / "scripts"
    if not src_root.is_dir():
        raise RuntimeError(f"missing source directory {src_root}")
    files = sorted(
        p for p in src_root.rglob("*")
        if p.is_file()
        and "__pycache__" not in p.parts
        and p.name != ".DS_Store"
        and p.suffix != ".pyc"
    )
    with tarfile.open(dest, "w", format=tarfile.USTAR_FORMAT) as tar:
        for p in files:
            info = tarfile.TarInfo(
                name="rnr_scripts/" + p.relative_to(src_root).as_posix())
            info.size = p.stat().st_size
            info.mtime = 0
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            info.mode = 0o644
            with open(p, "rb") as f:
                tar.addfile(info, f)


def build_sqlite(entry: dict, dest: Path) -> None:
    """Synthetic relational database with seeded rows (class 3.2)."""
    params = entry["params"]
    rnd = random.Random(params["seed"] ^ 0x5EED)
    dest.unlink(missing_ok=True)
    con = sqlite3.connect(dest)
    try:
        con.execute(f"PRAGMA page_size={params['page_size']}")
        con.execute("PRAGMA journal_mode=MEMORY")
        con.executescript("""
            CREATE TABLE sensors(
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                lat REAL, lon REAL,
                install_ts INTEGER
            );
            CREATE TABLE readings(
                id INTEGER PRIMARY KEY,
                sensor_id INTEGER NOT NULL REFERENCES sensors(id),
                ts INTEGER NOT NULL,
                value REAL NOT NULL,
                status TEXT NOT NULL
            );
            CREATE TABLE logs(
                id INTEGER PRIMARY KEY,
                ts INTEGER NOT NULL,
                level TEXT NOT NULL,
                message TEXT NOT NULL
            );
        """)
        statuses = ["ok", "ok", "ok", "ok", "warn", "drift", "recal"]
        levels = ["INFO", "INFO", "INFO", "WARN", "ERROR"]
        t0 = 1700000000
        sensors = []
        for i in range(1, params["sensors"] + 1):
            sensors.append((
                i,
                f"sensor-{i:04d}-{rnd.choice('ABCDEF')}",
                round(rnd.uniform(-90, 90), 6),
                round(rnd.uniform(-180, 180), 6),
                t0 - rnd.randint(0, 10**7),
            ))
        con.executemany("INSERT INTO sensors VALUES (?,?,?,?,?)", sensors)

        # Time-series with per-sensor random walk: realistic column
        # correlations for the Type-II positional model.
        walk = {i: rnd.uniform(10.0, 30.0) for i in range(1, params["sensors"] + 1)}
        rows = []
        for j in range(1, params["readings"] + 1):
            sid = rnd.randint(1, params["sensors"])
            walk[sid] += rnd.gauss(0.0, 0.05)
            rows.append((j, sid, t0 + j * 5 + rnd.randint(0, 4),
                         round(walk[sid], 4), rnd.choice(statuses)))
            if len(rows) >= 10000:
                con.executemany("INSERT INTO readings VALUES (?,?,?,?,?)", rows)
                rows = []
        if rows:
            con.executemany("INSERT INTO readings VALUES (?,?,?,?,?)", rows)

        logrows = []
        for j in range(1, params["logs"] + 1):
            nwords = rnd.randint(4, 12)
            msg = " ".join(rnd.choice(_STANDIN_WORDS) for _ in range(nwords))
            logrows.append((j, t0 + j * 60, rnd.choice(levels), msg))
        con.executemany("INSERT INTO logs VALUES (?,?,?,?)", logrows)

        con.execute("CREATE INDEX idx_readings_sensor_ts "
                    "ON readings(sensor_id, ts)")
        con.commit()
        con.execute("VACUUM")
    finally:
        con.close()


def build_telemetry(entry: dict, dest: Path) -> None:
    """Spatially correlated float32 sensor field (class 3.3)."""
    import numpy as np
    params = entry["params"]
    n0, n1 = params["shape"]
    ell = params["corr_length_px"]
    rng = np.random.default_rng(params["seed"])
    white = rng.standard_normal((n0, n1))
    # Gaussian low-pass in the Fourier domain: correlation length ell px.
    fx = np.fft.fftfreq(n0)[:, None]
    fy = np.fft.rfftfreq(n1)[None, :]
    transfer = np.exp(-2.0 * (np.pi ** 2) * (ell ** 2) * (fx ** 2 + fy ** 2))
    field = np.fft.irfft2(np.fft.rfft2(white) * transfer, s=(n0, n1))
    field /= field.std()
    # Smooth deterministic trend on top of the stochastic texture.
    yy, xx = np.mgrid[0:n0, 0:n1]
    field += 0.5 * np.sin(2 * np.pi * 3 * xx / n1) * np.cos(2 * np.pi * 2 * yy / n0)
    with open(dest, "wb") as f:
        f.write(np.ascontiguousarray(field, dtype="<f4").tobytes(order="C"))


def build_urandom(entry: dict, dest: Path) -> None:
    total = entry["params"]["bytes"]
    with open(dest, "wb") as f:
        written = 0
        while written < total:
            n = min(CHUNK, total - written)
            f.write(os.urandom(n))
            written += n


def build_zstd_blob(entry: dict, dest: Path, manifest: dict) -> None:
    import zstandard
    params = entry["params"]
    src_entry = manifest["entries"].get(params["input"])
    if not src_entry:
        raise RuntimeError(f"input corpus '{params['input']}' not fetched yet")
    src = DATA_DIR / src_entry["file"]
    with open(src, "rb") as f:
        raw = f.read(params["input_prefix_bytes"])
    cctx = zstandard.ZstdCompressor(level=params["zstd_level"])
    with open(dest, "wb") as f:
        f.write(cctx.compress(raw))


def build_base64(entry: dict, dest: Path) -> None:
    params = entry["params"]
    rnd = random.Random(params["seed"])
    raw = rnd.randbytes(params["raw_bytes"])
    encoded = base64.b64encode(raw)
    line = params["line_length"]
    with open(dest, "wb") as f:
        for i in range(0, len(encoded), line):
            f.write(encoded[i:i + line])
            f.write(b"\n")


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def materialize(entry: dict, dest: Path, manifest: dict) -> dict:
    """Build/download one payload; return manifest-entry status fields."""
    status = {"standin": False}
    name = entry["name"]
    if entry["provenance"] == "download":
        try:
            fetch_download(entry, dest)
        except (urllib.error.URLError, OSError, zipfile.BadZipFile) as exc:
            print(f"  WARNING: download failed for {name} ({exc}); "
                  f"synthesizing deterministic stand-in")
            synthesize_standin(entry, dest)
            status["standin"] = True
            status["standin_reason"] = str(exc)
    elif name == "scripts_src":
        build_scripts_tar(entry, dest)
    elif name == "sqlite_synth":
        build_sqlite(entry, dest)
        status["sqlite_library_version"] = sqlite3.sqlite_version
    elif name == "telemetry_grid":
        build_telemetry(entry, dest)
    elif name == "ctrl_urandom":
        build_urandom(entry, dest)
    elif name == "ctrl_zstd":
        build_zstd_blob(entry, dest, manifest)
    elif name == "ctrl_base64":
        build_base64(entry, dest)
    else:
        raise RuntimeError(f"no builder for entry {name}")
    return status


def process_entry(entry: dict, manifest: dict, include_deferred: bool,
                  force: bool) -> None:
    name = entry["name"]
    dest = DATA_DIR / entry["file"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    recorded = manifest["entries"].get(name, {})

    if entry.get("deferred") and not include_deferred:
        print(f"[{name}] deferred ({entry['deferred']}); manifest entry only "
              f"(use --include-deferred to fetch)")
        merged = {**entry, **{k: recorded[k] for k in
                              ("bytes", "sha256", "fetched_utc", "standin")
                              if k in recorded}}
        manifest["entries"][name] = merged
        return

    # Idempotence: skip when the payload matches the recorded checksum.
    if dest.exists() and not force:
        if recorded.get("sha256") and recorded.get("bytes") == dest.stat().st_size:
            if sha256_file(dest) == recorded["sha256"]:
                print(f"[{name}] up to date ({recorded['bytes']} bytes); skipping")
                manifest["entries"][name] = {**entry, **{
                    k: recorded[k] for k in
                    ("bytes", "sha256", "fetched_utc", "standin",
                     "standin_reason", "sqlite_library_version")
                    if k in recorded}}
                return
            print(f"[{name}] checksum mismatch with manifest; rebuilding")
        elif not recorded.get("sha256"):
            print(f"[{name}] present but unrecorded; adopting existing payload")
            record = dict(entry)
            record["bytes"] = dest.stat().st_size
            record["sha256"] = sha256_file(dest)
            record["fetched_utc"] = utcnow()
            record["standin"] = recorded.get("standin", False)
            manifest["entries"][name] = record
            return

    print(f"[{name}] materializing ({entry['provenance']}) ...")
    status = materialize(entry, dest, manifest)
    record = dict(entry)
    record.update(status)
    record["bytes"] = dest.stat().st_size
    record["sha256"] = sha256_file(dest)
    record["fetched_utc"] = utcnow()
    print(f"[{name}] done: {record['bytes']} bytes, "
          f"sha256={record['sha256'][:16]}...")
    manifest["entries"][name] = record


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", metavar="NAME", action="append",
                    help="process only the named entry (repeatable)")
    ap.add_argument("--include-deferred", action="store_true",
                    help="also fetch entries marked deferred (enwik9, large)")
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if payload matches the manifest")
    ap.add_argument("--list", action="store_true",
                    help="print the registry and exit")
    args = ap.parse_args()

    if args.list:
        for e in REGISTRY:
            flag = f" [deferred={e['deferred']}]" if e.get("deferred") else ""
            print(f"{e['name']:<16} {e['class']:<28} {e['provenance']}{flag}")
        return 0

    manifest = load_manifest()
    selected = [e for e in REGISTRY
                if not args.only or e["name"] in args.only]
    if args.only:
        missing = set(args.only) - {e["name"] for e in selected}
        if missing:
            print(f"ERROR: unknown entries: {sorted(missing)}")
            return 2

    # ctrl_zstd derives from enwik8; registry order already satisfies this.
    failures = 0
    for entry in selected:
        try:
            process_entry(entry, manifest, args.include_deferred, args.force)
        except Exception as exc:  # keep going; report at the end
            failures += 1
            print(f"[{entry['name']}] ERROR: {exc}")
        save_manifest(manifest)   # persist progress after each entry

    print(f"\nManifest written to {MANIFEST_PATH}")
    if failures:
        print(f"{failures} entr{'y' if failures == 1 else 'ies'} FAILED")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
