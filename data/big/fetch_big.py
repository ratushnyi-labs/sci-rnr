#!/usr/bin/env python3
"""fetch_big.py -- assemble the Type-I scale-benchmark corpus ladder (data/big/).

Size convention: decimal megabytes, 1 MB = 10^6 bytes.  Tier ladder is EXACTLY
500 / 1500 / 5000 MB = 500_000_000 / 1_500_000_000 / 5_000_000_000 bytes,
realised as PREFIX-NESTED masters: one master payload per content type; the
smaller tiers are byte-prefixes of the same master, so per-type results at
different scales are measured on nested data.  Exception: `text` is capped at
1000 MB (enwik9 is exactly 10^9 bytes; looping/repeating text would corrupt
compressibility statistics), so the text ladder is 500 / 1000 MB only and the
1500/5000 tiers are documented as capped.

Everything here is deterministic and idempotent:
  * downloads are resumable (curl -C -) and verified by size + sha256;
  * synthetic payloads (white_noise) are regenerable from a recorded seed;
  * derived payloads (video concat, JPEG stream, archive_mix tar,
    precompressed zstd) follow fixed, seeded, recorded recipes.

Usage:
  fetch_big.py --fetch [--jobs N] [--only NAME]   # download sources
  fetch_big.py --build TYPE                       # build one master
                                                  # (white_noise|text|video|
                                                  #  images|archive_mix|
                                                  #  precompressed|all)
  fetch_big.py --manifest                         # (re)write MANIFEST-BIG.json
  fetch_big.py --status                           # report what exists

Payloads live in downloads/ and masters/ (both git-ignored); only scripts and
MANIFEST-BIG.json are committed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BIG_DIR = Path(__file__).resolve().parent
DL_DIR = BIG_DIR / "downloads"
MASTER_DIR = BIG_DIR / "masters"
TMP_DIR = BIG_DIR / "tmp"
MANIFEST_PATH = BIG_DIR / "MANIFEST-BIG.json"
DATA_DIR = BIG_DIR.parent          # data/ (existing small-corpus payloads)

MB = 10**6                          # decimal megabyte convention
TIERS_MB = [500, 1500, 5000]
TEXT_TIERS_MB = [500, 1000]         # text is capped at enwik9 = 10^9 bytes
MASTER_BYTES = 5000 * MB
TEXT_BYTES = 1000 * MB
SEED = 20260710                     # global corpus seed (date of assembly)
CHUNK = 1 << 26                     # 64 MiB build/hash chunk

# --------------------------------------------------------------------------
# Download sources.  bytes = authoritative Content-Length observed 2026-07-10.
# sha256 is filled in after first successful fetch (cached in *.sha256 side
# files and recorded in MANIFEST-BIG.json); when present it is verified.
# --------------------------------------------------------------------------
SOURCES = [
    # -- text -------------------------------------------------------------
    dict(name="enwik9.zip",
         url="https://mattmahoney.net/dc/enwik9.zip",
         bytes=322592222,
         license="Wikipedia text (2006-03-03 dump); GFDL/CC-BY-SA",
         role="text"),
    # -- video: distinct open movies first (they fill the low tiers), then
    #    distinct encodes of already-used titles (tail of the 5000 MB tier).
    dict(name="big_buck_bunny_1080p_h264.mov",
         url="https://download.blender.org/peach/bigbuckbunny_movies/big_buck_bunny_1080p_h264.mov",
         bytes=725106140,
         license="(c) Blender Foundation 2008, CC-BY 3.0",
         role="video"),
    dict(name="Sintel.2010.1080p.mkv",
         url="https://download.blender.org/durian/movies/Sintel.2010.1080p.mkv",
         bytes=1180090590,
         license="(c) Blender Foundation 2010, CC-BY 3.0",
         role="video"),
    dict(name="tears_of_steel_1080p.mov",
         url="https://download.blender.org/demo/movies/ToS/tears_of_steel_1080p.mov",
         bytes=583774083,
         license="(c) Blender Foundation 2012, CC-BY 3.0",
         role="video"),
    dict(name="ED_1024.avi",
         url="https://download.blender.org/ED/ED_1024.avi",
         bytes=445866736,
         license="(c) Blender Foundation 2006, CC-BY 2.5 (Elephants Dream)",
         role="video"),
    dict(name="caminandes_gran_dillama.mp4.zip",
         url="https://download.blender.org/demo/movies/caminandes_gran_dillama.mp4.zip",
         bytes=125632282,
         license="(c) Blender Foundation 2013, CC-BY 3.0 (Caminandes: Gran Dillama)",
         role="video-zip"),
    dict(name="ToS-4k-1920.mov",
         url="https://download.blender.org/demo/movies/ToS/ToS-4k-1920.mov",
         bytes=738876331,
         license="(c) Blender Foundation 2012, CC-BY 3.0 (Tears of Steel, 4K master encode)",
         role="video"),
    dict(name="Sintel.2010.720p.mkv",
         url="https://download.blender.org/durian/movies/Sintel.2010.720p.mkv",
         bytes=681285280,   # matches server Content-Length and official .md5

         license="(c) Blender Foundation 2010, CC-BY 3.0 (Sintel, 720p encode)",
         role="video"),
    dict(name="bbb_sunflower_2160p_30fps_normal.mp4.zip",
         url="https://download.blender.org/demo/movies/BBB/bbb_sunflower_2160p_30fps_normal.mp4.zip",
         bytes=632204510,
         license="(c) Blender Foundation 2008/2013, CC-BY 3.0 (BBB sunflower 2160p re-render)",
         role="video-zip"),
]

# Video master concatenation order: distinct titles first, then distinct
# encodes of already-used titles.  Low tiers therefore contain the most
# distinct material; no file appears twice.  For .zip sources the extracted
# .mp4 member is concatenated, not the zip container.
VIDEO_CONCAT_ORDER = [
    "big_buck_bunny_1080p_h264.mov",
    "Sintel.2010.1080p.mkv",
    "tears_of_steel_1080p.mov",
    "ED_1024.avi",
    "caminandes_gran_dillama.mp4",          # extracted from .zip
    "ToS-4k-1920.mov",
    "Sintel.2010.720p.mkv",
    "bbb_sunflower_2160p_30fps_normal.mp4", # extracted from .zip
]

# JPEG-stream recipe: (source video, keep every Nth frame, byte cap).
# Frames are emitted as one concatenated MJPEG stream (q:v 2) per source, in
# this fixed order, each source contributing at most `cap` bytes, until the
# 5000 MB master target is reached.  Provenance: video-derived JPEG frames.
IMAGES_RECIPE = [
    ("ToS-4k-1920.mov", 2, 1200 * MB),
    ("bbb_sunflower_2160p_30fps_normal.mp4", 2, 1200 * MB),
    ("Sintel.2010.1080p.mkv", 3, 900 * MB),
    ("big_buck_bunny_1080p_h264.mov", 3, 900 * MB),
    ("tears_of_steel_1080p.mov", 3, 900 * MB),
    ("ED_1024.avi", 2, 900 * MB),
    ("Sintel.2010.720p.mkv", 2, 900 * MB),
]

MASTERS = {
    "white_noise": "masters/white_noise.bin",
    "text": "masters/text_enwik9.bin",
    "video": "masters/video_concat.bin",
    "images": "masters/images_jpegstream.bin",
    "archive_mix": "masters/archive_mix.tar",
    "precompressed": "masters/precompressed_zstd19.bin",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256_file(path: Path, limit: int | None = None) -> str:
    h = hashlib.sha256()
    remaining = limit if limit is not None else float("inf")
    with open(path, "rb") as f:
        while remaining > 0:
            block = f.read(int(min(CHUNK, remaining)))
            if not block:
                break
            h.update(block)
            remaining -= len(block)
    return h.hexdigest()


def prefix_hashes(path: Path, cut_points: list[int]) -> dict[int, str]:
    """sha256 of each prefix length in cut_points, in ONE pass over the file."""
    cuts = sorted(cut_points)
    out: dict[int, str] = {}
    h = hashlib.sha256()
    pos = 0
    ci = 0
    with open(path, "rb") as f:
        while ci < len(cuts):
            need = cuts[ci] - pos
            block = f.read(int(min(CHUNK, need)))
            if not block:
                break
            h.update(block)
            pos += len(block)
            while ci < len(cuts) and pos == cuts[ci]:
                out[cuts[ci]] = h.copy().hexdigest()
                ci += 1
    return out


# --------------------------------------------------------------------------
# fetch
# --------------------------------------------------------------------------

def fetch_one(src: dict) -> tuple[str, bool, str]:
    dest = DL_DIR / src["name"]
    side = dest.with_suffix(dest.suffix + ".sha256")
    if dest.exists() and dest.stat().st_size == src["bytes"]:
        if not side.exists():
            side.write_text(sha256_file(dest) + "\n")
        return src["name"], True, "already fetched"
    DL_DIR.mkdir(exist_ok=True)
    # resumable download via curl (falls back to fresh start on -C failure)
    cmd = ["curl", "-sSL", "--retry", "3", "--retry-delay", "2",
           "-C", "-", "-o", str(dest), src["url"]]
    t0 = time.time()
    rc = subprocess.run(cmd, capture_output=True, text=True)
    if rc.returncode != 0 and dest.exists():
        # some servers reject resume of a complete file; retry fresh if short
        if dest.stat().st_size != src["bytes"]:
            dest.unlink()
            rc = subprocess.run(cmd, capture_output=True, text=True)
    if rc.returncode != 0:
        return src["name"], False, f"curl failed: {rc.stderr.strip()[:200]}"
    size = dest.stat().st_size
    if size != src["bytes"]:
        return src["name"], False, f"size mismatch: got {size}, want {src['bytes']}"
    digest = sha256_file(dest)
    side.write_text(digest + "\n")
    dt = time.time() - t0
    return src["name"], True, f"{size/1e6:.1f} MB in {dt:.0f}s, sha256={digest[:16]}..."


def do_fetch(jobs: int, only: str | None) -> bool:
    todo = [s for s in SOURCES if only is None or s["name"] == only]
    ok = True
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(fetch_one, s): s for s in todo}
        for fut in as_completed(futs):
            name, good, msg = fut.result()
            log(f"fetch {name}: {'OK' if good else 'FAIL'} ({msg})")
            ok &= good
    return ok


def extract_zips() -> None:
    """Extract the single .mp4 member from each *-zip source (idempotent)."""
    for src in SOURCES:
        if src["role"] != "video-zip":
            continue
        zpath = DL_DIR / src["name"]
        if not zpath.exists():
            continue
        with zipfile.ZipFile(zpath) as zf:
            members = [m for m in zf.namelist() if m.lower().endswith(".mp4")]
            assert len(members) == 1, members
            out = DL_DIR / Path(members[0]).name
            if out.exists() and out.stat().st_size == zf.getinfo(members[0]).file_size:
                continue
            log(f"extracting {members[0]} from {src['name']}")
            with zf.open(members[0]) as fin, open(out, "wb") as fout:
                shutil.copyfileobj(fin, fout, CHUNK)


# --------------------------------------------------------------------------
# builders
# --------------------------------------------------------------------------

def build_white_noise() -> Path:
    """Seeded, regenerable noise: SHAKE-256 XOF with per-chunk domain
    separation.  Chunk i (64 MiB) = SHAKE256("rnr-big:white_noise:v1:
    seed=<SEED>:chunk=<i>").  Independent of any third-party library."""
    out = MASTER_DIR / "white_noise.bin"
    if out.exists() and out.stat().st_size == MASTER_BYTES:
        log("white_noise: already built")
        return out
    MASTER_DIR.mkdir(exist_ok=True)
    tmp = out.with_suffix(".part")
    t0 = time.time()
    with open(tmp, "wb") as f:
        written = 0
        i = 0
        while written < MASTER_BYTES:
            n = min(CHUNK, MASTER_BYTES - written)
            block = white_noise_chunk(i, n)
            f.write(block)
            written += n
            i += 1
    tmp.rename(out)
    log(f"white_noise: {MASTER_BYTES/1e6:.0f} MB in {time.time()-t0:.0f}s")
    return out


def white_noise_chunk(i: int, n: int = CHUNK) -> bytes:
    dom = f"rnr-big:white_noise:v1:seed={SEED}:chunk={i}".encode()
    return hashlib.shake_256(dom).digest(n)


def build_text() -> Path:
    out = MASTER_DIR / "text_enwik9.bin"
    if out.exists() and out.stat().st_size == TEXT_BYTES:
        log("text: already built")
        return out
    zpath = DL_DIR / "enwik9.zip"
    if not zpath.exists():
        raise FileNotFoundError("enwik9.zip not fetched yet")
    MASTER_DIR.mkdir(exist_ok=True)
    tmp = out.with_suffix(".part")
    with zipfile.ZipFile(zpath) as zf:
        with zf.open("enwik9") as fin, open(tmp, "wb") as fout:
            shutil.copyfileobj(fin, fout, CHUNK)
    assert tmp.stat().st_size == TEXT_BYTES, tmp.stat().st_size
    tmp.rename(out)
    log(f"text: enwik9 extracted, {TEXT_BYTES} bytes")
    return out


def build_video() -> Path:
    out = MASTER_DIR / "video_concat.bin"
    if out.exists() and out.stat().st_size == MASTER_BYTES:
        log("video: already built")
        return out
    extract_zips()
    parts = [DL_DIR / n for n in VIDEO_CONCAT_ORDER]
    missing = [p.name for p in parts if not p.exists()]
    if missing:
        raise FileNotFoundError(f"video parts missing: {missing}")
    total = sum(p.stat().st_size for p in parts)
    if total < MASTER_BYTES:
        raise RuntimeError(f"video parts sum {total} < {MASTER_BYTES}")
    MASTER_DIR.mkdir(exist_ok=True)
    tmp = out.with_suffix(".part")
    written = 0
    with open(tmp, "wb") as fout:
        for p in parts:
            if written >= MASTER_BYTES:
                break
            with open(p, "rb") as fin:
                while written < MASTER_BYTES:
                    block = fin.read(min(CHUNK, MASTER_BYTES - written))
                    if not block:
                        break
                    fout.write(block)
                    written += len(block)
            log(f"video: after {p.name}: {written/1e6:.0f} MB")
    assert written == MASTER_BYTES
    tmp.rename(out)
    log("video: master built (last file truncated at the 5000 MB boundary)")
    return out


def ffmpeg_binary() -> str:
    """ffmpeg executable bundled by the pip package imageio-ffmpeg
    (venv-local; version recorded in the manifest)."""
    import imageio_ffmpeg
    return imageio_ffmpeg.get_ffmpeg_exe()


def build_images(target: int = MASTER_BYTES) -> Path:
    """Concatenated JPEG (MJPEG q:v 2) frame stream, extracted from the
    downloaded open movies per IMAGES_RECIPE.  Provenance: video-derived
    JPEG frames (not an independent photographic corpus)."""
    out = MASTER_DIR / "images_jpegstream.bin"
    if out.exists():
        if out.stat().st_size >= target:
            log("images: already built")
            return out
        log(f"images: partial build ({out.stat().st_size/1e6:.0f} MB < "
            f"{target/1e6:.0f} MB); rebuilding deterministically")
        out.unlink()
    ff = ffmpeg_binary()
    MASTER_DIR.mkdir(exist_ok=True)
    tmp = out.with_suffix(".part")
    written = 0
    with open(tmp, "wb") as fout:
        for src_name, every_n, cap in IMAGES_RECIPE:
            if written >= target:
                break
            src = DL_DIR / src_name
            if not src.exists():
                log(f"images: source {src_name} missing, skipping")
                continue
            budget = min(cap, target - written)
            log(f"images: {src_name} every {every_n}th frame, budget {budget/1e6:.0f} MB")
            cmd = [ff, "-v", "error", "-i", str(src),
                   "-vf", f"select=not(mod(n\\,{every_n}))", "-fps_mode", "vfr",
                   "-q:v", "2", "-f", "image2pipe", "-vcodec", "mjpeg", "-"]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL)
            got = 0
            while got < budget:
                block = proc.stdout.read(min(CHUNK, budget - got))
                if not block:
                    break
                fout.write(block)
                got += len(block)
            proc.kill()
            proc.wait()
            written += got
            log(f"images: +{got/1e6:.0f} MB from {src_name} (total {written/1e6:.0f} MB)")
    if written > target:
        with open(tmp, "r+b") as f:
            f.truncate(target)
        written = target
    tmp.rename(out)
    log(f"images: master built, {written/1e6:.0f} MB"
        + ("" if written >= target else " (CAPPED below 5000 MB target)"))
    return out


def build_archive_mix() -> Path:
    """Deterministic tar: existing data/payloads/* plus seeded slices of the
    big masters.  Entry metadata is normalised (mtime=0, uid=gid=0) so the
    byte stream is reproducible given the masters.  Truncated at exactly
    5000 MB (tar tail cut at the tier boundary; benchmark treats payloads as
    raw byte streams)."""
    import random
    out = MASTER_DIR / "archive_mix.tar"
    if out.exists() and out.stat().st_size == MASTER_BYTES:
        log("archive_mix: already built")
        return out
    rng = random.Random(SEED)
    small = sorted((DATA_DIR / "payloads").glob("*"))
    donors = [MASTER_DIR / "white_noise.bin", MASTER_DIR / "text_enwik9.bin",
              MASTER_DIR / "video_concat.bin", MASTER_DIR / "images_jpegstream.bin"]
    donors = [d for d in donors if d.exists()]
    if not donors:
        raise RuntimeError("archive_mix: no big masters exist yet")
    MASTER_DIR.mkdir(exist_ok=True)
    tmp = out.with_suffix(".part")
    raw_target = MASTER_BYTES + 64 * MB     # overshoot, then truncate
    with open(tmp, "wb") as fraw:
        with tarfile.open(fileobj=fraw, mode="w", format=tarfile.GNU_FORMAT) as tf:
            def add_bytes(arcname: str, payload_path: Path,
                          offset: int = 0, length: int | None = None):
                size = length if length is not None else payload_path.stat().st_size
                info = tarfile.TarInfo(name=arcname)
                info.size = size
                info.mtime = 0
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                with open(payload_path, "rb") as f:
                    f.seek(offset)
                    tf.addfile(info, _LimitedReader(f, size))
            for p in small:
                add_bytes(f"mix/small/{p.name}", p)
            i = 0
            while fraw.tell() < raw_target:
                d = donors[i % len(donors)]
                dsize = d.stat().st_size
                length = rng.randrange(32 * MB, 96 * MB)
                offset = rng.randrange(0, max(1, dsize - length))
                add_bytes(f"mix/slices/{d.stem}_{i:03d}.bin", d, offset, length)
                i += 1
    with open(tmp, "r+b") as f:
        f.truncate(MASTER_BYTES)
    tmp.rename(out)
    log(f"archive_mix: master built ({i} slices + {len(small)} small payloads, "
        f"truncated to {MASTER_BYTES} bytes)")
    return out


class _LimitedReader:
    def __init__(self, f, remaining):
        self.f, self.remaining = f, remaining

    def read(self, n=-1):
        if self.remaining <= 0:
            return b""
        n = self.remaining if n < 0 else min(n, self.remaining)
        block = self.f.read(n)
        self.remaining -= len(block)
        return block


def zstd_version() -> str:
    rc = subprocess.run(["zstd", "--version"], capture_output=True, text=True)
    return rc.stdout.strip() or rc.stderr.strip()


def build_precompressed() -> Path:
    """zstd -19 of the text master followed by zstd -19 of the video master,
    concatenated, truncated at exactly 5000 MB: structured incompressible
    input for the benchmark."""
    out = MASTER_DIR / "precompressed_zstd19.bin"
    if out.exists() and out.stat().st_size == MASTER_BYTES:
        log("precompressed: already built")
        return out
    text = MASTER_DIR / "text_enwik9.bin"
    video = MASTER_DIR / "video_concat.bin"
    for p in (text, video):
        if not p.exists():
            raise FileNotFoundError(f"precompressed: needs {p.name}")
    MASTER_DIR.mkdir(exist_ok=True)
    parts = []
    for p in (text, video):
        z = TMP_DIR / (p.stem + ".zst")
        TMP_DIR.mkdir(exist_ok=True)
        if not z.exists():
            log(f"precompressed: zstd -19 -T0 {p.name} ...")
            t0 = time.time()
            rc = subprocess.run(["zstd", "-19", "-T0", "-f", "-o", str(z), str(p)],
                                capture_output=True, text=True)
            if rc.returncode != 0:
                raise RuntimeError(f"zstd failed: {rc.stderr[:300]}")
            log(f"precompressed: {z.name} {z.stat().st_size/1e6:.0f} MB "
                f"in {time.time()-t0:.0f}s")
        parts.append(z)
    total = sum(p.stat().st_size for p in parts)
    if total < MASTER_BYTES:
        log(f"precompressed: WARNING concat {total} < target {MASTER_BYTES}; "
            "master will be CAPPED at what exists")
    tmp = out.with_suffix(".part")
    written = 0
    with open(tmp, "wb") as fout:
        for p in parts:
            if written >= MASTER_BYTES:
                break
            with open(p, "rb") as fin:
                while written < MASTER_BYTES:
                    block = fin.read(min(CHUNK, MASTER_BYTES - written))
                    if not block:
                        break
                    fout.write(block)
                    written += len(block)
    tmp.rename(out)
    log(f"precompressed: master built, {written/1e6:.0f} MB")
    return out


BUILDERS = {
    "white_noise": build_white_noise,
    "text": build_text,
    "video": build_video,
    "images": build_images,
    "archive_mix": build_archive_mix,
    "precompressed": build_precompressed,
}


# --------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------

def source_record(src: dict) -> dict:
    dest = DL_DIR / src["name"]
    side = dest.with_suffix(dest.suffix + ".sha256")
    rec = dict(file=src["name"], url=src["url"], bytes=src["bytes"],
               license=src["license"])
    if side.exists():
        rec["sha256"] = side.read_text().strip()
    rec["fetched"] = dest.exists() and dest.stat().st_size == src["bytes"]
    if src["role"] == "video-zip":
        member = DL_DIR / src["name"].removesuffix(".zip")
        if member.exists():
            mside = member.with_suffix(member.suffix + ".sha256")
            if not mside.exists():
                mside.write_text(sha256_file(member) + "\n")
            rec["extracted_member"] = dict(
                file=member.name, bytes=member.stat().st_size,
                sha256=mside.read_text().strip())
    return rec


def type_entry(tname: str) -> dict:
    rel = MASTERS[tname]
    path = BIG_DIR / rel
    tiers_mb = TEXT_TIERS_MB if tname == "text" else TIERS_MB
    e: dict = {
        "type": tname,
        "file": rel,
        "tiers_mb": tiers_mb,
        "size_convention": "decimal MB = 10^6 bytes",
        "prefix_nested": True,
    }
    notes = {
        "white_noise": (
            "Seeded regenerable noise. Chunk i (64 MiB) = SHAKE-256 XOF of "
            f"'rnr-big:white_noise:v1:seed={SEED}:chunk=i'; concatenated and "
            "truncated to 5000 MB. Regenerate with fetch_big.py --build white_noise."),
        "text": (
            "enwik9 (first 10^9 bytes of English Wikipedia XML dump "
            "2006-03-03; Hutter-prize corpus). Ladder capped at 1000 MB: "
            "text is NOT looped/repeated because repetition corrupts "
            "compressibility statistics; 1500/5000 tiers do not exist for "
            "this type by design."),
        "video": (
            "Concatenation of distinct open-licensed movie files (Blender "
            "open movies), in the recorded order; no file repeated; distinct "
            "titles first, additional distinct encodes of already-used "
            "titles last; final file truncated at the exact 5000 MB "
            "boundary."),
        "images": (
            "Video-derived JPEG frames: single concatenated MJPEG (q:v 2) "
            "stream extracted from the downloaded open movies with the "
            "bundled ffmpeg of the pip package imageio-ffmpeg, per the "
            "fixed IMAGES_RECIPE (source order, every-Nth-frame, per-source "
            "byte caps). Honest provenance: derived from the video corpus, "
            "not an independent photographic corpus."),
        "archive_mix": (
            "Deterministic GNU tar (mtime=0, uid=gid=0): all existing "
            "data/payloads/* small corpora plus seeded slices "
            f"(random.Random({SEED}), 32-96 MB each) cycling over the "
            "white_noise/text/video/images masters; truncated at exactly "
            "5000 MB."),
        "precompressed": (
            "zstd -19 -T0 of text master followed by zstd -19 -T0 of video "
            "master, concatenated, truncated at exactly 5000 MB; structured "
            "incompressible input. zstd version recorded below."),
    }
    e["recipe"] = notes[tname]
    if tname == "video":
        e["concat_order"] = VIDEO_CONCAT_ORDER
    if tname == "images":
        e["extraction_recipe"] = [
            dict(source=s, every_nth_frame=n, byte_cap=c)
            for s, n, c in IMAGES_RECIPE]
        try:
            import imageio_ffmpeg
            e["ffmpeg"] = dict(
                package="imageio-ffmpeg " + imageio_ffmpeg.__version__,
                ffmpeg_version=imageio_ffmpeg.get_ffmpeg_version())
        except Exception:
            pass
    if tname == "precompressed":
        try:
            e["zstd_version"] = zstd_version()
        except Exception:
            pass
    if path.exists():
        size = path.stat().st_size
        e["bytes"] = size
        cuts = [mb * MB for mb in tiers_mb if mb * MB <= size]
        if size not in cuts:
            cuts.append(size)
        ph = prefix_hashes(path, cuts)
        e["sha256"] = ph.get(size)
        e["tiers"] = [
            {"mb": mb, "bytes": mb * MB, "sha256": ph.get(mb * MB),
             "materialized": mb * MB <= size}
            for mb in tiers_mb]
        e["status"] = ("complete" if size >= tiers_mb[-1] * MB
                       else f"partial ({size/1e6:.0f} MB of {tiers_mb[-1]} MB)")
    else:
        e["bytes"] = None
        e["tiers"] = [{"mb": mb, "bytes": mb * MB, "sha256": None,
                       "materialized": False} for mb in tiers_mb]
        e["status"] = "deferred (master not built)"
    return e


def write_manifest() -> None:
    man = {
        "suite": "rnr-big",
        "version": 1,
        "size_convention": "decimal MB = 10^6 bytes; tiers 500/1500/5000 MB "
                           "(text: 500/1000 MB, capped at enwik9)",
        "prefix_nesting": "each tier is the byte-prefix of the per-type "
                          "master; tier sha256 = sha256 of that prefix",
        "seed": SEED,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "sources": [source_record(s) for s in SOURCES],
        "entries": {t: type_entry(t) for t in MASTERS},
    }
    MANIFEST_PATH.write_text(json.dumps(man, indent=2) + "\n")
    log(f"manifest written: {MANIFEST_PATH}")


def status() -> None:
    for s in SOURCES:
        p = DL_DIR / s["name"]
        got = p.stat().st_size if p.exists() else 0
        print(f"  src {s['name']:<44} {got:>12}/{s['bytes']} "
              f"{'OK' if got == s['bytes'] else 'partial' if got else 'missing'}")
    for t, rel in MASTERS.items():
        p = BIG_DIR / rel
        print(f"  master {t:<18} "
              f"{p.stat().st_size if p.exists() else 0:>12} bytes")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--only")
    ap.add_argument("--build", choices=list(BUILDERS) + ["all"])
    ap.add_argument("--manifest", action="store_true")
    ap.add_argument("--status", action="store_true")
    args = ap.parse_args()
    ok = True
    if args.fetch:
        ok &= do_fetch(args.jobs, args.only)
    if args.build:
        names = list(BUILDERS) if args.build == "all" else [args.build]
        for n in names:
            try:
                BUILDERS[n]()
            except Exception as exc:
                log(f"build {n}: FAILED ({exc})")
                ok = False
    if args.manifest:
        write_manifest()
    if args.status:
        status()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
