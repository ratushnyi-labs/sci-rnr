#!/usr/bin/env python3
"""Deterministic corpus materialization for the Type-I scale benchmark.

One corpus file per (tier, content-type) cell, built exactly to the
pre-registered byte count (protocol.md section 2) and pinned by sha256 in
a per-tier manifest (corpora/<tier>/manifest.json).  Materialization is
idempotent: an existing file whose size and sha256 match the manifest is
never rebuilt.

Recipes (protocol.md section 3; summary):

  white_noise   numpy Philox counter-based PRG keyed by the METHODS.md
                seed schedule -- reproducible at any size.
  text          byte slice/cycle of enwik8 (smoke) or enwik9 (tiers).
                Tiers REQUIRE enwik9 (1 GB period keeps every coder
                window below the repetition period; see protocol 3.2).
  video         synthetic raw-video generator v1: 640x360 8-bit luma
                frames, AR(1) temporal correlation over spatially
                box-smoothed innovations (seeded, self-contained).
  images        smoke: the pinned telemetry_grid payload; tiers:
                sequence of seeded correlated 2048x2048 float32 rasters
                (same class as data/ telemetry generator).
  archive_mix   cycled concatenation calgary.tar | canterbury.tar |
                rnr_scripts_src.tar | sqlite_synth.db, truncated.
  precompressed zstd -3 -T1 stream (the zstd CLI default level -- the
                most common precompressed payload in practice) over the
                text recipe's byte stream, truncated to the cell size.

All generators stream in bounded chunks: memory stays O(chunk), disk
holds only the corpus file being written.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
DATA = REPO / "data" / "payloads"
CORPORA_DIR = HERE / "corpora"

# Real tiers (500/1500/5000 MB) are byte-prefixes of the pinned data/big
# masters (real video/images/text corpora; data/big/MANIFEST-BIG.json).
# The synthetic generators below serve ONLY the smoke tiers.
BIG_DIR = REPO / "data" / "big"
BIG_MANIFEST = BIG_DIR / "MANIFEST-BIG.json"
REAL_TIERS = ("500", "1500", "5000")

sys.path.insert(0, str(REPO / "bench"))
import metrics as _metrics  # noqa: E402  (derive_seed; frozen seed schedule)

MB = 1_000_000
MIB = 1_048_576

# Pre-registered cell sizes in bytes (protocol.md section 2).
TIER_BYTES = {
    "smoke": 16 * MIB,          # machinery-validation tier
    "smoke48": 48 * MIB,        # second smoke size: exercises the
                                # multi-size figure paths before real tiers
    "500": 500 * MB,
    "1500": 1500 * MB,
    "5000": 5000 * MB,
}

CONTENT_TYPES = (
    "white_noise", "text", "video", "images", "archive_mix", "precompressed",
)

_CHUNK = 8 * MIB

ZSTD_BIN = "/opt/homebrew/bin/zstd"


def _seed(tier: str, ctype: str) -> int:
    return _metrics.derive_seed(f"type1_scale:{tier}:{ctype}", 0, tag="corpus")


# ---------------------------------------------------------------------------
# Generators: each writes exactly `nbytes` to the open file object.
# ---------------------------------------------------------------------------

def _gen_white_noise(out, nbytes: int, tier: str) -> None:
    rng = np.random.Generator(np.random.Philox(key=_seed(tier, "white_noise")))
    left = nbytes
    while left > 0:
        n = min(_CHUNK, left)
        out.write(rng.bytes(n))
        left -= n


def _payload(name: str) -> Path:
    p = DATA / name
    if not p.exists():
        raise FileNotFoundError(
            f"required data/ payload missing: {p} (run data/fetch.py"
            + (" --include-deferred" if name == "enwik9" else "")
            + ")")
    return p


def _cycle_files(out, paths: list, nbytes: int) -> None:
    """Concatenate the files cyclically until exactly nbytes are written."""
    left = nbytes
    while left > 0:
        for p in paths:
            with open(p, "rb") as f:
                while left > 0:
                    chunk = f.read(min(_CHUNK, left))
                    if not chunk:
                        break
                    out.write(chunk)
                    left -= len(chunk)
            if left == 0:
                return


def _text_sources(tier: str) -> list:
    if tier in ("smoke", "smoke48"):
        return [_payload("enwik8")]
    # Tiers require enwik9: 1 GB repetition period exceeds every coder
    # window in the matrix (max 128 MiB), so cycling stays artifact-free.
    return [_payload("enwik9")]


def _gen_text(out, nbytes: int, tier: str) -> None:
    _cycle_files(out, _text_sources(tier), nbytes)


def _gen_archive_mix(out, nbytes: int, tier: str) -> None:
    paths = [_payload("calgary.tar"), _payload("canterbury.tar"),
             _payload("rnr_scripts_src.tar"), _payload("sqlite_synth.db")]
    _cycle_files(out, paths, nbytes)


def _box_smooth(a: np.ndarray, k: int) -> np.ndarray:
    """Separable 2D box filter (window k) via cumulative sums.

    Edge windows are truncated and count-normalized; float32 throughout.
    """
    for axis in (0, 1):
        n = a.shape[axis]
        c = np.cumsum(a, axis=axis, dtype=np.float32)
        zshape = list(a.shape)
        zshape[axis] = 1
        c = np.concatenate([np.zeros(zshape, np.float32), c], axis=axis)
        pad = k // 2
        hi_idx = np.minimum(np.arange(n) + pad + 1, n)
        lo_idx = np.maximum(np.arange(n) - pad, 0)
        hi = np.take(c, hi_idx, axis=axis)
        lo = np.take(c, lo_idx, axis=axis)
        cnt_shape = [1, 1]
        cnt_shape[axis] = n
        cnt = (hi_idx - lo_idx).astype(np.float32).reshape(cnt_shape)
        a = (hi - lo) / cnt
    return a


def _gen_video(out, nbytes: int, tier: str) -> None:
    """Synthetic raw-video model v1 (protocol 3.4): 640x360 luma frames,
    frame_t = rho*frame_{t-1} + (1-rho)*smooth(noise_t), quantized uint8."""
    H, W, rho, k = 360, 640, np.float32(0.94), 15
    rng = np.random.Generator(np.random.Philox(key=_seed(tier, "video")))
    state = _box_smooth(
        rng.standard_normal((H, W), dtype=np.float32), k) * np.float32(k)
    left = nbytes
    while left > 0:
        innov = _box_smooth(
            rng.standard_normal((H, W), dtype=np.float32), k) * np.float32(k)
        state = rho * state + (np.float32(1) - rho) * innov
        frame = np.clip(state * np.float32(48) + np.float32(128),
                        0, 255).astype(np.uint8)
        raw = frame.tobytes()
        out.write(raw[: min(len(raw), left)])
        left -= min(len(raw), left)


def _gen_images(out, nbytes: int, tier: str) -> None:
    if tier == "smoke":
        # Existing pinned payload; exactly 16 MiB.
        with open(_payload("telemetry_grid.f32"), "rb") as f:
            left = nbytes
            while left > 0:
                chunk = f.read(min(_CHUNK, left))
                if not chunk:
                    raise ValueError("telemetry_grid shorter than smoke size")
                out.write(chunk)
                left -= len(chunk)
        return
    # Tier recipe: independent seeded correlated float32 rasters,
    # 2048x2048 (16 MiB each), same class as the data/ telemetry field.
    N, k = 2048, 49
    rng = np.random.Generator(np.random.Philox(key=_seed(tier, "images")))
    yy, xx = np.mgrid[0:N, 0:N].astype(np.float32) / np.float32(N)
    left = nbytes
    i = 0
    while left > 0:
        field = _box_smooth(
            rng.standard_normal((N, N), dtype=np.float32), k) * np.float32(k)
        phase = np.float32(0.37 * (i + 1))
        trend = np.sin(2 * np.pi * (xx + phase)) * np.cos(
            2 * np.pi * (yy - phase))
        raw = (field + trend.astype(np.float32)).astype("<f4").tobytes()
        out.write(raw[: min(len(raw), left)])
        left -= min(len(raw), left)
        i += 1


INNER_ZSTD_LEVEL = 3  # zstd CLI default; pre-registered (protocol 3.6)


def _gen_precompressed(out, nbytes: int, tier: str) -> None:
    """zstd -3 -T1 over the text recipe's stream, truncated to size.
    The inner compressor is pinned by version in the manifest record."""
    proc = subprocess.Popen(
        [ZSTD_BIN, "-q", f"-{INNER_ZSTD_LEVEL}", "-T1", "-c"],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL)

    class _Feeder:
        def __init__(self, pipe):
            self.pipe = pipe

        def write(self, b):
            self.pipe.write(b)

    import threading

    def feed():
        try:
            # Feed up to 8x the target size of text; zstd -19 on enwik*
            # yields ~0.26x, so 8x input always over-fills the target.
            _gen_text(_Feeder(proc.stdin), 8 * nbytes, tier)
        except BrokenPipeError:
            pass
        finally:
            try:
                proc.stdin.close()
            except BrokenPipeError:
                pass

    t = threading.Thread(target=feed, daemon=True)
    t.start()
    left = nbytes
    while left > 0:
        chunk = proc.stdout.read(min(_CHUNK, left))
        if not chunk:
            raise ValueError("precompressed stream ended before target size")
        out.write(chunk)
        left -= len(chunk)
    proc.stdout.close()
    proc.kill()
    proc.wait()
    t.join(timeout=60)


_GENERATORS = {
    "white_noise": _gen_white_noise,
    "text": _gen_text,
    "video": _gen_video,
    "images": _gen_images,
    "archive_mix": _gen_archive_mix,
    "precompressed": _gen_precompressed,
}


# ---------------------------------------------------------------------------
# Materialization with manifest pinning
# ---------------------------------------------------------------------------

def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(_CHUNK)
            if not chunk:
                return h.hexdigest()
            h.update(chunk)


def _manifest_path(tier: str) -> Path:
    return CORPORA_DIR / tier / "manifest.json"


def _load_manifest(tier: str) -> dict:
    p = _manifest_path(tier)
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {"tier": tier, "entries": {}}


def _save_manifest(tier: str, manifest: dict) -> None:
    p = _manifest_path(tier)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    tmp.replace(p)


def corpus_path(tier: str, ctype: str) -> Path:
    return CORPORA_DIR / tier / f"{ctype}.bin"


def _big_entries() -> dict:
    with open(BIG_MANIFEST) as f:
        return json.load(f)["entries"]


def available(tier: str, ctype: str) -> bool:
    """Is this (tier, ctype) cell defined?  Real tiers follow the data/big
    manifest (text is capped at 1000 MB, so text@1500/5000 do not exist by
    design); smoke tiers always exist."""
    if ctype not in CONTENT_TYPES:
        return False
    if tier not in REAL_TIERS:
        return True
    try:
        ent = _big_entries()[ctype]
    except (FileNotFoundError, KeyError):
        return False
    return int(tier) in ent.get("tiers_mb", [])


def _materialize_from_big(tier: str, ctype: str, log=print) -> dict:
    """Real tiers: stream the exact byte-prefix of the pinned data/big
    master into the corpus file and verify it against the manifest's
    per-tier sha256."""
    ent = _big_entries()[ctype]
    tier_mb = int(tier)
    tinfo = next(t for t in ent["tiers"] if t["mb"] == tier_mb)
    nbytes = tinfo["bytes"]
    master = BIG_DIR / ent["file"]
    if not master.exists():
        raise RuntimeError(
            f"data/big master missing for {ctype}: {master} -- run "
            f"data/big/fetch_big.py first")
    path = corpus_path(tier, ctype)
    manifest = _load_manifest(tier)
    entry = manifest["entries"].get(ctype)
    if entry and path.exists() and path.stat().st_size == entry["bytes"]:
        return entry
    path.parent.mkdir(parents=True, exist_ok=True)
    log(f"[corpora] slicing {tier}/{ctype} = first {nbytes} bytes of "
        f"{master.name} ...")
    t0 = time.perf_counter()
    tmp = path.with_suffix(".building")
    h = hashlib.sha256()
    left = nbytes
    with open(master, "rb") as src, open(tmp, "wb") as out:
        while left > 0:
            chunk = src.read(min(_CHUNK, left))
            if not chunk:
                raise RuntimeError(f"{master} shorter than {nbytes} bytes")
            out.write(chunk)
            h.update(chunk)
            left -= len(chunk)
    live = h.hexdigest()
    if live != tinfo["sha256"]:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"{tier}/{ctype}: prefix sha256 {live[:16]} != manifest "
            f"{tinfo['sha256'][:16]} -- data/big master corrupted?")
    tmp.replace(path)
    entry = {
        "type": ctype,
        "tier": tier,
        "file": path.name,
        "bytes": nbytes,
        "sha256": live,
        "source": "data/big prefix",
        "master_file": ent["file"],
        "master_sha256": ent["sha256"],
        "provenance": ent.get("recipe", ""),
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "build_s": round(time.perf_counter() - t0, 3),
    }
    manifest["entries"][ctype] = entry
    _save_manifest(tier, manifest)
    return entry


def materialize(tier: str, ctype: str, log=print) -> dict:
    """Ensure the (tier, ctype) corpus exists; return its manifest entry."""
    if tier not in TIER_BYTES:
        raise ValueError(f"unknown tier {tier!r}; known {sorted(TIER_BYTES)}")
    if ctype not in _GENERATORS:
        raise ValueError(f"unknown type {ctype!r}; known {CONTENT_TYPES}")
    if not available(tier, ctype):
        raise ValueError(
            f"cell {tier}/{ctype} does not exist (see data/big manifest: "
            f"text is capped at 1000 MB)")
    if tier in REAL_TIERS:
        return _materialize_from_big(tier, ctype, log=log)
    nbytes = TIER_BYTES[tier]
    path = corpus_path(tier, ctype)
    manifest = _load_manifest(tier)
    entry = manifest["entries"].get(ctype)
    if entry and path.exists() and path.stat().st_size == entry["bytes"]:
        return entry  # trusted: pinned at build time; campaign re-hashes once
    path.parent.mkdir(parents=True, exist_ok=True)
    log(f"[corpora] building {tier}/{ctype} ({nbytes} bytes) ...")
    t0 = time.perf_counter()
    tmp = path.with_suffix(".building")
    with open(tmp, "wb") as f:
        _GENERATORS[ctype](f, nbytes, tier)
    if tmp.stat().st_size != nbytes:
        raise RuntimeError(
            f"{tier}/{ctype}: generator wrote {tmp.stat().st_size} bytes, "
            f"expected {nbytes}")
    tmp.replace(path)
    entry = {
        "type": ctype,
        "tier": tier,
        "file": path.name,
        "bytes": nbytes,
        "sha256": _sha256_file(path),
        "seed": _seed(tier, ctype),
        "numpy": np.__version__,
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "build_s": round(time.perf_counter() - t0, 3),
    }
    if ctype == "precompressed":
        out = subprocess.run([ZSTD_BIN, "--version"],
                             capture_output=True, text=True)
        entry["inner_zstd"] = (out.stdout or out.stderr).strip()
    manifest["entries"][ctype] = entry
    _save_manifest(tier, manifest)
    log(f"[corpora] {tier}/{ctype} done in {entry['build_s']}s "
        f"sha256={entry['sha256'][:16]}...")
    return entry


if __name__ == "__main__":
    tier = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    types = sys.argv[2].split(",") if len(sys.argv) > 2 else CONTENT_TYPES
    for ct in types:
        materialize(tier, ct)
