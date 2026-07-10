#!/usr/bin/env python3
"""Baseline lossless-compressor runners for the RNR empirical program.

Implements the Python-library subset of the baseline matrix in the
experimental-design companion, section 4.1:

    gzip    (stdlib, zlib)          level 9            low baseline
    bz2     (stdlib)                level 9            block-sorting baseline
    lzma    (stdlib, .xz format)    preset 9 | EXTREME dictionary ceiling (xz -9e)
    zstd    (zstandard package)     level 22           modern high-throughput
    brotli  (brotli package)        quality 11         second modern baseline

zpaq, cmix, and NNCP require external binaries and are out of scope for
this in-process harness; at integration they plug in through the same
BaselineResult interface via subprocess wrappers.

Every runner performs a full round-trip (compress, decompress, compare
against the original; a mismatch raises RoundTripError) and returns a
BaselineResult carrying (compressed_size, enc_s, dec_s,
peak_rss_estimate).

Peak-RSS measurement modes
--------------------------
in-process (default): peak_rss_estimate is the increase of the process
    resident-set high-watermark (ru_maxrss) across the encode+decode
    calls. Because ru_maxrss is a monotone high-watermark, this is a
    LOWER BOUND and is frequently 0 after the first call in a process.
    Cheap; adequate for smoke tests and relative sanity checks only.
isolated (run_baseline(..., isolated=True)): the measurement runs in a
    fresh child interpreter; peak_rss_estimate is the child's own
    ru_maxrss at exit. This includes the interpreter baseline
    (~15-25 MB) but is a true per-measurement high-watermark. Use this
    mode for reported numbers.

This module deliberately imports neither numpy nor scipy so the
isolated-mode child stays lightweight.
"""

from __future__ import annotations

import bz2
import dataclasses
import gzip
import json
import lzma
import resource
import subprocess
import sys
import tempfile
import time
import zlib
from pathlib import Path

try:
    import zstandard

    _HAVE_ZSTD = True
except ImportError:  # pragma: no cover - environment without the package
    _HAVE_ZSTD = False

try:
    import brotli

    _HAVE_BROTLI = True
except ImportError:  # pragma: no cover - environment without the package
    _HAVE_BROTLI = False


class RoundTripError(RuntimeError):
    """Raised when decompress(compress(x)) != x for a baseline coder."""


@dataclasses.dataclass(frozen=True)
class BaselineResult:
    coder: str
    config: dict
    original_size: int
    compressed_size: int
    enc_s: float
    dec_s: float
    peak_rss_estimate: int  # bytes; see module docstring for semantics
    rss_mode: str  # "in-process-delta" or "isolated-child"

    @property
    def bits_per_byte(self) -> float:
        return 8.0 * self.compressed_size / self.original_size


def _ru_maxrss_bytes() -> int:
    """Process resident-set high-watermark in bytes (portable units)."""
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB.
    return raw if sys.platform == "darwin" else raw * 1024


# ---------------------------------------------------------------------------
# Codec table: name -> (default_config, compress_fn, decompress_fn)
# ---------------------------------------------------------------------------

def _gzip_c(data: bytes, cfg: dict) -> bytes:
    return gzip.compress(data, compresslevel=cfg["level"])


def _gzip_d(blob: bytes, cfg: dict) -> bytes:
    return gzip.decompress(blob)


def _bz2_c(data: bytes, cfg: dict) -> bytes:
    return bz2.compress(data, compresslevel=cfg["level"])


def _bz2_d(blob: bytes, cfg: dict) -> bytes:
    return bz2.decompress(blob)


def _lzma_c(data: bytes, cfg: dict) -> bytes:
    preset = cfg["preset"] | (lzma.PRESET_EXTREME if cfg["extreme"] else 0)
    return lzma.compress(data, preset=preset)


def _lzma_d(blob: bytes, cfg: dict) -> bytes:
    return lzma.decompress(blob)


def _zstd_c(data: bytes, cfg: dict) -> bytes:
    return zstandard.ZstdCompressor(level=cfg["level"]).compress(data)


def _zstd_d(blob: bytes, cfg: dict) -> bytes:
    return zstandard.ZstdDecompressor().decompress(blob)


def _brotli_c(data: bytes, cfg: dict) -> bytes:
    return brotli.compress(data, quality=cfg["quality"])


def _brotli_d(blob: bytes, cfg: dict) -> bytes:
    return brotli.decompress(blob)


_CODECS: dict = {
    "gzip": ({"level": 9}, _gzip_c, _gzip_d),
    "bz2": ({"level": 9}, _bz2_c, _bz2_d),
    "lzma": ({"preset": 9, "extreme": True}, _lzma_c, _lzma_d),
}
if _HAVE_ZSTD:
    _CODECS["zstd"] = ({"level": 22}, _zstd_c, _zstd_d)
if _HAVE_BROTLI:
    _CODECS["brotli"] = ({"quality": 11}, _brotli_c, _brotli_d)


def available_codecs() -> list:
    """Names of baseline coders usable in this environment."""
    return sorted(_CODECS)


def codec_versions() -> dict:
    """Pinned versions of every baseline coder, for the results record."""
    versions = {
        "python": sys.version.split()[0],
        "gzip": f"zlib {zlib.ZLIB_RUNTIME_VERSION} (stdlib gzip)",
        "bz2": "stdlib bz2 (bundled libbz2)",
        "lzma": "stdlib lzma (bundled liblzma, .xz container)",
    }
    if _HAVE_ZSTD:
        lib = ".".join(str(v) for v in zstandard.ZSTD_VERSION)
        versions["zstd"] = f"zstandard {zstandard.__version__} (libzstd {lib})"
    if _HAVE_BROTLI:
        versions["brotli"] = f"brotli {getattr(brotli, 'version', 'unknown')}"
    return versions


def _run_in_process(coder: str, data: bytes, config: dict) -> BaselineResult:
    default_cfg, comp, decomp = _CODECS[coder]
    cfg = dict(default_cfg)
    if config:
        cfg.update(config)

    rss_before = _ru_maxrss_bytes()

    t0 = time.perf_counter()
    blob = comp(data, cfg)
    enc_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    restored = decomp(blob, cfg)
    dec_s = time.perf_counter() - t0

    rss_after = _ru_maxrss_bytes()

    if restored != data:
        raise RoundTripError(
            f"{coder}: decompressed output differs from input "
            f"({len(restored)} vs {len(data)} bytes)"
        )

    return BaselineResult(
        coder=coder,
        config=cfg,
        original_size=len(data),
        compressed_size=len(blob),
        enc_s=enc_s,
        dec_s=dec_s,
        peak_rss_estimate=max(rss_after - rss_before, 0),
        rss_mode="in-process-delta",
    )


def _run_isolated(coder: str, data: bytes, config: dict) -> BaselineResult:
    """Run one measurement in a fresh child interpreter for a clean RSS
    high-watermark. See module docstring."""
    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as f:
        f.write(data)
        in_path = f.name
    try:
        proc = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--child",
                in_path,
                coder,
                json.dumps(config or {}),
            ],
            capture_output=True,
            text=True,
            timeout=1200,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"isolated {coder} child failed: {proc.stderr.strip()}"
            )
        payload = json.loads(proc.stdout)
    finally:
        Path(in_path).unlink(missing_ok=True)

    return BaselineResult(
        coder=coder,
        config=payload["config"],
        original_size=payload["original_size"],
        compressed_size=payload["compressed_size"],
        enc_s=payload["enc_s"],
        dec_s=payload["dec_s"],
        peak_rss_estimate=payload["child_ru_maxrss_bytes"],
        rss_mode="isolated-child",
    )


def run_baseline(
    coder: str,
    data: bytes,
    config: dict | None = None,
    *,
    isolated: bool = False,
) -> BaselineResult:
    """Compress+decompress `data` with a named baseline coder.

    Returns a BaselineResult with (compressed_size, enc_s, dec_s,
    peak_rss_estimate). Raises KeyError for unknown coders and
    RoundTripError if the round-trip fails.
    """
    if coder not in _CODECS:
        raise KeyError(
            f"unknown or unavailable coder {coder!r}; "
            f"available: {available_codecs()}"
        )
    if isolated:
        return _run_isolated(coder, data, config or {})
    return _run_in_process(coder, data, config or {})


def _child_main(argv: list) -> int:
    """Entry point for isolated-mode children: measure one coder run and
    print a JSON payload on stdout."""
    in_path, coder, config_json = argv
    data = Path(in_path).read_bytes()
    result = _run_in_process(coder, data, json.loads(config_json))
    print(
        json.dumps(
            {
                "config": result.config,
                "original_size": result.original_size,
                "compressed_size": result.compressed_size,
                "enc_s": result.enc_s,
                "dec_s": result.dec_s,
                "child_ru_maxrss_bytes": _ru_maxrss_bytes(),
            }
        )
    )
    return 0


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "--child":
        sys.exit(_child_main(sys.argv[2:]))
    # Bare invocation: print environment summary.
    print("available baseline coders:", ", ".join(available_codecs()))
    for name, ver in codec_versions().items():
        print(f"  {name}: {ver}")
