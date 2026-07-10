#!/usr/bin/env python3
"""Coder registry + single-run measurement harness for the Type-I scale
benchmark.

Every coder is exercised as an external process so that encode/decode
wall time and peak RSS are measured the same way for every entry in the
matrix (protocol.md section 4):

  wall time  parent time.perf_counter() around the child's lifetime;
  peak RSS   /usr/bin/time -l (macOS) "maximum resident set size",
             reported in bytes, one child per phase (rss_mode
             "time-l-child");
  round trip compressed tempfile is size-measured, then decode output is
             STREAMED into sha256 (never written to disk for stdout
             coders) and compared to the corpus hash; the tempfile is
             deleted before the function returns, so disk usage stays
             bounded by corpus + one compressed file.

The RNR coders decode to a tempfile (their CLI writes files, not
stdout); the file is hashed and deleted immediately.  This difference
is recorded per record in `decode_io` and discussed in protocol 4.3.
"""

from __future__ import annotations

import hashlib
import os
import signal
import subprocess
import sys
import time
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PYTHON = sys.executable

TIME_BIN = "/usr/bin/time"
ZSTD = "/opt/homebrew/bin/zstd"
XZ = "/opt/homebrew/bin/xz"
BROTLI = "/opt/homebrew/bin/brotli"
GZIP = "/usr/bin/gzip"
BZIP2 = "/usr/bin/bzip2"

RNR_REF = REPO / "impl" / "rnr1.py"
IMPL_FAST_DIR = REPO / "impl_fast"
IMPL_FAST_ADAPTER = IMPL_FAST_DIR / "BENCH_ADAPTER.json"

_CHUNK = 8 * 1024 * 1024


class CoderUnavailable(RuntimeError):
    pass


class CellTimeout(RuntimeError):
    def __init__(self, phase: str, elapsed_s: float):
        super().__init__(f"time cap hit during {phase} after {elapsed_s:.0f}s")
        self.phase = phase
        self.elapsed_s = elapsed_s


# ---------------------------------------------------------------------------
# Registry.  Each entry:
#   enc(in_path)  -> (argv, stdout_mode)  stdout_mode: True => compressed
#                    stream on stdout; False => argv contains {out}
#   dec(enc_path) -> (argv, stdout_mode)  stdout_mode: True => decoded
#                    stream on stdout; False => argv contains {out}
#   tiers: set of tiers the coder runs at (None = all)
# ---------------------------------------------------------------------------

def _stdout_coder(enc_argv_fn, dec_argv_fn):
    return {
        "enc": lambda ip: (enc_argv_fn(str(ip)), True),
        "dec": lambda ep: (dec_argv_fn(str(ep)), True),
    }


REGISTRY: dict = {
    "zstd-19": {
        "family": "zstd",
        "config": {"level": 19, "threads": 1},
        **_stdout_coder(lambda i: [ZSTD, "-q", "-19", "-T1", "-c", i],
                        lambda e: [ZSTD, "-d", "-q", "-c", e]),
        "version": [ZSTD, "--version"],
        "tiers": None,
    },
    "zstd-22u": {
        "family": "zstd",
        "config": {"level": 22, "ultra": True, "threads": 1},
        **_stdout_coder(
            lambda i: [ZSTD, "-q", "--ultra", "-22", "-T1", "-c", i],
            lambda e: [ZSTD, "-d", "-q", "--memory=2048MB", "-c", e]),
        "version": [ZSTD, "--version"],
        "tiers": None,
    },
    "xz-6": {
        "family": "xz",
        "config": {"preset": 6, "threads": 1},
        **_stdout_coder(lambda i: [XZ, "-6", "-T1", "-k", "-c", i],
                        lambda e: [XZ, "-d", "-T1", "-c", e]),
        "version": [XZ, "--version"],
        "tiers": None,
    },
    "xz-9e": {
        "family": "xz",
        "config": {"preset": "9e", "threads": 1},
        **_stdout_coder(lambda i: [XZ, "-9", "-e", "-T1", "-k", "-c", i],
                        lambda e: [XZ, "-d", "-T1", "-c", e]),
        "version": [XZ, "--version"],
        "tiers": None,
    },
    "brotli-q11": {
        "family": "brotli",
        "config": {"quality": 11, "lgwin": 24},
        **_stdout_coder(
            lambda i: [BROTLI, "-q", "11", "--lgwin=24", "-c", i],
            lambda e: [BROTLI, "-d", "-c", e]),
        "version": [BROTLI, "--version"],
        # q11 pre-registered as infeasible at 5000 MB (protocol 4.2).
        "tiers": {"smoke", "500", "1500"},
    },
    "brotli-q9": {
        "family": "brotli",
        "config": {"quality": 9, "lgwin": 24},
        **_stdout_coder(
            lambda i: [BROTLI, "-q", "9", "--lgwin=24", "-c", i],
            lambda e: [BROTLI, "-d", "-c", e]),
        "version": [BROTLI, "--version"],
        # The 5000 MB stand-in for brotli (and smoke, to validate it).
        "tiers": {"smoke", "5000"},
    },
    "gzip-9": {
        "family": "gzip",
        # -n omits the input name+mtime from the header: without it gzip
        # embeds the source file's timestamp and archives are NOT
        # repeatable across regenerated inputs (caught by the
        # compressed_sha256 determinism witness).
        "config": {"level": 9, "no_name": True},
        **_stdout_coder(lambda i: [GZIP, "-9", "-n", "-c", i],
                        lambda e: [GZIP, "-d", "-c", e]),
        "version": [GZIP, "--version"],
        "tiers": None,
    },
    "bzip2-9": {
        "family": "bzip2",
        "config": {"level": 9},
        **_stdout_coder(lambda i: [BZIP2, "-9", "-c", i],
                        lambda e: [BZIP2, "-d", "-c", e]),
        "version": [BZIP2, "--help"],  # bzip2 prints version in --help header
        "tiers": None,
    },
    # Reference RNR Type-I coder: normative but ~0.05-0.07 MB/s, so it is
    # smoke-only (a 500 MB cell would need ~2.5 h > the 90 min cap;
    # protocol 4.1).  The scale tiers use rnr1-fast from impl_fast/.
    "rnr1-ref": {
        "family": "rnr1",
        "config": {"W": 3, "K": 65536, "impl": "reference"},
        "enc": lambda ip: ([PYTHON, "-u", str(RNR_REF), "pack", str(ip),
                            "{out}", "--W", "3", "--K", "65536"], False),
        "dec": lambda ep: ([PYTHON, "-u", str(RNR_REF), "unpack", str(ep),
                            "{out}"], False),
        "version": None,
        "tiers": {"smoke"},
    },
}


def _load_impl_fast() -> None:
    """Register rnr1-fast / rnr1-cm-fast when impl_fast/ provides them.

    Contract (protocol 4.1): impl_fast/BENCH_ADAPTER.json maps coder
    names to {"pack": argv, "unpack": argv, "config": {...}} templates
    with {in}/{out} placeholders; relative argv[0] paths resolve inside
    impl_fast/.  Fallback: impl_fast/rnr1_fast.py with an
    impl/rnr1.py-compatible pack/unpack CLI registers rnr1-fast only.
    """
    import json

    def _tpl(argv_tpl, inp, out):
        argv = [a.replace("{in}", str(inp)).replace("{out}", str(out))
                for a in argv_tpl]
        if argv and not Path(argv[0]).is_absolute() and (
                IMPL_FAST_DIR / argv[0]).exists():
            argv[0] = str(IMPL_FAST_DIR / argv[0])
        if argv and argv[0].endswith(".py"):
            argv = [PYTHON, "-u"] + argv
        return argv

    if IMPL_FAST_ADAPTER.exists():
        with open(IMPL_FAST_ADAPTER) as f:
            spec = json.load(f)
        for name, c in spec.get("coders", {}).items():
            pack_t, unpack_t = list(c["pack"]), list(c["unpack"])
            REGISTRY[name] = {
                "family": "rnr1",
                "config": dict(c.get("config", {})),
                "enc": (lambda ip, t=pack_t: (_tpl(t, ip, "{out}"), False)),
                "dec": (lambda ep, t=unpack_t: (_tpl(t, ep, "{out}"), False)),
                "version": None,
                "tiers": None,
            }
        return
    fast_py = IMPL_FAST_DIR / "rnr1_fast.py"
    if fast_py.exists():
        REGISTRY["rnr1-fast"] = {
            "family": "rnr1",
            "config": {"W": 3, "K": 65536, "impl": "fast"},
            "enc": lambda ip: ([PYTHON, "-u", str(fast_py), "pack", str(ip),
                                "{out}", "--W", "3", "--K", "65536"], False),
            "dec": lambda ep: ([PYTHON, "-u", str(fast_py), "unpack",
                                str(ep), "{out}"], False),
            "version": None,
            "tiers": None,
        }


_load_impl_fast()


def available_coders(tier: str | None = None) -> list:
    names = []
    for name, spec in REGISTRY.items():
        if tier is not None and spec["tiers"] is not None \
                and tier not in spec["tiers"]:
            continue
        names.append(name)
    return sorted(names)


def coder_version(name: str) -> str:
    spec = REGISTRY[name]
    if spec["version"] is None:
        return f"repo:{REPO.name} (python {sys.version.split()[0]})"
    try:
        out = subprocess.run(spec["version"], capture_output=True,
                             text=True, timeout=30)
        text = (out.stdout or out.stderr).strip().splitlines()
        return text[0].strip() if text else "unknown"
    except OSError:
        return "unavailable"


# ---------------------------------------------------------------------------
# Measured subprocess execution
# ---------------------------------------------------------------------------

def _parse_max_rss(stderr_text: str) -> int | None:
    for line in stderr_text.splitlines():
        parts = line.split()
        if "maximum" in parts and "resident" in parts:
            try:
                return int(parts[0])  # bytes on darwin
            except ValueError:
                return None
    return None


def _run_timed(argv: list, *, stdout_path: Path | None,
               stream_hash: bool, timeout_s: float, phase: str):
    """Run argv under /usr/bin/time -l.

    stdout_path: file that receives the child's stdout (compressed
        output or decoded output), or None when stream_hash is set.
    stream_hash: drain child stdout into sha256 without touching disk.
    Returns (wall_s, max_rss_bytes, sha256_hex_or_None, stdout_bytes).
    Raises CellTimeout on cap expiry (the process group is killed) and
    RuntimeError on a nonzero exit.
    """
    full = [TIME_BIN, "-l"] + argv
    t0 = time.perf_counter()
    out_fh = open(stdout_path, "wb") if stdout_path is not None else None
    try:
        proc = subprocess.Popen(
            full,
            stdout=(subprocess.PIPE if stream_hash else out_fh),
            stderr=subprocess.PIPE,
            start_new_session=True,
        )
        hasher = hashlib.sha256() if stream_hash else None
        nbytes = 0
        try:
            if stream_hash:
                deadline = t0 + timeout_s
                while True:
                    if time.perf_counter() > deadline:
                        raise subprocess.TimeoutExpired(full, timeout_s)
                    chunk = proc.stdout.read(_CHUNK)
                    if not chunk:
                        break
                    hasher.update(chunk)
                    nbytes += len(chunk)
                proc.stdout.close()
                remaining = max(1.0, deadline - time.perf_counter())
                stderr = proc.communicate(timeout=remaining)[1]
            else:
                stderr = proc.communicate(timeout=timeout_s)[1]
        except subprocess.TimeoutExpired:
            elapsed = time.perf_counter() - t0
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()
            raise CellTimeout(phase, elapsed) from None
        wall = time.perf_counter() - t0
    finally:
        if out_fh is not None:
            out_fh.close()
    if proc.returncode != 0:
        tail = stderr.decode(errors="replace").strip()[-400:]
        raise RuntimeError(f"{phase} failed rc={proc.returncode}: {tail}")
    rss = _parse_max_rss(stderr.decode(errors="replace"))
    return wall, rss, (hasher.hexdigest() if hasher else None), nbytes


def _hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(_CHUNK)
            if not chunk:
                return h.hexdigest()
            h.update(chunk)


def run_cell(coder: str, corpus_file: Path, corpus_sha256: str,
             scratch_dir: Path, cap_s: float) -> dict:
    """One full encode -> size -> decode-verify -> cleanup measurement.

    Returns a measurement dict; raises CellTimeout (with .phase) when
    the per-cell cap expires, RuntimeError on coder failure.  All
    tempfiles are deleted on every path.
    """
    spec = REGISTRY[coder]
    scratch_dir.mkdir(parents=True, exist_ok=True)
    tag = uuid.uuid4().hex[:12]
    enc_path = scratch_dir / f"{coder}-{tag}.enc"
    dec_path = scratch_dir / f"{coder}-{tag}.dec"
    t_cell = time.perf_counter()

    def _left() -> float:
        return cap_s - (time.perf_counter() - t_cell)

    try:
        # --- encode ---
        argv, enc_stdout = spec["enc"](corpus_file)
        if enc_stdout:
            enc_s, enc_rss, _, _ = _run_timed(
                argv, stdout_path=enc_path, stream_hash=False,
                timeout_s=_left(), phase="encode")
        else:
            argv = [a.replace("{out}", str(enc_path)) for a in argv]
            enc_s, enc_rss, _, _ = _run_timed(
                argv, stdout_path=None, stream_hash=False,
                timeout_s=_left(), phase="encode")
        compressed_size = enc_path.stat().st_size
        # Repeatability witness: every repetition of a cell must produce a
        # byte-identical archive (all registered coders are deterministic
        # single-config invocations).  Hashed OUTSIDE the timed sections.
        compressed_sha256 = _hash_file(enc_path)

        # --- decode + verify ---
        argv, dec_stdout = spec["dec"](enc_path)
        if dec_stdout:
            dec_s, dec_rss, digest, ndec = _run_timed(
                argv, stdout_path=None, stream_hash=True,
                timeout_s=_left(), phase="decode")
            decode_io = "stdout-stream"
        else:
            argv = [a.replace("{out}", str(dec_path)) for a in argv]
            dec_s, dec_rss, _, _ = _run_timed(
                argv, stdout_path=None, stream_hash=False,
                timeout_s=_left(), phase="decode")
            digest = _hash_file(dec_path)
            ndec = dec_path.stat().st_size
            decode_io = "tempfile"

        original_size = corpus_file.stat().st_size
        return {
            "original_size": original_size,
            "compressed_size": compressed_size,
            "compressed_sha256": compressed_sha256,
            "bpb": 8.0 * compressed_size / original_size,
            "ratio": compressed_size / original_size,
            "enc_s": enc_s,
            "dec_s": dec_s,
            "enc_peak_rss": enc_rss,
            "dec_peak_rss": dec_rss,
            # bench/results.py summarize() compatibility field:
            "peak_rss_estimate": max(enc_rss or 0, dec_rss or 0) or None,
            "roundtrip_ok": (digest == corpus_sha256
                             and ndec == original_size),
            "decoded_sha256": digest,
            "decode_io": decode_io,
        }
    finally:
        enc_path.unlink(missing_ok=True)
        dec_path.unlink(missing_ok=True)


if __name__ == "__main__":
    tier = sys.argv[1] if len(sys.argv) > 1 else None
    print(f"coders available{' at tier ' + tier if tier else ''}:")
    for name in available_coders(tier):
        print(f"  {name:<12} {coder_version(name)}")
