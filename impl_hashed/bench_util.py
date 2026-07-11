#!/usr/bin/env python3
"""bench_util.py -- shared measurement helpers for the impl_hashed benchmarks.

  * corpus slices        text (enwik8 prefix) and code (assembled source);
  * external compressors  gzip/bzip2/xz/zstd/brotli via /usr/bin/time -l so
                          every run reports wall time AND peak RSS;
  * model cross-entropy   ideal code length of a predictor's distribution.

All external runs are single-threaded (-T1 where supported) for a fair
comparison against the single-process reference coder, and peak RSS is read
from the macOS `/usr/bin/time -l` "maximum resident set size" line (bytes).
"""

import json
import os
import tempfile
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_RNR = _HERE.parent
sys.path.insert(0, str(_RNR / "data"))

MB = 1 << 20

# Code corpus: assembled deterministically from Python/C source trees (this
# host's stdlib + venv site-packages), path-sorted, truncated to 64 MiB.
# Not checked in (measurement input only); regenerated on demand so the
# scripts are self-contained.  Override with $HASHED_CODE_CORPUS.
_SCRATCH = os.environ.get(
    "HASHED_SCRATCH",
    os.path.join(tempfile.gettempdir(), "rnr_hashed_scratch"))
os.makedirs(_SCRATCH, exist_ok=True)
CODE_CORPUS = Path(os.environ.get("HASHED_CODE_CORPUS",
                                  os.path.join(_SCRATCH, "code_corpus.bin")))


def ensure_code_corpus(target=64 * MB):
    """Build the code corpus if missing; return its path."""
    if CODE_CORPUS.exists() and CODE_CORPUS.stat().st_size >= target:
        return CODE_CORPUS
    import sysconfig
    import site
    roots = [sysconfig.get_path("stdlib"), site.getsitepackages()[0]]
    files = []
    for r in roots:
        for dp, _, fns in os.walk(r):
            for fn in fns:
                if fn.endswith((".py", ".c", ".h", ".pyi")):
                    files.append(os.path.join(dp, fn))
    files.sort()
    CODE_CORPUS.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(CODE_CORPUS, "wb") as w:
        for f in files:
            try:
                b = open(f, "rb").read()
            except OSError:
                continue
            w.write(b)
            n += len(b)
            if n >= target:
                break
    return CODE_CORPUS


def load_slice(corpus, size):
    """Return the first `size` bytes of corpus in {'text','code'}."""
    if corpus == "text":
        import loader
        return loader.read_range("enwik8", 0, size)
    if corpus == "code":
        ensure_code_corpus(max(size, 64 * MB))
        with open(CODE_CORPUS, "rb") as f:
            return f.read(size)
    raise ValueError("unknown corpus %r" % corpus)


# ---------------------------------------------------------------------------
# External compressors
# ---------------------------------------------------------------------------

def _tool(name):
    return shutil.which(name)


EXTERNAL = {
    "gzip-9":   (["gzip", "-9"],                 ["gzip", "-d"]),
    "bz2-9":    (["bzip2", "-9"],                ["bzip2", "-d"]),
    "xz-6":     (["xz", "-6", "-T1"],            ["xz", "-d", "-T1"]),
    "xz-9e":    (["xz", "-9", "-e", "-T1"],      ["xz", "-d", "-T1"]),
    "zstd-19":  (["zstd", "-19", "-T1", "-q"],   ["zstd", "-d", "-T1", "-q"]),
    "zstd-22u": (["zstd", "--ultra", "-22", "-T1", "-q"],
                 ["zstd", "-d", "-T1", "-q"]),
    "brotli-11": (["brotli", "-q", "11"],        ["brotli", "-d"]),
}


def _timed(argv, out_path):
    """Run argv (a tool + flags + input file) with stdout -> out_path under
    /usr/bin/time -l.  Return (wall_s, peak_rss_mb, returncode)."""
    tool = shutil.which(argv[0])
    if tool is None:
        return None
    full = ["/usr/bin/time", "-l", tool] + argv[1:]
    t0 = time.perf_counter()
    with open(out_path, "wb") as fo:
        p = subprocess.run(full, stdout=fo, stderr=subprocess.PIPE)
    wall = time.perf_counter() - t0
    m = re.search(rb"(\d+)\s+maximum resident set size", p.stderr)
    peak = int(m.group(1)) / 1e6 if m else float("nan")
    return wall, peak, p.returncode


def run_external(name, in_path, tmpdir):
    """Compress + decompress in_path with the named external tool.

    Returns dict{comp_bytes,bpb,enc_s,enc_mb_s,dec_s,dec_mb_s,
    enc_peak_mb,dec_peak_mb,ok} or None if the tool is unavailable.
    """
    if name not in EXTERNAL:
        raise KeyError(name)
    encflags, decflags = EXTERNAL[name]
    if shutil.which(encflags[0]) is None:
        return None
    in_path = Path(in_path)
    n = in_path.stat().st_size
    comp = Path(tmpdir) / (name + ".comp")
    dec = Path(tmpdir) / (name + ".dec")
    enc = _timed(encflags + ["-c", str(in_path)], comp)
    if enc is None or enc[2] != 0:
        return None
    dec_r = _timed(decflags + ["-c", str(comp)], dec)
    if dec_r is None or dec_r[2] != 0:
        return None
    comp_bytes = comp.stat().st_size
    ok = (dec.stat().st_size == n and
          _same_file(dec, in_path))
    res = {
        "comp_bytes": comp_bytes,
        "bpb": 8.0 * comp_bytes / max(1, n),
        "enc_s": enc[0], "enc_mb_s": n / 1e6 / max(enc[0], 1e-9),
        "enc_peak_mb": enc[1],
        "dec_s": dec_r[0], "dec_mb_s": n / 1e6 / max(dec_r[0], 1e-9),
        "dec_peak_mb": dec_r[1],
        "ok": ok,
    }
    comp.unlink(missing_ok=True)
    dec.unlink(missing_ok=True)
    return res


def _same_file(a, b, chunk=1 << 20):
    with open(a, "rb") as fa, open(b, "rb") as fb:
        while True:
            ba, bb = fa.read(chunk), fb.read(chunk)
            if ba != bb:
                return False
            if not ba:
                return True


# ---------------------------------------------------------------------------
# Model cross-entropy (rate lower bound; the coder redundancy is < 1e-4 bpb)
# ---------------------------------------------------------------------------

def model_ce_bpb(pred, data):
    """Ideal code length (bpb) of the predictor's own distribution."""
    import math
    bits = 0.0
    for x in data:
        freq, total = pred.dist()
        bits += -math.log2(freq[x] / total)
        pred.update(x)
    return bits / max(1, len(data))


def csv_write(path, rows, fields):
    import csv
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


if __name__ == "__main__":
    ensure_code_corpus()
    print("code corpus:", CODE_CORPUS, CODE_CORPUS.stat().st_size, "bytes")
    print("external tools available:",
          [k for k in EXTERNAL if shutil.which(EXTERNAL[k][0][0])])
