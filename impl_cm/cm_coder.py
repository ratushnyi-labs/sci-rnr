#!/usr/bin/env python3
"""cm_coder.py -- RNR Type-I container with the context-mixing predictor.

Composition, not duplication: this module IMPORTS the reference coder
(impl/rnr1.py: arithmetic coder, Type-I class/payload partition,
sub-block container, seek index, CLI machinery) and injects the
CMPredictor of impl_cm/cm_predictor.py through the module's predictor
seam (the NGramPredictor / model_description_hash globals used by
encode_subblock, decode_subblock, pack and Archive).  The only rnr1
logic transcribed here is the ~15-line per-sub-block encode loop, which
is re-stated verbatim so an optional cross-entropy trace (offline float
analysis, R-3.6) can observe the predictor's integer distribution at
each position without touching the coded bitstream.

Usage:
    from cm_coder import pack_cm, unpack_cm, read_cm, TraceRecorder
    raw = pack_cm(data, W=3, K=65536)          # bit-identical archive
    assert unpack_cm(raw) == data
    rec = TraceRecorder()
    raw = pack_cm(data, W=3, K=65536, trace=rec)   # + per-block CE

CLI (two-process determinism harness):
    python -u cm_coder.py pack <in> <out> [--W 3] [--K 65536]
    python -u cm_coder.py unpack <in> <out>
    python -u cm_coder.py info <archive>
"""

import argparse
import contextlib
import hashlib
import math
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "impl"))

import rnr1  # noqa: E402
import cm_predictor  # noqa: E402
from cm_predictor import CMPredictor  # noqa: E402


def cm_model_description_hash(W):
    """R-13.3 model description hash for the CM predictor variant."""
    desc = (
        "rnr-type1-code12a;predictor=cm-logistic;W=%d;pscale=%d;"
        "kt=1/2;lr=%d;wfix=%d;wclamp=%d;acc=%d;cap=%d;cand=%d;floor=1;"
        "classes=c0,c1,c2,cand,esc;squash=%s;e12=%s"
        % (
            W,
            cm_predictor.PSCALE,
            cm_predictor.LR_NUM,
            cm_predictor.WEIGHT_FIX_BITS,
            cm_predictor.WEIGHT_CLAMP,
            cm_predictor.ACC_INIT,
            cm_predictor.COUNT_CAP,
            rnr1.NUM_CANDIDATES,
            hashlib.sha256(
                cm_predictor.SQUASH.astype("<i8").tobytes()
            ).hexdigest()[:16],
            hashlib.sha256(rnr1.E12.astype("<i8").tobytes()).hexdigest()[:16],
        )
    ).encode("ascii")
    return hashlib.sha256(desc).digest()[:8]


class TraceRecorder:
    """Per-sub-block offline analysis trace (floats permitted, R-3.6).

    blocks: list of dicts {n, ce_bits, blob_bytes}; ce_bits is the
    model cross-entropy sum of -log2(freq[x]/total) over the block --
    the ideal code length of the predictor's own integer distribution.
    """

    def __init__(self):
        self.blocks = []


_TRACE = None  # active TraceRecorder or None


def _traced_encode_subblock(block, W):
    """rnr1.encode_subblock with an optional cross-entropy observer.

    The integer coding path is identical statement-for-statement to
    rnr1.encode_subblock (predictor resolved through the patched
    rnr1.NGramPredictor seam); the trace lines are analysis-only and do
    not influence any coded value.
    """
    pred = rnr1.NGramPredictor(W)
    enc = rnr1.ArithmeticEncoder()
    ce_bits = 0.0
    for x in block:
        freq, total = pred.dist()
        bstar, cand, classid, class_freq = rnr1._partition(freq)
        cls = int(classid[x])
        ccum = np.concatenate(([0], np.cumsum(class_freq)))
        enc.encode(int(ccum[cls]), int(ccum[cls + 1]), total)
        if cls != rnr1.CLS_0:
            members = rnr1._class_members(cls, cand, classid)
            mcum = np.concatenate(([0], np.cumsum(freq[members])))
            idx = int(np.nonzero(members == x)[0][0])
            enc.encode(int(mcum[idx]), int(mcum[idx + 1]), int(mcum[-1]))
        if _TRACE is not None:
            ce_bits += -math.log2(freq[x] / total)
        pred.update(x)
    blob = enc.finish()
    if _TRACE is not None:
        _TRACE.blocks.append(
            {"n": len(block), "ce_bits": ce_bits, "blob_bytes": len(blob)}
        )
    return blob


@contextlib.contextmanager
def cm_patched(trace=None):
    """Swap the CM predictor (and its model hash) into the rnr1 seam."""
    global _TRACE
    saved = (
        rnr1.NGramPredictor,
        rnr1.model_description_hash,
        rnr1.encode_subblock,
        _TRACE,
    )
    rnr1.NGramPredictor = CMPredictor
    rnr1.model_description_hash = cm_model_description_hash
    rnr1.encode_subblock = _traced_encode_subblock
    _TRACE = trace
    try:
        yield
    finally:
        (
            rnr1.NGramPredictor,
            rnr1.model_description_hash,
            rnr1.encode_subblock,
            _TRACE,
        ) = saved


@contextlib.contextmanager
def ngram_traced(trace=None):
    """Keep the reference n-gram predictor but enable the CE trace."""
    global _TRACE
    saved = (rnr1.encode_subblock, _TRACE)
    rnr1.encode_subblock = _traced_encode_subblock
    _TRACE = trace
    try:
        yield
    finally:
        rnr1.encode_subblock, _TRACE = saved


def pack_cm(data, W=rnr1.DEFAULT_W, K=rnr1.DEFAULT_K, trace=None):
    with cm_patched(trace=trace):
        return rnr1.pack(data, W=W, K=K)


def unpack_cm(raw, verify=True):
    with cm_patched():
        return rnr1.unpack(raw, verify=verify)


def read_cm(raw, p, k, verify=False):
    with cm_patched():
        return rnr1.Archive(raw).read(p, k, verify=verify)


def pack_ngram(data, W=rnr1.DEFAULT_W, K=rnr1.DEFAULT_K, trace=None):
    """Reference n-gram coder with the same trace instrumentation."""
    with ngram_traced(trace=trace):
        return rnr1.pack(data, W=W, K=K)


def archive_stats(raw):
    """Container decomposition of an archive (either predictor).

    Returns dict with n, archive_bytes, repair_bytes, container_bytes,
    per_block_coded_bytes (list), modes (list).
    """
    # Header fields are predictor-agnostic; bypass the model-hash check
    # by parsing with whichever hash matches.
    try:
        arc = rnr1.Archive(raw)
    except ValueError:
        with cm_patched():
            arc = rnr1.Archive(raw)
    per_block = []
    for j in range(arc.m):
        start = arc.entries[j][1] // 8
        end = arc.entries[j + 1][1] // 8 if j + 1 < arc.m else len(arc.repair)
        per_block.append(end - start)
    return {
        "n": arc.n,
        "K": arc.K,
        "W": arc.W,
        "m": arc.m,
        "archive_bytes": len(raw),
        "repair_bytes": len(arc.repair),
        "container_bytes": len(raw) - len(arc.repair),
        "per_block_coded_bytes": per_block,
        "modes": [e[3] for e in arc.entries],
    }


def _cli(argv=None):
    ap = argparse.ArgumentParser(
        prog="cm_coder",
        description="RNR Type-I coder with the context-mixing predictor "
        "(container/CLI delegated to impl/rnr1.py).",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--W", type=int, default=rnr1.DEFAULT_W)
    p.add_argument("--K", type=int, default=rnr1.DEFAULT_K)
    p = sub.add_parser("unpack")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--no-verify", action="store_true")
    p = sub.add_parser("info")
    p.add_argument("archive")
    args = ap.parse_args(argv)

    if args.cmd == "pack":
        data = Path(args.input).read_bytes()
        raw = pack_cm(data, W=args.W, K=args.K)
        Path(args.output).write_bytes(raw)
        print(
            "packed %d bytes -> %d bytes (%.4f bpb)"
            % (len(data), len(raw), 8.0 * len(raw) / max(1, len(data)))
        )
    elif args.cmd == "unpack":
        raw = Path(args.input).read_bytes()
        data = unpack_cm(raw, verify=not args.no_verify)
        Path(args.output).write_bytes(data)
        print("unpacked %d bytes -> %d bytes" % (len(raw), len(data)))
    else:
        st = archive_stats(Path(args.archive).read_bytes())
        for k in ("n", "W", "K", "m", "archive_bytes", "repair_bytes",
                  "container_bytes"):
            print("  %-16s %s" % (k, st[k]))


if __name__ == "__main__":
    _cli()
