#!/usr/bin/env python3
"""hashed_coder.py -- RNR Type-I container with the hashed CM predictor.

Composition, not duplication: this module IMPORTS the reference coder
(impl/rnr1.py: arithmetic coder, Type-I class/payload partition, sub-block
container, seek index) and injects HashedCMPredictor through rnr1's
predictor seam (rnr1.NGramPredictor / rnr1.model_description_hash), exactly
as impl_cm/cm_coder.py injects the exact CM predictor.

Two differences from cm_coder, both required by the memory-bounded model:

  1.  W may exceed the reference cap of 8 (the hashed table decouples the
      resident set from W, so high orders are the whole point).  The header
      W field is uint16, and Archive imposes no W range, so the only place
      that rejects W>8 is rnr1.pack's guard; pack_hashed re-states the
      ~25-line container loop verbatim without that guard.

  2.  WARM (whole-block) default.  The hashed tables must fill before they
      predict, and a fresh 2^HBITS-bucket table is almost empty after a
      64 KiB sub-block -- so the reference K=65536 would leave the model
      cold and the rate near 8 bpb.  pack_hashed therefore defaults to a
      SINGLE warm sub-block (K >= len(data)); pass an explicit smaller K to
      get the reset-per-sync random-access variant (each sync point zeroes
      the buckets, trading warm-up for seekability).  This is the
      "whole-block warm variant" the task permits.

Usage:
    from hashed_coder import pack_hashed, unpack_hashed
    raw = pack_hashed(data, W=12, HBITS=22)     # bit-identical archive
    assert unpack_hashed(raw) == data

CLI:
    python -u hashed_coder.py pack   <in> <out> [--W 12] [--HBITS 22] [--K N]
    python -u hashed_coder.py unpack <in> <out>
    python -u hashed_coder.py info   <archive>
"""

import argparse
import contextlib
import hashlib
import struct
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "impl_cm"))
sys.path.insert(0, str(_HERE.parent / "impl"))

import rnr1  # noqa: E402
import cm_predictor  # noqa: E402
import hashed_predictor  # noqa: E402
from hashed_predictor import HashedCMPredictor  # noqa: E402


# Bound HBITS for the current pack/unpack (rnr1's seam calls the predictor
# factory with W only; HBITS + any semi-static weights travel via globals).
_HBITS = HashedCMPredictor.__init__.__defaults__[0]
_FIXED_WEIGHTS = None      # (W+1,256) int64 shipped weights, or None (online)
W_MAX = 24


def _predictor_factory(W):
    p = HashedCMPredictor(W, HBITS=_HBITS)
    if _FIXED_WEIGHTS is not None:      # semi-static (multi-pass) coded pass
        p.load_weights(_FIXED_WEIGHTS)
        p.freeze = True
    return p


def hashed_model_description_hash(W):
    """R-13.3 model description hash for the hashed CM predictor variant."""
    desc = (
        "rnr-type1-code12a;predictor=hashed-cm-logistic;W=%d;hbits=%d;"
        "pscale=%d;kt=1/2;lr=%d;wfix=%d;wclamp=%d;acc=%d;bitcap=%d;cand=%d;"
        "floor=1;classes=c0,c1,c2,cand,esc;squash=%s;e12=%s"
        % (
            W,
            _HBITS,
            cm_predictor.PSCALE,
            cm_predictor.LR_NUM,
            cm_predictor.WEIGHT_FIX_BITS,
            cm_predictor.WEIGHT_CLAMP,
            cm_predictor.ACC_INIT,
            hashed_predictor.BIT_COUNT_CAP,
            rnr1.NUM_CANDIDATES,
            hashlib.sha256(
                cm_predictor.SQUASH.astype("<i8").tobytes()
            ).hexdigest()[:16],
            hashlib.sha256(rnr1.E12.astype("<i8").tobytes()).hexdigest()[:16],
        )
    ).encode("ascii")
    return hashlib.sha256(desc).digest()[:8]


@contextlib.contextmanager
def _hashed_patched(HBITS, fixed_weights=None):
    """Swap the hashed predictor + its model hash into the rnr1 seam."""
    global _HBITS, _FIXED_WEIGHTS
    saved = (rnr1.NGramPredictor, rnr1.model_description_hash,
             _HBITS, _FIXED_WEIGHTS)
    _HBITS = int(HBITS)
    _FIXED_WEIGHTS = fixed_weights
    rnr1.NGramPredictor = _predictor_factory
    rnr1.model_description_hash = hashed_model_description_hash
    try:
        yield
    finally:
        (rnr1.NGramPredictor, rnr1.model_description_hash,
         _HBITS, _FIXED_WEIGHTS) = saved


def _pack_container(data, W, K):
    """rnr1.pack, verbatim, minus the W in [1,8] guard (W up to W_MAX here)."""
    if W < 1 or W > W_MAX:
        raise ValueError("W must be in [1, %d]" % W_MAX)
    if K < 256:
        raise ValueError("K must be >= 256")
    n = len(data)
    m = (n + K - 1) // K
    blobs, entries = [], []
    repair_bit_offset = 0
    for j in range(m):
        block = data[j * K: (j + 1) * K]
        blob = rnr1.encode_subblock(block, W)
        mode = rnr1.MODE_CODED
        if len(blob) >= len(block):
            blob, mode = block, rnr1.MODE_RAW
        entries.append(struct.pack(
            rnr1.INDEX_ENTRY_FMT, j * K, repair_bit_offset, 0, mode,
            rnr1.h_v(block)))
        blobs.append(blob)
        repair_bit_offset += 8 * len(blob)
    header = struct.pack(
        rnr1.HEADER_FMT, rnr1.MAGIC, bytes(rnr1.VERSION), 0,
        rnr1.FLAG_SUBBLOCK_HASHES, W, K, n,
        rnr1.model_description_hash(W), rnr1.h_v(data), m)
    return header + b"".join(entries) + b"".join(blobs)


def pack_hashed(data, W=12, HBITS=22, K=None):
    """Encode into an RNR1 archive with the hashed CM predictor.

    K defaults to a single warm sub-block (K = len(data)); pass an explicit
    K for the reset-per-sync random-access variant.
    """
    if K is None:
        K = max(256, len(data))
    with _hashed_patched(HBITS):
        return _pack_container(data, W=W, K=K)


def unpack_hashed(raw, HBITS=22, verify=True):
    with _hashed_patched(HBITS):
        return rnr1.unpack(raw, verify=verify)


def read_hashed(raw, p, k, HBITS=22, verify=False):
    with _hashed_patched(HBITS):
        return rnr1.Archive(raw).read(p, k, verify=verify)


def info_hashed(raw, HBITS=22):
    with _hashed_patched(HBITS):
        arc = rnr1.Archive(raw)
        return {"n": arc.n, "W": arc.W, "K": arc.K, "m": arc.m,
                "archive_bytes": len(raw), "repair_bytes": len(arc.repair),
                "bpb": 8.0 * len(raw) / max(1, arc.n)}


# ---------------------------------------------------------------------------
# Multi-pass (semi-static) variant: epoch-train the mixer weights, ship them
# ---------------------------------------------------------------------------
#
# A memory-bounded, self-decoding codec cannot ship its hashed COUNTS (they
# are 2^HBITS buckets) and cannot look ahead, so the counts stay causal and
# single-pass (encoder and decoder rebuild them identically left-to-right).
# The MIXER WEIGHTS, by contrast, are only (W+1) x 256 integers -- small
# enough to ship.  The multi-pass variant epoch-trains the weights over the
# data P times (counts re-warmed each epoch, weights carried over) so they
# converge to this file's statistics, then codes ONE causal pass with the
# converged weights FROZEN.  The decoder loads the shipped weights and
# replays the same causal count adaptation, so the archive is exactly
# reproducible.  The weight blob is tiny (e.g. 17*256*4 = 17 KB), i.e.
# ~2e-6 bpb on a 64 MB file and ~0.13 bpb on a 1 MB file -- accounted for in
# the reported archive size, never hidden.

_H2_MAGIC = b"RNRH"
_H2_FMT = "<4s H H H H I"          # magic, W, HBITS, passes, reserved, wlen
_H2_SIZE = struct.calcsize(_H2_FMT)


def train_weights(data, W, HBITS, passes, log=None):
    """Epoch-train the mixer weights; return (weights int64[(W+1),256],
    per_pass_time_s list).  Counts re-warm from empty each epoch."""
    import time
    pred = HashedCMPredictor(W, HBITS=HBITS)
    per_pass = []
    for ep in range(passes):
        t0 = time.perf_counter()
        pred.reset_counts()                 # keep weights, re-warm counts
        for x in data:
            pred.dist()
            pred.update(x)
        dt = time.perf_counter() - t0
        per_pass.append(dt)
        if log:
            log("    epoch %d/%d  %.1fs" % (ep + 1, passes, dt))
    return pred.weights.copy(), per_pass


def pack_hashed_mp(data, W=12, HBITS=22, passes=1, K=None, log=None):
    """Multi-pass encode.  passes==1 is the plain online archive (no wrapper,
    nothing shipped); passes>1 ships the epoch-trained frozen weights.

    Returns (archive_bytes, meta) where meta has per_pass_time_s / n_passes.
    """
    import time
    if passes <= 1:
        t0 = time.perf_counter()
        raw = pack_hashed(data, W=W, HBITS=HBITS, K=K)
        return raw, {"n_passes": 1, "per_pass_time_s": [],
                     "train_time_s": 0.0, "encode_time_s": time.perf_counter() - t0}
    w, per_pass = train_weights(data, W, HBITS, passes, log=log)
    t0 = time.perf_counter()
    if K is None:
        K = max(256, len(data))
    with _hashed_patched(HBITS, fixed_weights=w):
        inner = _pack_container(data, W=W, K=K)
    wblob = w.astype("<i4").tobytes()
    head = struct.pack(_H2_FMT, _H2_MAGIC, W, HBITS, passes, 0, len(wblob))
    raw = head + wblob + inner
    enc_dt = time.perf_counter() - t0
    return raw, {"n_passes": passes, "per_pass_time_s": per_pass,
                 "train_time_s": float(sum(per_pass)), "encode_time_s": enc_dt}


def unpack_hashed_mp(raw, verify=True):
    """Decode either a plain online archive (RNR1) or a multi-pass wrapper
    (RNRH: shipped weights + inner RNR1)."""
    if raw[:4] == _H2_MAGIC:
        magic, W, HBITS, passes, _r, wlen = struct.unpack_from(_H2_FMT, raw, 0)
        wblob = raw[_H2_SIZE:_H2_SIZE + wlen]
        w = np.frombuffer(wblob, dtype="<i4").astype(np.int64).reshape(W + 1, 256)
        inner = raw[_H2_SIZE + wlen:]
        with _hashed_patched(HBITS, fixed_weights=w):
            return rnr1.unpack(inner, verify=verify)
    # plain online archive: HBITS is not stored; the model-hash check will
    # reject a wrong HBITS.  Recover HBITS by trying the header's W's hash.
    return _unpack_online_autohbits(raw, verify=verify)


def _unpack_online_autohbits(raw, verify=True, hbits_candidates=range(12, 27)):
    """Online (passes==1) archives don't carry HBITS; find the HBITS whose
    model-description hash matches the header (cheap: a few sha256s)."""
    W = struct.unpack_from(rnr1.HEADER_FMT, raw, 0)[4]
    want = struct.unpack_from(rnr1.HEADER_FMT, raw, 0)[7]
    global _HBITS
    for hb in hbits_candidates:
        _HBITS = hb
        if hashed_model_description_hash(W) == want:
            return unpack_hashed(raw, HBITS=hb, verify=verify)
    raise ValueError("no HBITS in %s matches archive model hash"
                     % list(hbits_candidates))


def _cli(argv=None):
    ap = argparse.ArgumentParser(prog="hashed_coder")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack")
    p.add_argument("input"); p.add_argument("output")
    p.add_argument("--W", type=int, default=12)
    p.add_argument("--HBITS", type=int, default=22)
    p.add_argument("--K", type=int, default=0, help="0 = warm whole-block")
    p = sub.add_parser("unpack")
    p.add_argument("input"); p.add_argument("output")
    p.add_argument("--HBITS", type=int, default=22)
    p.add_argument("--no-verify", action="store_true")
    p = sub.add_parser("info")
    p.add_argument("archive")
    p.add_argument("--HBITS", type=int, default=22)
    args = ap.parse_args(argv)

    if args.cmd == "pack":
        data = Path(args.input).read_bytes()
        K = None if args.K <= 0 else args.K
        raw = pack_hashed(data, W=args.W, HBITS=args.HBITS, K=K)
        Path(args.output).write_bytes(raw)
        print("packed %d -> %d bytes (%.4f bpb) W=%d HBITS=%d"
              % (len(data), len(raw), 8.0 * len(raw) / max(1, len(data)),
                 args.W, args.HBITS))
    elif args.cmd == "unpack":
        raw = Path(args.input).read_bytes()
        data = unpack_hashed(raw, HBITS=args.HBITS, verify=not args.no_verify)
        Path(args.output).write_bytes(data)
        print("unpacked %d -> %d bytes" % (len(raw), len(data)))
    else:
        st = info_hashed(Path(args.archive).read_bytes(), HBITS=args.HBITS)
        for k in ("n", "W", "K", "m", "archive_bytes", "repair_bytes", "bpb"):
            print("  %-14s %s" % (k, st[k]))


if __name__ == "__main__":
    _cli()
