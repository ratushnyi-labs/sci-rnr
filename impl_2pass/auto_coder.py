#!/usr/bin/env python3
"""Adaptive two-pass coder: per-block MDL order (--W auto) + model compression.

A working round-tripping coder (not a projection) that realizes the two
measured wins of the W/block study:

  * --W auto : pass 1 builds the frozen model at each candidate order W and
    keeps the one with the smallest projected TOTAL bits (stream coded size +
    compressed model section).  The chosen W is recorded in the header, so
    the order is data/block-adaptive rather than a fixed guess.
  * model compression : the serialized model section is stored zlib-compressed
    iff that is smaller (never-lose; flagged in the header).  This is exact --
    the decoder inflates it before use.

Container (self-contained lab format; correctness-first, reuses the frozen
encode/decode of two_pass.py):
    magic 'RNRA'(4) | ver(1) | W(1) | Klog2(1) | n(8) | model_mode(1:0 raw/1 zlib)
    | model_len(4) | model_blob | [ sub_len(4) | frozen_subblock ]*
Every field little-endian.  Round-trip is bit-exact; the archive-sha
determinism witness holds (integer training, deterministic prune).

Usage:
    python auto_coder.py pack   IN OUT [--W auto|N] [--K 65536] [--no-modelz]
    python auto_coder.py unpack IN OUT
    python auto_coder.py demo   [--mb 1]        # W-auto vs fixed-W=3, verified
"""
from __future__ import annotations

import argparse
import struct
import sys
import time
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import two_pass  # noqa: E402

MAGIC = b"RNRA"
VER = 1
W_CANDIDATES = (2, 3, 4, 5, 6)


def _projected_total_bits(data, W):
    """(chosen model, blob, stream_bytes, model_store_bytes, mode)."""
    model, blob, _ = two_pass.build_frozen_model(data, W)
    stream = two_pass.encode_subblock_frozen(data, model)
    z = zlib.compress(blob, 9)
    if len(z) < len(blob):
        mstore, mode = z, 1
    else:
        mstore, mode = blob, 0
    total_bits = 8 * (len(stream) + len(mstore))
    return model, blob, mstore, mode, len(stream), total_bits


def pack(data: bytes, W="auto", K: int = 65536, model_compress=True):
    n = len(data)
    if W == "auto":
        best = None
        for cand in W_CANDIDATES:
            model, blob, mstore, mode, slen, tb = _projected_total_bits(
                data, cand)
            if best is None or tb < best[0]:
                best = (tb, cand, model, blob, mstore, mode)
        _, W, model, blob, mstore, mode = best
    else:
        W = int(W)
        model, blob, _ = two_pass.build_frozen_model(data, W)
        z = zlib.compress(blob, 9)
        mstore, mode = ((z, 1) if (model_compress and len(z) < len(blob))
                        else (blob, 0))
    Klog2 = K.bit_length() - 1
    assert (1 << Klog2) == K, "K must be a power of two"
    out = bytearray()
    out += MAGIC + bytes([VER, W, Klog2])
    out += struct.pack("<Q", n)
    out += bytes([mode]) + struct.pack("<I", len(mstore)) + mstore
    for off in range(0, n, K):
        block = data[off:off + K]
        sub = two_pass.encode_subblock_frozen(block, model)
        out += struct.pack("<I", len(sub)) + sub
    return bytes(out), W, mode


def unpack(arc: bytes) -> bytes:
    assert arc[:4] == MAGIC, "bad magic"
    ver, W, Klog2 = arc[4], arc[5], arc[6]
    assert ver == VER
    K = 1 << Klog2
    n = struct.unpack_from("<Q", arc, 7)[0]
    mode = arc[15]
    mlen = struct.unpack_from("<I", arc, 16)[0]
    off = 20
    mblob = arc[off:off + mlen]
    off += mlen
    if mode == 1:
        mblob = zlib.decompress(mblob)
    model = two_pass.FrozenModel(two_pass.deserialize_model(mblob, W), W)
    out = bytearray()
    while len(out) < n:
        slen = struct.unpack_from("<I", arc, off)[0]
        off += 4
        sub = arc[off:off + slen]
        off += slen
        nblk = min(K, n - len(out))
        out += two_pass.decode_subblock_frozen(sub, nblk, model)
    return bytes(out)


def demo(mb: int) -> int:
    import hashlib
    for corpus, src in (("text", HERE.parent / "data/payloads/enwik8"),
                        ("code", HERE.parent / "data/payloads/rnr_scripts_src.tar")):
        data = src.read_bytes()
        data = data[:mb * 1_000_000] if corpus == "text" else data
        n = len(data)
        sha = hashlib.sha256(data).hexdigest()
        # fixed W=3, no model compression (baseline two-pass)
        a3, _, _ = pack(data, W=3, K=65536, model_compress=False)
        # W=auto + model compression
        t0 = time.time()
        aa, wsel, mode = pack(data, W="auto", K=65536)
        dt = time.time() - t0
        # round-trip both
        ok3 = hashlib.sha256(unpack(a3)).hexdigest() == sha
        oka = hashlib.sha256(unpack(aa)).hexdigest() == sha
        bpb3 = 8 * len(a3) / n
        bpba = 8 * len(aa) / n
        print(f"{corpus:>5} ({n/1e6:.2f} MB): "
              f"fixed W=3 {bpb3:.4f} bpb  ->  W-auto(W={wsel},"
              f"{'zlib' if mode else 'raw'} model) {bpba:.4f} bpb  "
              f"gain {bpb3-bpba:+.4f}  roundtrip {ok3 and oka}  "
              f"(select {dt:.1f}s)")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack")
    p.add_argument("inp"); p.add_argument("out")
    p.add_argument("--W", default="auto")
    p.add_argument("--K", type=int, default=65536)
    p.add_argument("--no-modelz", action="store_true")
    u = sub.add_parser("unpack")
    u.add_argument("inp"); u.add_argument("out")
    d = sub.add_parser("demo")
    d.add_argument("--mb", type=int, default=1)
    args = ap.parse_args(argv)
    if args.cmd == "pack":
        data = Path(args.inp).read_bytes()
        W = args.W if args.W == "auto" else int(args.W)
        arc, wsel, mode = pack(data, W=W, K=args.K,
                               model_compress=not args.no_modelz)
        Path(args.out).write_bytes(arc)
        print(f"packed W={wsel} model={'zlib' if mode else 'raw'} "
              f"{len(arc)} bytes ({8*len(arc)/len(data):.4f} bpb)")
    elif args.cmd == "unpack":
        Path(args.out).write_bytes(unpack(Path(args.inp).read_bytes()))
    elif args.cmd == "demo":
        return demo(args.mb)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
