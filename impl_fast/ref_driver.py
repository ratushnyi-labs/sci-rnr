#!/usr/bin/env python3
"""Parallel driver around the NORMATIVE reference coder impl/rnr1.py.

The reference pack() encodes sub-blocks strictly sequentially; each sub-block
encode is a pure function of (block bytes, W) with fresh predictor and coder
state (impl/rnr1.py encode_subblock).  This driver farms the per-sub-block
calls of rnr1.encode_subblock / rnr1.decode_subblock out to a process pool
and assembles the container with the very same code paths (rnr1.h_v,
rnr1.model_description_hash, rnr1 struct formats), so its output is the same
bytes rnr1.pack() would produce -- this is additionally verified at runtime
by run_checks.py (driver output == rnr1.pack output on multi-megabyte
inputs, check F1) before any parallel result is trusted.

The reference implementation itself is NOT modified and stays normative.

Usage:
  ref_driver.py pack   <input> <output> --W w --K k [--jobs J]
  ref_driver.py unpack <archive> <output> [--jobs J]   # full verify
"""

import argparse
import os
import struct
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "impl"))
import rnr1  # noqa: E402


def _enc_one(args):
    j, block, W = args
    return j, rnr1.encode_subblock(block, W)


def parallel_pack(data, W, K, jobs):
    """Byte-identical reimplementation of rnr1.pack with a process pool."""
    if W < 1 or W > 8:
        raise ValueError("W must be in [1, 8]")
    if K < 256:
        raise ValueError("K must be >= 256")
    n = len(data)
    m = (n + K - 1) // K
    tasks = ((j, data[j * K: (j + 1) * K], W) for j in range(m))
    blobs = [None] * m
    if m:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for j, blob in ex.map(_enc_one, tasks, chunksize=4):
                blobs[j] = blob
    entries = []
    repair_bit_offset = 0
    out_blobs = []
    for j in range(m):
        block = data[j * K: (j + 1) * K]
        blob = blobs[j]
        mode = rnr1.MODE_CODED
        if len(blob) >= len(block):
            blob = block
            mode = rnr1.MODE_RAW
        entries.append(struct.pack(
            rnr1.INDEX_ENTRY_FMT, j * K, repair_bit_offset, 0, mode,
            rnr1.h_v(block)))
        out_blobs.append(blob)
        repair_bit_offset += 8 * len(blob)
    header = struct.pack(
        rnr1.HEADER_FMT, rnr1.MAGIC, bytes(rnr1.VERSION), 0,
        rnr1.FLAG_SUBBLOCK_HASHES, W, K, n,
        rnr1.model_description_hash(W), rnr1.h_v(data), m)
    return header + b"".join(entries) + b"".join(out_blobs)


def _dec_one(args):
    j, blob, nbytes, W, mode, want_hash = args
    if mode == rnr1.MODE_RAW:
        out = blob[:nbytes]
    else:
        out = rnr1.decode_subblock(blob, nbytes, W)
    if rnr1.h_v(out) != want_hash:
        raise ValueError("sub-block %d hash mismatch" % j)
    return j, out


def parallel_unpack(raw, jobs):
    """Byte-identical reimplementation of rnr1.unpack (verify=True) with a
    process pool; all header/index validation is rnr1.Archive's own."""
    arc = rnr1.Archive(raw)
    tasks = ((j, arc._block_blob(j), arc._block_len(j), arc.W,
              arc.entries[j][3], arc.entries[j][4]) for j in range(arc.m))
    outs = [None] * arc.m
    if arc.m:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for j, out in ex.map(_dec_one, tasks, chunksize=4):
                outs[j] = out
    data = b"".join(outs)
    if rnr1.h_v(data) != arc.v:
        raise ValueError("archive verification failed: H_v mismatch")
    return data


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack")
    p.add_argument("input"); p.add_argument("output")
    p.add_argument("--W", type=int, required=True)
    p.add_argument("--K", type=int, default=65536)
    p.add_argument("--jobs", type=int, default=os.cpu_count())
    p = sub.add_parser("unpack")
    p.add_argument("input"); p.add_argument("output")
    p.add_argument("--jobs", type=int, default=os.cpu_count())
    args = ap.parse_args()
    with open(args.input, "rb") as f:
        data = f.read()
    if args.cmd == "pack":
        out = parallel_pack(data, args.W, args.K, args.jobs)
    else:
        out = parallel_unpack(data, args.jobs)
    with open(args.output, "wb") as f:
        f.write(out)
    print("%s %s: %d -> %d bytes" % (args.cmd, args.input, len(data), len(out)))


if __name__ == "__main__":
    main()
