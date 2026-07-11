#!/usr/bin/env python3
"""Context-order (W) and block-size study for the two-pass coder.

Motivating question: a fixed order-W = 3 context is weak; the best order
depends on the data and on how much of it a block sees.  Two-pass coding can
choose W per block in pass 1 (measure the frozen-model coded size at several
orders, keep the smallest, record it in the header).  This script measures:

  (1) rate vs context order W, per corpus -- the frozen-model coded rate
      (pure prediction quality, no sync overhead) plus the model-section
      share, so the total and the best W are visible;
  (2) model-section compression: raw vs zlib/lzma, and the resulting total;
  (3) block-size (outer B) effect via the model share = model_bits / B.

Rate is bit-deterministic.  Coded size per W is measured with
encode_subblock_frozen over the whole input as one block (K = whole file),
i.e. the no-sync two-pass limit -- isolating what the order buys, since the
two-pass rate is already nearly K-independent (measured separately).

Usage: python block_study.py [--text-mb 4] [--out out/block_study.csv]
"""
from __future__ import annotations

import argparse
import bz2
import csv
import lzma
import sys
import time
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import two_pass  # noqa: E402

MIB = 1024 * 1024
WS = [2, 3, 4, 5, 6]
OUTER_B = [1 * MIB, 4 * MIB, 16 * MIB, 64 * MIB, 256 * MIB]


def load(name: str, nbytes: int) -> bytes:
    src = {
        "text": REPO / "data" / "payloads" / "enwik8",
        "code": REPO / "data" / "payloads" / "rnr_scripts_src.tar",
    }[name]
    b = src.read_bytes()
    return b[:nbytes] if nbytes else b


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--text-mb", type=int, default=4)
    ap.add_argument("--out", default=str(HERE / "out" / "block_study.csv"))
    args = ap.parse_args(argv)

    rows = []
    for corpus in ("text", "code"):
        data = load(corpus, args.text_mb * 1_000_000 if corpus == "text" else 0)
        n = len(data)
        print(f"\n=== {corpus}  ({n/1e6:.2f} MB) ===")
        print(f"{'W':>3} {'stream bpb':>11} {'model KB':>9} {'model bpb':>10} "
              f"{'total bpb':>10} {'total+mz bpb':>13} {'t(s)':>6}")
        best = None
        for W in WS:
            t0 = time.time()
            model, blob, stats = two_pass.build_frozen_model(data, W)
            coded = two_pass.encode_subblock_frozen(data, model)
            dt = time.time() - t0
            stream_bits = 8 * len(coded)
            m_raw = len(blob)
            m_zlib = len(zlib.compress(blob, 9))
            m_lzma = len(lzma.compress(blob, preset=9))
            m_best = min(m_zlib, m_lzma, len(bz2.compress(blob, 9)))
            stream_bpb = stream_bits / n
            model_bpb = 8 * m_raw / n
            total_bpb = (stream_bits + 8 * m_raw) / n
            total_mz = (stream_bits + 8 * m_zlib) / n
            print(f"{W:>3} {stream_bpb:>11.4f} {m_raw/1024:>9.1f} "
                  f"{model_bpb:>10.4f} {total_bpb:>10.4f} {total_mz:>13.4f} "
                  f"{dt:>6.1f}")
            rows.append({
                "corpus": corpus, "n_bytes": n, "W": W,
                "stream_bytes": len(coded), "stream_bpb": round(stream_bpb, 5),
                "model_raw_bytes": m_raw, "model_zlib_bytes": m_zlib,
                "model_lzma_bytes": m_lzma, "model_best_bytes": m_best,
                "model_bpb": round(model_bpb, 5),
                "total_bpb": round(total_bpb, 5),
                "total_zlibmodel_bpb": round(total_mz, 5),
                "encode_s": round(dt, 2),
            })
            if best is None or total_mz < best[1]:
                best = (W, total_mz)
        wbest, bbest = best
        w3 = next(r for r in rows
                  if r["corpus"] == corpus and r["W"] == 3)
        gain = w3["total_zlibmodel_bpb"] - bbest
        print(f"  best W = {wbest} (total+mz {bbest:.4f} bpb); "
              f"vs fixed W=3 {w3['total_zlibmodel_bpb']:.4f} -> "
              f"gain {gain:+.4f} bpb by letting W depend on the block")
        # model share vs outer block B, at the best W's model
        mb = next(r for r in rows
                  if r["corpus"] == corpus and r["W"] == wbest)
        print(f"{'outer B':>10} {'model share raw':>16} {'model share best':>17}")
        for B in OUTER_B:
            sr = 8 * mb["model_raw_bytes"] / B
            sb = 8 * mb["model_best_bytes"] / B
            print(f"{B//MIB:>9}M {sr:>16.4f} {sb:>17.4f}")
            rows.append({
                "corpus": corpus, "n_bytes": n, "W": wbest,
                "outer_block_B": B,
                "model_share_raw_bpb": round(sr, 5),
                "model_share_best_bpb": round(sb, 5),
                "model_raw_bytes": mb["model_raw_bytes"],
                "model_best_bytes": mb["model_best_bytes"],
            })

    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    cols = sorted({k for r in rows for k in r})
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"\n[block_study] {len(rows)} rows -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
