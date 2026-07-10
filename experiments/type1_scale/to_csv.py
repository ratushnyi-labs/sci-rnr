#!/usr/bin/env python3
"""Flatten all campaign results.jsonl files into one results.csv.

The JSONL files under out/<tier>/results.jsonl are the source of truth
(append-only, written by run_campaign.py).  This exporter produces the
flat CSV used for plotting: one row per completed cell run, stable
column order, deterministic row order (tier, type, coder, run_idx).

Usage:
    python to_csv.py [--out out/results.csv] [--tiers 500,1500,5000,...]
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"

COLUMNS = [
    "tier", "size_bytes", "size_mb", "content_type", "coder", "family",
    "coder_version", "run_idx", "compressed_bytes", "bpb", "ratio",
    "enc_s", "dec_s", "enc_mb_s", "dec_mb_s",
    "enc_peak_rss_mb", "dec_peak_rss_mb",
    "roundtrip_ok", "censored", "status", "threads",
    "corpus_sha256", "seed", "timestamp",
]

# smoke tiers are machinery calibration, never benchmark evidence
DEFAULT_TIERS = ["500", "1500", "5000"]


def rows_from(path: Path):
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("status") not in ("done", "censored"):
                continue
            m = r.get("measured", {})
            osize = m.get("original_size")
            enc_s, dec_s = m.get("enc_s"), m.get("dec_s")
            cfg = r.get("config", {}) or {}
            yield {
                "tier": r.get("tier"),
                "size_bytes": osize,
                "size_mb": round(osize / 1e6, 3) if osize else None,
                "content_type": r.get("ctype"),
                "coder": r.get("coder"),
                "family": cfg.get("impl") or r.get("coder", "").split("-")[0],
                "coder_version": r.get("coder_version"),
                "run_idx": r.get("run_idx"),
                "compressed_bytes": m.get("compressed_size"),
                "bpb": m.get("bpb"),
                "ratio": m.get("ratio"),
                "enc_s": enc_s,
                "dec_s": dec_s,
                "enc_mb_s": (round(osize / 1e6 / enc_s, 3)
                             if osize and enc_s else None),
                "dec_mb_s": (round(osize / 1e6 / dec_s, 3)
                             if osize and dec_s else None),
                "enc_peak_rss_mb": (round(m["enc_peak_rss"] / 1e6, 1)
                                    if m.get("enc_peak_rss") else None),
                "dec_peak_rss_mb": (round(m["dec_peak_rss"] / 1e6, 1)
                                    if m.get("dec_peak_rss") else None),
                "roundtrip_ok": m.get("roundtrip_ok"),
                "censored": r.get("censored", False),
                "status": r.get("status"),
                "threads": cfg.get("threads", 1),
                "corpus_sha256": r.get("corpus_sha256"),
                "seed": r.get("seed"),
                "timestamp": r.get("timestamp"),
            }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT / "results.csv"))
    ap.add_argument("--tiers", default=",".join(DEFAULT_TIERS),
                    help="comma-separated tier dirs under out/ "
                         "(default: the paper tiers; add smoke explicitly "
                         "if you want calibration rows)")
    args = ap.parse_args(argv)

    tiers = [t.strip() for t in args.tiers.split(",") if t.strip()]
    rows = []
    seen = []
    for tier in tiers:
        p = OUT / tier / "results.jsonl"
        if p.exists():
            rows.extend(rows_from(p))
            seen.append(tier)
    rows.sort(key=lambda r: (str(r["tier"]), str(r["content_type"]),
                             str(r["coder"]), int(r["run_idx"] or 0)))
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"[to_csv] {len(rows)} rows from tiers {seen} -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
