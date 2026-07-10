#!/usr/bin/env python3
"""Campaign driver for the H5/H6/H7 Type-III experiments.

Usage:
    python run_experiment.py [--scale reduced|full]

Writes out/summary_<scale>.json (full detail, minus bulky per-symbol
arrays) and out/results_<scale>.jsonl (bench-compatible append-only
records via bench/results.py).  All randomness derives from the frozen
master seed (bench/METHODS.md section 4).
"""

from __future__ import annotations

import argparse
import json
import time

import common
import h5_modes
import h6_cmi
import h7_dict
import results

SCALES = {
    "reduced": {
        "h5_compositions": 1, "h5_tiles": 32, "h5_tile_bytes": 8192,
        "h6_slice": 2 * 1024 * 1024, "h6_cm_slice": 65536, "h6_boot": 200,
        "h7_ladder": {"enwik8": [32768, 65536, 131072, 262144],
                      "scripts_src": [65536, 131072]},
        "h7_bitsback_slice": 16384,
    },
    "full": {
        "h5_compositions": 2, "h5_tiles": 64, "h5_tile_bytes": 16384,
        "h6_slice": 8 * 1024 * 1024, "h6_cm_slice": 262144, "h6_boot": 1000,
        "h7_ladder": {"enwik8": [65536, 131072, 262144, 524288, 1048576],
                      "scripts_src": [131072, 262144, 524288]},
        "h7_bitsback_slice": 32768,
    },
}


def run(scale: str) -> dict:
    cfg = SCALES[scale]
    t0 = time.time()
    jsonl = common.OUT / f"results_{scale}.jsonl"
    if jsonl.exists():
        jsonl.unlink()      # regenerate this campaign's store from scratch

    # ---------------- H5 ----------------
    h5 = []
    for i in range(cfg["h5_compositions"]):
        r = h5_modes.run_h5(i, cfg["h5_tiles"], cfg["h5_tile_bytes"])
        h5.append(r)
        corpus = f"hetero_composite_s{i}"
        base_cfg = {"tile_bytes": r["tile_bytes"], "tiles": r["tiles"],
                    "scale": scale}
        for coder, bits in (
            [("oracle_pertile", r["oracle_bits"]),
             ("modesel_prefix", r["selector_prefix"]["bits"]),
             ("modesel_causal", r["selector_causal"]["bits"])]
            + [("fixed_" + m, b) for m, b in r["fixed_mode_bits"].items()]
        ):
            rec = results.make_record(
                corpus, coder, base_cfg, r["seed"],
                {"bpb": bits / r["n_bytes"],
                 "compressed_size": bits // 8,
                 "original_size": r["n_bytes"]},
                hypothesis="H5", composition=r["segments"],
                label="secondary")
            results.append(jsonl, rec)

    # ---------------- H6 ----------------
    h6 = h6_cmi.run_h6(cfg["h6_slice"], cfg["h6_cm_slice"], cfg["h6_boot"])
    for name, r in h6.items():
        for inst, v in r["instruments"].items():
            if inst == "per_plane":
                continue
            rec = results.make_record(
                name, "cmi_" + inst,
                {"slice_bytes": r["slice_bytes"], "scale": scale,
                 "B": cfg["h6_boot"]},
                common.derive_seed("h6-boot", 0),
                {"cmi_debiased_bpb": v["debiased"], "cmi_mm_bpb": v["mm"],
                 "cmi_raw_bpb": v["raw"], "null_mm_bpb": v["null_mm"],
                 "ci_lo": v["debiased_lo"], "ci_hi": v["debiased_hi"]},
                hypothesis="H6", corpus_sha256=r["corpus_sha256"],
                ci_method=v["method"], label="secondary")
            results.append(jsonl, rec)

    # ---------------- H7 ----------------
    h7 = h7_dict.run_h7(cfg["h7_ladder"], cfg["h7_bitsback_slice"])
    for corpus, rows in h7["ladder"].items():
        for row in rows:
            rec = results.make_record(
                corpus, "dict3c_greedy_mdl",
                {"N": row["N"], "scale": scale},
                common.derive_seed("h7", 0),
                {"bpb": row["rate_bpb"],
                 "compressed_size": int(row["rate_bpb"] * row["N"] / 8),
                 "original_size": row["N"],
                 "overhead_fraction": row["overhead_fraction"],
                 "h_hat_bpb": row["h_hat_bpb"],
                 "rate_over_Nh": row["rate_over_Nh"],
                 "K": row["K"]},
                hypothesis="H7", label="secondary")
            results.append(jsonl, rec)

    # strip bulky arrays before persisting the summary
    for r in h5:
        r.pop("cost_matrix_bits", None)
    for rows in h7["ladder"].values():
        for row in rows:
            row.pop("admission_log", None)
            neut = row.get("index_neutrality", {})
            neut.pop("symbols", None)

    summary = {
        "scale": scale, "config": cfg,
        "elapsed_s": time.time() - t0,
        "h5": h5, "h6": h6, "h7": h7,
    }
    with open(common.OUT / f"summary_{scale}.json", "w") as f:
        json.dump(summary, f, indent=1, default=str)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", choices=list(SCALES), default="reduced")
    args = ap.parse_args()
    s = run(args.scale)
    print(f"[{args.scale}] done in {s['elapsed_s']:.1f}s -> "
          f"out/summary_{args.scale}.json")


if __name__ == "__main__":
    main()
