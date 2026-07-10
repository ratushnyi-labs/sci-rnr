#!/usr/bin/env python3
"""H2/H3/H4 Type-II separation campaign (exp-design sections 2.2-2.3).

Per corpus (see corpora.py for the five record/field decompositions):

  (a) factorized-vs-joint coding gap with the reference Type-I coder at
      matched W and per-block cold starts, with paired per-block
      bootstrap CIs (frozen methodology: BCa, B=10000, METHODS.md) and
      gzip/lzma/zstd baselines on the same byte streams as controls;
  (b) plug-in multi-information TC(D_N) on the same field decomposition
      (Miller-Madow + multinomial bootstrap; analytic TC for the
      synthetic generators) and the gap-vs-TC comparison of Theorem 5.5;
  (c) H4 permutation control: per-field independent permutation across
      records (marginals preserved exactly, joint structure destroyed),
      under which TC -> 0 and the gap must vanish.

Rates are deterministic given the derived seeds (no unlogged
randomness); run-to-run spread comes from the seeded record shuffle /
permutation (seed indices 0..n_seeds-1, tag "run" per METHODS.md).

Outputs (all under this experiment directory):
  results/h2_h4_typeii_<scale>.jsonl   bench/results.py-compatible store
  results/summary_<scale>.json         aggregated numbers for the gate
  logs/run_<scale>.log                 full stdout copy
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in ("bench", "data", "impl"):
    sys.path.insert(0, str(REPO / p))
sys.path.insert(0, str(HERE))

import baselines  # noqa: E402
import coding  # noqa: E402
import corpora  # noqa: E402
import metrics  # noqa: E402
import results as store  # noqa: E402
import tc as tcmod  # noqa: E402

W = 3               # matched context order for BOTH arms (coder default)
N_BLOCKS = 16       # per-corpus block partition for paired CIs
BASELINE_CODERS = ("gzip", "lzma", "zstd")

SCALES = {
    "reduced": {"M": corpora.REDUCED_M, "n_boot": 1000, "n_seeds": 1},
    "full": {"M": corpora.FULL_M, "n_boot": 10000, "n_seeds": 3},
}


class Tee:
    def __init__(self, path: Path):
        self.f = open(path, "w", encoding="utf-8")
        self.stdout = sys.stdout

    def write(self, s):
        self.f.write(s)
        self.stdout.write(s)

    def flush(self):
        self.f.flush()
        self.stdout.flush()


def run_corpus(name: str, scale: str, cfg: dict, jsonl: Path) -> dict:
    t_corpus = time.perf_counter()
    M_arg = cfg["M"][name]
    built = corpora.BUILDERS[name]() if M_arg is None \
        else corpora.BUILDERS[name](M_arg)
    records = built["records"]
    M, R = records.shape
    nbytes = M * R
    print(f"\n=== {name}: M={M} R={R} ({nbytes} bytes/arm, scale={scale}) ===")

    exp_id = f"{corpora.EXPERIMENT_ID}:{name}"
    common_cfg = {"W": W, "R": R, "M": M, "n_blocks": N_BLOCKS,
                  "scale": scale, "fields": built["field_names"]}
    blk_bytes = coding.block_sizes_bytes(M, R, N_BLOCKS)

    # ---- (b) TC estimates: structured (shuffle-invariant, computed once)
    boot_seed = metrics.derive_seed(exp_id, 0, tag="boot")
    tc_struct = tcmod.estimate_tc(records, n_boot=cfg["n_boot"],
                                  boot_seed=boot_seed)
    cmi_struct = tcmod.chain_mi(records)
    print(f"  TC(structured): {tc_struct['tc_bpb']:.4f} bpb "
          f"[{tc_struct['tc_boot_lo_bpb']:.4f}, "
          f"{tc_struct['tc_boot_hi_bpb']:.4f}] "
          f"(distinct={tc_struct['distinct_records']}, "
          f"coverage={tc_struct['coverage']:.3f}, "
          f"estimable={tc_struct['estimable']})")
    print(f"  chain-MI(structured) [rigorous TC lower bound]: "
          f"{cmi_struct['chain_mi_bpb']:.4f} bpb")
    if built["analytic_tc_bpb"] is not None:
        print(f"  TC(analytic model): {built['analytic_tc_bpb']:.4f} bpb")
    store.append(jsonl, store.make_record(
        name, "tc-plugin", {**common_cfg, "variant": "structured"}, 0,
        {"tc_bpb": tc_struct["tc_bpb"],
         "tc_lo": tc_struct["tc_boot_lo_bpb"],
         "tc_hi": tc_struct["tc_boot_hi_bpb"],
         "chain_mi_bpb": cmi_struct["chain_mi_bpb"]},
        corpus_sha256=built["corpus_sha256"], experiment=exp_id,
        tc_detail=tc_struct, chain_mi_detail=cmi_struct,
        analytic_tc_bpb=built["analytic_tc_bpb"]))

    prefix_curve = None
    if R >= 8:
        prefixes = sorted({max(2, R // 4), R // 2, 3 * R // 4, R})
        prefix_curve = tcmod.prefix_tc_curve(records, prefixes)
        pts = ", ".join(f"r={p['prefix_bytes']}: {p['tc_bits']:.2f}b"
                        for p in prefix_curve)
        print(f"  TC prefix curve (Theta(N) proxy): {pts}")

    # ---- (a)+(c) coding arms per seed and variant
    per_seed = []
    tc_perm_list = []
    for si in range(cfg["n_seeds"]):
        struct = corpora.shuffle_records(records, name, si)
        perm = corpora.permute_fields(struct, name, si)
        seed_val = metrics.derive_seed(exp_id, si, tag="run")
        row = {"seed_index": si}
        for variant, recs in (("structured", struct), ("permuted", perm)):
            t0 = time.perf_counter()
            jb = coding.blockwise_bits(coding.joint_streams(recs),
                                       N_BLOCKS, W)
            fb = coding.blockwise_bits(coding.factorized_streams(recs),
                                       N_BLOCKS, W)
            enc_s = time.perf_counter() - t0
            jr, fr = jb / blk_bytes, fb / blk_bytes  # per-block bpb
            gap = fr - jr
            ci = metrics.paired_bootstrap_diff(
                fr, jr, n_resamples=cfg["n_boot"],
                seed=metrics.derive_seed(f"{exp_id}:{variant}", si,
                                         tag="boot"))
            d_eff = metrics.cohens_d_paired(fr, jr)
            row[variant] = {
                "joint_bpb": float(jb.sum() / nbytes),
                "fact_bpb": float(fb.sum() / nbytes),
                "gap_bpb": float(gap.mean()),
                "gap_lo": ci.lo, "gap_hi": ci.hi, "gap_ci_method": ci.method,
                "cohens_d": d_eff,
                "joint_blocks_bpb": jr.tolist(),
                "fact_blocks_bpb": fr.tolist(),
            }
            print(f"  seed {si} {variant:>10}: joint {row[variant]['joint_bpb']:.4f} "
                  f"fact {row[variant]['fact_bpb']:.4f} "
                  f"gap {row[variant]['gap_bpb']:.4f} "
                  f"[{ci.lo:.4f}, {ci.hi:.4f}] bpb ({enc_s:.1f}s)")
            for arm, arr, tot in (("joint", jr, jb), ("factorized", fr, fb)):
                store.append(jsonl, store.make_record(
                    name, f"rnr1-replay-{arm}",
                    {**common_cfg, "variant": variant}, seed_val,
                    {"bpb": float(tot.sum() / nbytes), "enc_s": enc_s / 2},
                    corpus_sha256=built["corpus_sha256"], experiment=exp_id,
                    blocks_bpb=arr.tolist(), rate_deterministic=True))
        tcp = tcmod.estimate_tc(
            perm, n_boot=max(1000, cfg["n_boot"] // 10),
            boot_seed=metrics.derive_seed(f"{exp_id}:perm", si, tag="boot"))
        cmip = tcmod.chain_mi(perm)
        tc_perm_list.append({"tc_bpb": tcp["tc_bpb"],
                             "estimable": tcp["estimable"],
                             "coverage": tcp["coverage"],
                             "chain_mi_bpb": cmip["chain_mi_bpb"]})
        store.append(jsonl, store.make_record(
            name, "tc-plugin", {**common_cfg, "variant": "permuted"},
            seed_val, {"tc_bpb": tcp["tc_bpb"],
                       "chain_mi_bpb": cmip["chain_mi_bpb"]},
            corpus_sha256=built["corpus_sha256"], experiment=exp_id,
            tc_detail={k: tcp[k] for k in
                       ("distinct_records", "coverage", "estimable")}))
        per_seed.append(row)
    perm_est = [t["tc_bpb"] for t in tc_perm_list if t["estimable"]]
    print(f"  permuted control (true TC = 0): chain-MI "
          f"{np.mean([t['chain_mi_bpb'] for t in tc_perm_list]):.4f} bpb; "
          f"full TC "
          + (f"{np.mean(perm_est):.4f} bpb" if perm_est
             else "not estimable (support ~ sample size)"))

    # ---- natural-order joint rate (secondary; cross-record structure)
    jb_nat = coding.blockwise_bits(coding.joint_streams(records), N_BLOCKS, W)
    nat_bpb = float(jb_nat.sum() / nbytes)
    print(f"  natural-order joint (secondary): {nat_bpb:.4f} bpb")

    # ---- replay-vs-pack cross-validation (structured joint, seed 0)
    val = coding.validate_replay(corpora.shuffle_records(records, name, 0), W)
    print(f"  replay vs pack: rel_err={val['rel_err']:.4%} "
          f"raw_blocks={val['raw_mode_blocks']} "
          f"roundtrip={'ok' if val['roundtrip_ok'] else 'FAIL'}")

    # ---- baselines as controls (structured joint / planar / permuted)
    struct0 = corpora.shuffle_records(records, name, 0)
    perm0 = corpora.permute_fields(struct0, name, 0)
    streams = {
        "joint": np.ascontiguousarray(struct0).tobytes(),
        "planar": b"".join(np.ascontiguousarray(struct0[:, f]).tobytes()
                           for f in range(R)),
        "permuted-joint": np.ascontiguousarray(perm0).tobytes(),
    }
    base = {}
    for coder in BASELINE_CODERS:
        if coder not in baselines.available_codecs():
            continue
        base[coder] = {}
        for layout, data in streams.items():
            r = baselines.run_baseline(coder, data)
            b = 8.0 * r.compressed_size / len(data)
            base[coder][layout] = b
            store.append(jsonl, store.make_record(
                name, coder, {**common_cfg, "layout": layout}, 0,
                {"bpb": b, "compressed_size": r.compressed_size,
                 "original_size": len(data), "enc_s": r.enc_s,
                 "dec_s": r.dec_s},
                corpus_sha256=built["corpus_sha256"], experiment=exp_id,
                rss_mode="in-process"))
        print(f"  baseline {coder}: " + " ".join(
            f"{k}={v:.3f}" for k, v in base[coder].items()))

    return {
        "corpus": name, "M": M, "R": R, "W": W,
        "analytic_tc_bpb": built["analytic_tc_bpb"],
        "tc_struct": {k: tc_struct[k] for k in
                      ("tc_bpb", "tc_boot_lo_bpb", "tc_boot_hi_bpb",
                       "tc_mm_bpb", "tc_plug_bpb",
                       "distinct_records", "coverage", "estimable")},
        "h_joint_grass_bpb":
            tc_struct["estimators_bits"]["grass"]["h_joint_bits"] / R,
        "h_marg_grass_bpb":
            tc_struct["estimators_bits"]["grass"]["h_marg_sum_bits"] / R,
        "chain_mi_struct_bpb": cmi_struct["chain_mi_bpb"],
        "tc_perm": tc_perm_list,
        "prefix_tc_curve": prefix_curve,
        "per_seed": per_seed,
        "natural_order_joint_bpb": nat_bpb,
        "replay_validation": val,
        "baselines_bpb": base,
        "elapsed_s": time.perf_counter() - t_corpus,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scale", choices=("reduced", "full"), default="reduced")
    ap.add_argument("--corpus", action="append",
                    help="restrict to named corpora (repeatable)")
    args = ap.parse_args()
    cfg = SCALES[args.scale]

    log_path = HERE / "logs" / f"run_{args.scale}.log"
    sys.stdout = Tee(log_path)
    jsonl = HERE / "results" / f"h2_h4_typeii_{args.scale}.jsonl"
    if jsonl.exists():
        jsonl.unlink()  # deterministic re-run replaces the derived store

    names = args.corpus or list(corpora.BUILDERS)
    t0 = time.perf_counter()
    print(f"H2/H3/H4 Type-II separation campaign -- scale={args.scale}, "
          f"W={W}, n_blocks={N_BLOCKS}, n_boot={cfg['n_boot']}, "
          f"n_seeds={cfg['n_seeds']}")
    summary = {"scale": args.scale, "W": W, "n_blocks": N_BLOCKS,
               "n_boot": cfg["n_boot"], "n_seeds": cfg["n_seeds"],
               "corpora": {}}
    for name in names:
        summary["corpora"][name] = run_corpus(name, args.scale, cfg, jsonl)
    summary["total_elapsed_s"] = time.perf_counter() - t0

    out = HERE / "results" / f"summary_{args.scale}.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1, sort_keys=True)
    print(f"\nwrote {out} and {jsonl} in {summary['total_elapsed_s']:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
