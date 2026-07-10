#!/usr/bin/env python3
"""Precondition-statistics campaign over the local corpus suite.

Resolves, at dev-slice scale, the measurements that backlog items
BUG-001-D, BUG-003-D, and BUG-010-D name as external:

  BUG-001-D : the conditional-advantage ratio of the Abstract's scoped
              claim (Definition 2.1 / inequality (10.1)), measured with
              the reference Type-I coder's order-W predictor as the
              side-information instrument at W in {2,3(,4)};
  BUG-003-D : the sync-point ratio-overhead datum at typical (K, W)
              (exp-design H10(c)), on named corpora;
  BUG-010-D : which corpora satisfy the section-10.2.1 precondition --
              the qualification table (corpus x precondition ->
              qualifies / marginal / fails) and the per-hypothesis
              scoping conclusion, emitted as a MANIFEST patch file.

SCALE IS DECLARED HONESTLY: coding measurements run on the first
1 MiB (full) / 256 KiB (reduced) of each corpus, indicator statistics
on the first 8 MiB / 2 MiB. This is the local dev-slice version of the
external campaign, not the full-corpus run those items ultimately ask
for; enwik9 (deferred payload) is not measured.

Usage:
  python -u measure_preconditions.py            # reduced (< 10 min)
  python -u measure_preconditions.py --full     # full dev-slice campaign
  python -u measure_preconditions.py --corpus enwik8 ...   # subset

Outputs (out/): results JSONL (bench/results.py schema),
qualification/advantage/indicators/ladder CSVs, summary_<scale>.json;
on --full additionally the top-level deliverables findings.md,
qualification.csv, out/manifest_qualifies_patch.json.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import precond_lib as P

import loader  # noqa: E402  (path set up by precond_lib)
import metrics  # noqa: E402
import results as results_store  # noqa: E402  (bench/results.py)

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"

SCALES = {
    "reduced": dict(
        label="reduced",
        code_slice=256 * 1024,
        w_list=(2, 3),
        ind_slice=2 * 1024 * 1024,
        sync_window=128 * 1024,
    ),
    "full": dict(
        label="full",
        code_slice=1024 * 1024,
        w_list=(2, 3, 4),
        ind_slice=8 * 1024 * 1024,
        sync_window=128 * 1024,
    ),
}

BLOCK_BASELINES = ("zstd", "lzma")   # matched-granularity (seekable proxy)
STREAM_BASELINES = ("zstd", "lzma", "bz2", "gzip", "brotli")

# Pre-declared thresholds (fixed before the campaign ran; see findings).
TH = dict(
    stat_q=0.05, stat_m=0.20,      # rel. spread of per-third ladder levels
    mix_q=0.05, mix_m=0.15,        # max |acf| at lags >= 256
    pred_q=0.15, pred_m=0.05,      # relative held-out ladder drop
    h3_q=0.8, h3_m=0.1,            # bpb drop => TC-rate lower bound
    h6_q=0.10, h6_m=0.05,          # I(H;L|C) bits/byte (H6 falsif. = 0.1)
    h16_bpb=7.9,                   # incompressibility floor for H16 role
    periodic=0.05,                 # |acf| flag for positional periodicity
    base64_lo=5.90, base64_hi=6.10,  # closed-form H1-calibration anchor
)

STRUCTURAL_LAGS = (77, 4096, 8192)
FAR_LAGS = (256, 1024, 4096, 8192)


# ---------------------------------------------------------------------------
# Per-corpus measurement
# ---------------------------------------------------------------------------


def measure_corpus(name: str, cfg: dict) -> dict:
    entry = loader.entry(name)
    n_code = min(cfg["code_slice"], loader.corpus_bytes(name))
    data = P.slice_bytes_of(name, n_code)
    blocks = [data[j: j + P.K_SYNC] for j in range(0, len(data), P.K_SYNC)]

    print(f"[{name}] slice={len(data)} bytes, {len(blocks)} blocks @64KiB",
          flush=True)

    # --- (a) reference-coder measurement at each W --------------------------
    rnr = {}
    for W in cfg["w_list"]:
        t0 = time.perf_counter()
        rnr[W] = P.rnr_archive_measure(data, W)
        rnr[W]["enc_s"] = time.perf_counter() - t0
        print(f"[{name}]   rnr1 W={W}: total {rnr[W]['total_bpb']:.4f} bpb "
              f"(repair {rnr[W]['repair_bpb']:.4f}, coded "
              f"{rnr[W]['coded_fraction']:.2f}) in {rnr[W]['enc_s']:.1f}s",
              flush=True)
    best_W = min(rnr, key=lambda w: rnr[w]["total_bpb"])

    # --- baselines ----------------------------------------------------------
    import baselines as B  # noqa: PLC0415

    avail = set(B.available_codecs())
    base_block, base_stream = {}, {}
    for coder in BLOCK_BASELINES:
        if coder in avail:
            bpbs = P.baseline_block_bpbs(blocks, coder)
            base_block[coder] = {
                "bpbs": bpbs,
                "mean_bpb": float(np.mean(bpbs)),
                "total_bytes": int(
                    sum(round(b * len(blk) / 8) for b, blk in zip(bpbs, blocks))
                ),
            }
    for coder in STREAM_BASELINES:
        if coder in avail:
            base_stream[coder] = P.baseline_stream_bpb(data, coder)
    stream_best = min(base_stream, key=base_stream.get)
    print(f"[{name}]   baselines: block "
          + ", ".join(f"{c}={v['mean_bpb']:.4f}" for c, v in base_block.items())
          + " | stream best " f"{stream_best}={base_stream[stream_best]:.4f}",
          flush=True)

    # --- paired advantage statistics (frozen methodology) -------------------
    rnr_block_bpbs = [b["bpb"] for b in rnr[best_W]["per_block"]]
    primary = "zstd" if "zstd" in base_block else sorted(base_block)[0]
    base_bpbs = base_block[primary]["bpbs"]
    seed_boot = metrics.derive_seed(f"precond:{name}:adv", 0, tag="boot")
    ci = metrics.paired_bootstrap_diff(rnr_block_bpbs, base_bpbs, seed=seed_boot)
    pval = P.bootstrap_p_two_sided(
        np.asarray(rnr_block_bpbs) - np.asarray(base_bpbs), seed=seed_boot
    )
    dval = metrics.cohens_d_paired(rnr_block_bpbs, base_bpbs)
    base_block_archive = base_block[primary]["total_bytes"]
    adv = {
        "best_W": best_W,
        "primary_baseline": f"{primary}-block64k",
        "diff_mean_bpb": float(np.mean(rnr_block_bpbs) - np.mean(base_bpbs)),
        "ci_lo": ci.lo,
        "ci_hi": ci.hi,
        "ci_method": ci.method,
        "p_boot_secondary": pval,
        "cohens_d": dval,
        # Def-2.1 ratio, conservative: full RNR archive (all container
        # overhead) over bare concatenated baseline block blobs.
        "ratio_block": rnr[best_W]["archive_bytes"] / base_block_archive,
        "ratio_stream_best": rnr[best_W]["total_bpb"] / base_stream[stream_best],
        "stream_best_coder": stream_best,
    }
    if ci.hi < 0:
        adv["verdict"] = "qualifies"
    elif ci.lo > 0:
        adv["verdict"] = "fails"
    else:
        adv["verdict"] = "marginal"
    # Definition 2.1 compares FULL encoder output against the baseline:
    # a per-block repair-stream win that vanishes once container overhead
    # is charged (ratio_block >= 1) is not an advantage. This guards the
    # store-mode controls, where raw fallback (8.0000 bpb/block) "beats"
    # the baseline's per-block framing overhead trivially.
    if adv["verdict"] == "qualifies" and adv["ratio_block"] >= 1.0:
        adv["verdict"] = "marginal"
        adv["note"] = ("per-block repair rate beats the baseline, but the "
                       "full-archive Def-2.1 ratio is >= 1 (container "
                       "overhead); downgraded")

    # --- sync-point overhead (BUG-003-D / H10(c) datum) ---------------------
    sync = None
    win = min(cfg["sync_window"], len(data))
    if win >= 2 * P.K_SYNC:
        wdata = data[:win]
        with_sync = P.rnr_archive_measure(wdata, 3, K=P.K_SYNC)
        if with_sync["coded_fraction"] == 1.0:
            no_sync = P.rnr_archive_measure(wdata, 3, K=win)
            sync = {
                "window_bytes": win,
                "W": 3,
                "K": P.K_SYNC,
                "overhead_total": with_sync["archive_bytes"]
                / no_sync["archive_bytes"] - 1.0,
                "overhead_repair": with_sync["repair_bytes"]
                / no_sync["repair_bytes"] - 1.0,
                "note": "coded; cold-context + index cost per sync",
            }
        else:
            import rnr1  # noqa: PLC0415

            m_sync = (win + P.K_SYNC - 1) // P.K_SYNC
            nosync_size = rnr1.HEADER_SIZE + rnr1.INDEX_ENTRY_SIZE + win
            sync = {
                "window_bytes": win,
                "W": 3,
                "K": P.K_SYNC,
                "overhead_total": (m_sync - 1) * rnr1.INDEX_ENTRY_SIZE
                / nosync_size,
                "overhead_repair": 0.0,
                "note": "store-mode data; overhead = index entries only "
                "(no cold-context term)",
            }
        print(f"[{name}]   sync overhead @K=64KiB: "
              f"{sync['overhead_total']*100:.3f}% ({sync['note']})", flush=True)

    # --- (b) stationarity / mixing indicators -------------------------------
    n_ind = min(cfg["ind_slice"], loader.corpus_bytes(name))
    x = np.frombuffer(P.slice_bytes_of(name, n_ind), dtype=np.uint8)
    thirds = P.thirds_ladder(x)
    spread_h0 = P.rel_spread([t["h0"] for t in thirds])
    spread_h2 = P.rel_spread([t["h2"] for t in thirds])
    stat_metric = max(spread_h0, spread_h2)
    stat = {
        "thirds": thirds,
        "spread_h0": spread_h0,
        "spread_h2": spread_h2,
        "metric": stat_metric,
        "verdict": ("qualifies" if stat_metric <= TH["stat_q"]
                    else "marginal" if stat_metric <= TH["stat_m"] else "fails"),
    }

    acf = P.byte_acf(x)
    far_vals = [abs(acf[l]) for l in FAR_LAGS if acf.get(l) is not None]
    far_max = max(far_vals) if far_vals else 0.0
    periodic = [l for l in STRUCTURAL_LAGS
                if acf.get(l) is not None and abs(acf[l]) >= TH["periodic"]]
    mix = {
        "acf": {str(k): v for k, v in acf.items()},
        "far_max": far_max,
        "periodic_lags": periodic,
        "verdict": ("qualifies" if far_max <= TH["mix_q"]
                    else "marginal" if far_max <= TH["mix_m"] else "fails"),
    }

    # --- (c) entropy-rate ladder --------------------------------------------
    holdout = {k: P.holdout_cross_entropy(x, k) for k in range(0, 5)}
    h0 = holdout[0]["ce"]
    min_order = min(holdout, key=lambda k: holdout[k]["ce"])
    min_ce = holdout[min_order]["ce"]
    drop_bpb = h0 - min_ce
    drop_rel = drop_bpb / h0 if h0 > 0 else 0.0
    pred = {
        "holdout": {str(k): v for k, v in holdout.items()},
        "h0": h0,
        "min_ce": min_ce,
        "min_order": min_order,
        "drop_bpb": drop_bpb,
        "drop_rel": drop_rel,
        "verdict": ("qualifies" if drop_rel >= TH["pred_q"]
                    else "marginal" if drop_rel >= TH["pred_m"] else "fails"),
    }
    cmi = P.nibble_cmi(x)
    print(f"[{name}]   indicators: H0={h0:.3f} minCE={min_ce:.3f}@k={min_order}"
          f" drop={drop_bpb:.3f} | spread={stat_metric:.4f} far|acf|="
          f"{far_max:.4f} periodic={periodic} | I(H;L|C)={cmi:.3f}", flush=True)

    return {
        "corpus": name,
        "class": entry["class"],
        "prior_qualifies_for": list(entry.get("qualifies_for", [])),
        "slice_bytes": len(data),
        "ind_bytes": int(x.size),
        "rnr": rnr,
        "baseline_block": base_block,
        "baseline_stream": base_stream,
        "adv": adv,
        "sync": sync,
        "stat": stat,
        "mix": mix,
        "pred": pred,
        "cmi": cmi,
    }


# ---------------------------------------------------------------------------
# Per-hypothesis verdict engine (role-aware; rules documented in findings)
# ---------------------------------------------------------------------------


def h_verdicts(r: dict) -> dict:
    name, cls = r["corpus"], r["class"]
    pred_v = r["pred"]["verdict"]
    stat_v = r["stat"]["verdict"]
    out = {}

    def put(h, verdict, reason):
        out[h] = {"verdict": verdict, "reason": reason}

    candidates = set(r["prior_qualifies_for"])
    # measured H6 addition candidate for ASCII-text corpora (class 3.1)
    if cls.startswith("3.1") and "H6" not in candidates:
        candidates.add("H6")

    for h in sorted(candidates):
        if h == "H1":
            if name == "ctrl_base64":
                ok = TH["base64_lo"] <= r["pred"]["h0"] <= TH["base64_hi"]
                put("H1", "qualifies" if ok else "fails",
                    f"calibration control: measured H0={r['pred']['h0']:.3f} "
                    f"vs closed-form ~5.99 bpb "
                    f"({'inside' if ok else 'outside'} "
                    f"[{TH['base64_lo']},{TH['base64_hi']}])")
            elif pred_v == "qualifies" and stat_v in ("qualifies", "marginal"):
                put("H1", "qualifies",
                    f"held-out drop {r['pred']['drop_rel']*100:.0f}% and "
                    f"per-third spread {r['stat']['metric']*100:.1f}%")
            elif pred_v in ("qualifies", "marginal"):
                put("H1", "marginal",
                    f"predictive structure {pred_v}, stationarity {stat_v}")
            else:
                put("H1", "fails", "no exploitable predictive structure")
        elif h == "H2":
            periodic = r["mix"]["periodic_lags"]
            if cls.startswith("3.2") and (periodic or pred_v == "qualifies"):
                put("H2", "qualifies",
                    f"positional periodicity at lags {periodic}, "
                    f"held-out drop {r['pred']['drop_bpb']:.2f} bpb")
            elif cls.startswith("3.2"):
                put("H2", "marginal", "structured-binary class but no "
                    "periodicity/prediction signal at this scale")
            else:
                put("H2", "fails", "not a structured-binary corpus")
        elif h == "H3":
            d = r["pred"]["drop_bpb"]
            v = ("qualifies" if d >= TH["h3_q"]
                 else "marginal" if d >= TH["h3_m"] else "fails")
            put("H3", v,
                f"TC-rate lower bound H0-minCE = {d:.2f} bpb "
                f"(TC(D_N) >= {d:.2f}*N bits; H3 needs Theta(N))")
        elif h == "H4":
            if pred_v == "fails":
                put("H4", "fails", "no low-dimensional predictive structure "
                    "signal (screen); LB estimator not run here")
            else:
                put("H4", "screen-pass",
                    "intrinsic-dimension estimator not run here; "
                    f"predictive-structure screen {pred_v}")
        elif h == "H5":
            if cls.startswith("3.4") and stat_v in ("fails", "marginal") \
                    and pred_v != "fails":
                put("H5", "qualifies",
                    f"heterogeneity confirmed: per-third spread "
                    f"{r['stat']['metric']*100:.1f}% (scheduler target)")
            elif cls.startswith("3.4"):
                put("H5", "marginal",
                    f"heterogeneous class but spread only "
                    f"{r['stat']['metric']*100:.1f}%")
            else:
                put("H5", "fails", "not a heterogeneous-archive corpus")
        elif h == "H6":
            i = r["cmi"]
            v = ("qualifies" if i >= TH["h6_q"]
                 else "marginal" if i >= TH["h6_m"] else "fails")
            band = ("0.3-0.7 (ASCII)" if cls.startswith("3.1")
                    else "1.0-2.0 (float32)" if cls.startswith("3.3")
                    else "n/a")
            put("H6", v, f"I(H;L|C)={i:.3f} bpb, predicted band {band}, "
                f"falsification floor {TH['h6_q']}")
        elif h == "H8":
            put("H8", "qualifies",
                "byte-identity conformance corpus role is data-agnostic")
        elif h == "H15":
            put("H15", "not-measured", "enwik9 payload deferred; not fetched")
        elif h == "H16":
            best = min(r["baseline_stream"].values())
            raw_frac = 1.0 - r["rnr"][r["adv"]["best_W"]]["coded_fraction"]
            if best >= TH["h16_bpb"] and r["pred"]["drop_bpb"] < 0.05:
                put("H16", "qualifies",
                    f"incompressible confirmed: best baseline {best:.3f} bpb, "
                    f"ladder drop {r['pred']['drop_bpb']:.3f}; RNR store-mode "
                    f"fraction {raw_frac:.2f}")
            elif best >= 7.5:
                put("H16", "marginal", f"best baseline {best:.3f} bpb")
            else:
                put("H16", "fails",
                    f"compressible (best baseline {best:.3f} bpb): "
                    "not a worst-case control")
        else:
            put(h, "not-instrumented", "no precondition instrument maps "
                "to this hypothesis in this campaign")
    return out


def proposed_qualifies(r: dict) -> dict:
    """Manifest patch entry: keep prior minus clear fails, add measured
    qualifiers not previously listed (conservative)."""
    prior = list(r["prior_qualifies_for"])
    verdicts = r["h_verdicts"]
    remove = [h for h in prior
              if verdicts.get(h, {}).get("verdict") == "fails"]
    add = [h for h, v in verdicts.items()
           if h not in prior and v["verdict"] == "qualifies"]
    proposed = [h for h in prior if h not in remove] + sorted(add)
    return {
        "prior": prior,
        "proposed": proposed,
        "remove": remove,
        "add": sorted(add),
        "verdicts": {h: v["verdict"] for h, v in verdicts.items()},
        "justification": "; ".join(
            f"{h}: {v['verdict']} ({v['reason']})"
            for h, v in sorted(verdicts.items())
        ),
    }


# ---------------------------------------------------------------------------
# Persistence: JSONL (bench schema), CSVs, findings, patch
# ---------------------------------------------------------------------------


def append_records(res: dict, scale: str, jsonl_path: Path) -> None:
    name = res["corpus"]
    sha = P.corpus_sha(name)
    common = dict(
        corpus_sha256=sha,
        protocol="experiments/preconditions",
        methods="bench/METHODS.md v1.0 (secondary, labeled)",
        secondary=True,
        deterministic_coder=True,
        scale=scale,
    )

    def rec(coder, config, measured):
        seed = metrics.derive_seed(f"precond:{name}:{coder}", 0, tag="run")
        results_store.append(
            jsonl_path,
            results_store.make_record(name, coder, config, seed, measured,
                                      **common),
        )

    for W, m in res["rnr"].items():
        rec("rnr1-ref",
            {"W": W, "K": m["K"], "slice_bytes": m["n"], "offset": 0},
            {"bpb": m["total_bpb"], "repair_bpb": m["repair_bpb"],
             "a_ideal_bpb": m["a_ideal_bpb"], "b_model_bpb": m["b_model_bpb"],
             "c_rho_bpb": m["c_rho_bpb"], "d_nu_bpb": m["d_nu_bpb"],
             "e_mu_bpb": m["e_mu_bpb"], "coded_fraction": m["coded_fraction"],
             "compressed_size": m["archive_bytes"], "original_size": m["n"],
             "enc_s": m["enc_s"]})
    for coder, v in res["baseline_block"].items():
        rec(f"{coder}-block64k",
            {"block_bytes": P.K_SYNC, "slice_bytes": res["slice_bytes"]},
            {"bpb": v["mean_bpb"], "compressed_size": v["total_bytes"],
             "original_size": res["slice_bytes"]})
    for coder, bpb in res["baseline_stream"].items():
        rec(f"{coder}-stream", {"slice_bytes": res["slice_bytes"]},
            {"bpb": bpb, "original_size": res["slice_bytes"]})
    a = res["adv"]
    rec("precond-advantage",
        {"best_W": a["best_W"], "vs": a["primary_baseline"],
         "slice_bytes": res["slice_bytes"]},
        {k: a[k] for k in ("diff_mean_bpb", "ci_lo", "ci_hi",
                           "p_boot_secondary", "cohens_d", "ratio_block",
                           "ratio_stream_best")})
    if res["sync"]:
        s = res["sync"]
        rec("rnr1-sync-overhead",
            {"W": s["W"], "K": s["K"], "window_bytes": s["window_bytes"]},
            {"overhead_total": s["overhead_total"],
             "overhead_repair": s["overhead_repair"], "note": s["note"]})
    acf = res["mix"]["acf"]
    rec("precond-indicators", {"ind_bytes": res["ind_bytes"]},
        {"h0_holdout": res["pred"]["h0"], "min_ce": res["pred"]["min_ce"],
         "min_order": res["pred"]["min_order"],
         "drop_bpb": res["pred"]["drop_bpb"],
         "drop_rel": res["pred"]["drop_rel"],
         "stat_spread": res["stat"]["metric"],
         "acf_far_max": res["mix"]["far_max"],
         "acf_lag1": acf.get("1"), "acf_lag256": acf.get("256"),
         "acf_lag4096": acf.get("4096"), "acf_lag8192": acf.get("8192"),
         "nibble_cmi": res["cmi"]})


def write_csvs(all_res: list, scale: str, top_level: bool) -> None:
    qual_rows = []
    for r in all_res:
        row = {
            "corpus": r["corpus"],
            "class": r["class"],
            "P_ADV": r["adv"]["verdict"],
            "adv_diff_bpb": f"{r['adv']['diff_mean_bpb']:.4f}",
            "adv_ci": f"[{r['adv']['ci_lo']:.4f},{r['adv']['ci_hi']:.4f}]",
            "adv_ratio_block": f"{r['adv']['ratio_block']:.4f}",
            "adv_ratio_stream": f"{r['adv']['ratio_stream_best']:.4f}",
            "P_STAT": r["stat"]["verdict"],
            "stat_spread": f"{r['stat']['metric']:.4f}",
            "P_MIX": r["mix"]["verdict"],
            "acf_far_max": f"{r['mix']['far_max']:.4f}",
            "periodic_lags": "|".join(map(str, r["mix"]["periodic_lags"])),
            "P_PRED": r["pred"]["verdict"],
            "drop_bpb": f"{r['pred']['drop_bpb']:.4f}",
            "drop_rel": f"{r['pred']['drop_rel']:.4f}",
            "nibble_cmi": f"{r['cmi']:.4f}",
            "sync_overhead": (f"{r['sync']['overhead_total']:.5f}"
                              if r["sync"] else ""),
            "h_verdicts": "|".join(
                f"{h}:{v['verdict']}" for h, v in sorted(r["h_verdicts"].items())
            ),
        }
        qual_rows.append(row)
    paths = [OUT / f"qualification_{scale}.csv"]
    if top_level:
        paths.append(HERE / "qualification.csv")
    for path in paths:
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(qual_rows[0].keys()))
            w.writeheader()
            w.writerows(qual_rows)

    with open(OUT / f"advantage_{scale}.csv", "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["corpus", "W", "total_bpb", "repair_bpb", "a_ideal_bpb",
                    "c_rho_bpb", "b_model_bpb", "d_nu_bpb", "e_mu_bpb",
                    "coded_fraction", "zstd_block_bpb", "lzma_block_bpb",
                    "stream_best", "stream_best_bpb"])
        for r in all_res:
            for W, m in sorted(r["rnr"].items()):
                w.writerow([
                    r["corpus"], W, f"{m['total_bpb']:.4f}",
                    f"{m['repair_bpb']:.4f}",
                    "" if m["a_ideal_bpb"] is None else f"{m['a_ideal_bpb']:.4f}",
                    "" if m["c_rho_bpb"] is None else f"{m['c_rho_bpb']:.5f}",
                    f"{m['b_model_bpb']:.6f}", f"{m['d_nu_bpb']:.5f}",
                    f"{m['e_mu_bpb']:.5f}", f"{m['coded_fraction']:.3f}",
                    f"{r['baseline_block'].get('zstd', {}).get('mean_bpb', float('nan')):.4f}",
                    f"{r['baseline_block'].get('lzma', {}).get('mean_bpb', float('nan')):.4f}",
                    r["adv"]["stream_best_coder"],
                    f"{r['baseline_stream'][r['adv']['stream_best_coder']]:.4f}",
                ])

    with open(OUT / f"ladder_{scale}.csv", "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["corpus", "order", "holdout_ce_bpb", "holdout_coverage",
                    "third1_h", "third2_h", "third3_h"])
        for r in all_res:
            for k in range(0, 5):
                hk = r["pred"]["holdout"][str(k)]
                thirds = [t.get(f"h{k}") for t in r["stat"]["thirds"]] \
                    if k <= 3 else [None, None, None]
                w.writerow([r["corpus"], k, f"{hk['ce']:.4f}",
                            f"{hk['coverage']:.4f}"]
                           + ["" if t is None else f"{t:.4f}" for t in thirds])

    with open(OUT / f"indicators_{scale}.csv", "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        lags = [str(l) for l in P.ACF_LAGS]
        w.writerow(["corpus"] + [f"acf_{l}" for l in lags]
                   + ["stat_spread_h0", "stat_spread_h2", "nibble_cmi"])
        for r in all_res:
            w.writerow([r["corpus"]]
                       + [("" if r["mix"]["acf"].get(l) is None
                           else f"{r['mix']['acf'][l]:.5f}") for l in lags]
                       + [f"{r['stat']['spread_h0']:.5f}",
                          f"{r['stat']['spread_h2']:.5f}",
                          f"{r['cmi']:.4f}"])


def write_patch(all_res: list, scale: str) -> None:
    entries = {r["corpus"]: proposed_qualifies(r) for r in all_res}
    patch = {
        "patch_target": "data/MANIFEST.json",
        "field": "entries.<name>.qualifies_for",
        "generated_by": "experiments/preconditions/measure_preconditions.py",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "scale": scale,
        "scale_caveat": (
            "Dev-slice measurement (coding <= 1 MiB, indicators <= 8 MiB "
            "per corpus, front-of-file slices) with the reference order-W "
            "counting predictor as the side-information instrument. "
            "Verdicts are precondition screens, not full-corpus "
            "hypothesis outcomes; BUG-001-D/BUG-002-D full-scale external "
            "runs remain open. enwik9 (deferred) not measured."
        ),
        "apply_rule": (
            "For each entry, replace qualifies_for with 'proposed' after "
            "human review; 'remove' lists prior hypotheses whose "
            "precondition measurably FAILED, 'add' lists measured "
            "qualifiers not previously listed."
        ),
        "entries": entries,
    }
    with open(OUT / "manifest_qualifies_patch.json", "w", encoding="utf-8") as f:
        json.dump(patch, f, indent=2, sort_keys=True)


def render_findings(all_res: list, scale: str, cfg: dict, holm: dict,
                    path: Path) -> None:
    L = []
    add = L.append
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    add("# Precondition statistics — qualification of the corpus suite")
    add("")
    add(f"Campaign scale: **{scale}** (coding slice "
        f"{cfg['code_slice']//1024} KiB, W in {list(cfg['w_list'])}, "
        f"K = 64 KiB; indicator slice {cfg['ind_slice']//1024//1024} MiB; "
        f"front-of-file slices). Generated {now} by "
        "`experiments/preconditions/measure_preconditions.py`.")
    add("")
    add("Resolves, at dev-slice scale, the measurements named external by "
        "`docs/backlog/blocked/BUG-001-D` (conditional-advantage ratio of "
        "the Abstract's scoped claim, Definition 2.1 / inequality (10.1)), "
        "`BUG-003-D` (sync-point ratio overhead at typical (K, W), "
        "exp-design H10(c)), and `BUG-010-D` (which corpora satisfy the "
        "section-10.2.1 precondition). Statistics follow the frozen "
        "methodology (`bench/METHODS.md` v1.0): BCa bootstrap, B = 10000, "
        "derived seeds; this screen is **secondary/labeled** relative to "
        "the frozen H1–H18 primary families, so no amendment is required.")
    add("")

    add("## 1. Conditional-advantage ratio (BUG-001-D instrument)")
    add("")
    add("Instrument: reference Type-I coder (`impl/rnr1.py`), order-W "
        "adaptive counting predictor as side information, sync K = 64 KiB. "
        "Primary comparison: paired per-block repair bpb vs per-block "
        "zstd-22 on the same 64 KiB blocks (both coders restart cold at "
        "each block — the in-process proxy for seekable formats at matched "
        "granularity). `ratio_block` charges the FULL RNR archive "
        "(header+index+hashes) against bare baseline blobs (conservative); "
        "`ratio_stream` compares against the best whole-slice streaming "
        "baseline (harder). Def. 2.1 advantage needs ratio < 1.")
    add("")
    add("| corpus | best W | RNR total bpb | zstd-block bpb | lzma-block bpb "
        "| paired diff bpb [95% BCa] | p (boot, secondary) | Holm reject | "
        "ratio_block | stream best | ratio_stream |")
    add("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in all_res:
        a = r["adv"]
        bw = a["best_W"]
        zb = r["baseline_block"].get("zstd", {}).get("mean_bpb")
        lb = r["baseline_block"].get("lzma", {}).get("mean_bpb")
        add(f"| {r['corpus']} | {bw} | {r['rnr'][bw]['total_bpb']:.4f} | "
            f"{zb:.4f} | {lb:.4f} | "
            f"{a['diff_mean_bpb']:+.4f} [{a['ci_lo']:+.4f}, {a['ci_hi']:+.4f}]"
            f" | {a['p_boot_secondary']:.4f} | "
            f"{'yes' if holm[r['corpus']] else 'no'} | "
            f"{a['ratio_block']:.4f} | {a['stream_best_coder']} "
            f"{r['baseline_stream'][a['stream_best_coder']]:.4f} | "
            f"{a['ratio_stream_best']:.4f} |")
    add("")

    add("### 1a. Inequality (10.1) term breakdown at best W (bits/byte)")
    add("")
    add("| corpus | (a) H(X|Y) proxy | (b) model | (c) rho | (d) nu | "
        "(e) mu | total | coded frac |")
    add("|---|---|---|---|---|---|---|---|")
    for r in all_res:
        m = r["rnr"][r["adv"]["best_W"]]
        fa = "n/a" if m["a_ideal_bpb"] is None else f"{m['a_ideal_bpb']:.4f}"
        fc = "n/a" if m["c_rho_bpb"] is None else f"{m['c_rho_bpb']:.5f}"
        add(f"| {r['corpus']} | {fa} | {m['b_model_bpb']:.6f} | {fc} | "
            f"{m['d_nu_bpb']:.5f} | {m['e_mu_bpb']:.5f} | "
            f"{m['total_bpb']:.4f} | {m['coded_fraction']:.2f} |")
    add("")
    add("Term (b) is only the 8-byte model-description hash: the counting "
        "predictor carries no trained parameters, so the neural L(M)/mN "
        "term of (10.1) is NOT exercised by this instrument (it would only "
        "make the inequality harder). (a)+(c) sum to the repair rate over "
        "coded blocks; store-mode (raw) blocks pay 8 bpb.")
    add("")

    add("## 2. Sync-point overhead at K = 64 KiB, W = 3 (BUG-003-D datum)")
    add("")
    add("| corpus | window | overhead (total archive) | overhead "
        "(repair only) | note |")
    add("|---|---|---|---|---|")
    for r in all_res:
        s = r["sync"]
        if not s:
            add(f"| {r['corpus']} | — | — | — | window smaller than 2K |")
            continue
        add(f"| {r['corpus']} | {s['window_bytes']//1024} KiB | "
            f"{s['overhead_total']*100:.3f}% | {s['overhead_repair']*100:.3f}%"
            f" | {s['note']} |")
    add("")
    add("H10(c) target: < 2% at K = 64 KiB; falsification threshold: > 10%. "
        "Measured on a 2-block window (the per-sync cold-context cost is "
        "deterministic and scale-free per sync point); the adaptive "
        "counting predictor restarts from empty counts at each sync, which "
        "is the worst-case cold-context regime for this instrument.")
    add("")

    add("## 3. Stationarity / mixing indicators")
    add("")
    add("| corpus | per-third spread (max of H0,H2) | verdict | "
        "max abs acf (lags >= 256) | periodic lags | verdict |")
    add("|---|---|---|---|---|---|")
    for r in all_res:
        add(f"| {r['corpus']} | {r['stat']['metric']*100:.2f}% | "
            f"{r['stat']['verdict']} | {r['mix']['far_max']:.4f} | "
            f"{r['mix']['periodic_lags'] or '—'} | {r['mix']['verdict']} |")
    add("")
    add("A large per-third spread is DISqualifying for stationary-source "
        "hypotheses (H1/H3) but QUALIFYING for the heterogeneous-scheduler "
        "hypothesis (H5). Likewise, strong far-lag autocorrelation "
        "(positional periodicity: lag 4096 = SQLite page, 8192 = telemetry "
        "row stride, 77 = base64 line) is qualifying evidence for Type-II "
        "positional structure (H2). Full ACF table: "
        "`out/indicators_" + scale + ".csv`.")
    add("")

    add("## 4. Entropy-rate ladder (held-out order-0..4) and I(H;L|C)")
    add("")
    add("| corpus | H0 | min CE @ order | drop bpb | drop % | coverage@min |"
        " I(H;L|C) |")
    add("|---|---|---|---|---|---|---|")
    for r in all_res:
        p = r["pred"]
        cov = p["holdout"][str(p["min_order"])]["coverage"]
        add(f"| {r['corpus']} | {p['h0']:.4f} | {p['min_ce']:.4f} @ "
            f"{p['min_order']} | {p['drop_bpb']:.4f} | "
            f"{p['drop_rel']*100:.1f}% | {cov:.3f} | {r['cmi']:.4f} |")
    add("")
    add("Held-out cross-entropy (counts from the first half, Laplace-"
        "smoothed evaluation on the second half) is a real code length, "
        "hence an honest UPPER bound on the order-k entropy rate — no "
        "plug-in sparsity under-bias. drop = H0 − min_k CE_k lower-bounds "
        "the total-correlation rate: TC(D_N)/N >= H0 − h_k − o(1), the H3 "
        "screening quantity. Per-order values: `out/ladder_" + scale
        + ".csv`.")
    add("")

    add("## 5. QUALIFICATION TABLE (corpus x precondition)")
    add("")
    add("Thresholds (pre-declared): P-ADV qualifies iff the paired-diff "
        "95% BCa CI lies below 0 (fails iff above 0); P-STAT spread "
        f"<= {TH['stat_q']} / <= {TH['stat_m']}; P-MIX far-|acf| "
        f"<= {TH['mix_q']} / <= {TH['mix_m']}; P-PRED relative drop "
        f">= {TH['pred_q']} / >= {TH['pred_m']}.")
    add("")
    add("| corpus | class | P-ADV | P-STAT | P-MIX | P-PRED | "
        "per-hypothesis verdicts |")
    add("|---|---|---|---|---|---|---|")
    for r in all_res:
        hv = ", ".join(f"{h}: **{v['verdict']}**"
                       for h, v in sorted(r["h_verdicts"].items()))
        add(f"| {r['corpus']} | {r['class']} | {r['adv']['verdict']} | "
            f"{r['stat']['verdict']} | {r['mix']['verdict']} | "
            f"{r['pred']['verdict']} | {hv} |")
    add("")

    add("## 6. Per-hypothesis scoping conclusions")
    add("")
    for r in all_res:
        add(f"**{r['corpus']}** — prior `qualifies_for`: "
            f"{r['prior_qualifies_for']}")
        for h, v in sorted(r["h_verdicts"].items()):
            add(f"- {h}: {v['verdict']} — {v['reason']}")
        pq = proposed_qualifies(r)
        if pq["remove"] or pq["add"]:
            add(f"- MANIFEST patch: remove {pq['remove'] or 'nothing'}, "
                f"add {pq['add'] or 'nothing'} -> proposed {pq['proposed']}")
        else:
            add("- MANIFEST patch: no change (all prior scopes survive the "
                "screen)")
        add("")
    add("Patch file: `out/manifest_qualifies_patch.json` (data/ is NOT "
        "edited directly).")
    add("")

    add("## 7. Honest caveats")
    add("")
    add("1. **Scale**: dev slices (coding <= "
        f"{cfg['code_slice']//1024} KiB, indicators <= "
        f"{cfg['ind_slice']//1024//1024} MiB, front-of-file). "
        "BUG-001-D/BUG-010-D name FULL-corpus campaigns with a neural "
        "predictor in the amortized regime; this campaign is the local, "
        "instrumented version, and the external full-scale items remain "
        "open. Front-of-file slices bias heterogeneous tars (calgary: the "
        "coding slice covers only the first member(s); indicator slices "
        "cover the whole file for corpora <= the indicator budget).")
    add("2. **Instrument strength**: the order-W counting predictor is a "
        "deliberately weak side-information instrument. A P-ADV failure "
        "does NOT falsify the paper's conditional claim (the claim is "
        "scoped to decoder-reproducible NEURAL side information with "
        "amortized L(M)); it means the advantage is not demonstrated by "
        "this instrument at this scale. A P-ADV pass, conversely, is "
        "strong evidence the precondition is satisfiable.")
    add("3. **Determinism**: coder and indicators are deterministic; "
        "coding-rate numbers are single runs (the METHODS 5-seed "
        "repetition governs TIMING hypotheses, none reported here). All "
        "bootstrap seeds derive from the master seed per METHODS section "
        "4 and are recorded in the JSONL.")
    add("4. **Secondary status**: the precondition screen is not one of "
        "the frozen H1–H18 primary tests; all p-values are labeled "
        "secondary, with Holm-Bonferroni applied within this screening "
        "family (9 corpora) for transparency only.")
    add("5. **Plug-in bias**: per-third ladder values are plug-in "
        "estimates (biased low at high orders/sparse contexts); "
        "qualification uses held-out cross-entropy where bias direction "
        "is safe, and per-third plug-in values only through the "
        "RELATIVE spread.")
    add("6. **enwik9**: deferred payload, not fetched, not measured; the "
        "H15 scope is untouched.")
    add("7. **Store-mode corpora**: for incompressible controls the "
        "sync-overhead datum reduces to index bytes (no cold-context "
        "term); reported as such, not blended with coded corpora.")
    add("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------


def run_campaign(scale: str, only: list | None = None) -> dict:
    cfg = SCALES[scale]
    OUT.mkdir(exist_ok=True)
    names = [n for n in loader.corpus_names() if not only or n in only]
    t0 = time.perf_counter()
    all_res = []
    jsonl_path = OUT / "results.jsonl"
    for name in names:
        res = measure_corpus(name, cfg)
        res["h_verdicts"] = h_verdicts(res)
        all_res.append(res)

    # Holm within the ADV screening family (secondary, labeled).
    pvals = [r["adv"]["p_boot_secondary"] for r in all_res]
    rejected = metrics.holm_bonferroni(pvals)
    holm = {r["corpus"]: rej for r, rej in zip(all_res, rejected)}

    full = scale == "full"
    for res in all_res:
        append_records(res, scale, jsonl_path)
    write_csvs(all_res, scale, top_level=full)
    if full:
        write_patch(all_res, scale)
    findings_path = HERE / ("findings.md" if full
                            else "out/findings_reduced.md")
    render_findings(all_res, scale, cfg, holm, findings_path)

    summary = {
        "scale": scale,
        "config": {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in cfg.items()},
        "corpora": [r["corpus"] for r in all_res],
        "elapsed_s": time.perf_counter() - t0,
        "holm_rejected": holm,
        "results": [
            {k: r[k] for k in ("corpus", "class", "slice_bytes", "ind_bytes",
                               "adv", "sync", "stat", "mix", "pred", "cmi",
                               "h_verdicts", "prior_qualifies_for")}
            | {"rnr": {str(W): {kk: vv for kk, vv in m.items()
                                if kk != "per_block"}
                       for W, m in r["rnr"].items()},
               "rnr_per_block_bpb": {
                   str(W): [b["bpb"] for b in m["per_block"]]
                   for W, m in r["rnr"].items()},
               "baseline_block": r["baseline_block"],
               "baseline_stream": r["baseline_stream"]}
            for r in all_res
        ],
    }
    with open(OUT / f"summary_{scale}.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1, sort_keys=True)
    print(f"\ncampaign [{scale}] done in {summary['elapsed_s']:.0f}s; "
          f"findings -> {findings_path}", flush=True)
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--full", action="store_true",
                    help="full dev-slice campaign (1 MiB coding slices, "
                         "W in {2,3,4}, 8 MiB indicators)")
    ap.add_argument("--corpus", action="append", default=None,
                    help="restrict to named corpora (repeatable)")
    args = ap.parse_args(argv)
    run_campaign("full" if args.full else "reduced", only=args.corpus)
    return 0


if __name__ == "__main__":
    sys.exit(main())
