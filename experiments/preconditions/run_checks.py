#!/usr/bin/env python3
"""PASS/FAIL gate for the precondition-statistics campaign.

Default: unit checks + the REDUCED campaign (< 10 min total).
--full : unit checks + the FULL dev-slice campaign (the real numbers;
         writes findings.md / qualification.csv / manifest patch).

Every check prints [PASS]/[FAIL]; exit code 0 iff all pass.

  python -u run_checks.py [--full]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

import precond_lib as P
import measure_preconditions as M

import loader  # noqa: E402
import metrics  # noqa: E402
import results as results_store  # noqa: E402
import rnr1  # noqa: E402

HERE = Path(__file__).resolve().parent
CHECKS = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    CHECKS.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""),
          flush=True)
    return ok


# ---------------------------------------------------------------------------
# 1. Corpus integrity (manifest sha256 of every measured payload)
# ---------------------------------------------------------------------------


def check_corpus_integrity() -> None:
    for name in loader.corpus_names():
        path = loader.payload_path(name)
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        ok = h.hexdigest() == loader.entry(name)["sha256"]
        check(f"corpus sha256: {name}", ok)


# ---------------------------------------------------------------------------
# 2. Reference-coder instrumentation consistency
# ---------------------------------------------------------------------------


def check_coder_instrumentation() -> None:
    probe = P.slice_bytes_of("scripts_src", 128 * 1024)
    W = 3

    # (i) instrumented encoder emits byte-identical blobs
    blob_a, ideal = P.encode_block_stats(probe[: P.K_SYNC], W)
    blob_ref = rnr1.encode_subblock(probe[: P.K_SYNC], W)
    check("instrumented encoder blob == reference blob", blob_a == blob_ref)

    # (ii) determinism: re-encode identical
    blob_b, ideal_b = P.encode_block_stats(probe[: P.K_SYNC], W)
    check("encoder determinism (blob + ideal bits)",
          blob_b == blob_a and ideal_b == ideal)

    # (iii) realized >= ideal (AC cannot beat its own model), small excess
    realized = 8.0 * len(blob_a)
    ok = ideal <= realized <= ideal * 1.002 + 64.0
    check("ideal <= realized <= ideal*(1+2e-3)+64 bits", ok,
          f"ideal={ideal:.1f} realized={realized:.1f}")

    # (iv) analytic archive size == len(rnr1.pack)
    meas = P.rnr_archive_measure(probe, W)
    packed = rnr1.pack(probe, W=W, K=P.K_SYNC)
    check("analytic archive size == len(pack)",
          meas["archive_bytes"] == len(packed),
          f"{meas['archive_bytes']} vs {len(packed)}")

    # (v) round-trip
    check("round-trip unpack(pack(x)) == x", rnr1.unpack(packed) == probe)

    # (vi) term identity: total = repair + b + d + e exactly
    ident = abs(meas["total_bpb"] - (meas["repair_bpb"] + meas["b_model_bpb"]
                                     + meas["d_nu_bpb"] + meas["e_mu_bpb"]))
    check("term identity total = repair+(b)+(d)+(e)", ident < 1e-9,
          f"residual {ident:.2e}")


# ---------------------------------------------------------------------------
# 3. Indicator calibration (closed-form anchors)
# ---------------------------------------------------------------------------


def check_indicator_calibration() -> None:
    # base64 control: order-0 entropy has a closed form (~5.99 bpb)
    x64 = np.frombuffer(P.slice_bytes_of("ctrl_base64", 1 << 21), np.uint8)
    h0 = P.plugin_cond_entropy(x64, 0)["h"]
    check("ctrl_base64 H0 in closed-form band [5.90, 6.10]",
          5.90 <= h0 <= 6.10, f"H0={h0:.4f}")

    # urandom control: H0 ~ 8, held-out drop ~ 0
    xr = np.frombuffer(P.slice_bytes_of("ctrl_urandom", 1 << 21), np.uint8)
    h0r = P.plugin_cond_entropy(xr, 0)["h"]
    check("ctrl_urandom H0 in [7.98, 8.0]", 7.98 <= h0r <= 8.0,
          f"H0={h0r:.5f}")
    ce0 = P.holdout_cross_entropy(xr, 0)["ce"]
    ce2 = P.holdout_cross_entropy(xr, 2)["ce"]
    check("ctrl_urandom held-out CE does not drop (order 2 >= order 0 - 0.02)",
          ce2 >= ce0 - 0.02, f"ce0={ce0:.4f} ce2={ce2:.4f}")

    # ACF sanity: iid bytes have |acf| < 0.01 at all measured lags
    acf = P.byte_acf(xr)
    worst = max(abs(v) for v in acf.values() if v is not None)
    check("ctrl_urandom |acf| < 0.01 at all lags", worst < 0.01,
          f"max|acf|={worst:.5f}")

    # nibble CMI sanity: iid uniform bytes -> I(H;L|C) ~ 0
    cmi = P.nibble_cmi(xr)
    check("ctrl_urandom I(H;L|C) < 0.01", cmi < 0.01, f"I={cmi:.5f}")


# ---------------------------------------------------------------------------
# 4. Frozen-methodology statistics behave deterministically
# ---------------------------------------------------------------------------


def check_statistics() -> None:
    rng = np.random.default_rng(7)
    xs = rng.normal(1.0, 0.5, size=24)
    seed = metrics.derive_seed("precond:gate", 0, tag="boot")
    ci1 = metrics.bootstrap_ci(xs, seed=seed)
    ci2 = metrics.bootstrap_ci(xs, seed=seed)
    check("BCa bootstrap deterministic under derived seed",
          (ci1.lo, ci1.hi) == (ci2.lo, ci2.hi),
          f"[{ci1.lo:.4f},{ci1.hi:.4f}] method={ci1.method}")
    p_shift = P.bootstrap_p_two_sided(xs, seed=seed)
    p_null = P.bootstrap_p_two_sided(xs - xs.mean(), seed=seed)
    check("bootstrap p: shifted sample significant, centered not",
          p_shift < 0.01 < p_null, f"p_shift={p_shift:.4f} p_null={p_null:.3f}")
    rej = metrics.holm_bonferroni([0.001, 0.04, 0.9])
    check("Holm-Bonferroni known vector", rej == [True, False, False],
          str(rej))


# ---------------------------------------------------------------------------
# 5. Campaign + output validation
# ---------------------------------------------------------------------------


def check_campaign(scale: str) -> None:
    summary = M.run_campaign(scale)
    n_corpora = len(loader.corpus_names())
    check("campaign covered every non-deferred corpus",
          len(summary["corpora"]) == n_corpora,
          f"{len(summary['corpora'])}/{n_corpora}")

    # JSONL store loads, and every record is bench-schema complete
    records = results_store.load(M.OUT / "results.jsonl")
    ok = all(all(k in r for k in results_store.REQUIRED_KEYS) for r in records)
    ok = ok and all("corpus_sha256" in r for r in records)
    check("results.jsonl bench-schema complete (incl. corpus_sha256)", ok,
          f"{len(records)} records")

    # term identity + verdict consistency on the summary
    ident_ok, verdict_ok, ac_ok = True, True, True
    for r in summary["results"]:
        for W, m in r["rnr"].items():
            resid = abs(m["total_bpb"] - (m["repair_bpb"] + m["b_model_bpb"]
                                          + m["d_nu_bpb"] + m["e_mu_bpb"]))
            ident_ok &= resid < 1e-9
            if m["a_ideal_bpb"] is not None:
                ac_ok &= m["c_rho_bpb"] >= 0.0
        a = r["adv"]
        v = ("qualifies" if a["ci_hi"] < 0
             else "fails" if a["ci_lo"] > 0 else "marginal")
        if v == "qualifies" and a["ratio_block"] >= 1.0:
            v = "marginal"  # Def-2.1 full-cost guard (see measure_corpus)
        verdict_ok &= v == a["verdict"]
    check("term identity holds for every (corpus, W)", ident_ok)
    check("coder residual rho >= 0 on all coded measurements", ac_ok)
    check("P-ADV verdicts consistent with BCa CIs", verdict_ok)

    # negative controls behave: urandom must fail P-PRED and engage
    # store mode; base64 must anchor the H1 calibration
    by = {r["corpus"]: r for r in summary["results"]}
    ur = by["ctrl_urandom"]
    check("ctrl_urandom fails P-PRED (no fake structure)",
          ur["pred"]["verdict"] == "fails",
          f"drop={ur['pred']['drop_bpb']:.4f}")
    raw_ok = all(m["coded_fraction"] == 0.0 for m in ur["rnr"].values())
    check("ctrl_urandom fully store-mode (fail-safe engaged)", raw_ok)
    b64 = by["ctrl_base64"]
    check("ctrl_base64 H1 calibration verdict recorded",
          "H1" in b64["h_verdicts"],
          b64["h_verdicts"].get("H1", {}).get("verdict", "missing"))

    # artifacts exist
    want = [M.OUT / f"qualification_{scale}.csv",
            M.OUT / f"advantage_{scale}.csv",
            M.OUT / f"ladder_{scale}.csv",
            M.OUT / f"indicators_{scale}.csv",
            M.OUT / f"summary_{scale}.json"]
    if scale == "full":
        want += [HERE / "findings.md", HERE / "qualification.csv",
                 M.OUT / "manifest_qualifies_patch.json"]
    else:
        want += [M.OUT / "findings_reduced.md"]
    missing = [str(p) for p in want if not p.exists()]
    check("all campaign artifacts written", not missing,
          f"missing: {missing}" if missing else "")

    if scale == "full":
        patch = json.loads((M.OUT / "manifest_qualifies_patch.json")
                           .read_text(encoding="utf-8"))
        known = set(loader.corpus_names())
        ok = set(patch["entries"]) == known
        check("manifest patch covers exactly the measured corpora", ok)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args(argv)
    scale = "full" if args.full else "reduced"

    print(f"== precondition-statistics gate ({scale}) ==", flush=True)
    check_corpus_integrity()
    check_coder_instrumentation()
    check_indicator_calibration()
    check_statistics()
    check_campaign(scale)

    ok = all(CHECKS)
    print(f"\nOVERALL: {'PASS' if ok else 'FAIL'} "
          f"({sum(CHECKS)}/{len(CHECKS)} checks)", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
