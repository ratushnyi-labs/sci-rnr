#!/usr/bin/env python3
"""PASS/FAIL gate for the H5/H6/H7 Type-III experiment.

Default (no flags): runs the REDUCED-scale campaign end to end
(< 10 minutes) and evaluates the gates below.  --full runs the full
campaign; --no-rerun reuses out/summary_<scale>.json if present.

Gates (declared here before the first campaign read):

  G1 h5-instrument  Theorem 6.1's per-tile bound L_sel <= L_oracle +
                    2*eps_t holds on EVERY tile for BOTH selectors with
                    the exactly-instrumented eps_t (a violation means
                    the harness, not the theorem, is broken), sampled
                    rnr sub-block round-trips are byte-exact, and the
                    prefix selector's aggregate regret respects
                    sum(2 eps_t) + flag bits.
  G2 h5-selector    the prefix selector beats or matches the BEST FIXED
                    single mode on every composition (the operative
                    Type-III-A promise: composition without loss), and
                    its regret-vs-oracle ratio is < 20% (H5's
                    falsification threshold; the 5% target is a
                    finding, not a harness gate).
  G3 h6-anchor      estimator anchors, computed here at gate time:
                    (a) i.i.d.-uniform control: |MM-CMI| <= 0.02 bpb;
                    (b) synthetic coupled-nibble source with analytic
                    I(H;L|C) = 4 - H(L|H): match within 0.05 bpb.
  G4 h6-controls    campaign negative controls ctrl_urandom and
                    ctrl_zstd read |shuffle-null-debiased MM-CMI| <=
                    0.03 bpb under BOTH the prev1 and cm instruments
                    (the raw cm-instrument MM value on thin supports
                    retains ~0.36 bpb of residual plug-in inflation at
                    the 64 KiB CM slice -- measured, then removed by
                    the paired shuffle null; both are reported).
  G5 h7-coder       byte-exact round-trip at every ladder point; every
                    admitted phrase had Delta_L < 0 (Theorem 6.5 rule
                    applied consistently); rate < 8 bpb everywhere.
  G6 h7-oN          dictionary overhead fraction strictly decreases
                    along each corpus's N-ladder (Theorem 6.6b o(N)
                    shape at dev scale).
  G7 h7-index       Lemma 6.6a: ideal two-stage == ideal flat within
                    1e-8 relative; realized adaptive two-stage within
                    2% + 512 bits of realized flat; realized flat
                    within the declared adaptive-redundancy envelope
                    0.75*K*log2(2 + n/K) + 64 bits of the substream's
                    empirical entropy ("index cost == its entropy").
  G8 bitsback       Theorem 6.6 identity (6.4) pointwise-exact
                    (<= 1e-9) on the enumerable micro model for both
                    q = posterior and q = point mass; bits-back saving
                    >= 0; real-coder KL term >= 0.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time

import numpy as np

import common
import run_experiment

FAILURES = []


def gate(name: str, ok: bool, detail: str):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


def check_h5(h5):
    for r in h5:
        sp, sc = r["selector_prefix"], r["selector_causal"]
        gate("G1 h5-instrument",
             sp["bound_holds_all_tiles"] and sc["bound_holds_all_tiles"]
             and r["roundtrip_sampled_ok"]
             and sp["regret_le_2eps_aggregate"],
             f"comp s{r['seed_index']}: 2eps bound tile-exact both "
             f"selectors, round-trip ok, aggregate regret within "
             f"sum(2eps)+flags")
        best_fixed = min(r["fixed_mode_bpb"].values())
        reg = sp["regret_ratio"]["point"]
        gate("G2 h5-selector",
             sp["bpb"] <= best_fixed + 1e-12 and reg < 0.20,
             f"comp s{r['seed_index']}: selector {sp['bpb']:.4f} bpb vs "
             f"best fixed {best_fixed:.4f}; regret {100*reg:.2f}% "
             f"[{100*sp['regret_ratio']['lo']:.2f}, "
             f"{100*sp['regret_ratio']['hi']:.2f}]% (<20%)")


def check_h6_anchors():
    rng = np.random.default_rng(common.derive_seed("h6-anchor", 0))
    n = 1 << 18
    # (a) i.i.d. uniform: true CMI = 0
    x = rng.integers(0, 256, size=n).astype(np.int64)
    a = common.cmi_mm(x[:-1], x[1:], 256)["mm"]
    # (b) coupled nibbles independent of C: L = H w.p. 1/2 else uniform
    h = rng.integers(0, 16, size=n)
    coin = rng.random(n) < 0.5
    l = np.where(coin, h, rng.integers(0, 16, size=n))
    s = (h * 16 + l).astype(np.int64)
    c = rng.integers(0, 256, size=n).astype(np.int64)
    b = common.cmi_mm(c, s, 256)["mm"]
    p_eq = 0.5 + 0.5 / 16
    h_l_given_h = -(p_eq * math.log2(p_eq)
                    + 15 * (0.5 / 16) * math.log2(0.5 / 16))
    analytic = 4.0 - h_l_given_h
    # (c) thin-support regime (the CM-instrument failure mode): 64 Ki
    # samples, independent 256-valued context, WITH the coupled H-L
    # signal -- raw MM inflates; the conditional-permutation debias
    # must kill the bias while RETAINING the analytic signal
    n2 = 1 << 16
    c2 = rng.integers(0, 256, size=n2).astype(np.int64)
    h2 = rng.integers(0, 16, size=n2)
    l2 = np.where(rng.random(n2) < 0.5, h2, rng.integers(0, 16, size=n2))
    s2 = (h2 * 16 + l2).astype(np.int64)
    mm2 = common.cmi_mm(c2, s2, 256)["mm"]
    null2 = common.cmi_perm_null(c2, s2, 256, 12345)
    deb2 = mm2 - null2
    # and on a signal-free thin-support control the debias must read ~0
    s3 = rng.integers(0, 256, size=n2).astype(np.int64)
    deb3 = (common.cmi_mm(c2, s3, 256)["mm"]
            - common.cmi_perm_null(c2, s3, 256, 12346))
    gate("G3 h6-anchor",
         abs(a) <= 0.02 and abs(b - analytic) <= 0.05
         and abs(deb2 - analytic) <= 0.06 and abs(deb3) <= 0.03,
         f"uniform control {a:+.4f} (|.|<=0.02); coupled {b:.4f} vs "
         f"analytic {analytic:.4f} (tol 0.05); thin-support coupled raw "
         f"{mm2:+.4f} -> debiased {deb2:+.4f} vs {analytic:.4f} "
         f"(tol 0.06); thin-support null-free {deb3:+.4f} (|.|<=0.03)")


def check_h6(h6):
    for name in ("ctrl_urandom", "ctrl_zstd"):
        inst = h6[name]["instruments"]
        vals = {k: inst[k]["debiased"] for k in ("prev1", "cm")}
        raw = {k: inst[k]["mm"] for k in ("prev1", "cm")}
        gate("G4 h6-controls",
             all(abs(v) <= 0.03 for v in vals.values()),
             f"{name}: debiased " + ", ".join(f"{k}={v:+.4f}"
                                              for k, v in vals.items())
             + " (raw mm " + ", ".join(f"{k}={v:+.4f}"
                                       for k, v in raw.items()) + ")")


def check_h7(h7):
    for corpus, rows in h7["ladder"].items():
        ok5 = all(r["roundtrip_ok"] and r["admitted_all_negative_delta"]
                  and r["rate_bpb"] < 8.0 for r in rows)
        gate("G5 h7-coder", ok5,
             f"{corpus}: round-trip + admission rule + rate<8 at "
             f"N={[r['N'] for r in rows]}")
        fr = [r["overhead_fraction"] for r in rows]
        gate("G6 h7-oN", all(fr[i + 1] < fr[i] for i in range(len(fr) - 1)),
             f"{corpus}: overhead fractions "
             + " > ".join(f"{v:.4f}" for v in fr))
        for r in rows:
            neut = r["index_neutrality"]
            if neut.get("n_ids", 0) == 0:
                continue
            ideal_ok = (neut["ideal_identity_gap"]
                        <= 1e-8 * max(1.0, neut["ideal_flat_bits"]))
            two_ok = (abs(neut["realized_two_stage_bits"]
                          - neut["realized_flat_bits"])
                      <= 0.02 * neut["realized_flat_bits"] + 512)
            K = r["K"]
            envelope = 0.75 * K * math.log2(2 + neut["n_ids"] / max(1, K)) + 64
            emp_ok = (neut["realized_flat_bits"]
                      <= neut["empirical_entropy_bits"] + envelope)
            gate("G7 h7-index", ideal_ok and two_ok and emp_ok,
                 f"{corpus} N={r['N']}: ideal gap "
                 f"{neut['ideal_identity_gap']:.2e}; flat "
                 f"{neut['realized_flat_bits']:.0f}b vs two-stage "
                 f"{neut['realized_two_stage_bits']:.0f}b vs entropy "
                 f"{neut['empirical_entropy_bits']:.0f}b "
                 f"(envelope {envelope:.0f}b)")
    m = h7["bitsback_micro"]
    rc = h7["bitsback_real"]
    gate("G8 bitsback",
         m["max_pointwise_identity_gap"] <= 1e-9
         and m["bitsback_saving_bits"] >= -1e-12
         and m["data_mass_covered"] > 0.999
         and rc["kl_term_bits"] >= -1e-6,
         f"micro identity gap {m['max_pointwise_identity_gap']:.2e}, "
         f"saving {m['bitsback_saving_bits']:.4f} bits/string; real-coder "
         f"KL term {rc['kl_term_bpb']:.4f} bpb (>=0)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--no-rerun", action="store_true",
                    help="reuse out/summary_<scale>.json if present")
    args = ap.parse_args()
    scale = "full" if args.full else "reduced"
    t0 = time.time()
    path = common.OUT / f"summary_{scale}.json"
    if args.no_rerun and path.exists():
        with open(path) as f:
            summary = json.load(f)
        print(f"reusing {path}")
    else:
        print(f"running {scale} campaign ...")
        summary = run_experiment.run(scale)
    print(f"campaign ready ({summary['elapsed_s']:.1f}s); gates:")
    check_h5(summary["h5"])
    check_h6_anchors()
    check_h6(summary["h6"])
    check_h7(summary["h7"])
    verdict = "PASS" if not FAILURES else "FAIL"
    print(f"OVERALL {verdict} ({len(FAILURES)} failing gate(s); "
          f"{time.time() - t0:.1f}s wall)")
    if FAILURES:
        print("failing:", sorted(set(FAILURES)))
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    sys.exit(main())
