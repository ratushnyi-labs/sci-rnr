#!/usr/bin/env python3
"""Study (c): moderate deviations and the RA-penalty fade (Theorem 7.33).

Setting: i.i.d. source, ideal predictor; budget N h + a_N sqrt(N V) with
a_N = N^{1/4} (-> oo, = o(sqrt N)).  The sub-block random-access overhead
(N/K) delta at K = c sqrt(N) equals a CONSTANT standardized shift
beta = delta/(c sqrt V):

  P(L >= B)      = exp(-(1+o(1)) a_N^2 / 2)                (MDP)
  P(L + RA >= B) = exp(-(1+o(1)) (a_N - beta)^2 / 2)        (shifted MDP)

The FADE: the empirical exponent ratio  ln P_RA / ln P -> 1 (RA invisible at
the exponent scale, Thm 7.31), while (1/a_N) ln (P_RA / P) -> beta > 0 (the
overflow probability stays inflated by exp((1+o(1)) beta a_N)); at the CLT end
(a_N = O(1)) the exact ratio is Q(a_N - beta)/Q(a_N) (Thm 7.29).

Empirical program (tilted importance sampling for the deep tails; direct MC
cross-check where reachable):
  C1  IS estimates agree with direct MC at the shallow end (N=256).
  C2  MDP rate: -ln P / (a_N^2/2) -> 1 along a_N = N^{1/4}.
  C3  the fade: empirical exponent ratio ln P_RA / ln P increases
      monotonically toward 1 along the ladder.
  C4  the penalty stays in the probability: empirical (1/a_N) ln(P_RA/P)
      increases toward beta (never above), the Cramer-Petrov O(N^{-1/4})
      approach.
  C5  CLT end: at fixed a = 1.5 the direct-MC tail ratio matches
      Q(a - beta)/Q(a).

Outputs: out/study_c_fade.csv, out/study_c.png
"""
import argparse
import math
import os
import sys

import numpy as np
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (IIDSource, LN2, OUT_DIR, PALETTE, report, save_csv,
                    save_fig, style_axes)


def run(full=False, outdir=OUT_DIR):
    results = []
    rng = np.random.default_rng(20260712)
    print("=" * 72)
    print("STUDY C -- moderate deviations + RA fade (Thm 7.33)")
    print("=" * 72)

    src = IIDSource([0.6, 0.3, 0.1])
    h, V = src.h, src.V
    delta, c = 0.25, 1.0
    beta = delta / (c * math.sqrt(V))
    print(f"P0=(0.6,0.3,0.1)  h={h:.4f}  V={V:.4f}  delta={delta}  c={c}  "
          f"beta=delta/(c sqrt V)={beta:.4f}")

    ladder = [256, 1024, 4096, 16384]
    if full:
        ladder.append(65536)
    R_is = 400_000 if full else 200_000

    # ---------------------------------------------------------------- C1
    N0 = 256
    a0 = N0 ** 0.25
    B0 = N0 * h + a0 * math.sqrt(N0 * V)
    ra0 = math.sqrt(N0) * delta / c
    R_mc = 8_000_000 if full else 4_000_000
    hits1 = 0
    hits2 = 0
    for _ in range(4):  # chunked to bound memory
        L = src.sample_lengths(N0, R_mc // 4, rng)
        hits1 += int((L >= B0 - 1e-9).sum())
        hits2 += int((L >= B0 - ra0 - 1e-9).sum())
    p_mc = hits1 / R_mc
    p_mc_ra = hits2 / R_mc
    lp_is, _, _ = src.is_tail_log2(N0, B0, R_is, rng)
    lp_is_ra, _, _ = src.is_tail_log2(N0, B0 - ra0, R_is, rng)
    d1 = abs(math.log2(p_mc) - lp_is)
    d2 = abs(math.log2(p_mc_ra) - lp_is_ra)
    # MC noise floor: SE of log2 p_hat ~ 1/(sqrt(hits) ln 2); gate at ~2.5 SE
    tol1 = max(0.12, 2.5 / (math.sqrt(max(hits1, 1)) * LN2))
    tol2 = max(0.12, 2.5 / (math.sqrt(max(hits2, 1)) * LN2))
    print(f"  N={N0}: log2P  MC={math.log2(p_mc):.3f} IS={lp_is:.3f} |d|={d1:.3f} "
          f"(tol {tol1:.3f}, hits={hits1}); log2P_RA MC={math.log2(p_mc_ra):.3f} "
          f"IS={lp_is_ra:.3f} |d|={d2:.3f} (tol {tol2:.3f}, hits={hits2})")
    report("C1", d1 < tol1 and d2 < tol2 and hits1 > 50,
           "tilted-IS tail estimates match direct MC at the shallow end "
           "(within 2.5 MC standard errors)", results)

    # ------------------------------------------------------- ladder estimates
    rows = []
    mdp_rates, exp_ratios, lr_over_a = [], [], []
    for N in ladder:
        aN = N ** 0.25
        B = N * h + aN * math.sqrt(N * V)
        ra = math.sqrt(N) * delta / c
        lp, hits, rse = src.is_tail_log2(N, B, R_is, rng)
        lp_ra, hits2, rse2 = src.is_tail_log2(N, B - ra, R_is, rng)
        lnP = lp * LN2
        lnP_ra = lp_ra * LN2
        mdp = -lnP / (aN ** 2 / 2)
        ratio = lnP_ra / lnP
        lr = (lnP_ra - lnP) / aN
        mdp_rates.append((N, mdp))
        exp_ratios.append((N, ratio))
        lr_over_a.append((N, lr))
        pred_ratio = (aN - beta) ** 2 / aN ** 2
        rows.append([N, aN, lnP, lnP_ra, mdp, ratio, pred_ratio, lr, beta,
                     rse, rse2])
        print(f"  N={N:6d} a_N={aN:6.2f}: -lnP/(a^2/2)={mdp:.4f}  "
              f"lnP_RA/lnP={ratio:.4f} (pred {pred_ratio:.4f})  "
              f"(1/a)ln(P_RA/P)={lr:.4f} (->beta={beta:.4f})")

    # ---------------------------------------------------------------- C2
    m = [r for _, r in mdp_rates]
    c2 = all(m[i + 1] < m[i] for i in range(len(m) - 1)) and abs(m[-1] - 1) < 0.12
    report("C2", c2,
           f"MDP rate -lnP/(a_N^2/2) decreases toward 1 along a_N=N^(1/4) "
           f"({m[0]:.3f}->{m[-1]:.3f})", results)

    # ---------------------------------------------------------------- C3
    r = [x for _, x in exp_ratios]
    c3 = (all(r[i + 1] > r[i] for i in range(len(r) - 1))
          and r[-1] > 0.92 and r[0] < 0.90)
    report("C3", c3,
           f"exponent ratio lnP_RA/lnP rises {r[0]:.3f} -> {r[-1]:.3f} toward 1: "
           f"the RA penalty FADES from the exponent (-> Thm 7.31)", results)

    # ---------------------------------------------------------------- C4
    g = [x for _, x in lr_over_a]
    c4 = (all(g[i + 1] > g[i] for i in range(len(g) - 1))
          and all(x < beta + 0.01 for x in g)
          and (beta - g[-1]) < (beta - g[0]))
    report("C4", c4,
           f"(1/a_N) ln(P_RA/P) rises {g[0]:.3f} -> {g[-1]:.3f} toward "
           f"beta={beta:.3f} (from below, Cramer-Petrov O(N^-1/4)): the penalty "
           f"STAYS in the probability", results)

    # ---------------------------------------------------------------- C5
    a_fix = 1.5
    N5 = 4096
    B5 = N5 * h + a_fix * math.sqrt(N5 * V)
    ra5 = math.sqrt(N5) * delta / c
    L5 = src.sample_lengths(N5, R_mc, rng)
    p5 = float((L5 >= B5 - 1e-9).mean())
    p5r = float((L5 >= B5 - ra5 - 1e-9).mean())
    emp_ratio = p5r / p5
    pred = norm.sf(a_fix - beta) / norm.sf(a_fix)
    print(f"  CLT end: a={a_fix}, N={N5}: P_RA/P = {emp_ratio:.4f}  "
          f"Q(a-beta)/Q(a) = {pred:.4f}")
    c5 = abs(emp_ratio / pred - 1) < 0.10
    report("C5", c5,
           f"CLT-end tail ratio matches Q(a-beta)/Q(a)={pred:.3f} within 10% "
           f"(-> Thm 7.29: RA fully visible at a_N=O(1))", results)

    save_csv(os.path.join(outdir, "study_c_fade.csv"),
             ["N", "a_N", "lnP", "lnP_RA", "mdp_rate", "exp_ratio",
              "exp_ratio_pred", "lr_over_aN", "beta", "rse_IS", "rse_IS_RA"],
             rows)

    # ---------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    Ns = [n for n, _ in mdp_rates]
    ax = axes[0]
    ax.plot(Ns, m, "o-", color=PALETTE[0], label="empirical")
    ax.axhline(1.0, color=PALETTE[3], ls="--", label="MDP limit 1")
    ax.set_xscale("log")
    ax.set_xlabel("N  ($a_N = N^{1/4}$)")
    ax.set_ylabel(r"$-\ln P \,/\, (a_N^2/2)$")
    ax.set_title("MDP rate")
    ax.legend(fontsize=8)
    style_axes(ax)

    ax = axes[1]
    ax.plot(Ns, r, "o-", color=PALETTE[0], label="empirical")
    ax.plot(Ns, [(n ** 0.25 - beta) ** 2 / math.sqrt(n) for n in Ns], "--",
            color=PALETTE[1], label=r"$(a_N-\beta)^2/a_N^2$")
    ax.axhline(1.0, color=PALETTE[3], ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$\ln P_{RA} / \ln P$")
    ax.set_title("exponent ratio -> 1: RA fades\nfrom the exponent")
    ax.legend(fontsize=8)
    style_axes(ax)

    ax = axes[2]
    ax.plot(Ns, g, "o-", color=PALETTE[0], label="empirical")
    ax.axhline(beta, color=PALETTE[3], ls="--", label=f"$\\beta$={beta:.3f}")
    ax.set_xscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$(1/a_N)\,\ln(P_{RA}/P)$")
    ax.set_title("...but stays in the probability:\n" r"$\to \beta$ from below")
    ax.legend(fontsize=8)
    style_axes(ax)
    fig.suptitle("Study C -- moderate-deviation RA fade (Thm 7.33)", y=1.02)
    save_fig(fig, os.path.join(outdir, "study_c.png"))
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    res = run(full=args.full)
    ok = all(x[1] for x in res)
    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
