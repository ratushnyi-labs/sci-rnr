#!/usr/bin/env python3
"""Study (a): archive-overflow LDP (Theorem 7.31).

Law under test: for i.i.d. / finite-state Markov sources under a (possibly
mismatched) model, -(1/N) log2 P(L_N >= N(h+delta)) -> E(h+delta),
E(R) = sup_t [tR - Lambda(t)], with Lambda the base-2 scaled CGF (memoryless:
log2 sum P M^{-t}; Markov: log2 of the Perron root of the tilted transfer
matrix).  Empirical program:

  A1  i.i.d.: direct Monte Carlo tail estimates at 2 delta values x N-ladder;
      Bahadur-Rao-corrected slope fit of -log2 P vs N recovers E(R).
  A2  i.i.d.: tilted importance sampling agrees with direct MC where both are
      feasible and extends the same law N-fold deeper into the tail (the
      matching LOWER bound: the probability really is 2^{-NE(R)(1+o(1))},
      not merely bounded above by it).
  A3  Markov: same 2 delta x N-ladder program with the tilted-transfer-matrix
      rate function (direct MC + tilted-kernel importance sampling).
  A4  Markov ENDPOINT = MAX CYCLE MEAN: on a mismatched 2-state pair with
      r_+ = max cycle mean < ell_max = max single-step cost:
        (i)   lim_t Lambda'(t) equals the exhaustive-cycle max cycle mean;
        (ii)  the DP maximum path mean converges to r_+ (never to ell_max);
        (iii) overflow at R in (r_+, ell_max) is IMPOSSIBLE (DP-certified,
              0 hits in MC), while at R < r_+ the empirical rate is finite
              and tracks E(R).

Outputs: out/study_a_rates.csv, out/study_a_endpoint.csv, out/study_a.png
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (IIDSource, MarkovSource, OUT_DIR, PALETTE, report, save_csv,
                    save_fig, style_axes)


def bahadur_offset(N):
    # Bahadur-Rao prefactor: -log2 P = N E + (1/2) log2 N + O(1)
    return 0.5 * math.log2(N)


def run(full=False, outdir=OUT_DIR):
    results = []
    rng = np.random.default_rng(20260710)
    print("=" * 72)
    print("STUDY A -- archive-overflow LDP (Thm 7.31): empirical tail vs E(R)")
    print("=" * 72)

    R_mc = 800_000 if full else 400_000
    R_is = 100_000 if full else 40_000

    rows = []

    # ---------------------------------------------------------------- A1+A2 iid
    src = IIDSource([0.6, 0.3, 0.1])
    print(f"[iid] P0=(0.6,0.3,0.1)  h={src.h:.4f} bits/sym  V={src.V:.4f}")
    deltas_iid = [0.10, 0.18]
    ladders_mc = {0.10: [150, 300, 600, 1200], 0.18: [80, 160, 320]}
    ladders_is = {0.10: [150, 300, 600, 1200, 2400], 0.18: [80, 160, 320, 640, 1280]}
    if full:
        ladders_mc[0.10].append(2400)
        ladders_is[0.10].append(4800)
        ladders_is[0.18].append(2560)

    slopes_ok, agree_ok = True, True
    iid_curves = {}
    for delta in deltas_iid:
        R = src.h + delta
        E = src.legendre_E(R)
        mc_pts, is_pts = [], []
        for N in ladders_mc[delta]:
            L = src.sample_lengths(N, R_mc, rng)
            hits = int((L >= N * R - 1e-12).sum())
            if hits >= 20:
                lp = math.log2(hits / R_mc)
                mc_pts.append((N, lp, hits))
                rows.append(["iid", "MC", delta, N, lp, -lp / N, E, hits])
        for N in ladders_is[delta]:
            lp, hits, rse = src.is_tail_log2(N, N * R, R_is, rng)
            is_pts.append((N, lp, hits))
            rows.append(["iid", "IS", delta, N, lp, -lp / N, E, hits])
        # Bahadur-Rao-corrected slope of -log2 P vs N on the IS curve
        Ns = np.array([p[0] for p in is_pts], float)
        ys = np.array([-p[1] for p in is_pts], float)
        Xoff = ys - 0.5 * np.log2(Ns)
        slope = float(np.linalg.lstsq(np.vstack([Ns, np.ones_like(Ns)]).T, Xoff,
                                      rcond=None)[0][0])
        rel = abs(slope - E) / E
        print(f"  [iid] delta={delta}: E(R)={E:.5f}  BR-corrected slope={slope:.5f} "
              f"(rel err {rel:.3f})")
        slopes_ok &= rel < 0.10
        # MC vs IS agreement on common N
        for (N, lp_mc, _) in mc_pts:
            match = [p for p in is_pts if p[0] == N]
            if match:
                d = abs(lp_mc - match[0][1])
                print(f"  [iid] delta={delta} N={N}: log2P  MC={lp_mc:.3f}  "
                      f"IS={match[0][1]:.3f}  |diff|={d:.3f}")
                agree_ok &= d < 0.25
        iid_curves[delta] = (is_pts, E)
    report("A1", slopes_ok,
           "iid overflow: BR-corrected slope of -log2 P vs N matches E(R) within "
           "10% at both delta values", results)
    report("A2", agree_ok,
           "iid: tilted importance sampling agrees with direct MC (<0.25 bits in "
           "log2 P) and extends the law deeper into the tail", results)

    # ---------------------------------------------------------------- A3 Markov
    mk = MarkovSource([[0.85, 0.15], [0.15, 0.85]])
    print(f"[markov] 2-state p_stay=0.85  h={mk.h:.4f}  V={mk.V:.4f}")
    deltas_mk = [0.10, 0.18]
    lad_mc = {0.10: [100, 200, 400], 0.18: [60, 120, 240]}
    lad_is = {0.10: [100, 200, 400, 800], 0.18: [60, 120, 240, 480]}
    if full:
        lad_is[0.10].append(1600)
        lad_is[0.18].append(960)
    mk_slopes_ok, mk_agree_ok = True, True
    mk_curves = {}
    R_mc_mk = R_mc // 2
    for delta in deltas_mk:
        R = mk.h + delta
        E = mk.legendre_E(R)
        mc_pts, is_pts = [], []
        for N in lad_mc[delta]:
            L = mk.sample_lengths(N, R_mc_mk, rng)
            hits = int((L >= N * R - 1e-12).sum())
            if hits >= 20:
                lp = math.log2(hits / R_mc_mk)
                mc_pts.append((N, lp, hits))
                rows.append(["markov", "MC", delta, N, lp, -lp / N, E, hits])
        for N in lad_is[delta]:
            lp, hits, rse = mk.is_tail_log2(N, N * R, R_is // 2, rng)
            is_pts.append((N, lp, hits))
            rows.append(["markov", "IS", delta, N, lp, -lp / N, E, hits])
        Ns = np.array([p[0] for p in is_pts], float)
        ys = np.array([-p[1] for p in is_pts], float)
        Xoff = ys - 0.5 * np.log2(Ns)
        slope = float(np.linalg.lstsq(np.vstack([Ns, np.ones_like(Ns)]).T, Xoff,
                                      rcond=None)[0][0])
        rel = abs(slope - E) / E
        print(f"  [markov] delta={delta}: E(R)={E:.5f}  BR-corrected slope={slope:.5f} "
              f"(rel err {rel:.3f})")
        mk_slopes_ok &= rel < 0.12
        for (N, lp_mc, _) in mc_pts:
            match = [p for p in is_pts if p[0] == N]
            if match:
                d = abs(lp_mc - match[0][1])
                print(f"  [markov] delta={delta} N={N}: log2P  MC={lp_mc:.3f}  "
                      f"IS={match[0][1]:.3f}  |diff|={d:.3f}")
                mk_agree_ok &= d < 0.30
        mk_curves[delta] = (is_pts, E)
    report("A3", mk_slopes_ok and mk_agree_ok,
           "Markov overflow: tilted-transfer-matrix E(R) recovered from the "
           "empirical tail (slope within 12%; MC/IS agree)", results)

    # ---------------------------------------------------------------- A4 endpoint
    print("-" * 72)
    print("A4  Markov endpoint: r_+ = MAX CYCLE MEAN (not the single-step max)")
    # mismatched pair with a WIDE interior (h_cross, r_+) and r_+ < ell_max:
    # the expensive edge 1->2 (cost 2.737) lies on no sustainable high-mean cycle.
    ep = MarkovSource([[0.7, 0.3], [0.6, 0.4]], model=[[0.85, 0.15], [0.3, 0.7]])
    r_plus = ep.max_cycle_mean()
    ell_max = ep.edge_max()
    lam_inf = ep.lambda_prime(200.0)
    print(f"  cycles: max cycle mean r_+={r_plus:.4f}; single-step max "
          f"ell_max={ell_max:.4f}; lim Lambda'(t) (t=200) = {lam_inf:.4f}")
    ok_i = abs(lam_inf - r_plus) < 5e-3 and (ell_max - r_plus) > 0.3

    ep_rows = []
    dp_means = []
    for N in [25, 50, 100, 200, 400]:
        mean = ep.max_path_total(N) / N
        dp_means.append((N, mean))
        ep_rows.append([N, mean, r_plus, ell_max])
        print(f"  DP max path mean  N={N:4d}: {mean:.4f}  (r_+={r_plus:.4f})")
    gaps = [abs(m - r_plus) for _, m in dp_means]
    ok_ii = all(gaps[i + 1] <= gaps[i] + 1e-12 for i in range(len(gaps) - 1)) \
        and gaps[-1] < 0.01 and dp_means[-1][1] < ell_max - 0.3

    # (iii) impossibility above r_+, finite rate below
    N_imp = 200
    R_mid = 2.5  # in (r_+=2.237, ell_max=2.737)
    dp_max_mean = ep.max_path_total(N_imp) / N_imp
    L = ep.sample_lengths(N_imp, 200_000, rng)
    hits_mid = int((L >= N_imp * R_mid).sum())
    print(f"  R={R_mid} in (r_+, ell_max): DP max mean at N={N_imp} is "
          f"{dp_max_mean:.4f} < {R_mid} (event impossible); MC hits={hits_mid}/200000")
    R_below = 1.4  # in the interior (h_cross=1.073, r_+=2.237)
    E_below = ep.legendre_E(R_below, t_hi=200.0, n_grid=20001)
    emp_rates = []
    for N in ([40, 60, 90] if not full else [40, 60, 90, 130]):
        Lb = ep.sample_lengths(N, R_mc, rng)
        h_b = int((Lb >= N * R_below).sum())
        if h_b > 0:
            emp_rates.append((N, -math.log2(h_b / R_mc) / N))
            ep_rows.append([N, None, r_plus, ell_max])
    for N, r in emp_rates:
        print(f"  R={R_below} (< r_+): N={N:3d} empirical rate={r:.4f}  "
              f"(E({R_below})={E_below:.4f})")
    # empirical rate is finite, above E (prefactor), decreasing toward E
    rates = [r for _, r in emp_rates]
    ok_iii = (hits_mid == 0 and dp_max_mean < R_mid
              and len(rates) >= 3 and all(np.diff(rates) < 0)
              and rates[-1] > E_below * 0.9 and rates[-1] < E_below * 2.2)
    report("A4", ok_i and ok_ii and ok_iii,
           f"endpoint: lim Lambda' = max cycle mean {r_plus:.3f} != ell_max "
           f"{ell_max:.3f}; DP path mean -> r_+; overflow impossible on "
           f"(r_+, ell_max), finite rate below r_+", results)

    save_csv(os.path.join(outdir, "study_a_rates.csv"),
             ["source", "estimator", "delta", "N", "log2P", "rate", "E_theory",
              "hits"], rows)
    save_csv(os.path.join(outdir, "study_a_endpoint.csv"),
             ["N", "dp_max_path_mean", "r_plus", "ell_max"],
             [r for r in ep_rows if r[1] is not None])

    # ---------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    ax = axes[0]
    for i, delta in enumerate(deltas_iid):
        pts, E = iid_curves[delta]
        Ns = [p[0] for p in pts]
        ax.plot(Ns, [-p[1] / p[0] for p in pts], "o-", color=PALETTE[i],
                label=f"empirical rate, $\\delta$={delta}")
        ax.axhline(E, color=PALETTE[i], ls="--", lw=1,
                   label=f"E(h+{delta})={E:.4f}")
    ax.set_xscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$-\frac{1}{N}\log_2 \hat P(L_N \geq N(h+\delta))$")
    ax.set_title("iid overflow rate vs E(R)")
    ax.legend(fontsize=7)
    style_axes(ax)

    ax = axes[1]
    for i, delta in enumerate(deltas_mk):
        pts, E = mk_curves[delta]
        Ns = [p[0] for p in pts]
        ax.plot(Ns, [-p[1] / p[0] for p in pts], "s-", color=PALETTE[i],
                label=f"empirical rate, $\\delta$={delta}")
        ax.axhline(E, color=PALETTE[i], ls="--", lw=1,
                   label=f"E(h+{delta})={E:.4f}")
    ax.set_xscale("log")
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("N")
    ax.set_ylabel("empirical rate (bits/sym)")
    ax.set_title("Markov overflow rate vs tilted-matrix E(R)")
    ax.legend(fontsize=7)
    style_axes(ax)

    ax = axes[2]
    Ns = [n for n, _ in dp_means]
    ax.plot(Ns, [m for _, m in dp_means], "o-", color=PALETTE[0],
            label="DP max path mean")
    ax.axhline(r_plus, color=PALETTE[2], ls="--", label=f"$r_+$ = max cycle mean = {r_plus:.3f}")
    ax.axhline(ell_max, color=PALETTE[3], ls=":", label=f"$\\ell_{{max}}$ = {ell_max:.3f}")
    ax.set_xlabel("N")
    ax.set_ylabel("max achievable mean cost")
    ax.set_title("LDP endpoint = max cycle mean, not $\\ell_{max}$")
    ax.legend(fontsize=8)
    style_axes(ax)
    fig.suptitle("Study A -- archive-overflow LDP (Thm 7.31)", y=1.02)
    save_fig(fig, os.path.join(outdir, "study_a.png"))
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    res = run(full=args.full)
    ok = all(r[1] for r in res)
    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
