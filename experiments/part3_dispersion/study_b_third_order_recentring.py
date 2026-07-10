#!/usr/bin/env python3
"""
study_b_third_order_recentring.py
============================================================================
STUDY (b): third-order recentring on the Gray grid -- the 1/2 log n term.

Laws under test (Remarks 7.34n'/7.34b):
  * on Gray the d-tilted information is a deterministic shift of the
    surprisal, j_n = i_n - n c, c = lambda* D log2(e) + log2(gamma) = h2(D):
    the entire distribution recentres, so P(j_n <= nR) is governed by the
    lossless CLT (empirical CDF check);
  * the third-order 1/2 log n correction: the operational (ball-covering)
    achievability quantity is the ball log-probability
        G_n(x) = -log2 P_{Y*}(B_{nD}(x)) = j_n(x) - log2 S_n(x),
    whose excess over j_n is a LOCAL-limit scale factor: the ball mass is a
    windowed local mass of order n^{-1/2}, so
        E[G_n - j_n] = (1/2) log2 n + O(1),
    and the same drift appears at the (1-eps) quantile -- exactly the
    + (1/2) log n of L_lossy(n,eps,D) = nR(D) + sqrt(nV) Q^{-1}(eps)
    + (1/2) log n + O(1) (Remark 7.34n').

DESIGN (two Gray series with integer nD, all laws EXACT over the full 2^n
block law -- no Monte Carlo noise in the CDFs):
  series 1: p = 0.48, D = 0.25  (D_c = 0.3077), n = 8,12,16,20 (+24 full);
  series 2: p = 0.45, D = 0.20  (D_c = 0.2125), n = 10,15,20   (+25 full).

CHECKS:
  B0  Gray validity (P_Y >= 0) and the recentring identity j_n = i_n - nc
      deterministic to 1e-8 at every ladder point (7.34n' N6-style).
  B1  empirical CDF P(j_n <= nR(D)) within the lattice+CLT band of 1/2:
      |P - 1/2| <= (span + 1 - h2(p)) / sqrt(2 pi n V)  (i_n is lattice with
      span log2((1-p)/p); the band is the max lattice jump + mean offset).
  B2  mean recentred drift: fit E[G_n - j_n] = a log2 n + b over the ladder
      gives a in [0.30, 0.70] (predicted 1/2; sign + magnitude), drift
      increasing in n, on both series.
  B3  quantile drift at eps = 0.1: same fit and gate for
      q_{1-eps}(G_n) - q_{1-eps}(j_n).
  B4  bounded remainder: |E[G_n - j_n] - (1/2) log2(2 pi n D(1-D))| < 1 bit
      at every ladder point (the O(1) term stays O(1)).

Outputs: out/study_b_recentring.csv, out/study_b_recentring.png.

Caveats (honest): the ladder slopes sit slightly below 1/2 (0.40-0.44,
rising with n) -- the finite-n approach to the asymptotic 1/2 from below;
i_n is lattice for the BSMS, so the CDF checks carry the lattice band
(7.34n' notes the lattice caveat for the clean 1/2 log n statement).
"""

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def quantile_exact(vals, weights, q):
    order = np.argsort(vals)
    csum = np.cumsum(weights[order])
    k = np.searchsorted(csum, q)
    return float(vals[order][min(k, len(vals) - 1)])


def series_measurements(p, D, ns, eps=0.1):
    rows = []
    for n in ns:
        t = round(n * D)
        assert abs(t - n * D) < 1e-9, "ladder must keep nD integer"
        PY = cm.deconv_output(n, p, D)
        P, i_bits, _ = cm.bsms_block(n, p)
        j = cm.dtilted_bits(n, p, D, PY)
        G = cm.ball_logmass_bits(n, p, D, t, PY)
        resid = float(np.abs(j - (i_bits - n * cm.h2(D))).max())
        gap_mean = float(np.dot(P, G - j))
        qj = quantile_exact(j, P, 1 - eps)
        qG = quantile_exact(G, P, 1 - eps)
        R = cm.rate_gray_bits(p, D)
        Pleq = float(P[j <= n * R].sum())
        rows.append(dict(n=n, t=t, minPY=float(PY.min()), resid=resid,
                         gap_mean=gap_mean, gap_q=qG - qj, Pleq=Pleq))
    return rows


def fit_log2(ns, ys):
    xs = [math.log2(n) for n in ns]
    slope, b = np.polyfit(xs, ys, 1)
    return float(slope), float(b)


def run(full=False):
    ck = cm.Checker()
    print("=" * 78)
    print("STUDY (b): third-order recentring -- P(j_n <= nR) and the "
          "1/2 log n drift")
    print("=" * 78)

    series = [
        dict(p=0.48, D=0.25, ns=[8, 12, 16, 20] + ([24] if full else [])),
        dict(p=0.45, D=0.20, ns=[10, 15, 20] + ([25] if full else [])),
    ]
    eps = 0.1
    all_rows = []
    okB0 = okB1 = okB2 = okB3 = okB4 = True
    for s in series:
        p, D, ns = s["p"], s["D"], s["ns"]
        V = cm.V_lossless_bits(p)
        span = math.log2((1 - p) / p)
        print(f"  series p={p}, D={D} (D_c={cm.D_c(p):.4f}), "
              f"R(D)={cm.rate_gray_bits(p, D):.4f} b/sym, ladder n={ns}")
        rows = series_measurements(p, D, ns, eps)
        s["rows"] = rows
        for r in rows:
            band = (span + 1 - cm.h2(p)) / math.sqrt(2 * math.pi * r["n"] * V)
            okB0 &= r["minPY"] >= -1e-12 and r["resid"] < 1e-8
            okB1 &= abs(r["Pleq"] - 0.5) <= band
            rem = r["gap_mean"] - 0.5 * math.log2(
                2 * math.pi * r["n"] * D * (1 - D))
            okB4 &= abs(rem) < 1.0
            print(f"     n={r['n']:3d} t={r['t']}  resid={r['resid']:.1e}  "
                  f"P(j<=nR)={r['Pleq']:.4f} (band +/-{band:.4f})  "
                  f"E[G-j]={r['gap_mean']:.4f}  q-drift={r['gap_q']:.4f}  "
                  f"O(1)-rem={rem:+.3f}")
            all_rows.append([p, D, r["n"], r["t"], r["resid"], r["Pleq"],
                             band, r["gap_mean"], r["gap_q"], rem])
        a_m, b_m = fit_log2(ns, [r["gap_mean"] for r in rows])
        a_q, b_q = fit_log2(ns, [r["gap_q"] for r in rows])
        s["fit"] = (a_m, b_m)
        inc = all(rows[i + 1]["gap_mean"] > rows[i]["gap_mean"]
                  for i in range(len(rows) - 1))
        print(f"     fit E[G-j] = {a_m:.3f} log2(n) + {b_m:+.3f}   "
              f"[predicted slope 0.5]   monotone: {inc}")
        print(f"     fit q-drift = {a_q:.3f} log2(n) + {b_q:+.3f}")
        okB2 &= 0.30 <= a_m <= 0.70 and inc
        okB3 &= 0.30 <= a_q <= 0.70

    ck.rep("B0 Gray validity + recentring identity j_n = i_n - nc (<1e-8)", okB0)
    ck.rep("B1 P(j_n <= nR) within the lattice+CLT band of 1/2", okB1)
    ck.rep("B2 mean drift E[G_n - j_n]: slope in [0.30,0.70] of log2 n, "
           "increasing", okB2)
    ck.rep(f"B3 (1-{eps}) quantile drift: slope in [0.30,0.70] of log2 n", okB3)
    ck.rep("B4 remainder E[G-j] - (1/2)log2(2 pi n D(1-D)) bounded (<1 bit)",
           okB4)

    os.makedirs(OUT, exist_ok=True)
    cm.write_csv(os.path.join(OUT, "study_b_recentring.csv"),
                 ["p", "D", "n", "t", "identity_resid", "P_j_le_nR",
                  "clt_band", "gap_mean_bits", "gap_quantile_bits",
                  "O1_remainder_bits"], all_rows)
    _plot(series)
    print("  artifacts: out/study_b_recentring.csv, out/study_b_recentring.png")
    return ck


def _plot(series):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for k, s in enumerate(series):
        rows = s["rows"]
        ns = [r["n"] for r in rows]
        gm = [r["gap_mean"] for r in rows]
        gq = [r["gap_q"] for r in rows]
        a, b = s["fit"]
        lbl = f"p={s['p']}, D={s['D']}"
        ax1.plot(ns, gm, "o-", color=f"C{k}", label=f"E[G-j], {lbl}")
        ax1.plot(ns, gq, "s--", color=f"C{k}", alpha=0.6,
                 label=f"q-drift, {lbl}")
        nn = np.linspace(min(ns), max(ns), 100)
        ax1.plot(nn, 0.5 * np.log2(nn) + (gm[-1] - 0.5 * math.log2(ns[-1])),
                 ":", color=f"C{k}", lw=1,
                 label=f"(1/2)log2 n + const")
        ax2.plot(ns, [r["Pleq"] for r in rows], "o-", color=f"C{k}",
                 label=f"P(j_n<=nR), {lbl}")
        V = cm.V_lossless_bits(s["p"])
        span = math.log2((1 - s["p"]) / s["p"])
        band = [(span + 1 - cm.h2(s["p"])) / math.sqrt(2 * math.pi * n * V)
                for n in ns]
        ax2.fill_between(ns, [0.5 - b_ for b_ in band],
                         [0.5 + b_ for b_ in band], color=f"C{k}", alpha=0.12)
    ax1.set_xscale("log")
    ax1.set_xlabel("n")
    ax1.set_ylabel("bits")
    ax1.set_title("third-order drift: ball log-mass minus d-tilted info")
    ax1.legend(fontsize=7)
    ax2.axhline(0.5, color="k", lw=1, ls="--")
    ax2.set_xlabel("n")
    ax2.set_ylabel("P(j_n <= nR(D))")
    ax2.set_title("recentred CDF at the rate point (lattice+CLT band)")
    ax2.legend(fontsize=7)
    fig.suptitle("Study (b): the 1/2 log n recentring on the Gray grid")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "study_b_recentring.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    ck = run(full=args.full)
    print("=" * 78)
    print(f"OVERALL (study b) -> {'PASS' if ck.ok else 'FAIL'}")
    sys.exit(0 if ck.ok else 1)
