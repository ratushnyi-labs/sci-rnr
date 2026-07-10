#!/usr/bin/env python3
"""
study_c_beyond_gray_law.py
============================================================================
STUDY (c): the beyond-Gray finite-n threshold law D_c^(n) - D_c = K(p)^2/n^2.

Laws under test (Remarks 7.34o/7.34o'):
  * D_c^(n) -- the smallest D at which the n-block deconvolution
    (K^{-1})^{(x)n} P_X acquires a negative entry, on the ALTERNATING
    binding word -- decreases monotonically to D_c = inf_n D_c^(n);
  * the finite-n law:  D_c^(n) - D_c = K(p)^2/n^2 (1+o(1)),
        K(p) = 2 pi / c(p) = pi p (1-2p)^{1/4} / (2 (1-p)^{3/2}),
    from the eigenvalue-collision phase slope phi_2(D) ~ c(p) sqrt(D-D_c):
    the n-block mass first goes negative when n phi_2 ~ 2 pi.

METHOD:
  * main measurement: the alternating-word mass P_Y(0101...) via exact
    2x2 transfer products in 40-digit arithmetic (no float64 roundoff bias
    -- the full-WHT float64 bisection is 2-3% biased at n >= 12 for small
    p, per the committed probe's own control note); first-sign-change scan
    + 80-step bisection per (p, n);
  * cross-check: the full 2^n WHT deconvolution min-entry bisection at
    n = 10 (verifies the binding word IS the alternating one).

CHECKS (p in {0.1, 0.2, 0.3}; even ladder n = 8..24 reduced, ..32 full):
  C1  binding word: transfer threshold == full-WHT min-entry threshold at
      n = 10 (diff < 2e-6); the n = 12 comparison is reported descriptively
      (float64 bias grows at small p -- honest caveat, matches the
      committed probe's 80-bit control note).
  C2  monotone approach: D_c^(n) strictly decreasing in n toward D_c.
  C3  rate: log-log fit of D_c^(n) - D_c vs n has slope -2 (+/- 0.10).
  C4  constant: n^2 (D_c^(n) - D_c) at the last rung within 10% of the
      closed-form K(p)^2 at each p (fitted-vs-predicted table printed,
      including a 1/n-extrapolated intercept).

Outputs: out/study_c_beyond_gray.csv, out/study_c_beyond_gray.png.
"""

import argparse
import math
import os
import sys

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
mp.mp.dps = 40


def alt_mass(n, p, D):
    """P_Y(alternating word 0101...) via exact 2x2 transfer products."""
    D = mp.mpf(D)
    p = mp.mpf(p)
    a0 = 1 - 2 * D
    Ki = [[(1 - D) / a0, -D / a0], [-D / a0, (1 - D) / a0]]
    T = [[1 - p, p], [p, 1 - p]]
    y = [i % 2 for i in range(n)]
    v = [mp.mpf('0.5') * Ki[y[0]][0], mp.mpf('0.5') * Ki[y[0]][1]]
    for i in range(1, n):
        w0 = (v[0] * T[0][0] + v[1] * T[1][0]) * Ki[y[i]][0]
        w1 = (v[0] * T[0][1] + v[1] * T[1][1]) * Ki[y[i]][1]
        v = [w0, w1]
    return v[0] + v[1]


def Dc_n_transfer(n, p):
    """Per-n threshold via first sign change of the alternating-word mass."""
    Dc = cm.D_c(p)
    K2 = cm.K_beyond_gray(p) ** 2
    step = mp.mpf(K2) / n ** 2 / 8
    lo = mp.mpf(Dc)
    hi = None
    d = lo + step
    for _ in range(4000):
        if alt_mass(n, p, d) < 0:
            hi = d
            break
        lo = d
        d += step
    if hi is None:
        raise RuntimeError(f"no sign change found for p={p}, n={n}")
    for _ in range(80):
        mid = (lo + hi) / 2
        if alt_mass(n, p, mid) > 0:
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2)


def Dc_n_wht(n, p):
    """Per-n threshold via full 2^n deconvolution min-entry (float64)."""
    lo, hi = cm.D_c(p), 0.5 - 1e-9
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if cm.deconv_output(n, p, mid).min() >= -1e-14:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def run(full=False):
    ck = cm.Checker()
    print("=" * 78)
    print("STUDY (c): beyond-Gray finite-n threshold law "
          "D_c^(n) - D_c = K(p)^2 / n^2")
    print("=" * 78)

    ps = [0.1, 0.2, 0.3] + ([0.15, 0.25] if full else [])
    ladder = list(range(8, 26, 4)) + ([28, 32] if full else [])
    okC1 = okC2 = okC3 = okC4 = True
    rows = []
    results = {}
    for p in sorted(ps):
        Dc = cm.D_c(p)
        K2 = cm.K_beyond_gray(p) ** 2
        # cross-check vs full-WHT at n = 10 (+ descriptive n = 12)
        t10, w10 = Dc_n_transfer(10, p), Dc_n_wht(10, p)
        t12, w12 = Dc_n_transfer(12, p), Dc_n_wht(12, p)
        okC1 &= abs(t10 - w10) < 2e-6
        print(f"  p={p}: D_c={Dc:.7f}, K(p)^2={K2:.6f}")
        print(f"     binding-word check n=10: transfer={t10:.9f} "
              f"WHT={w10:.9f} diff={abs(t10 - w10):.1e}")
        print(f"     (n=12 descriptive:      transfer={t12:.9f} "
              f"WHT={w12:.9f} diff={abs(w12 - t12):.1e} -- float64 bias)")
        deltas = []
        for n in ladder:
            dn = Dc_n_transfer(n, p)
            delta = dn - Dc
            deltas.append((n, dn, delta))
            print(f"     n={n:2d}: D_c^(n)={dn:.9f}   n^2 (D_c^(n)-D_c) = "
                  f"{delta * n * n:.6f}")
            rows.append([p, n, dn, delta, delta * n * n, K2])
        results[p] = deltas
        okC2 &= all(deltas[i + 1][1] < deltas[i][1]
                    for i in range(len(deltas) - 1))
        # rate fit
        xs = [math.log(n) for n, _, _ in deltas]
        ys = [math.log(d) for _, _, d in deltas]
        slope, _ = np.polyfit(xs, ys, 1)
        # constant: last rung + 1/n-extrapolated intercept
        yn = [d * n * n for n, _, d in deltas]
        inv = [1.0 / n for n, _, _ in deltas]
        beta, K2fit = np.polyfit(inv, yn, 1)
        last = yn[-1]
        print(f"     rate fit: slope={slope:.4f} (predicted -2);  constant: "
              f"last rung {last:.6f}, 1/n-extrapolated {K2fit:.6f}, "
              f"predicted {K2:.6f}  (last/pred = {last / K2:.4f})")
        okC3 &= abs(slope + 2.0) < 0.10
        okC4 &= abs(last / K2 - 1.0) < 0.10

    ck.rep("C1 binding word is the alternating word (transfer == WHT, n=10)",
           okC1)
    ck.rep("C2 D_c^(n) decreases monotonically toward D_c", okC2)
    ck.rep("C3 log-log rate: slope = -2 (+/- 0.10) at every p", okC3)
    ck.rep("C4 constant: n^2 (D_c^(n)-D_c) within 10% of K(p)^2 at every p",
           okC4)

    os.makedirs(OUT, exist_ok=True)
    cm.write_csv(os.path.join(OUT, "study_c_beyond_gray.csv"),
                 ["p", "n", "Dc_n", "Dc_n_minus_Dc", "n2_scaled",
                  "K2_predicted"], rows)
    _plot(results)
    print("  artifacts: out/study_c_beyond_gray.csv, "
          "out/study_c_beyond_gray.png")
    return ck


def _plot(results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for k, (p, deltas) in enumerate(sorted(results.items())):
        K2 = cm.K_beyond_gray(p) ** 2
        ns = np.array([n for n, _, _ in deltas], dtype=float)
        ds = np.array([d for _, _, d in deltas])
        ax1.loglog(ns, ds, "o-", color=f"C{k}", label=f"p={p}")
        ax1.loglog(ns, K2 / ns ** 2, "--", color=f"C{k}", lw=1,
                   label=f"K({p})^2/n^2")
        ax2.plot(ns, ds * ns ** 2, "o-", color=f"C{k}", label=f"p={p}")
        ax2.axhline(K2, color=f"C{k}", ls="--", lw=1)
    ax1.set_xlabel("n")
    ax1.set_ylabel("D_c^(n) - D_c")
    ax1.set_title("finite-n threshold excess vs the 1/n^2 law")
    ax1.legend(fontsize=7)
    ax2.set_xlabel("n")
    ax2.set_ylabel("n^2 (D_c^(n) - D_c)")
    ax2.set_title("scaled excess vs closed-form K(p)^2 (dashed)")
    ax2.legend(fontsize=7)
    fig.suptitle("Study (c): beyond-Gray finite-n law, "
                 "K(p) = pi p (1-2p)^{1/4} / (2(1-p)^{3/2})")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "study_c_beyond_gray.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    ck = run(full=args.full)
    print("=" * 78)
    print(f"OVERALL (study c) -> {'PASS' if ck.ok else 'FAIL'}")
    sys.exit(0 if ck.ok else 1)
