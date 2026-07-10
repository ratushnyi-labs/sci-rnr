#!/usr/bin/env python3
"""Study (d): worst-case query cost EVT / read-buffer provisioning (Theorem 7.35).

Law under test: with m = N/K blocks of K symbols, L_i the per-block codelength,
h = entropy rate, V = varentropy rate:
  (a) leading order:  max_i L_i = K h + sqrt(2 K V ln m) (1 + o_P(1));
  (b) refined Gumbel (valid iff ln m = o(K^{1/3}), i.e. K >> (ln m)^3):
      a_m ((max - Kh)/sqrt(KV) - b_m) => Gumbel,  E[max] = Kh +
      sqrt(KV)(b_m + gamma/a_m) + o(.),  a_m = sqrt(2 ln m),
      b_m = a_m - (ln ln m + ln 4 pi)/(2 a_m).

Empirical program:
  D1  i.i.d.: standardized mean max tracks b_m + gamma/a_m across an m-ladder
      in the valid regime K >> (ln m)^3 (per-block multinomial sampling).
  D2  i.i.d.: leading-order ratio (max - Kh)/sqrt(2KV ln m) rises toward 1
      from below (the sqrt(2 ln m) provisioning law).
  D3  Gumbel SHAPE: the KS distance between a_m(Z-b_m) and the Gumbel law is
      small at the valid corner (distributional, not just the mean).
  D4  Markov source: same law with the fundamental-matrix (h, V) --
      cross-checked against empirical block moments.
  D5  VALIDITY BOUNDARY: on a skewed source at fixed m, the refined-Gumbel
      error grows monotonically as K falls through (ln m)^3, tracking the
      Cramer-Petrov correction x^3/sqrt(K) -- the ln m = o(K^{1/3}) boundary
      is real and simulation-visible.

Outputs: out/study_d_evt.csv, out/study_d_boundary.csv, out/study_d.png
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (EULER, IIDSource, MarkovSource, OUT_DIR, PALETTE,
                    gumbel_norm, report, save_csv, save_fig, style_axes)


def block_maxes_iid(src, K, m, R, rng):
    """R realizations of max over m iid blocks of K symbols.  Per-block
    codelength depends only on the block's counts -> sample multinomial counts
    directly (O(R m A) instead of O(R m K))."""
    counts = rng.multinomial(K, src.p, size=(R, m))
    L = counts @ src.l
    return L.max(axis=1)


def block_maxes_markov(mk, K, m, R, rng):
    """R independent contiguous chains of m*K steps; per-block sums of the
    per-step cost; max over the m (dependent) blocks.  Streams block sums
    (O(R) memory) instead of materializing the R x (m K) cost matrix."""
    cur = rng.choice(mk.n, size=R, p=mk.pi)
    best = np.full(R, -np.inf)
    blk = np.zeros(R)
    for t in range(m * K):
        u = rng.random(R)
        nxt = np.minimum((u[:, None] > mk.cdf[cur]).sum(axis=1), mk.n - 1)
        blk += mk.g[cur, nxt]
        cur = nxt
        if (t + 1) % K == 0:
            np.maximum(best, blk, out=best)
            blk[:] = 0.0
    return best


def gumbel_ks(z_std):
    """KS distance of standardized sample vs the standard Gumbel CDF."""
    x = np.sort(z_std)
    cdf = np.exp(-np.exp(-x))
    n = len(x)
    emp_hi = np.arange(1, n + 1) / n
    emp_lo = np.arange(0, n) / n
    return float(max(np.max(np.abs(emp_hi - cdf)), np.max(np.abs(emp_lo - cdf))))


def run(full=False, outdir=OUT_DIR):
    results = []
    rng = np.random.default_rng(20260713)
    print("=" * 72)
    print("STUDY D -- worst-case query EVT / read buffer (Thm 7.35)")
    print("=" * 72)

    src = IIDSource([0.5, 0.25, 0.15, 0.10])
    h, V = src.h, src.V
    print(f"iid p=(0.5,0.25,0.15,0.10)  h={h:.4f}  V={V:.4f}")

    R = 4000 if full else 1500
    K = 700
    ms = [16, 64, 256, 1024]
    rows = []

    # ---------------------------------------------------------------- D1 + D2
    refined_errs, lead_ratios = [], []
    z_samples_valid = None
    for m in ms:
        mx = block_maxes_iid(src, K, m, R, rng)
        z = (mx.mean() - K * h) / math.sqrt(K * V)
        a, b = gumbel_norm(m)
        pred = b + EULER / a
        lead = (mx.mean() - K * h) / math.sqrt(2 * K * V * math.log(m))
        rel = abs(z - pred) / pred
        refined_errs.append(rel)
        lead_ratios.append(lead)
        rows.append(["iid", K, m, mx.mean(), K * h, z, pred, rel, lead])
        print(f"  iid K={K} m={m:5d}: std max z={z:.4f} refined pred={pred:.4f} "
              f"rel err={rel:.3f}  leading ratio={lead:.4f}")
        if m == 1024:
            z_samples_valid = (mx - K * h) / math.sqrt(K * V)
    d1 = all(e < 0.06 for e in refined_errs)
    report("D1", d1,
           f"iid E[max] tracks Kh + sqrt(KV)(b_m + gamma/a_m) within 6% "
           f"(std units) across m={ms} at K={K} (K >> (ln m)^3)", results)
    d2 = (all(lead_ratios[i] < lead_ratios[i + 1] for i in range(len(ms) - 1))
          and lead_ratios[-1] > 0.87 and lead_ratios[-1] <= 1.02)
    report("D2", d2,
           f"leading-order ratio (max-Kh)/sqrt(2KV ln m) rises "
           f"{lead_ratios[0]:.3f}->{lead_ratios[-1]:.3f} toward 1 from below",
           results)

    # ---------------------------------------------------------------- D3 shape
    a, b = gumbel_norm(1024)
    zstd = a * (z_samples_valid / 1.0 - b)
    ks = gumbel_ks(zstd)
    # Gumbel convergence for Gaussian-type maxima is O(1/ln m): at m=1024 the
    # residual shape error is a few %, so gate the KS distance, not a p-value.
    d3 = ks < 0.06
    report("D3", d3,
           f"KS distance of a_m(Z-b_m) to the Gumbel law = {ks:.4f} < 0.06 at "
           f"m=1024, K={K} (distributional fit, O(1/ln m) residual)", results)

    # ---------------------------------------------------------------- D4 Markov
    mk = MarkovSource([[0.7, 0.2, 0.1], [0.1, 0.6, 0.3], [0.2, 0.3, 0.5]])
    print(f"Markov 3-state: h={mk.h:.4f}  V(fundamental matrix)={mk.V:.4f}")
    Km = 300
    Rm = 1200 if full else 500
    # cross-check the closed-form (h, V) against empirical block moments
    bs = mk.sample_cost_paths(400, 4000 if full else 2500, rng).sum(axis=1)
    h_emp, V_emp = float(bs.mean()) / 400, float(bs.var(ddof=1)) / 400
    print(f"  cross-check @K=400: h_emp={h_emp:.4f} (vs {mk.h:.4f}), "
          f"V_emp={V_emp:.4f} (vs {mk.V:.4f})")
    mk_ok = abs(h_emp / mk.h - 1) < 0.01 and abs(V_emp / mk.V - 1) < 0.10
    for m in ([64, 256, 1024] if full else [64, 256]):
        mx = block_maxes_markov(mk, Km, m, Rm, rng)
        z = (mx.mean() - Km * mk.h) / math.sqrt(Km * mk.V)
        a, b = gumbel_norm(m)
        pred = b + EULER / a
        rel = abs(z - pred) / pred
        rows.append(["markov", Km, m, mx.mean(), Km * mk.h, z, pred, rel,
                     (mx.mean() - Km * mk.h) / math.sqrt(2 * Km * mk.V * math.log(m))])
        print(f"  markov K={Km} m={m:5d}: std max z={z:.4f} pred={pred:.4f} "
              f"rel err={rel:.3f}")
        mk_ok &= rel < 0.10
    report("D4", mk_ok,
           "Markov blocks obey the same EVT law with the fundamental-matrix "
           "(h,V) (cross-checked against empirical block moments)", results)

    # ---------------------------------------------------------------- D5 boundary
    print("-" * 72)
    print("D5  validity boundary: refined Gumbel needs K >> (ln m)^3 "
          "(skewed source, fixed m)")
    src5 = IIDSource([0.85, 0.08, 0.05, 0.02])
    m5 = 1024
    R5 = 4000 if full else 2000
    a, b = gumbel_norm(m5)
    pred = b + EULER / a
    lnm3 = math.log(m5) ** 3
    brows = []
    errs = []
    for K5 in [1200, 333, 100, 30, 10]:
        mx = block_maxes_iid(src5, K5, m5, R5, rng)
        z = (mx.mean() - K5 * src5.h) / math.sqrt(K5 * src5.V)
        rel = abs(z - pred) / pred
        x3 = a ** 3 / math.sqrt(K5)
        errs.append(rel)
        brows.append([K5, K5 / lnm3, x3, z, pred, rel])
        print(f"  K={K5:5d}: K/(ln m)^3={K5/lnm3:5.2f}  x^3/sqrt K={x3:5.2f}  "
              f"z={z:.4f}  rel err={rel:.4f}")
    d5 = (all(errs[i] < errs[i + 1] for i in range(len(errs) - 1))
          and errs[0] < 0.04 and errs[-1] > 0.10 and errs[-1] > 4 * errs[0])
    report("D5", d5,
           f"refined-Gumbel error grows monotonically {errs[0]:.3f}->{errs[-1]:.3f} "
           f"as K falls through (ln m)^3: the ln m = o(K^{{1/3}}) boundary is "
           f"real (Cramer-Petrov x^3/sqrt K)", results)

    save_csv(os.path.join(outdir, "study_d_evt.csv"),
             ["source", "K", "m", "mean_max", "Kh", "z_std", "refined_pred",
              "rel_err", "leading_ratio"], rows)
    save_csv(os.path.join(outdir, "study_d_boundary.csv"),
             ["K", "K_over_lnm3", "x3_over_sqrtK", "z_std", "refined_pred",
              "rel_err"], brows)

    # ---------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    ax = axes[0]
    zs = [r[5] for r in rows if r[0] == "iid"]
    preds = [r[6] for r in rows if r[0] == "iid"]
    ax.plot(ms, zs, "o-", color=PALETTE[0], label="empirical std max")
    ax.plot(ms, preds, "--", color=PALETTE[3],
            label=r"$b_m + \gamma/a_m$ (refined Gumbel)")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("m (blocks)")
    ax.set_ylabel(r"$(E[\max_i L_i] - Kh)/\sqrt{KV}$")
    ax.set_title(f"read-buffer margin, iid, K={K}")
    ax.legend(fontsize=8)
    style_axes(ax)

    ax = axes[1]
    xs = np.sort(zstd)
    ax.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=PALETTE[0],
            label="empirical CDF")
    grid = np.linspace(xs[0], xs[-1], 300)
    ax.plot(grid, np.exp(-np.exp(-grid)), "--", color=PALETTE[3],
            label="Gumbel CDF")
    ax.set_xlabel(r"$a_m(Z - b_m)$")
    ax.set_ylabel("CDF")
    ax.set_title(f"Gumbel shape, m=1024 (KS={ks:.3f})")
    ax.legend(fontsize=8)
    style_axes(ax)

    ax = axes[2]
    ax.plot([r[1] for r in brows], [r[5] for r in brows], "o-", color=PALETTE[0],
            label="refined-Gumbel rel. error")
    ax.axvline(1.0, color=PALETTE[3], ls="--", label=r"$K = (\ln m)^3$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$K / (\ln m)^3$")
    ax.set_ylabel("relative error")
    ax.set_title("validity boundary $\\ln m = o(K^{1/3})$\n(skewed source, m=1024)")
    ax.legend(fontsize=8)
    style_axes(ax)
    fig.suptitle("Study D -- worst-case query EVT (Thm 7.35)", y=1.02)
    save_fig(fig, os.path.join(outdir, "study_d.png"))
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    res = run(full=args.full)
    ok = all(x[1] for x in res)
    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
