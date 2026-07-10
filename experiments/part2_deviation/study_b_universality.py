#!/usr/bin/env python3
"""Study (b): universality dividend (Theorem 7.32).

Laws under test, for a full-support i.i.d. source P0 on A symbols (d = A-1):
  * OPTIMAL overflow exponent  E*(R) = inf{ D(Q||P0) : H(Q) >= R }, attained by
    the universal (KT / Dirichlet(1/2)) mixture code;
  * a FIXED model M -- even the TRUE one -- is overflow-tail-suboptimal:
    its exponent is the Legendre rate E_M(R) of Thm 7.31 (matched M=P0 gives
    E_arith(R) = inf{D : H+D >= R} < E*(R) strictly for R > H(P0); a WRONG
    fixed model is worse still);
  * mixture cost: L_KT - L_ideal = (d/2) log2 N + O(1) (Clarke-Barron) --
    third-order, dispersion- and LDP-invisible, yet exponent-IMPROVING.

Empirical program (exact multinomial type sums as ground truth + Monte Carlo):
  B1  exponent deficit: exact overflow curves -log2 P(L >= NR) for the three
      codes (wrong fixed / true fixed / KT mixture) across an N-ladder; the
      BR-corrected slopes recover E_wrong < E_arith < E*; direct MC agrees
      with the exact tails at small N.
  B2  the ordering is the theory ordering: slope_wrong < slope_fix <
      slope_KT and each within tolerance of its theoretical exponent.
  B3  mixture cost fit: E[L_KT - L_fix] regressed on log2 N has slope
      d/2 (the Clarke-Barron coefficient), and the residual constant is O(1).

Outputs: out/study_b_exponents.csv, out/study_b_redundancy.csv, out/study_b.png
"""
import argparse
import math
import os
import sys

import numpy as np
from scipy.special import gammaln

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (IIDSource, LN2, OUT_DIR, PALETTE, report, save_csv,
                    save_fig, style_axes)

P0 = np.array([0.6, 0.3, 0.1])
M_WRONG = np.array([0.45, 0.45, 0.10])
A = 3
D_DIM = A - 1  # simplex dimension d


def H_of(q):
    q = np.asarray(q, float)
    qs = q[q > 0]
    return float(-(qs * np.log2(qs)).sum())


def D_of(q, p):
    q = np.asarray(q, float)
    mask = q > 0
    return float((q[mask] * np.log2(q[mask] / p[mask])).sum())


def tilt(s):
    w = P0 ** s
    return w / w.sum()


def Estar(R):
    """E*(R) = D(Q_s||P0) at the entropy-flattened tilt with H(Q_s)=R."""
    lo, hi = 1e-9, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if H_of(tilt(mid)) > R:
            lo = mid
        else:
            hi = mid
    return D_of(tilt(0.5 * (lo + hi)), P0)


def all_types(N):
    """All count triples (n1,n2,n3) with sum N, plus exact log2 multinomial mass."""
    n1 = np.repeat(np.arange(N + 1), np.arange(N + 1, 0, -1))
    n2 = np.concatenate([np.arange(N - k + 1) for k in range(N + 1)])
    n3 = N - n1 - n2
    counts = np.stack([n1, n2, n3], axis=1)
    logmass = (gammaln(N + 1) - gammaln(counts + 1).sum(axis=1)
               + (counts * np.log(P0)[None, :]).sum(axis=1)) / LN2
    return counts, logmass


def L_kt(counts):
    """KT (Dirichlet(1/2)) mixture codelength in bits, vectorized over rows."""
    N = counts.sum(axis=1)
    lp = (gammaln(A / 2.0) - gammaln(N + A / 2.0)
          + (gammaln(counts + 0.5) - gammaln(0.5)).sum(axis=1))
    return -lp / LN2


def exact_tail_log2(counts, logmass, L, thresh):
    sel = L >= thresh - 1e-9
    if not sel.any():
        return -np.inf
    m = logmass[sel].max()
    return float(m + math.log2(np.exp((logmass[sel] - m) * LN2).sum()))


def run(full=False, outdir=OUT_DIR):
    results = []
    rng = np.random.default_rng(20260711)
    print("=" * 72)
    print("STUDY B -- universality dividend (Thm 7.32): exponent deficit + "
          "(d/2)log N cost")
    print("=" * 72)

    src_true = IIDSource(P0)
    src_wrong = IIDSource(P0, model=M_WRONG)
    l_true = -np.log2(P0)
    l_wrong = -np.log2(M_WRONG)
    h = src_true.h
    R_thr = h + 0.16                      # overflow threshold rate, in (H, log2 A)
    E_star = Estar(R_thr)
    E_arith = src_true.legendre_E(R_thr)  # matched fixed model (Thm 7.31)
    E_wrong = src_wrong.legendre_E(R_thr)  # mismatched fixed model
    print(f"P0={P0.tolist()}  h={h:.4f}  R=h+0.16={R_thr:.4f}")
    print(f"theory: E_wrong={E_wrong:.5f} < E_arith={E_arith:.5f} < E*={E_star:.5f}")

    # ------------------------------------------------------------ exact curves
    ladder = [100, 200, 300, 450, 600] if not full else [100, 200, 300, 450, 600, 900]
    rows = []
    curves = {"wrong": [], "fixed": [], "KT": []}
    for N in ladder:
        counts, logmass = all_types(N)
        Lf = counts @ l_true
        Lw = counts @ l_wrong
        Lk = L_kt(counts)
        for name, L in [("wrong", Lw), ("fixed", Lf), ("KT", Lk)]:
            lp = exact_tail_log2(counts, logmass, L, N * R_thr)
            curves[name].append((N, lp))
            rows.append([name, N, lp, -lp / N])
    # BR-corrected slope fits
    slopes = {}
    for name, pts in curves.items():
        Ns = np.array([p[0] for p in pts], float)
        ys = np.array([-p[1] for p in pts], float) - 0.5 * np.log2(Ns)
        slopes[name] = float(np.linalg.lstsq(
            np.vstack([Ns, np.ones_like(Ns)]).T, ys, rcond=None)[0][0])
    print(f"exact-tail BR-corrected slopes: wrong={slopes['wrong']:.5f} "
          f"fixed={slopes['fixed']:.5f} KT={slopes['KT']:.5f}")

    # MC agreement at N=150 (empirical check that the exact sums are the truth)
    N_mc = 150
    counts_mc = rng.multinomial(N_mc, P0, size=1_000_000 if full else 500_000)
    R_mc = len(counts_mc)
    counts_ex, logmass_ex = all_types(N_mc)
    mc_ok = True
    for name, lvec in [("wrong", l_wrong), ("fixed", l_true)]:
        L = counts_mc @ lvec
        hits = int((L >= N_mc * R_thr - 1e-9).sum())
        lp_mc = math.log2(hits / R_mc) if hits else -np.inf
        lp_ex = exact_tail_log2(counts_ex, logmass_ex,
                                counts_ex @ lvec, N_mc * R_thr)
        print(f"  MC check {name} N={N_mc}: log2P MC={lp_mc:.3f} exact={lp_ex:.3f} "
              f"(hits={hits})")
        mc_ok &= hits > 50 and abs(lp_mc - lp_ex) < 0.25
    Lk_mc = L_kt(counts_mc)
    hits = int((Lk_mc >= N_mc * R_thr - 1e-9).sum())
    lp_mc = math.log2(hits / R_mc) if hits else -np.inf
    lp_ex = exact_tail_log2(counts_ex, logmass_ex, L_kt(counts_ex), N_mc * R_thr)
    print(f"  MC check KT    N={N_mc}: log2P MC={lp_mc:.3f} exact={lp_ex:.3f} "
          f"(hits={hits})")
    mc_ok &= hits > 50 and abs(lp_mc - lp_ex) < 0.25
    report("B1", mc_ok,
           "direct MC overflow frequencies agree with the exact multinomial "
           "tails for all three codes (<0.25 bits)", results)

    # ------------------------------------------------------------ B2 ordering
    tol = {"wrong": 0.10, "fixed": 0.10, "KT": 0.20}  # KT has the (d/2)logN/N shift
    theory = {"wrong": E_wrong, "fixed": E_arith, "KT": E_star}
    each_ok = all(abs(slopes[k] - theory[k]) / theory[k] < tol[k] for k in slopes)
    order_ok = slopes["wrong"] < slopes["fixed"] < slopes["KT"]
    gap_ok = slopes["KT"] > slopes["fixed"] * 1.05  # the dividend is strict
    for k in ("wrong", "fixed", "KT"):
        print(f"  {k:5s}: slope={slopes[k]:.5f}  theory={theory[k]:.5f}  "
              f"rel err={abs(slopes[k]-theory[k])/theory[k]:.3f}")
    report("B2", each_ok and order_ok and gap_ok,
           f"empirical exponents ordered wrong<fixed<KT and track "
           f"E_wrong={E_wrong:.4f} < E_arith={E_arith:.4f} < E*={E_star:.4f} "
           f"(strict dividend)", results)

    # ------------------------------------------------------------ B3 (d/2)log N
    lad_red = [64, 128, 256, 512, 1024, 2048, 4096]
    if full:
        lad_red.append(8192)
    R_red = 40_000 if full else 20_000
    red_rows = []
    red_means = []
    for N in lad_red:
        counts = rng.multinomial(N, P0, size=R_red)
        red = L_kt(counts) - counts @ l_true
        red_means.append((N, float(red.mean()), float(red.std(ddof=1))))
        red_rows.append([N, float(red.mean()), float(red.std(ddof=1))])
    Ns = np.array([r[0] for r in red_means], float)
    ys = np.array([r[1] for r in red_means], float)
    X = np.vstack([np.log2(Ns), np.ones_like(Ns)]).T
    coef, const = np.linalg.lstsq(X, ys, rcond=None)[0]
    pred_coef = D_DIM / 2.0
    print(f"  redundancy fit: E[L_KT - L_fix] = {coef:.4f} * log2 N + {const:.3f}  "
          f"(Clarke-Barron d/2 = {pred_coef})")
    resid = ys - (pred_coef * np.log2(Ns) + const)
    b3_ok = abs(coef - pred_coef) / pred_coef < 0.06 and np.max(np.abs(resid)) < 0.25
    report("B3", b3_ok,
           f"mixture cost fits (d/2) log2 N: slope {coef:.3f} vs d/2={pred_coef} "
           f"(within 6%), O(1) residual", results)

    save_csv(os.path.join(outdir, "study_b_exponents.csv"),
             ["code", "N", "log2P_exact", "rate"], rows)
    save_csv(os.path.join(outdir, "study_b_redundancy.csv"),
             ["N", "mean_redundancy_bits", "std"], red_rows)

    # ------------------------------------------------------------ figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ax = axes[0]
    labels = {"wrong": "fixed wrong model", "fixed": "fixed true model",
              "KT": "universal mixture (KT)"}
    for i, name in enumerate(("wrong", "fixed", "KT")):
        pts = curves[name]
        ax.plot([p[0] for p in pts], [-p[1] for p in pts], "o-",
                color=PALETTE[i], label=f"{labels[name]}")
        ax.plot([p[0] for p in pts],
                [theory[name] * p[0] + 0.5 * math.log2(p[0]) for p in pts],
                ls="--", lw=1, color=PALETTE[i])
    ax.set_xlabel("N")
    ax.set_ylabel(r"$-\log_2 P(L \geq NR)$ (exact)")
    ax.set_title("overflow exponents: dashed = theory slope\n"
                 r"$E_{wrong} < E_{arith} < E^*$")
    ax.legend(fontsize=8)
    style_axes(ax)

    ax = axes[1]
    ax.plot(np.log2(Ns), ys, "o", color=PALETTE[0], label="mean redundancy")
    ax.plot(np.log2(Ns), pred_coef * np.log2(Ns) + const, "--",
            color=PALETTE[3], label=f"(d/2) log$_2$ N + c  (d/2={pred_coef})")
    ax.set_xlabel(r"$\log_2 N$")
    ax.set_ylabel(r"$E[L_{KT} - L_{fix}]$ (bits)")
    ax.set_title("universality price = (d/2) log$_2$ N (third order)")
    ax.legend(fontsize=8)
    style_axes(ax)
    fig.suptitle("Study B -- universality dividend (Thm 7.32)", y=1.02)
    save_fig(fig, os.path.join(outdir, "study_b.png"))
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    res = run(full=args.full)
    ok = all(r[1] for r in res)
    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
