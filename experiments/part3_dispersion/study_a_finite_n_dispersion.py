#!/usr/bin/env python3
"""
study_a_finite_n_dispersion.py
============================================================================
STUDY (a): finite-n operational dispersion of the BSMS inside Gray.

Law under test (Remark 7.34b): on the Gray region 0 < D <= D_c(p) the
n-block d-tilted information (with the exact optimal output) satisfies the
all-n SLB identity j_n = i_n - n h2(D), hence

    Var(j_n)/n = V_lossless (n-1)/n  ->  V_lossless = p(1-p) log2^2((1-p)/p),

a D-independent plateau with an exactly-1/n finite-blocklength deficit.

DESIGN (three (p,D) Gray points: interior, near-edge, endpoint D = D_c):
  A0  the deconvolved output is RD-optimal: dense Blahut-Arimoto reference
      at n = 6 reproduces (R_n, D) of the exact deconvolution solver.
  A1  Gray validity: deconvolution P_Y >= 0 at n = 12 at all three points.
  A2  the SLB identity j_n = i_n - n h2(D) at n = 12 to machine precision
      (this is the exact operational j_n from the 2^n-dim solver, no
      identity assumed).
  A3  exact-grade dispersion: Var(j_n)/n from the full 2^n law equals
      V_lossless (n-1)/n to 1e-10 for n = 6..14.
  A4  simulation-grade dispersion: Monte Carlo source paths on an n-ladder
      (n = 64..2048 reduced / ..8192 full), empirical Var(j_n)/n within
      z*SE of V_lossless (n-1)/n, and the unbiased rescaled estimate
      Var(j_n)/(n-1) covers V_lossless at every ladder point.
  A5  convergence-rate fit: log-log fit of the deficit
      V_lossless - Var(j_n)/n against n on the exact-grade ladder gives
      slope -1 (|slope+1| < 0.02) and constant V_lossless (within 5%).

Outputs: out/study_a_dispersion.csv, out/study_a_dispersion.png.

Caveat (honest): at Monte Carlo blocklengths (n > 14) j_n is evaluated
through the all-n SLB identity -- proven in Remark 7.34b and machine-checked
here at the exact-grade blocklengths -- since the 2^n solver is infeasible;
the exact-grade rungs (A0-A3) are identity-free.
"""

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def dense_ba_reference(n, p, s, iters=40000, tol=1e-14):
    """Dense O(N^2) Blahut-Arimoto at slope s; returns (D_per_letter, R_bits_per_block)."""
    N = 1 << n
    P, _, _ = cm.bsms_block(n, p)
    x = np.arange(N, dtype=np.uint64)
    dmat = cm.popcount(x[:, None] ^ x[None, :]).astype(np.float64)
    W = np.exp(-s * dmat)
    q = np.full(N, 1.0 / N)
    for _ in range(iters):
        Z = np.maximum(W @ q, 1e-300)
        qn = q * (W.T @ (P / Z))
        qn /= qn.sum()
        if np.abs(qn - q).max() < tol:
            q = qn
            break
        q = qn
    Z = np.maximum(W @ q, 1e-300)
    Pyx = (q[None, :] * W) / Z[:, None]
    Dv = float(np.sum(P[:, None] * Pyx * dmat)) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = np.log2(np.where(Pyx > 0, Pyx / np.maximum(q[None, :], 1e-300), 1.0))
    R_bits = float(np.sum(np.where(Pyx > 0, P[:, None] * Pyx * lr, 0.0)))
    return Dv, R_bits


def exact_point(n, p, D):
    """Exact operational (R_n, Var(j_n)/n, identity residual, min P_Y) at (n,p,D)."""
    PY = cm.deconv_output(n, p, D)
    P, i_bits, _ = cm.bsms_block(n, p)
    j = cm.dtilted_bits(n, p, D, PY)
    Ej = float(np.dot(P, j))
    Vj = float(np.dot(P, (j - Ej) ** 2))
    resid = float(np.abs(j - (i_bits - n * cm.h2(D))).max())
    return dict(minPY=float(PY.min()), Ej=Ej, Vj_n=Vj / n, resid=resid)


def mc_var(n, p, D, nsamp, rng):
    """Monte Carlo Var(j_n)/n with moment-based SE, sampling BSMS paths.

    j_n is evaluated via the (machine-verified) SLB identity
    j_n = i_n - n h2(D), i_n = 1 + N_sw log2(1/p) + (n-1-N_sw) log2(1/(1-p)).
    """
    sw = rng.binomial(n - 1, p, size=nsamp).astype(np.float64)
    i_bits = 1.0 + sw * (-math.log2(p)) + (n - 1 - sw) * (-math.log2(1 - p))
    j = i_bits - n * cm.h2(D)
    s2 = float(np.var(j, ddof=1))
    m = float(np.mean(j))
    m4 = float(np.mean((j - m) ** 4))
    se = math.sqrt(max(m4 - s2 ** 2, 0.0) / nsamp)
    return s2, se


def run(full=False, seed=20260710):
    ck = cm.Checker()
    print("=" * 78)
    print("STUDY (a): finite-n operational dispersion, BSMS inside Gray")
    print("=" * 78)

    points = []
    for p, frac, tag in ((0.30, 0.5, "interior"), (0.38, 0.8, "near-edge"),
                         (0.45, 1.0, "endpoint")):
        D = frac * cm.D_c(p)
        points.append((p, D, tag))
        print(f"  point ({tag}): p={p}, D={D:.6f} (={frac} D_c, D_c={cm.D_c(p):.6f}), "
              f"V_lossless={cm.V_lossless_bits(p):.6f}")

    # ---- A0: dense BA reference vs deconvolution solver at n = 6
    print("-" * 78)
    p0, D0, _ = points[0]
    s0 = math.log((1 - D0) / D0)
    Dba, Rba = dense_ba_reference(6, p0, s0)
    ex = exact_point(6, p0, Dba)   # deconvolution j at the BA-achieved D
    okA0 = abs(Dba - D0) < 5e-3 and abs(ex["Ej"] - Rba) < 1e-6
    print(f"     dense BA (n=6, slope lambda*): D={Dba:.8f} R={Rba:.10f} bits;"
          f" deconv E[j_n]={ex['Ej']:.10f}")
    ck.rep("A0 deconvolved output is RD-optimal (dense BA reference, n=6)", okA0)

    # ---- A1 + A2: Gray validity and the SLB identity at n = 12
    okA1, okA2 = True, True
    for p, D, tag in points:
        e = exact_point(12, p, D)
        print(f"     n=12 {tag:>9}: min P_Y = {e['minPY']:.3e}, "
              f"|j_n - (i_n - n h2(D))|_max = {e['resid']:.3e}")
        okA1 &= e["minPY"] >= -1e-12
        okA2 &= e["resid"] < 1e-9
    ck.rep("A1 deconvolution valid (P_Y >= 0) at all three Gray points", okA1)
    ck.rep("A2 SLB identity j_n = i_n - n h2(D) machine-exact (n=12)", okA2)

    # ---- A3: exact-grade Var(j_n)/n == V (n-1)/n
    print("-" * 78)
    rows = []
    okA3 = True
    n_exact = list(range(6, 15, 2))
    for p, D, tag in points:
        V = cm.V_lossless_bits(p)
        for n in n_exact:
            e = exact_point(n, p, D)
            pred = V * (n - 1) / n
            err = abs(e["Vj_n"] - pred)
            okA3 &= err < 1e-10
            rows.append([tag, p, D, n, "exact", e["Vj_n"], pred, 0.0])
        print(f"     {tag:>9}: exact Var(j_n)/n at n={n_exact} matches "
              f"V (n-1)/n to <1e-10: {okA3}")
    ck.rep("A3 exact operational Var(j_n)/n = V_lossless (n-1)/n (n=6..14)", okA3)

    # ---- A4: Monte Carlo ladder with CIs
    print("-" * 78)
    n_ladder = [64, 128, 256, 512, 1024, 2048]
    nsamp = 40000
    if full:
        n_ladder += [4096, 8192]
        nsamp = 200000
    rng = np.random.default_rng(seed)
    okA4 = True
    z = 4.0
    for p, D, tag in points:
        V = cm.V_lossless_bits(p)
        worst = 0.0
        for n in n_ladder:
            s2, se = mc_var(n, p, D, nsamp, rng)
            pred_tot = V * (n - 1)
            dev = abs(s2 - pred_tot) / se
            worst = max(worst, dev)
            cover_inf = abs(s2 / (n - 1) - V) < z * se / (n - 1)
            okA4 &= dev < z and cover_inf
            rows.append([tag, p, D, n, "mc", s2 / n, pred_tot / n, se / n])
        print(f"     {tag:>9}: MC ladder n={n_ladder}, worst |dev|/SE = {worst:.2f} "
              f"(z gate {z})")
    ck.rep(f"A4 MC Var(j_n)/n covers V (n-1)/n and V_lossless (z={z}, "
           f"{nsamp} paths)", okA4)

    # ---- A5: convergence-rate fit on the exact-grade deficit
    print("-" * 78)
    okA5 = True
    for p, D, tag in points:
        V = cm.V_lossless_bits(p)
        xs, ys = [], []
        for n in n_exact:
            e = exact_point(n, p, D)
            deficit = V - e["Vj_n"]
            if deficit > 0:
                xs.append(math.log(n))
                ys.append(math.log(deficit))
        slope, logc = np.polyfit(xs, ys, 1)
        c = math.exp(logc)
        print(f"     {tag:>9}: deficit fit  V - Var/n ~ {c:.6f} * n^({slope:.4f})"
              f"   [predicted {V:.6f} * n^-1]")
        okA5 &= abs(slope + 1.0) < 0.02 and abs(c / V - 1.0) < 0.05
    ck.rep("A5 deficit rate fit: slope = -1 (+/-0.02), constant = V (+/-5%)", okA5)

    # ---- artifacts
    os.makedirs(OUT, exist_ok=True)
    cm.write_csv(os.path.join(OUT, "study_a_dispersion.csv"),
                 ["point", "p", "D", "n", "grade", "Var_j_over_n",
                  "predicted", "se"], rows)
    _plot(points, rows)
    print(f"  artifacts: out/study_a_dispersion.csv, out/study_a_dispersion.png")
    return ck


def _plot(points, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=False)
    for ax, (p, D, tag) in zip(axes, points):
        V = cm.V_lossless_bits(p)
        ex = [(r[3], r[5]) for r in rows if r[0] == tag and r[4] == "exact"]
        mc = [(r[3], r[5], r[7]) for r in rows if r[0] == tag and r[4] == "mc"]
        ax.axhline(V, color="k", lw=1, ls="--", label="V_lossless")
        ns = [n for n, _ in ex]
        ax.plot(ns, [v for _, v in ex], "o-", ms=4, label="exact 2^n solver")
        ax.plot([n for n, _, _ in mc], [v for _, v, _ in mc], "s", ms=4,
                label="Monte Carlo")
        ax.errorbar([n for n, _, _ in mc], [v for _, v, _ in mc],
                    yerr=[4 * s for _, _, s in mc], fmt="none", capsize=2)
        nn = np.geomspace(min(ns), max(n for n, _, _ in mc), 100)
        ax.plot(nn, V * (nn - 1) / nn, color="C3", lw=1,
                label="V (n-1)/n")
        ax.set_xscale("log")
        ax.set_title(f"{tag}: p={p}, D={D:.4f}")
        ax.set_xlabel("n")
        ax.set_ylabel("Var(j_n)/n  [bits^2/symbol]")
        ax.legend(fontsize=7)
    fig.suptitle("Study (a): BSMS Gray-region dispersion plateau, "
                 "Var(j_n)/n -> V_lossless")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "study_a_dispersion.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    ck = run(full=args.full)
    print("=" * 78)
    print(f"OVERALL (study a) -> {'PASS' if ck.ok else 'FAIL'}")
    sys.exit(0 if ck.ok else 1)
