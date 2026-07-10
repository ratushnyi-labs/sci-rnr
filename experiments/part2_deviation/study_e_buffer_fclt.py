#!/usr/bin/env python3
"""Study (e): functional CLT for the repair process; buffer laws (Theorem 7.36).

Setting: per-symbol repair cost ell_i = -log2 P(x_i | x_{i-1}) for a finite-state
Markov source; centered partial sum D_t = sum_{i<=t}(ell_i - h) is the streaming
write-buffer occupancy (emission ell_i, drain h per symbol).  Laws:
  * FCLT (Donsker): W_N(t) = D_{floor(Nt)}/sqrt(NV) => standard BM on [0,1];
  * net-surplus high-water-mark max_t D_t / sqrt(NV) => sup W =_d |N(0,1)|
    (half-normal; mean sqrt(2/pi));
  * PHYSICAL empty-start Lindley queue: high-water-mark = max draw-up
    max_j<=t (D_t - D_j) => sup |W| (Levy; theta law, mean sqrt(pi/2));
  * ONE varentropy V governs the archive dispersion (7.29), the read buffer
    (7.35) and the write buffer (7.36).

Empirical program (one set of simulated paths feeds every check):
  E1  Brownian scaling: Var W_N(t) = t, Cov = min(s,t); KS test of the
      rescaled non-overlapping increments against N(0,1); successive
      increments uncorrelated.
  E2  write buffer: max_t D_t / sqrt(NV) is half-normal (mean + KS distance).
  E3  Lindley high-water-mark: max draw-up / sqrt(NV) follows the sup|W|
      theta law (mean sqrt(pi/2)), strictly above the net-surplus law.
  E4  one V, three buffers: the SAME fundamental-matrix V standardizes
      (i) the end-value CLT (7.29 dispersion), (ii) the block-max EVT margin
      (7.35 read buffer), (iii) the half-normal write buffer (7.36).

Outputs: out/study_e_buffers.csv, out/study_e.png
"""
import argparse
import math
import os
import sys

import numpy as np
from scipy.stats import halfnorm, kstest, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (EULER, MarkovSource, OUT_DIR, PALETTE, gumbel_norm, report,
                    save_csv, save_fig, style_axes)

HALF_NORMAL_MEAN = math.sqrt(2 / math.pi)   # E sup W  = E|N(0,1)|
SUP_ABS_MEAN = math.sqrt(math.pi / 2)       # E sup|W| (theta law)


def sup_abs_cdf(a):
    """P(sup_{[0,1]} |W| <= a), alternating theta series."""
    a = np.asarray(a, dtype=float)
    out = np.zeros_like(a)
    pos = a > 1e-9
    s = np.zeros(pos.sum())
    for n in range(60):
        s += ((-1) ** n / (2 * n + 1)) * np.exp(
            -((2 * n + 1) ** 2) * math.pi ** 2 / (8 * a[pos] ** 2))
    out[pos] = (4 / math.pi) * s
    return np.clip(out, 0.0, 1.0)


def ks_distance(sample, cdf_fn):
    x = np.sort(np.asarray(sample))
    c = cdf_fn(x)
    n = len(x)
    hi = np.arange(1, n + 1) / n
    lo = np.arange(0, n) / n
    return float(max(np.max(np.abs(hi - c)), np.max(np.abs(lo - c))))


def run(full=False, outdir=OUT_DIR):
    results = []
    rng = np.random.default_rng(20260714)
    print("=" * 72)
    print("STUDY E -- buffer functional CLT (Thm 7.36): one V, three buffers")
    print("=" * 72)

    mk = MarkovSource([[0.7, 0.2, 0.1], [0.1, 0.6, 0.3], [0.2, 0.3, 0.5]])
    h, V = mk.h, mk.V
    N = 4000
    R = 8000 if full else 4000
    J = 8                      # non-overlapping increment blocks
    K = N // J                 # block length, reused as the 7.35 read-buffer K
    print(f"Markov 3-state: h={h:.4f}  V={V:.4f}  N={N}  R={R} paths")

    cost = mk.sample_cost_paths(N, R, rng)      # R x N per-step costs
    D = np.cumsum(cost - h, axis=1)             # centered partial sums
    s = math.sqrt(N * V)

    # ---------------------------------------------------------------- E1
    idx = {t: int(t * N) - 1 for t in (0.3, 0.5, 0.7, 1.0)}
    W = {t: D[:, i] / s for t, i in idx.items()}
    var_ok = all(abs(float(W[t].var()) - t) < 0.05 for t in (0.3, 0.5, 1.0))
    cov37 = float(np.cov(W[0.3], W[0.7])[0, 1])
    print(f"  Var W(0.3)={W[0.3].var():.4f} (0.3)  Var W(0.5)={W[0.5].var():.4f} "
          f"(0.5)  Var W(1)={W[1.0].var():.4f} (1.0)  Cov(0.3,0.7)={cov37:.4f} (0.3)")
    incr = cost.reshape(R, J, K).sum(axis=2) - K * h    # block increments
    z_incr = (incr / math.sqrt(K * V)).ravel()
    ks_incr = float(kstest(z_incr, "norm").statistic)
    lag1 = float(np.corrcoef(incr[:, :-1].ravel(), incr[:, 1:].ravel())[0, 1])
    print(f"  increments: KS vs N(0,1) = {ks_incr:.4f} (n={z_incr.size}); "
          f"lag-1 corr = {lag1:+.4f}")
    e1 = (var_ok and abs(cov37 - 0.3) < 0.04 and ks_incr < 0.02
          and abs(lag1) < 0.03)
    report("E1", e1,
           f"Brownian scaling: Var W(t)=t, Cov=min(s,t); rescaled increments "
           f"KS={ks_incr:.4f} vs N(0,1), lag-1 corr {lag1:+.3f}", results)

    # ---------------------------------------------------------------- E2
    M_net = D.max(axis=1) / s
    m_mean = float(M_net.mean())
    ks_hn = ks_distance(M_net, lambda x: 2 * norm.cdf(np.maximum(x, 0)) - 1)
    print(f"  net-surplus max: mean={m_mean:.4f} (sqrt(2/pi)={HALF_NORMAL_MEAN:.4f}); "
          f"KS vs half-normal={ks_hn:.4f}")
    e2 = abs(m_mean - HALF_NORMAL_MEAN) < 0.03 and ks_hn < 0.05
    report("E2", e2,
           f"write-buffer net-surplus high-water-mark is half-normal: "
           f"mean {m_mean:.3f} vs sqrt(2/pi)={HALF_NORMAL_MEAN:.4f}, "
           f"KS={ks_hn:.3f}", results)

    # ---------------------------------------------------------------- E3
    runmin = np.minimum.accumulate(np.minimum(D, 0.0), axis=1)
    drawup = ((D - runmin).max(axis=1)) / s     # Lindley max draw-up
    q_mean = float(drawup.mean())
    ks_th = ks_distance(drawup, sup_abs_cdf)
    print(f"  Lindley draw-up: mean={q_mean:.4f} (sqrt(pi/2)={SUP_ABS_MEAN:.4f}); "
          f"KS vs sup|W| theta law={ks_th:.4f}")
    e3 = (abs(q_mean - SUP_ABS_MEAN) < 0.05 and ks_th < 0.05
          and q_mean > m_mean + 0.15)
    report("E3", e3,
           f"physical queue high-water-mark follows the sup|W| theta law "
           f"(mean {q_mean:.3f} vs sqrt(pi/2)={SUP_ABS_MEAN:.4f}), strictly "
           f"above the half-normal net surplus", results)

    # ---------------------------------------------------------------- E4
    # one V standardizes three buffers (same closed-form fundamental-matrix V):
    # (i) archive dispersion 7.29: Var(D_N)/(N V) -> 1
    disp = float(W[1.0].var())
    # (ii) read buffer 7.35: block-max EVT margin with the same V
    m_evt = 64
    fam = (incr.ravel()[: (R * J // m_evt) * m_evt]).reshape(-1, m_evt)
    z_max = (fam.max(axis=1)) / math.sqrt(K * V)
    a_m, b_m = gumbel_norm(m_evt)
    evt_pred = b_m + EULER / a_m
    evt_emp = float(z_max.mean())
    # (iii) write buffer 7.36: half-normal mean ratio
    wb_ratio = m_mean / HALF_NORMAL_MEAN
    # empirical V from block sums vs closed form (independence of the check)
    V_emp = float(incr.var(ddof=1)) / K
    print(f"  (i) dispersion: Var(D_N)/(NV)={disp:.4f} (->1)")
    print(f"  (ii) read buffer: std block-max={evt_emp:.4f} vs b_m+gamma/a_m="
          f"{evt_pred:.4f} (m={m_evt}, K={K})")
    print(f"  (iii) write buffer: mean ratio={wb_ratio:.4f} (->1)")
    print(f"  empirical V from blocks = {V_emp:.4f} vs fundamental-matrix "
          f"V={V:.4f} (ratio {V_emp/V:.4f})")
    e4 = (abs(disp - 1) < 0.06 and abs(evt_emp / evt_pred - 1) < 0.06
          and abs(wb_ratio - 1) < 0.04 and abs(V_emp / V - 1) < 0.08)
    report("E4", e4,
           f"ONE V governs three buffers: dispersion ratio {disp:.3f}, EVT "
           f"ratio {evt_emp/evt_pred:.3f}, write-buffer ratio {wb_ratio:.3f} "
           f"(all ->1 under the same fundamental-matrix V)", results)

    save_csv(os.path.join(outdir, "study_e_buffers.csv"),
             ["quantity", "empirical", "theory", "ratio"],
             [["VarW(0.3)", float(W[0.3].var()), 0.3, float(W[0.3].var()) / 0.3],
              ["VarW(0.5)", float(W[0.5].var()), 0.5, float(W[0.5].var()) / 0.5],
              ["VarW(1.0)", disp, 1.0, disp],
              ["Cov(0.3,0.7)", cov37, 0.3, cov37 / 0.3],
              ["KS_incr_vs_normal", ks_incr, 0.0, ""],
              ["net_surplus_mean", m_mean, HALF_NORMAL_MEAN,
               m_mean / HALF_NORMAL_MEAN],
              ["KS_vs_halfnormal", ks_hn, 0.0, ""],
              ["lindley_drawup_mean", q_mean, SUP_ABS_MEAN,
               q_mean / SUP_ABS_MEAN],
              ["KS_vs_supabs_theta", ks_th, 0.0, ""],
              ["evt_block_max_std", evt_emp, evt_pred, evt_emp / evt_pred],
              ["V_blocks_vs_formula", V_emp, V, V_emp / V]])

    # ---------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    ax = axes[0]
    ts = np.linspace(0.05, 1.0, 20)
    emp_var = [float((D[:, int(t * N) - 1] / s).var()) for t in ts]
    ax.plot(ts, emp_var, "o", ms=4, color=PALETTE[0], label="Var $W_N(t)$")
    ax.plot(ts, ts, "--", color=PALETTE[3], label="t (Brownian)")
    ax.set_xlabel("t")
    ax.set_ylabel("variance")
    ax.set_title("Brownian variance profile")
    ax.legend(fontsize=8)
    style_axes(ax)

    ax = axes[1]
    grid = np.linspace(0, 3.2, 300)
    xs = np.sort(M_net)
    ax.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=PALETTE[0],
            label="net-surplus max (emp. CDF)")
    ax.plot(grid, 2 * norm.cdf(grid) - 1, "--", color=PALETTE[0], lw=1,
            label="half-normal")
    xq = np.sort(drawup)
    ax.plot(xq, np.arange(1, len(xq) + 1) / len(xq), color=PALETTE[3],
            label="Lindley draw-up (emp. CDF)")
    ax.plot(grid, sup_abs_cdf(grid), "--", color=PALETTE[3], lw=1,
            label=r"sup$|W|$ theta law")
    ax.set_xlabel(r"high-water-mark $/\sqrt{NV}$")
    ax.set_ylabel("CDF")
    ax.set_title("two buffer conventions, two laws")
    ax.legend(fontsize=7)
    style_axes(ax)

    ax = axes[2]
    names = ["dispersion\n(7.29)", "read buffer\nEVT (7.35)",
             "write buffer\n(7.36)"]
    ratios = [disp, evt_emp / evt_pred, wb_ratio]
    ax.bar(names, ratios, color=[PALETTE[0], PALETTE[1], PALETTE[2]],
           width=0.55)
    ax.axhline(1.0, color=PALETTE[3], ls="--", lw=1)
    for i, v in enumerate(ratios):
        ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=8)
    ax.set_ylim(0, 1.25)
    ax.set_ylabel("empirical / theory (same V)")
    ax.set_title("one varentropy V, three buffers")
    style_axes(ax)
    fig.suptitle("Study E -- functional CLT and buffer laws (Thm 7.36)", y=1.02)
    save_fig(fig, os.path.join(outdir, "study_e.png"))
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    res = run(full=args.full)
    ok = all(x[1] for x in res)
    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
