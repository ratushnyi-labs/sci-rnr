#!/usr/bin/env python3
"""
study_d_replica_margins.py
============================================================================
STUDY (d): replica-margin sampling -- R_A(s;D) - 3/2 over dense Gray grids.

Law under test (Theorem 7.34m / the Route B assembly chain; Remark 7.34m'):
the pattern-quotient replica spectral ratio

    R_A(s;D) = (1 - g_A(s)) / (D(1-D)(1 - cos s)),   g_A = |C_s/C_0|^2 rho(Q),

with Q the pattern-quotient replica matrix on agreement patterns of the
triple (source, replica, conjugate replica) -- 4x4 at A=2, 5x5 at A>=3 --
satisfies R_A(s;D) > 3/2 on the whole Gray region

    p in (0, (A-1)/A),  0 < D <= D_c^(A)(p),  s in (0, pi],

with the small-p binding-corner law (asymptotically sharp)

    margin(p) := R_A(pi; D_c(p)) - 3/2 = (A-2)/(A-1) p + O(p^2).

METHOD: the committed pattern-quotient builder (numeric port, validated
against the committed A=3/4/5 assembly scripts to 1e-10); the A-ary Gray
threshold D_c^(A)(p) as the SMALLEST positive root of the relevant-cubic
discriminant of the multiplicity-weighted 3-class alternating transfer
(the collision flag is non-monotone in D for A >= 4 -- first-sign-change
scan, then bisection).

CHECKS:
  D0  cross-validation against the committed assembly scripts
      (scripts/verify/thm_7_34_route_b_a{3,4,5}_assembly.py): D_c^(A)(p)
      and R_A(s;D) agree at reference points (1e-9 / 1e-10).
  D1  dense Gray grid: margin = R_A - 3/2 > 0 at EVERY grid point,
      A = 2..5 (p-grid x D/D_c-grid x s-grid; min margin reported per A).
  D2  binding corner (D = D_c, s = pi): margin > 0 across the full p-grid.
  D3  small-p law: Richardson (linear-in-p) extrapolation of margin/p on a
      corner p-ladder equals (A-2)/(A-1) within 0.01 absolute (A=3,4,5)
      and 0.01 at A=2 (where the predicted linear coefficient is 0).

Outputs: out/study_d_margins.csv (full margin grid),
         out/study_d_margins.png, out/study_d_corner_law.png.
"""

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as cm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
VERIFY_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "scripts", "verify")


def _load_committed(A):
    import importlib.util
    path = os.path.join(VERIFY_DIR, f"thm_7_34_route_b_a{A}_assembly.py")
    spec = importlib.util.spec_from_file_location(f"assembly_a{A}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(full=False):
    ck = cm.Checker()
    print("=" * 78)
    print("STUDY (d): replica margins R_A(s;D) - 3/2 on dense Gray grids, "
          "A = 2..5")
    print("=" * 78)

    # ---- D0: cross-validation against the committed assembly scripts
    okD0 = True
    for A in (3, 4, 5):
        try:
            mod = _load_committed(A)
        except FileNotFoundError:
            print(f"     committed A={A} assembly script missing")
            okD0 = False
            continue
        R_ref = getattr(mod, f"R{A}")
        for pv in (0.1, 0.4):
            dc_m, dc_r = cm.Dc_aary(A, pv), float(mod.Dc_of(pv))
            r_m = cm.replica_ratio(A, pv, 0.8 * dc_r, 2.0)
            r_r = float(R_ref(pv, 0.8 * dc_r, 2.0))
            okD0 &= abs(dc_m - dc_r) < 1e-9 and abs(r_m - r_r) < 1e-10
        print(f"     A={A}: D_c and R_A match the committed assembly "
              f"at reference points: {okD0}")
    ck.rep("D0 threshold + ratio cross-validated vs committed assemblies",
           okD0)

    # ---- D1/D2: dense Gray grids
    n_p = 32 if full else 16
    n_s = 20 if full else 10
    fracD = ([0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 1.0] if full
             else [0.15, 0.3, 0.5, 0.7, 0.85, 1.0])
    rows = []
    okD1 = okD2 = True
    grid_stats = {}
    for A in (2, 3, 4, 5):
        pmax = (A - 1) / A
        pgrid = np.concatenate([
            np.geomspace(0.02, 0.3, n_p // 2, endpoint=False),
            np.linspace(0.3, 0.98, n_p - n_p // 2),
        ]) * pmax
        sgrid = np.linspace(0.15, math.pi, n_s)
        minm, argmin = np.inf, None
        corner_min = np.inf
        for p in pgrid:
            Dc = cm.Dc_aary(A, float(p))
            for fD in fracD:
                D = fD * Dc
                for s in sgrid:
                    m = cm.replica_ratio(A, float(p), float(D), float(s)) - 1.5
                    rows.append([A, float(p), fD, float(D), float(s), m])
                    if m < minm:
                        minm, argmin = m, (float(p), fD, float(s))
                    if fD == 1.0 and abs(s - math.pi) < 1e-12:
                        corner_min = min(corner_min, m)
        okD1 &= minm > 0
        okD2 &= corner_min > 0
        grid_stats[A] = (minm, argmin, corner_min)
        print(f"     A={A}: {len(pgrid) * len(fracD) * len(sgrid)} grid pts, "
              f"min margin = {minm:.6f} at (p={argmin[0]:.4f}, "
              f"D={argmin[1]}*D_c, s={argmin[2]:.3f}); "
              f"corner-line min = {corner_min:.6f}")
    ck.rep("D1 min margin > 0 over the whole Gray grid, A = 2..5", okD1)
    ck.rep("D2 binding corner (D = D_c, s = pi) margin > 0, all p", okD2)

    # ---- D3: small-p corner law
    ladder = [0.04, 0.02, 0.01] + ([0.005, 0.0025] if full else [])
    okD3 = True
    corner_rows = []
    corner_data = {}
    for A in (2, 3, 4, 5):
        pred = (A - 2) / (A - 1)
        vals = []
        for p in ladder:
            dc = cm.Dc_aary(A, p)
            m = cm.replica_ratio(A, p, dc, math.pi) - 1.5
            vals.append((p, m / p))
            corner_rows.append([A, p, dc, m, m / p, pred])
        (p1, v1), (p2, v2) = vals[-2], vals[-1]
        extrap = v2 + (v2 - v1) * p2 / (p1 - p2)
        corner_data[A] = (vals, extrap, pred)
        okD3 &= abs(extrap - pred) < 0.01
        print(f"     A={A}: margin/p ladder {['%.4f' % v for _, v in vals]}"
              f" -> extrapolated {extrap:+.5f}, predicted (A-2)/(A-1) = "
              f"{pred:.5f}")
    ck.rep("D3 corner law margin ~ (A-2)/(A-1) p (extrapolation +/-0.01)",
           okD3)

    os.makedirs(OUT, exist_ok=True)
    cm.write_csv(os.path.join(OUT, "study_d_margins.csv"),
                 ["A", "p", "D_over_Dc", "D", "s", "margin"], rows)
    cm.write_csv(os.path.join(OUT, "study_d_corner_law.csv"),
                 ["A", "p", "Dc", "margin", "margin_over_p",
                  "predicted_slope"], corner_rows)
    _plot(rows, corner_data)
    print("  artifacts: out/study_d_margins.csv, out/study_d_corner_law.csv,")
    print("             out/study_d_margins.png, out/study_d_corner_law.png")
    return ck


def _plot(rows, corner_data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # margin heatmaps: min-over-s margin on the (p, D/D_c) plane, per A
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
    for ax, A in zip(axes, (2, 3, 4, 5)):
        pts = [(r[1], r[2], r[5]) for r in rows if r[0] == A]
        ps = sorted({x[0] for x in pts})
        fs = sorted({x[1] for x in pts})
        Z = np.full((len(fs), len(ps)), np.inf)
        for p, f, m in pts:
            i, j = fs.index(f), ps.index(p)
            Z[i, j] = min(Z[i, j], m)
        im = ax.pcolormesh(ps, fs, np.log10(Z), shading="nearest",
                           cmap="viridis")
        fig.colorbar(im, ax=ax, label="log10 min-over-s margin")
        ax.set_title(f"A={A}")
        ax.set_xlabel("p")
        ax.set_ylabel("D / D_c")
        ax.set_xscale("log")
    fig.suptitle("Study (d): replica margin R_A(s;D) - 3/2 over the Gray "
                 "region (min over s; all > 0)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "study_d_margins.png"), dpi=140)
    plt.close(fig)

    # corner law
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for k, (A, (vals, extrap, pred)) in enumerate(sorted(corner_data.items())):
        ps = [p for p, _ in vals]
        vs = [v for _, v in vals]
        ax.plot(ps, vs, "o-", color=f"C{k}", label=f"A={A} (measured)")
        ax.axhline(pred, color=f"C{k}", ls="--", lw=1)
        ax.plot([0], [extrap], "*", color=f"C{k}", ms=10)
    ax.set_xlabel("p")
    ax.set_ylabel("margin / p at (D_c, s=pi)")
    ax.set_title("binding-corner law: margin/p -> (A-2)/(A-1) "
                 "(dashed; stars = extrapolation)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "study_d_corner_law.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    ck = run(full=args.full)
    print("=" * 78)
    print(f"OVERALL (study d) -> {'PASS' if ck.ok else 'FAIL'}")
    sys.exit(0 if ck.ok else 1)
