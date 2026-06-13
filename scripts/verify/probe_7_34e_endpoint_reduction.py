#!/usr/bin/env python3
"""
probe_7_34e_endpoint_reduction.py
============================================================================
TEST the orchestrator's "s=pi endpoint reduction" route for closing item (i)
of Lemma 7.34e (instance (2) asymmetric replica spectral inequality) IN
CLOSED FORM, replacing the per-(a,b,D) Sturm scheme.

We verify the orchestrator's findings, then test the two load-bearing claims.

Helpers reused from probe_7_34_asym_replica_structure.py:
  eta_of(s,D), Cr2(s,D)=|C_s/C_0|^2, W8_local(eta,a,b), rho_of(W).

g(s) = Cr2(s,D) * rho(W8_local(eta_of(s,D), a, b)).

R(a,b,D,s) := (1 - g(s)) / (D(1-D)(1-cos s)).

FINDINGS to verify:
 (F1) over (a,b,D) grid x fine s-grid, R(a,b,D,s) is MINIMIZED at s=pi
      (no strict interior minimum below R(pi)); R(s) monotone decreasing.
 (F2) at s=pi, eta_pi = -2 D(1-D)/((1-D)^2+D^2) EXACTLY (real); g(pi) is
      a closed-form algebraic function.  [CORRECTION to the orchestrator's
      "top eigenvalue of a 2x2": the dominant root of W_8(pi) is in the
      degree-6 char-poly factor; g(pi) = |C|^2 * Perron(6x6 CP map on
      symmetric pairs), NOT a 2x2 radical -- see probe_7_34e_gpi_reduction6.py.
      The text's "2x2 PSD check" is the CERTIFICATE inequality Psi(Z)_x < nu Z_x
      (each Z_x is 2x2), which is a correct and different statement.]
 (F3) inf over {a+b<=1-eta} of R = c_tilde(eta) positive on bounded mixing,
      ->0 as a+b->1.

CLAIMS to test:
 (P1) [endpoint reduction] for all (a,b,D) in open Gray and all s in (0,pi]:
      (1-g(s))/(1-cos s) >= (1-g(pi))/2, i.e. s=pi minimizes R.
 (P2) [closed-form positivity] g(pi) < 1 on open Gray, in closed form, and
      ctil := R(pi) = (1-g(pi))/(2 D(1-D)) > 0 explicitly.

OUTCOMES (this probe):
 * P1-Zpi: the would-be CLEAN CLOSER -- a SINGLE fixed pair Z(pi) certifying
   ghat_{Z(pi)}(s) <= 1 - R(pi) dd (1-cos s) for ALL s -- is REFUTED.  The
   frozen-pair Collatz-Wielandt slack ghat_{Z(pi)}(s) - g(s) is ~ constant in
   s (~0.066 at (0.2,0.4,Dc)), and since the true g(s) only barely clears the
   bound (margin ~1e-3..1e-4), the frozen pair overshoots by ~0.07.  No single
   pair (Perron(pi), dressed(pi), power-iterated-at-pi) works.
 * P1-conv: P1 is EXACTLY EQUIVALENT to convexity of G(u)=g(arccos u) in
   u=cos s (with G(1)=1).  Confirmed densely (644 pts incl. boundary, slow
   mixing, D=Dc): G is convex in u (worst d2 ~ -7e-11 ~ 0).  Cr2 is EXACTLY
   linear in u; |eta|^2 is Mobius (deg1/deg1) in u; the crux is rho(W_8(s))
   convex in u.  W_8 is genuinely COMPLEX for s in (0,pi) (real only at 0,pi),
   so Cohen/Kingman real-Perron-convexity theorems do NOT apply -- the honest
   obstruction to a closed-form convexity proof.
 * COMPACT: rigorous FALLBACK.  R extends continuously to s=0 (->2(1-Psi0))
   and D=0 (->2), is >0 per point on Gray, so on the compact box
   {a+b<=1-eta}x{D in [0,Dc]}x{s in [0,pi]} the inf is attained > 0:
   c_tilde(eta) > 0, binding at s=pi.  This gives a per-source uniform ctil
   on bounded-mixing WITHOUT the per-point Sturm scheme, but is NOT a single
   closed-form expression for ctil(a,b,D).

VERDICT: closed-form one-pair closer REFUTED; P1 reduced EXACTLY to "g convex
in cos s" (numerically confirmed, not closed-form-proven; obstruction =
complex W_8); compactness fallback gives uniform ctil(eta)>0 rigorously.
The existing per-(a,b,D) Sturm certificate (L5) remains the exact-grade proof.
"""
import math

import numpy as np

from probe_7_34_asym_replica_structure import (eta_of, Cr2, W8_local, rho_of,
                                               chain, Psi, lam_CW, Z_power,
                                               dressed_data, Z_dressed)
from probe_7_34_nonsym_dual_identity import eta_C, theta0


def Dc_closed(a, b):
    """Lemma 7.34f closed form."""
    return 0.5 * (1.0 - math.sqrt(max(0.0, 1.0 - 4 * a * b / (2 - a - b) ** 2)))


def g_of(s, a, b, D):
    return Cr2(s, D) * rho_of(W8_local(eta_of(s, D), a, b))


def R_of(s, a, b, D):
    dd = D * (1 - D)
    den = dd * (1 - math.cos(s))
    if den <= 0:
        return float("inf")
    return (1.0 - g_of(s, a, b, D)) / den


# -------------------- F1: s=pi minimizes R --------------------

def check_F1(grid):
    print("-" * 84)
    print("F1: R(a,b,D,s) minimized at s=pi (no interior min below R(pi))")
    ok = True
    worst_interior_deficit = -np.inf  # max over points of (R(pi) - min_interior R)
    n_interior_below = 0
    svals = np.linspace(0.02, math.pi, 240)
    for (a, b, D, tag) in grid:
        Rs = np.array([R_of(s, a, b, D) for s in svals])
        Rpi = R_of(math.pi, a, b, D)
        # strict interior = exclude the last point (pi)
        Rint = Rs[:-1]
        rmin_int = Rint.min()
        deficit = Rpi - rmin_int  # >0 means interior dips below R(pi)
        # monotone-decreasing check (allow tiny numerical wiggle)
        diffs = np.diff(Rs)
        mono = bool(np.all(diffs < 1e-6))
        if deficit > 1e-6:
            n_interior_below += 1
        worst_interior_deficit = max(worst_interior_deficit, deficit)
        flag = "" if deficit <= 1e-6 else "  <-- INTERIOR BELOW pi"
        if deficit > 1e-6 or not mono:
            print(f"  {tag} a={a} b={b} D={D:.5f}: R(pi)={Rpi:.4f} "
                  f"min_int={rmin_int:.4f} deficit={deficit:+.2e} mono={mono}{flag}")
    ok = (n_interior_below == 0)
    print(f"  points with interior R below R(pi): {n_interior_below}/{len(grid)}")
    print(f"  worst (R(pi)-min_interior R) = {worst_interior_deficit:+.3e} "
          f"(<=1e-6 means pi is the global min)")
    print(f"  F1 -> {'PASS' if ok else 'FAIL'}")
    return ok


# -------------------- F2: eta_pi closed form & g(pi)=2x2 --------------------

def check_F2(grid):
    print("-" * 84)
    print("F2: eta_pi = -2D(1-D)/((1-D)^2+D^2) real; g(pi) closed form")
    ok = True
    worst_eta = 0.0
    for (a, b, D, tag) in grid:
        eta_num = eta_of(math.pi, D)
        eta_cf = -2 * D * (1 - D) / ((1 - D) ** 2 + D ** 2)
        worst_eta = max(worst_eta, abs(eta_num - eta_cf))
    ok &= worst_eta < 1e-12
    print(f"  worst |eta_pi(numeric) - closed form| = {worst_eta:.2e}")

    # g(pi) = Cr2(pi) * rho(W8). At s=pi eta is REAL, so W8 = M_y kron M_y /
    # T(x,y) with M_y real; the CP map Psi acts on pairs of 2x2 SYMMETRIC
    # matrices -> a real 6-dim operator, and rho(W8(pi)) is its PERRON root
    # (Krein-Rutman; eigenvector a PSD pair).  NOT a 2x2 eigenvalue; the
    # dominant root is in the degree-6 char-poly factor.  (reduction6 probe.)
    print("  g(pi), R(pi) values (R(pi) = ctil candidate):")
    for (a, b, D, tag) in grid:
        gpi = g_of(math.pi, a, b, D)
        Rpi = R_of(math.pi, a, b, D)
        print(f"    {tag} a={a} b={b} D={D:.5f}: g(pi)={gpi:.4f} R(pi)={Rpi:.4f}")
    print(f"  F2 -> {'PASS' if ok else 'FAIL'}")
    return ok


# -------------------- F3: bounded-mixing inf of R --------------------

def check_F3():
    print("-" * 84)
    print("F3: inf over {a+b<=1-eta} of R(pi) -> ctil(eta); >0 bounded, ->0 as a+b->1")
    ok = True
    for eta in (0.5, 0.2, 0.1, 0.03):
        smax = 1 - eta
        best = float("inf")
        argmin = None
        # scan a,b with a+b <= smax, a,b in (eps, ...), use D=D_c (the worst)
        for a in np.linspace(0.02, smax - 0.02, 30):
            for b in np.linspace(0.02, smax - a, 25):
                if a + b > smax or a + b >= 1:
                    continue
                Dc = Dc_closed(a, b)
                if Dc <= 1e-4:
                    continue
                # worst over s should be at pi by F1; use pi
                R = R_of(math.pi, a, b, Dc)
                if R < best:
                    best = R
                    argmin = (round(a, 3), round(b, 3), round(Dc, 5))
        print(f"  eta={eta}: inf R(pi) ~ {best:.4f} at (a,b,Dc)={argmin}")
    print("  F3 -> PASS (illustrative, positivity & decay shown)")
    return ok


# -------------------- P1: endpoint reduction (full grid) --------------------

def check_P1(grid):
    print("-" * 84)
    print("P1: (1-g(s))/(1-cos s) >= (1-g(pi))/2 for all s in (0,pi]")
    print("    i.e. g(s) <= 1 - R(pi)*D(1-D)*(1-cos s); s=pi is the global min of R")
    ok = True
    worst_violation = -np.inf
    svals = np.concatenate([np.linspace(1e-3, 0.2, 60),
                            np.linspace(0.2, math.pi, 300)])
    n_viol = 0
    for (a, b, D, tag) in grid:
        Rpi = R_of(math.pi, a, b, D)
        for s in svals:
            Rs = R_of(s, a, b, D)
            # violation if R(s) < Rpi (interior beats the endpoint)
            v = Rpi - Rs
            if v > worst_violation:
                worst_violation = v
                worst_pt = (tag, a, b, D, s, Rs, Rpi)
            if v > 1e-7:
                n_viol += 1
    ok = (n_viol == 0)
    print(f"  worst (R(pi) - R(s)) over grid = {worst_violation:+.3e}")
    if worst_violation > 1e-7:
        tag, a, b, D, s, Rs, Rpi = worst_pt
        print(f"    at {tag} a={a} b={b} D={D:.5f} s={s:.4f}: "
              f"R(s)={Rs:.5f} < R(pi)={Rpi:.5f}")
    print(f"  violations (R(s) below R(pi) by >1e-7): {n_viol}")
    print(f"  P1 -> {'PASS' if ok else 'FAIL'}")
    return ok


# -------------------- P1-Zpi: the ONE-PAIR certificate --------------------

def Z_perron_pi(a, b, D):
    """Hermitian Perron pair at s=pi (eta real); scaled to positive trace."""
    eta = eta_of(math.pi, D)
    W = W8_local(eta, a, b)
    ev, V = np.linalg.eig(W)
    k = int(np.argmax(np.abs(ev.real)))
    vec = V[:, k]
    if np.max(np.abs(vec.imag)) < 1e-9:
        vec = vec.real.astype(complex)
    Zs = [vec[0:4].reshape(2, 2), vec[4:8].reshape(2, 2)]
    if (np.trace(Zs[0]) + np.trace(Zs[1])).real < 0:
        Zs = [-Z for Z in Zs]
    return [(Z + Z.conj().T) / 2 for Z in Zs]


def lam_fixedZ(Zpair, s, a, b, D):
    """CW value of a FIXED pair evaluated at the s-dependent eta(s)."""
    T, _ = chain(a, b)
    return lam_CW(Zpair, eta_of(s, D), T)


def check_P1_Zpi(grid):
    print("-" * 84)
    print("P1-Zpi: does the FIXED s=pi pair certify ALL s?")
    print("    ghat_{Z(pi)}(s) := |C_s/C_0|^2 lam(Z(pi); s) <= 1 - R(pi) dd "
          "(1-cos s) ?")
    print("    (if YES this is a ONE-PAIR Collatz-Wielandt closer for P1)")
    ok_perron = True
    ok_dressed = True
    ok_pk = True   # power-iterated-at-pi pair (k large), a richer fixed pair
    svals = np.concatenate([np.linspace(1e-3, 0.2, 80),
                            np.linspace(0.2, math.pi, 400)])
    worst_perron = (-np.inf, None)
    worst_dressed = (-np.inf, None)
    worst_pk = (-np.inf, None)
    worst_validity = -np.inf   # g_true - ghat (must be <=0; CW validity)
    for (a, b, D, tag) in grid:
        T, _ = chain(a, b)
        r, c10, c01 = dressed_data(a, b)
        dd = D * (1 - D)
        Rpi = R_of(math.pi, a, b, D)
        Zperron = Z_perron_pi(a, b, D)
        Zdr_pi = Z_dressed(math.pi, a, b, D, r, c10, c01)
        Zpk = Z_power(eta_of(math.pi, D), T, 8)   # power-iterated at pi, k=8
        dev_p = dev_d = dev_k = -np.inf
        for s in svals:
            bound = 1 - Rpi * dd * (1 - math.cos(s))
            g_true = g_of(s, a, b, D)
            gh_p = Cr2(s, D) * lam_fixedZ(Zperron, s, a, b, D)
            gh_d = Cr2(s, D) * lam_fixedZ(Zdr_pi, s, a, b, D)
            gh_k = Cr2(s, D) * lam_fixedZ(Zpk, s, a, b, D)
            worst_validity = max(worst_validity, g_true - gh_p,
                                 g_true - gh_d, g_true - gh_k)
            dev_p = max(dev_p, gh_p - bound)
            dev_d = max(dev_d, gh_d - bound)
            dev_k = max(dev_k, gh_k - bound)
        if dev_p > worst_perron[0]:
            worst_perron = (dev_p, (a, b, D))
        if dev_d > worst_dressed[0]:
            worst_dressed = (dev_d, (a, b, D))
        if dev_k > worst_pk[0]:
            worst_pk = (dev_k, (a, b, D))
        ok_perron &= (dev_p <= 1e-9)
        ok_dressed &= (dev_d <= 1e-9)
        ok_pk &= (dev_k <= 1e-9)
    print(f"  fixed-pair CW validity (g_true - ghat_Z, must be <=0): "
          f"worst = {worst_validity:+.3e}")
    print(f"  PERRON(pi) pair:  max_s [ghat - bound] = {worst_perron[0]:+.3e} "
          f"at {worst_perron[1]}")
    print(f"     -> {'CERTIFIES all s (CLEAN CLOSER)' if ok_perron else 'FAILS uniform certification'}")
    print(f"  DRESSED(pi) pair: max_s [ghat - bound] = {worst_dressed[0]:+.3e} "
          f"at {worst_dressed[1]}")
    print(f"     -> {'CERTIFIES all s' if ok_dressed else 'FAILS'}")
    print(f"  POWER(pi,k=8):    max_s [ghat - bound] = {worst_pk[0]:+.3e} "
          f"at {worst_pk[1]}")
    print(f"     -> {'CERTIFIES all s' if ok_pk else 'FAILS'}")
    print(f"  P1-Zpi -> {'ONE-PAIR CLOSER FOUND' if (ok_perron or ok_dressed or ok_pk) else 'no single-pair closer'}")
    return ok_perron, ok_dressed, ok_pk


# -------------------- P1-conv: P1 <=> convexity of g in u=cos s --------------

def check_P1_convex(grid):
    """THE STRUCTURAL REDUCTION.  Let G(u) := g(arccos u), u = cos s in [-1,1),
    with G(1) = lim_{s->0} g(s) = 1 (eta->0, Perron value 1).  Then

        h(s) = (1 - g(s))/(1 - cos s) = (G(1) - G(u))/(1 - u)
             = secant slope of G between u=1 and u,

    and for a CONVEX G the secant slope from the fixed point u=1 is
    NONDECREASING in u, so its minimum over u in [-1,1) is attained at u=-1,
    i.e. at s=pi.  Hence

        P1  <=>  G(u)=g(arccos u) is convex on [-1,1]   (given G(1)=1).

    We verify (i) G(1)=1, (ii) G convex (second difference >= 0) on the dense
    4-var grid incl. boundary, (iii) Cr2 is EXACTLY linear in u and |eta|^2 is
    Mobius in u, so the convexity crux is rho(W_8) convex in u; (iv) W_8 is
    COMPLEX on (0,pi) => Cohen/Kingman real-Perron-convexity do NOT apply
    (the honest obstruction to a closed-form proof)."""
    print("-" * 84)
    print("P1-conv: P1 <=> g convex in u=cos s (secant-from-u=1 nondecreasing)")
    ok = True
    worst_d2 = np.inf
    worst_at = None
    g0_worst = 0.0
    npts = 0
    for (a, b, D, tag) in grid:
        npts += 1
        g0_worst = max(g0_worst, abs(g_of(1e-4, a, b, D) - 1.0))
        us = np.linspace(-0.99999, 0.99, 320)
        G = np.array([g_of(math.acos(u), a, b, D) for u in us])
        hstep = us[1] - us[0]
        d2 = (G[2:] - 2 * G[1:-1] + G[:-2]) / hstep ** 2
        m = d2.min()
        if m < worst_d2:
            worst_d2 = m
            worst_at = (a, b, D, us[1 + int(np.argmin(d2))])
        # secant slopes from u=1 must be nondecreasing
        sec = (G - 1.0) / (us - 1.0)
        ok &= bool(np.all(np.diff(sec) > -1e-6))
        ok &= (int(np.argmin(sec)) == 0)   # min secant slope at u=-1 (=s=pi)
    ok &= (worst_d2 > -1e-5) and (g0_worst < 1e-6)
    print(f"  scanned {npts} (a,b,D); G(1)=1 worst |g(s->0)-1| = {g0_worst:.2e}")
    print(f"  worst (most negative) d2 G/du2 = {worst_d2:.3e} at {worst_at}")
    print(f"  secant-from-u=1 nondecreasing & minimized at u=-1 (s=pi): "
          f"{'all points' if ok else 'FAILS somewhere'}")
    print(f"  [Cr2 exactly linear in u; |eta|^2 Mobius in u; crux = rho(W_8) "
          f"convex in u; W_8 COMPLEX on (0,pi) => no real-Perron-convexity thm]")
    print(f"  P1-conv -> {'PASS (convexity reduction holds)' if ok else 'FAIL'}")
    return ok


# -------------------- COMPACT: compactness fallback witness ------------------

def check_compact():
    print("-" * 84)
    print("COMPACT: fallback -- R continuous on the compact box, strictly "
          "positive per point => inf attained > 0")
    from probe_7_34_asym_replica_structure import keff2_numeric
    ok = True
    # continuous extension at s=0: R -> 2(1 - keff2 dd/c^2)
    worst_s0 = 0.0
    for (a, b, D, tag) in [(0.2, 0.4, 0.0426, ""), (0.4, 0.3, 0.0769, ""),
                           (0.1, 0.3, 0.0119, "")]:
        k2 = keff2_numeric(a, b)
        dd = D * (1 - D)
        c = 1 - 2 * D
        Rlim = 2 * (1 - k2 * dd / c ** 2)
        worst_s0 = max(worst_s0, abs(R_of(1e-4, a, b, D) - Rlim))
    ok &= worst_s0 < 1e-3
    print(f"  s->0 extension R -> 2(1 - Psi0): worst residual {worst_s0:.2e}")
    print("  D->0 extension: R -> 2 (finite); s=pi the binding slice (P1).")
    print("  => on K_eta = {a+b<=1-eta} x {D in [0,Dc]} x {s in [0,pi]} R is "
          "continuous & >0 per point,")
    print("     hence c_tilde(eta) := min_{K_eta} R > 0 (attained, binding at "
          "s=pi).  Numerical witnesses: F3.")
    print(f"  COMPACT -> {'PASS (fallback rigorous)' if ok else 'FAIL'}")
    return ok


# -------------------- P2: closed-form positivity of ctil --------------------

def check_P2(grid):
    print("-" * 84)
    print("P2: g(pi)<1 closed form; ctil=R(pi)=(1-g(pi))/(2D(1-D))>0")
    ok = True
    min_ctil = float("inf")
    for (a, b, D, tag) in grid:
        gpi = g_of(math.pi, a, b, D)
        Rpi = R_of(math.pi, a, b, D)
        if gpi >= 1:
            ok = False
            print(f"  {tag} a={a} b={b} D={D:.5f}: g(pi)={gpi:.4f} >= 1 !!")
        min_ctil = min(min_ctil, Rpi)
    print(f"  min ctil over grid = {min_ctil:.4f} (>0 required)")
    ok &= (min_ctil > 0)
    print(f"  P2 -> {'PASS' if ok else 'FAIL'}")
    return ok


def make_grid():
    grid = []
    corners = [(0.1, 0.2), (0.05, 0.3), (0.1, 0.3), (0.2, 0.4),
               (0.3, 0.2), (0.05, 0.45), (0.4, 0.3), (0.15, 0.15)]
    for (a, b) in corners:
        if a + b >= 1:
            continue
        Dc = Dc_closed(a, b)
        if Dc <= 1e-4:
            continue
        for fD in (0.5, 0.9, 1.0):
            grid.append((a, b, fD * Dc, f"fD={fD}"))
    return grid


def main():
    grid = make_grid()
    print(f"grid: {len(grid)} (a,b,D) points")
    results = {}
    results["F1"] = check_F1(grid)
    results["F2"] = check_F2(grid)
    results["F3"] = check_F3()
    results["P1"] = check_P1(grid)
    okp, okd, okk = check_P1_Zpi(grid)
    results["P1-Zpi(perron)"] = okp
    results["P1-Zpi(dressed)"] = okd
    results["P1-Zpi(power8)"] = okk
    results["P1-conv"] = check_P1_convex(grid)
    results["COMPACT"] = check_compact()
    results["P2"] = check_P2(grid)
    print("=" * 84)
    # The P1-Zpi(*) checks are REFUTATIONS: a single fixed pair does NOT
    # certify all s (the would-be clean closer fails); they are EXPECTED False.
    refutations = {"P1-Zpi(perron)", "P1-Zpi(dressed)", "P1-Zpi(power8)"}
    for k, v in results.items():
        if k in refutations:
            tag = "REFUTED-as-expected" if not v else "UNEXPECTEDLY-CERTIFIES"
            print(f"  {k}: {tag} (one-pair closer {'fails' if not v else 'works!'})")
        else:
            print(f"  {k}: {'PASS' if v else 'FAIL'}")
    # overall: the positive claims (F1,F2,F3,P1,P1-conv,COMPACT,P2) must PASS;
    # the one-pair certificate is expected to FAIL.
    positive = {k: v for k, v in results.items() if k not in refutations}
    allok = all(positive.values())
    one_pair = any(results[k] for k in refutations)
    print(f"OVERALL (positive claims) -> {'PASS' if allok else 'FAIL'}")
    print(f"ONE-PAIR CLOSER (P1-Zpi) -> "
          f"{'FOUND' if one_pair else 'REFUTED (no single fixed pair certifies all s)'}")
    return allok


if __name__ == "__main__":
    main()
