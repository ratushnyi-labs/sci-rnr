#!/usr/bin/env python3
"""
probe_7_34_dc_endpoint.py
=============================================================================
The dedicated D = D_c(p) ENDPOINT verification round for the Remark 7.34b
achievability chain (BSMS Gray-region RD dispersion).  The committed theorem
is stated on the OPEN interior 0 < D < D_c conservatively, with the scope
note "the endpoint would need its own verification round".  THIS IS THAT
ROUND: every machine-checkable chain step is run AT D = D_c exactly,
p in {0.1, 0.25, 0.4}.

Checks (PASS/FAIL each):
  DC0  endpoint constants, closed form: A(D), tau_max, delta_c, eps_strip,
       vbar >= (3/4)D(1-D) > 0, Psi_0(D_c) = (1-2p)/(4(1-p)^2) < 1/4,
       Psi_2(D_c) = (1-2p)/(4(1-p)^2) + p^4/(16(1-p)^4) < 5/16, phi < 1,
       c_mu, r_lb, delta_KP_lb all finite/positive; the V9-(iii) orderings.
  DC1  the proven replica spectral margin AT D = D_c exactly:
       g(s) <= 1 - (11/16) D(1-D)(1-cos s) (dense grid + mpmath 60 dps).
  DC2  SLB deconvolution at the MARGINAL discriminant (disc = 0) -- the new
       endpoint content (the pointwise |log2 S_n| = O(n) floor needs
       min_y P_{Y*n}(y) >= c^n, which the strict-interior margin no longer
       supplies at disc = 0):
       (a) at D=D_c the per-pair constants are RATIONAL: c*alpha = (1-p)^2;
           the 2x2 word transfer recursion P_Y(y) = 2^{-n} u^T prod M(s_i) v,
           M(s) = [[1, s], [c a s, a]], a = 1-2p, validated against the brute
           Hadamard deconvolution (n <= 12).
       (b) EXACT exhaustive min_y P_Y(y) > 0 for all n <= 18 (integer DP).
       (c) the MOBIUS-ORBIT FLOOR LEMMA (the new argument): with
           psi(x) = (1 - a x)/(1 - c a x), the worst-case ratio envelope
           x_{i+1} = psi(x_i), x_1 = 1 is monotone up to the double root
           x* = 1/(1-p) (exact identity psi(x) - x = c a (x-x*)^2/(1-c a x)),
           which sits STRICTLY below the pole: c a x* = 1-p < 1.  Hence every
           conditional P(y_{i+1}|y^i) >= p/2 and
               min_y P_Y(y) >= 2^{-n} p^{n-1}   (provable floor AT D_c).
           Verified exactly (Fraction arithmetic): the identity, the orbit,
           the floor vs the exhaustive min, the envelope on sampled words.
       (d) per-site factor n-stability: exact (min)^{1/n} for n <= 18; float
           exhaustive n = 20, 22, 24; exact greedy adversarial word to n=64;
           limit = p/2 (p=0.1: the committed n=12 value ~0.0547 -> 0.05).
       (e) SLB exactly tight at D = D_c: Blahut-Arimoto at n in {8, 10},
           slope ln((1-Dc)/Dc): R_n/n = H_n/n - h(D), D = D_c, and
           Var(j_n) = Var(i_n) (the converse identity at the endpoint).
  DC3  the committed GAP-1 write-out suite (probe_7_34_gap1_writeout_v2)
       run AT (0.25, D_c) and (0.4, D_c) IN FULL (V0-V9, incl. the explicit
       n_0(p, D_c) crossings); (0.1, D_c) reduced: D_c = 0.00307 puts the
       committed floor nD >= 5 first at n = 1631, so V0/V2/V8 run in full
       and the contour identity + remainder run at n = 2048 (t = 6).
  DC4  the (R2) central-chain main loop AT D = D_c (run_main with fD = 1.0)
       for p in {0.25, 0.4} (n to 2048); reduced large-n loop for p = 0.1
       (n in {2048, 4096}, t in {6, 12}: contour q_t, tail, W -> 1/(1-phi)).

VERDICT: dc-endpoint-extends-modulo-floor-lemma / dc-endpoint-ISSUE.
Deps: numpy, mpmath.  Runtime: --fast (DC0-DC2) ~2-4 min; full ~15-25 min.
"""
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib

g1 = importlib.import_module("probe_7_34_gap1_writeout_v2")
r2 = importlib.import_module("probe_7_34_r2_central_chain")
nm = importlib.import_module("remark_7_34b_bsms_rd_dispersion_numeric")

mp.mp.dps = 60

P_LIST = (Fraction(1, 10), Fraction(1, 4), Fraction(2, 5))


def Dc_f(p):
    return 0.5 * (1.0 - math.sqrt(1.0 - 2.0 * p) / (1.0 - p))


def Dc_mp(p):
    p = mp.mpf(p)
    return (1 - mp.sqrt(1 - 2 * p) / (1 - p)) / 2


# =============================== DC0 ===============================

def run_DC0():
    print("-" * 88)
    print("DC0: endpoint constants (closed form) at D = D_c, all three p")
    ok = True
    for pf in P_LIST:
        p = float(pf)
        D = Dc_f(p)
        A = D * (1 - D) / (1 - 2 * D)
        r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
        tau_max = min(math.log(2.0), math.log(1 + (1 - 2 * D) / (2 * D * D)))
        delta_c = min(r_lb / (4 * A), tau_max, 0.5)
        eps_strip = 0.5 * delta_c
        Psi0 = (1 - 2 * p) / (4 * (1 - p) ** 2)
        Psi0_raw = g1.Psi0(p, D)
        vbar = D * (1 - D) * (1 - Psi0_raw)
        Psi2 = Psi0 + p ** 4 / (16 * (1 - p) ** 4)
        k2 = (1 - 2 * p) ** 2 / (p * p * (1 - p) ** 2)
        Psi2_raw = k2 * D * (1 - D) * ((1 - 2 * D) ** 2
                   + 4 * D * D * (1 - D) ** 2) / (1 - 2 * D) ** 4
        phi = D / (1 - D)
        rho = p / (1 - p)
        c_mu = 2 * A * (1 / rho) * (1 / rho - rho)
        dKP = r_lb / A
        checks = dict(
            Dc_lt_half=D < 0.5,
            phi_lt_1=phi < 1.0,
            tau_pos=tau_max > 0,
            delta_c_pos=delta_c > 0,
            strip_pos=eps_strip > 0,
            Psi0_form=abs(Psi0 - Psi0_raw) < 1e-12,
            Psi0_lt_quarter=Psi0 < 0.25,
            Psi2_form=abs(Psi2 - Psi2_raw) < 1e-12,
            Psi2_lt_5_16=Psi2 < 5.0 / 16.0,
            vbar_floor=vbar >= 0.75 * D * (1 - D) - 1e-12,
            c_mu_fin=0 < c_mu < math.inf,
            dKP_pos=dKP > 0,
        )
        ok &= all(checks.values())
        bad = [k for k, v in checks.items() if not v]
        print(f"  p={p}: D_c={D:.6f} A={A:.4f} tau_max={tau_max:.3f} "
              f"delta_c={delta_c:.4f} vbar={vbar:.6f} Psi0={Psi0:.4f} "
              f"Psi2={Psi2:.4f} phi={phi:.4f} c_mu={c_mu:.3f} "
              f"{'ALL OK' if not bad else 'FAIL: ' + ','.join(bad)}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


# =============================== DC1 ===============================

def g_mp(s, p, D):
    s = mp.mpf(s)
    p = mp.mpf(p)
    t = 1 - mp.cos(s)
    F = 1 - 2 * D * (1 - D) * t
    k2 = (1 - 2 * p) ** 2 / (p * p * (1 - p) ** 2)
    G = 8 * k2 * D * D * (1 - D) ** 2 * t * ((1 - 2 * D) ** 2
        + 2 * D * D * (1 - D) ** 2 * t) / (1 - 2 * D) ** 4
    return (F + mp.sqrt(F * F + G)) / 2


def run_DC1():
    print("-" * 88)
    print("DC1: replica margin g(s) <= 1-(11/16)D(1-D)(1-cos s) AT D=D_c exactly")
    ok = True
    for pf in P_LIST:
        p = float(pf)
        D = Dc_f(p)
        s = np.linspace(1e-4, math.pi, 6000)
        viol = float(np.max(g1.g_closed(s, p, D)
                            - (1 - (11.0 / 16.0) * D * (1 - D) * (1 - np.cos(s)))))
        grid_ok = viol <= 1e-12
        # mpmath spot checks at 60 dps with the EXACT algebraic D_c
        Dm = Dc_mp(p)
        worst_mp = mp.mpf(-1)
        for sv in (mp.mpf("1e-3"), mp.mpf("0.5"), mp.mpf("2.0"), mp.pi):
            gv = g_mp(sv, p, Dm)
            bnd = 1 - mp.mpf(11) / 16 * Dm * (1 - Dm) * (1 - mp.cos(sv))
            worst_mp = max(worst_mp, gv - bnd)
        mp_ok = worst_mp < mp.mpf("1e-50")
        ok &= grid_ok and mp_ok
        print(f"  p={p}: grid worst violation = {viol:.2e} ({grid_ok}); "
              f"mpmath60 worst = {mp.nstr(worst_mp, 3)} ({mp_ok})")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


# =============================== DC2 ===============================

def transfer_PY_float(sbits, alpha, ca):
    """P_Y(y)*2^n for one word given sign array s (+-1)."""
    a, b = 1.0, float(sbits[0])
    for s in sbits[1:]:
        a, b = a + ca * s * b, s * a + alpha * b
    return a


def fwht_inplace(v):
    n = len(v)
    h = 1
    while h < n:
        for i in range(0, n, 2 * h):
            for j in range(i, i + h):
                x, y = v[j], v[j + h]
                v[j], v[j + h] = x + y, x - y
        h *= 2
    return v


def hadamard_PY(n, p, D):
    """Brute deconvolution: P_Y over all 2^n words (float64 FWHT)."""
    N = 1 << n
    logP2, _, _ = nm.source_logP2(n, 'bsms', p)
    P = np.exp(logP2 * math.log(2))
    P /= P.sum()
    # FWHT
    Ph = P.copy()
    h = 1
    while h < N:
        for i in range(0, N, 2 * h):
            a = Ph[i:i + h].copy()
            b = Ph[i + h:i + 2 * h].copy()
            Ph[i:i + h] = a + b
            Ph[i + h:i + 2 * h] = a - b
        h *= 2
    w = nm.popcount_table(N).astype(np.float64)
    Ph = Ph / (1.0 - 2.0 * D) ** w
    h = 1
    while h < N:
        for i in range(0, N, 2 * h):
            a = Ph[i:i + h].copy()
            b = Ph[i + h:i + 2 * h].copy()
            Ph[i:i + h] = a + b
            Ph[i + h:i + 2 * h] = a - b
        h *= 2
    return Ph / N


def exhaustive_min_float(alpha, ca, n):
    """min over all words (s_1=+1 wlog) of P_Y(y); returns (min, argmin bits)."""
    m = n - 1
    M = 1 << m
    idx = np.arange(M, dtype=np.uint32)
    a = np.ones(M)
    b = np.ones(M)
    for i in range(m):
        s = 1.0 - 2.0 * ((idx >> np.uint32(i)) & np.uint32(1)).astype(np.float64)
        a, b = a + ca * s * b, s * a + alpha * b
    j = int(np.argmin(a))
    return float(a[j]) / (1 << n), j


def exhaustive_min_exact(alpha_fr, ca_fr, n):
    """EXACT min over all words via DFS with shared prefixes (Fractions
    cleared to integers).  Returns (min_PY as Fraction, argmin sign tuple)."""
    den = (alpha_fr.denominator * ca_fr.denominator
           // math.gcd(alpha_fr.denominator, ca_fr.denominator))
    A_ = alpha_fr.numerator * (den // alpha_fr.denominator)
    C_ = ca_fr.numerator * (den // ca_fr.denominator)
    # state: (depth, a, b); P_Y = 2^{-n} a_n / den^{n-1}
    best = [None, None]
    stack = [(1, 1, 1, (1,))]
    while stack:
        d, a, b, sig = stack.pop()
        if d == n:
            if best[0] is None or a < best[0]:
                best[0] = a
                best[1] = sig
            continue
        for s in (1, -1):
            a2 = den * a + C_ * s * b
            b2 = den * s * a + A_ * b
            stack.append((d + 1, a2, b2, sig + (s,)))
    minPY = Fraction(best[0], (1 << n) * den ** (n - 1))
    return minPY, best[1]


def greedy_word_exact(alpha_fr, ca_fr, n):
    """The adversarial envelope word: s_{i+1} = -sign(r_i).  Exact.
    Returns (P_Y as Fraction, factor list as floats)."""
    a, b = Fraction(1), Fraction(1)
    facs = []
    for _ in range(n - 1):
        s = -1 if b > 0 else 1
        a2 = a + ca_fr * s * b
        b2 = s * a + alpha_fr * b
        facs.append(float(a2 / a))
        a, b = a2, b2
    return a / Fraction(1 << n), facs


def run_DC2():
    print("-" * 88)
    print("DC2: deconvolution at the marginal discriminant (D=D_c): exact "
          "positivity, the Mobius-orbit floor lemma, per-site stability, SLB")
    ok = True
    for pf in P_LIST:
        p = float(pf)
        alpha_fr = 1 - 2 * pf
        ca_fr = (1 - pf) ** 2          # c*alpha at D=D_c -- RATIONAL
        alpha, ca = float(alpha_fr), float(ca_fr)
        D = Dc_f(p)
        print(f"  --- p={pf} (alpha={alpha_fr}, c*alpha=(1-p)^2={ca_fr}, "
              f"D_c={D:.6f}) ---")
        # sanity: c*alpha == (1+alpha)^2/4 == alpha/(1-2Dc)^2
        s1 = abs(ca - (1 + alpha) ** 2 / 4) < 1e-15
        s2 = abs(ca - alpha / (1 - 2 * D) ** 2) < 1e-12
        # (a) transfer == Hadamard, n = 10, 12 (all words, float64)
        a_ok = True
        for n in (10, 12):
            PY = hadamard_PY(n, p, D)
            N = 1 << n
            rngloc = np.random.default_rng(1234 + n)
            sample = list(rngloc.integers(0, N, 200)) + [0, N - 1,
                          int("10" * (n // 2), 2)]
            err = 0.0
            for y in sample:
                sb = [1 - 2 * ((y >> k) & 1) for k in range(n)]
                tv = transfer_PY_float(sb, alpha, ca) / N
                err = max(err, abs(tv - PY[y]))
            a_ok &= err < 1e-12
            print(f"    (a) n={n}: transfer vs Hadamard max|diff| = {err:.2e} "
                  f"on {len(sample)} words (min_y P_Y[FWHT] = {PY.min():.3e})")
        # (b)+(c) exact exhaustive min, exact floor, orbit lemma
        xstar = 1 / (1 - pf)
        pole = 1 / ca_fr
        lem1 = ca_fr * xstar == 1 - pf and 1 - pf < 1      # strict pole margin
        # identity psi(x)-x = ca (x-x*)^2/(1-ca x) at 25 random rationals
        rngloc = np.random.default_rng(7)
        lem2 = True
        for _ in range(25):
            x = Fraction(int(rngloc.integers(0, 10 ** 6)), 10 ** 6) * xstar
            lhs = (1 - alpha_fr * x) / (1 - ca_fr * x) - x
            rhs = ca_fr * (x - xstar) ** 2 / (1 - ca_fr * x)
            lem2 &= lhs == rhs
        # orbit: monotone up, < x*, factors > p
        x = Fraction(1)
        lem3 = True
        orbit = [x]
        for _ in range(80):
            x2 = (1 - alpha_fr * x) / (1 - ca_fr * x)
            lem3 &= x2 >= x and x2 < xstar and 1 - ca_fr * x2 > pf
            x = x2
            orbit.append(x)
        print(f"    (c) orbit lemma: c*a*x* = 1-p = {float(ca_fr * xstar):.3f} < 1 "
              f"({lem1}); identity exact at 25 pts ({lem2}); orbit 1 -> "
              f"{float(orbit[-1]):.6f} < x* = {float(xstar):.6f} monotone, "
              f"factors > p ({lem3})")
        b_ok = True
        fl_ok = True
        per_site = {}
        for n in (8, 12, 16, 18):
            mPY, sig = exhaustive_min_exact(alpha_fr, ca_fr, n)
            floor = Fraction(pf ** (n - 1), 1 << n)
            pos = mPY > 0
            fl = mPY >= floor
            b_ok &= pos
            fl_ok &= fl
            per_site[n] = float(mPY) ** (1.0 / n)
            gP, _ = greedy_word_exact(alpha_fr, ca_fr, n)
            print(f"    (b) n={n}: EXACT min_y P_Y = {float(mPY):.4e} > 0 ({pos}); "
                  f">= 2^-n p^(n-1) = {float(floor):.2e} ({fl}); "
                  f"per-site = {per_site[n]:.5f}; greedy/min = {float(gP / mPY):.4f}")
        # envelope on sampled words, exact: |r_i| <= x_i and factor >= p
        env_ok = True
        n = 18
        for w in range(400):
            rb = np.random.default_rng(w).integers(0, 2, n - 1)
            a, b = Fraction(1), Fraction(1)
            for i, bit in enumerate(rb):
                s = 1 - 2 * int(bit)
                a2 = a + ca_fr * s * b
                b2 = s * a + alpha_fr * b
                env_ok &= abs(b) * 1 <= orbit[i] * a       # |r_i| <= x_i
                env_ok &= a2 * 1 >= pf * a                 # factor >= p
                a, b = a2, b2
        print(f"    (c) envelope |r_i| <= x_i and factor >= p exact on 400 "
              f"random words (n=18): {env_ok}")
        # (d) per-site factor stability: float exhaustive to n=24, greedy to 64
        d_rows = []
        for n in (20, 22, 24):
            mn, _ = exhaustive_min_float(alpha, ca, n)
            d_rows.append((n, mn ** (1.0 / n)))
        g_rows = []
        for n in (32, 48, 64):
            gP, facs = greedy_word_exact(alpha_fr, ca_fr, n)
            g_rows.append((n, float(gP) ** (1.0 / n), min(facs)))
        mono = all(d_rows[i][1] >= d_rows[i + 1][1] - 1e-12
                   for i in range(len(d_rows) - 1))
        toward = abs(g_rows[-1][1] - p / 2) < 0.25 * (per_site[8] - p / 2) + 1e-9
        above = all(v >= p / 2 - 1e-12 for _, v in d_rows) and \
            all(v >= p / 2 - 1e-15 for _, v, _ in g_rows)
        print(f"    (d) per-site (min)^(1/n): exact "
              + " ".join(f"n={n}:{v:.5f}" for n, v in per_site.items())
              + " | float " + " ".join(f"n={n}:{v:.5f}" for n, v in d_rows)
              + " | greedy " + " ".join(f"n={n}:{v:.5f}" for n, v, _ in g_rows)
              + f" -> p/2 = {p / 2:.4f} (decreasing {mono}, converging {toward}, "
              f"all >= p/2 {above})")
        # (e) SLB exactly tight at D_c (Blahut-Arimoto)
        e_ok = True
        for n in (8, 10):
            logP2, _, _ = nm.source_logP2(n, 'bsms', p)
            P = np.exp(logP2 * math.log(2))
            P /= P.sum()
            N = 1 << n
            pc = nm.popcount_table(N)
            xs = np.arange(N, dtype=np.int64)
            Dm = pc[(xs[:, None] ^ xs[None, :])]
            s_slope = math.log((1 - D) / D)
            q, Rn, Dtot, Z = nm.blahut_arimoto(P, Dm, s_slope, iters=4000,
                                               tol=1e-15)
            isur = -logP2
            Hn = float(np.sum(P * isur))
            slb = Hn - n * nm.h2(Dtot / n)
            j = (-s_slope * Dtot - np.log(np.maximum(Z, 1e-300))) / math.log(2)
            Ej = float(np.sum(P * j))
            Vj = float(np.sum(P * (j - Ej) ** 2))
            Ei = float(np.sum(P * isur))
            Vi = float(np.sum(P * (isur - Ei) ** 2))
            tight = abs(Rn - slb) < 1e-6 * n
            dmatch = abs(Dtot / n - D) < 1e-7
            vmatch = abs(Vj / Vi - 1) < 1e-6
            e_ok &= tight and dmatch and vmatch
            print(f"    (e) n={n}: BA at slope ln((1-Dc)/Dc): D/n = {Dtot / n:.8f} "
                  f"(=D_c {dmatch}); R_n - SLB = {Rn - slb:+.2e} (tight {tight}); "
                  f"Var(j)/Var(i) = {Vj / Vi:.10f} ({vmatch})")
        all_p = (s1 and s2 and a_ok and lem1 and lem2 and lem3 and b_ok
                 and fl_ok and env_ok and mono and toward and above and e_ok)
        ok &= all_p
        print(f"    => p={pf}: {'OK' if all_p else 'FAIL'}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


# =============================== DC3 ===============================

def run_DC3(rng):
    print("-" * 88)
    print("DC3: the committed GAP-1 write-out suite AT D = D_c")
    res = {}
    packs_raw = [g1.make_pack(0.25, 1.0), g1.make_pack(0.4, 1.0)]
    res["V0"] = g1.run_V0(rng, packs_raw)
    res["V1"] = g1.run_V1(rng)
    res["V2"], packs = g1.run_V2(rng, packs_raw)
    for P in packs:
        print(f"  [constants@Dc] p={P.p} D={P.D:.5f}: eps_strip={P.eps_strip:.5f} "
              f"eps2={P.eps2:.5f} eps_split={P.eps_split:.5f} C3={P.C3:.4f} "
              f"vbar={P.vbar:.5f} c_tail={P.c_tail:.2e} B3={P.B3:.2f} C={P.C:.1f}")
    res["V3"], res["V4"] = g1.run_V3_V4(rng, packs)
    res["V5"], res["V6"], res["V7"] = g1.run_V5_V6_V7(rng, packs)
    res["V9"] = g1.run_V9(rng, packs)
    # p = 0.1 reduced (nD >= 5 first at n = ceil(5/Dc) = 1631)
    print("  --- p=0.1 at D_c (reduced: in-regime blocklengths only) ---")
    P01_raw = g1.make_pack(0.1, 1.0)
    C3m = g1.measure_C3(rng, P01_raw, verbose=False)
    P01 = g1.make_pack_from(P01_raw, 1.05 * C3m)
    print(f"  [constants@Dc] p=0.1 D={P01.D:.6f}: nD>=5 first at "
          f"n={math.ceil(5 / P01.D)}; C3_meas={C3m:.4f} eps_strip="
          f"{P01.eps_strip:.5f} eps2={P01.eps2:.5f} eps_split={P01.eps_split:.5f} "
          f"vbar={P01.vbar:.6f} B3={P01.B3:.2f} C={P01.C:.1f}")
    n = 2048
    cap = min(P01.eps2, P01.eps_strip / 2)
    words, t = g1.small_beta_words(rng, P01, n, 2, cap)
    red_ok = len(words) > 0
    for (x, qt, mu, v, beta, y0, phi, s) in words:
        central, Vm, Vp, realz, _ = g1.deformed_pieces(x, P01, n, t, mu, v, y0)
        total = (central + Vm + Vp + realz).real / (2 * math.pi)
        rel = abs(total - qt) / abs(qt)
        Rn = qt / ((1 / math.sqrt(2 * math.pi * v)) * math.exp(-beta * beta / 2)) - 1
        bound = P01.C * (1 + abs(beta)) ** 3 / math.sqrt(n)
        c_ok = rel < 2e-3
        r_ok = abs(Rn) <= bound
        red_ok &= c_ok and r_ok
        print(f"  p=0.1 n={n} t={t}: contour identity rel = {rel:.1e} ({c_ok}); "
              f"beta={beta:+.3f} R_n={Rn:+.4f} <= C(1+|beta|)^3/sqrt(n)={bound:.1f} "
              f"({r_ok})")
    res["red01"] = red_ok
    # explicit n0 at the endpoint, all three p
    res["V8"] = g1.run_V8(packs + [P01])
    ok = all(res.values())
    for k, v in res.items():
        print(f"  DC3.{k}: {'PASS' if v else 'FAIL'}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


# =============================== DC4 ===============================

def run_DC4(rng):
    print("-" * 88)
    print("DC4: the (R2) central-chain main loop AT D = D_c")
    ok = True
    ok &= r2.run_main(rng, 0.25, 1.0, [256, 512, 1024, 2048], [120, 100, 70, 45])
    ok &= r2.run_main(rng, 0.4, 1.0, [256, 512, 1024, 2048], [120, 100, 70, 45])
    # p = 0.1 reduced large-n loop (t = 6, 12)
    p = 0.1
    D = Dc_f(p)
    phi_ratio = D / (1 - D)
    Winf = 1 / (1 - phi_ratio)
    print(f"  --- p=0.1 at D_c reduced loop (n=2048, 4096; t=6, 12) ---")
    red_ok = True
    for n, nw in ((2048, 24), (4096, 12)):
        t = int(math.floor(n * D))
        st = [r2.word_stats(x, p, D, t, phi_ratio)
              for x in r2.shell_words(rng, n, p, nw)]
        qpos = all(d["q"][t] > 0 for d in st)
        Ws = np.array([d["W"] for d in st])
        Zs = np.array([d["Z"] for d in st])
        sstar = min(d["sstar"] for d in st)
        eps = 0.9 * sstar
        gbar = 1 - (11 / 16) * D * (1 - D) * (1 - math.cos(eps))
        c4 = 0.25 * math.log(1 / gbar)
        tails = [float(np.max(np.abs(d["phi"][(d["s"] >= eps)
                 & (d["s"] <= math.pi)]))) for d in st]
        t_ok = max(tails) <= 2 * math.exp(-c4 * n)
        W_ok = float(np.max(np.abs(Ws / Winf - 1))) < 0.35
        v_ok = float(np.var(Zs, ddof=1)) < 3.0
        red_ok &= qpos and t_ok and W_ok and v_ok
        print(f"  n={n} t={t}: q_t>0 {qpos}; s*_min={sstar:.3f}; "
              f"tail sup|Phi|={max(tails):.2e} <= 2e^-cn={2 * math.exp(-c4 * n):.2e} "
              f"({t_ok}); max|W/Winf-1|={float(np.max(np.abs(Ws / Winf - 1))):.4f} "
              f"({W_ok}); Var(-log2 S_n)={float(np.var(Zs, ddof=1)):.4f} ({v_ok})")
    ok &= red_ok
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


# =============================== main ===============================

if __name__ == "__main__":
    fast = "--fast" in sys.argv
    seed = 42
    t0 = time.time()
    print("=" * 88)
    print("D = D_c ENDPOINT verification round (Remark 7.34b achievability "
          "chain), p in {0.1, 0.25, 0.4}")
    print("=" * 88)
    rng = np.random.default_rng(seed)
    res = {}
    res["DC0"] = run_DC0()
    res["DC1"] = run_DC1()
    res["DC2"] = run_DC2()
    if not fast:
        res["DC3"] = run_DC3(rng)
        res["DC4"] = run_DC4(rng)
    print("=" * 88)
    for k, v in res.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    allok = all(res.values())
    print(f"[{time.time() - t0:.0f}s] VERDICT: "
          + ("dc-endpoint-extends-modulo-floor-lemma (all endpoint checks pass; "
             "the one NEW ingredient is the Mobius-orbit floor lemma in DC2, "
             "supplying the pointwise |log2 S_n|=O(n) floor at disc=0)"
             if allok else "dc-endpoint-ISSUE"))
