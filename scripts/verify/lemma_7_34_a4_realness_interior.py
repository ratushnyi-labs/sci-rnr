#!/usr/bin/env python3
r"""
lemma_7_34_a4_realness_interior.py
============================================================================
BUG-009-D / Route B (A=4): REALNESS OF THE REPLICA SPECTRUM AT ALL
INTERIOR ANGLES 0 < s < pi, proven on the CURVE-RESTRICTED Gray domain.

THEOREM (certified here, exact integer/rational arithmetic end to end).
  For A = 4, all sg in (0,1), all D with 0 < D <= thetatilde4(sg) *
  Dbar(sg), and ALL spectral angles s in (0, pi), the 5x5 pattern-
  quotient replica matrix Q(p, D, w = e^{is}) has 5 REAL, DISTINCT
  eigenvalues, where
      p = (3/4)(1 - sg^2)                  (all p in (0, 3/4)),
      Dbar(sg)        = (3/4)(1-sg)^2/(1+sg^2),
      thetatilde4(sg) = 1 - sg^2/2  (the committed A=4 separating curve
                                     of lemma_7_34_a4_cert5_on_gray.py).
  Since D_c(p) < thetatilde4 * Dbar for all p (certificate (B) of the
  committed cert5-on-gray lemma), this covers the WHOLE A=4 Gray region
  at every interior angle.  Together with the committed s = pi lemma
  (lemma_7_34_a4_realness_pi.py), realness holds on the curve-restricted
  domain for ALL angles s in (0, pi].

WHY CURVE-RESTRICTED (vs the A=3 full-strip interior lemma
lemma_7_34_realness_interior.py): full-strip realness is FALSE at A=4
ALREADY AT INTERIOR ANGLES, not only at s = pi.  This script EXHIBITS
the complexification island at an interior angle exactly (gate R8b): at
(sg, th, t) = (13/20, 39/40, 39/20) -- strictly inside the open
superdomain strip (th < 1), strictly ABOVE the curve (39/40 > 631/800 =
thetatilde4(13/20)), and at t = 39/20 < 2 (an interior angle s < pi) --
the discriminant DISC below is < 0 (exact rational): two eigenvalues of
Q are genuinely complex there.  The curve restriction is necessary and
sharp-in-kind; the island clears the curve by ~0.15 in th.

REDUCTION TO A POLYNOMIAL DISCRIMINANT (each step gated in-script):
  (1) M := c Q, c = p(1-p) G Gt, G = D^2 w + (D-3)(1-D), Gt = D^2 +
      (D-3)(1-D) w, is a POLYNOMIAL matrix (R1a); chiM = charpoly(M) =
      sum_k (-1)^k E_k X^{5-k} with E_k w-symmetric at the central power
      (-w)^k:  E_k = (-w)^k Etld_k(p,D,t)|_{t=1-(w+1/w)/2} -- an EXACT
      Laurent identity per k (R1c).  On |w| = 1 (t = 1 - cos s):
        psi(Y) := Y^5 + sum_k 18^k Etld_k(p,D,t) Y^{5-k}
      has INTEGER polynomial coefficients (R3; the A=4 denominators are
      powers of 3, so the A=3 scale 2^k becomes 18^k = (2*9)^k) and
      satisfies psi(Y) = (18/w)^5 chiM(wY/18): its roots are exactly
        18 p(1-p) B(D,t) * lambda_i(Q),
      B(D,t) = (4D-3)^2 + 2tD^2(1-D)(3-D) = G Gt / w  (R2).  B > 0 on
      the curve domain: with N := 2(1+sg^2) - u(2-sg^2)(1-sg)^2,
        3 - 4D = (3/2) N/(1+sg^2)   and
        N = (1-u) * 2(1+sg^2) + u * sg * (sg(1-sg)^2 + 4)
      (exact identities, R2b), so N > 0 unless (sg,u) = (0,1); hence
      4D - 3 != 0 and B >= (4D-3)^2 > 0 for all sg in (0,1) -- and the
      multiplier 18 p(1-p) B is strictly positive on the whole claimed
      domain.  So realness/simplicity of Q's spectrum is EQUIVALENT to
      realness/simplicity of psi's roots.
  (2) DISC(p,D,t) := disc_Y(psi): ONE integer polynomial (34173 terms,
      deg (54,78,20) -- two more terms than the A=3 object of the same
      degrees), assembled through the universal quintic discriminant
      (59 terms, isobaric weight 20, R4/R5) and validated twice BEFORE
      any box run: exact rational equality with sympy's univariate
      discriminant at 16 random points (R6), and an mpmath eigenvalue-
      based sign+ratio cross-check at >= 25 points, dps 160 (R7);
      plus a targeted dps-40 adversarial scan of the tight-margin zone
      sg in [0.55,0.85], u in [0.9,1], t in [1.7,2] -- the region where
      the A=4 collision surface approaches the curve to ~1.005 x Dbar
      (R7b; numeric cross-check, not load-bearing).
  (3) t=2 face structure (R8): DISC(p,D,2) == 18^20 disc_X(f4) *
      f4(lam1M)^2, with chiM|_{w=-1} = f1 f4 the committed w=-1 split,
      f1 = X - lam1M, lam1M = (2/3) D(D-1)(4p-3)^2 (2D^2-4D+3) --
      rederived in-script (the committed A=4 s=pi lemma proves realness
      ON this face; here the identity is a structural gate only).
  (4) Rationalized CURVE domain (the committed move, with th = u *
      thetatilde4):  p = (3/4)(1-sg^2),
      D = (3/8) u (2-sg^2)(1-sg)^2/(1+sg^2), (sg,u,t) in
      [0,1]x[0,1]x[0,2].  target := 4^54 8^78 (1+sg^2)^78 DISC(subst)
      is an integer polynomial (R9a, exact Fraction gates), and by
      EXACT division (R9b/R9c)
        target = content * u^14 sg^48 t^7 (1-sg)^34 (1+sg^2)^10
                 (2-sg^2)^14 * core,
      every stripped factor > 0 on
        Omega := {0<sg<1} x {0<u<=1} x {0<t<2}
      (sg=0 <=> p=3/4, sg=1 <=> p=0, u=0 <=> D=0, t=0 <=> s=0 are all
      excluded limits).  Hence sign(disc psi) = sign(core) on Omega and
      the theorem reduces to:  core > 0 on Omega.  core: 143302 terms,
      degs (290, 64, 13).
  (5) CLOSED-BOX ZERO GEOGRAPHY (gated R9d; the exact analogue of the
      A=3 shape): the zero set of core on [0,1]^3-box consists of
        * the identically-zero edge  E := (sg=0, u=1, t in [0,2])
          (R9d-i; the A=3 edge (sg=0, th=1) survives the curve
          substitution since thetatilde4(0) = 1), and
        * the boundary pinch  P1 := (sg,u,t) = (0, 2/3, 2)
          (R9d-ii, the pinch-location identity: core(0,u,2) ==
          (1-u)^4 (3u-2)^4 E54(u) with E54 > 0 on [0,1] -- exact
          division + Sturm; u0 = 2/3 is the A=4 analogue of A=3's
          th0 = 3/4, i.e. the double root of the W-factor of
          disc(f4) on the sg=0 edge).
      Both are OUTSIDE Omega, but a zero at an interior point of an
      edge forces negative Bernstein coefficients at every box whose
      closure contains it: P1 needs a blow-up, exactly as A=3's Z3.
  (6) core > 0 ON Omega -- Bernstein decomposition (all certificates
      exact integer arithmetic; "no-negatives" = all scaled Bernstein
      coefficients >= 0, not all zero => f >= 0 on the closed box AND
      f > 0 on the OPEN box; "closed-strict" = min > 0 AND no exact-
      zero coefficients => f > 0 on the CLOSED box; closed-strict
      pieces may be bisected, that gluing is sound):
        C1   [0,1] x [0,1]       x [0,15/8]        no-negatives
        F1   the u=1 FACE, (sg,t) in [0,1]x[0,2]   no-negatives
             (covers the INCLUDED boundary u=1, all t, via open-face
             positivity -- stronger than the A=3 per-band face gates)
        C2   [0,1] x [0,5/8]     x [7/4,2]         no-negatives
        C3   [0,1] x [17/24,1]   x [7/4,2]         no-negatives
        C4   [1/64,1] x [5/8,17/24] x [7/4,2]      CLOSED-STRICT
        W3a  [0,1/64] x [5/8,17/24] x [7/4,127/64] CLOSED-STRICT
        W3b  [0,1/64] x [5/8,31/48] x [127/64,2]   CLOSED-STRICT
        W3c  [0,1/64] x [11/16,17/24] x [127/64,2] CLOSED-STRICT
        W3d  the blow-up window [0,1/64]x[31/48,11/16]x[127/64,2]
             around the boundary pinch P1 = (0, 2/3, 2).
      (In the recorded run every closed-strict box certifies AT ITS
      ROOT BOX -- zero bisections; the A=4 difficulty was the domain
      geography, not the polynomial, exactly as at A=3.)
  (7) W3d, quasi-homogeneous order-4 blow-up (R13).  Local coordinates
      s = sg, m = u - 2/3, v = 2 - t (window |m| <= 1/48, s <= 1/64,
      v <= 1/64); core_loc := 3^64 * core(s, 2/3+m, 2-v) / g in
      Z[s,m,v] (g = integer content; spot-gated against core).
      VERIFIED: minimal (1,1,2)-weighted degree of core_loc is EXACTLY
      4 (R13b), and the weighted-leading form is
        L4 = c002 v^2 + c021 (m - (4/3)s)^2 v
             + lam4 [ (3m - 4s)^4 + 3072 s^4 ],
      c002, c021, lam4 > 0 -- EXACT identities (R13c; the middle
      coefficient is a perfect square c111^2 = 4 c021 c201 with
      c111/c021 = -8/3, and the v=0 quartic is 81m^4 - 432m^3s +
      864m^2s^2 - 768ms^3 + 3328s^4 == (3m-4s)^4 + 3072 s^4, an
      integer identity).  L4 > 0 for v >= 0, (s,m,v) != 0.  Charts
      (monomial maps, bijective, exact inverse reconstruction gates
      R13d, Schwartz-Zippel identity gates mod 2^61-1 R13e):
        HS :  core_loc(s, s mh, s^2 vh)      = s^4  HS(s, mh, vh)
        HM+/-: core_loc(mu sh, +-mu, mu^2 vh) = mu^4 HM+-(mu, sh, vh)
        HU :  core_loc(nu sh, nu mh, nu^2)   = nu^4 HU(nu, sh, mh)
      certified CLOSED-STRICT on  [0,1/64]x[-1,1]x[0,1],
      [0,1/48]x[0,1]x[0,1] (both signs), [0,1/8]x[0,1]x[-1,1] (R13f).
      SOUNDNESS: for (s,m,v) != 0 in the window let rho = max(s, |m|,
      sqrt(v)); the corresponding chart point lies in its CLOSED box
      and the identity gives core = rho^4 * H(...) > 0.  Only P1
      itself (t = 2, outside Omega) is excluded.
  (8) COVERAGE (R14, exact rational comparisons).  u = 1: face F1
      (open-face > 0 for sg in (0,1), t in (0,2)).  u in (0,1),
      t in (0,15/8): C1 (open box).  u in (0,1), t in [15/8,2):
      u-union (0,5/8] (C2, closed) u [5/8,17/24] u (17/24,1) (C3,
      open box, t-range (7/4,2) contains [15/8,2)), where
      [5/8,17/24] is covered for sg >= 1/64 by C4 (closed) and for
      sg <= 1/64 by W3a (t <= 127/64), W3b/W3c (u outside the
      [31/48,11/16] window), and W3d-charts (the window, sg > 0).
      t-overlap: 7/4 < 15/8 < 127/64; window arithmetic 31/48 =
      2/3 - 1/48, 11/16 = 2/3 + 1/48, 127/64 = 2 - 1/64,
      sqrt(1/64) = 1/8.   =>  core > 0 on all of Omega.
  (9) ANCHOR + CONNECTEDNESS (R15): at (sg,u,t) = (1/2,3/4,1), i.e.
      (p,D) = (9/16, 63/640), s = pi/2: psi has 5 DISTINCT REAL roots
      (exact Sturm) and DISC > 0 (exact).  Omega is path-connected;
      disc psi > 0 on Omega; along any path the real-root count of a
      monic real quintic can change only through a collision
      (disc = 0), which never happens on Omega.  Hence 5 distinct real
      roots everywhere on Omega, i.e. the replica spectrum is REAL and
      SIMPLE.                                                     QED

COMBINED SPECTRAL STATE for A=4 after this certificate: realness holds
on the curve-restricted domain {D <= thetatilde4 * Dbar} -- a superset
of the Gray region -- for all angles s in (0, pi] (this script + the
committed s=pi lemma lemma_7_34_a4_realness_pi.py).  With the committed
A=4 positive-side certificates (chi', ..., chi'''' full-superdomain,
lemma_7_34_a4_chi_certificates.py; cert5-on-Gray,
lemma_7_34_a4_cert5_on_gray.py) and the committed A=4 negative-side
certificates (lemma_7_34_a4_negside_certificates.py), this closes the
realness leg of the A=4 replication of the A=3 spectral program.

ENGINE (packed-key integer tensors, key = a_sg<<20 | b_u<<10 | c_t):
scaled Bernstein coefficients b^_beta = sum_alpha C(d-alpha, beta-alpha)
a_alpha = C(d,beta) b_beta (positive rescale); affine box maps
v -> lo + (hi-lo)v are exact-integer (global positive scale q^deg).
Two implementations of the Bernstein axis map -- the committed A=3 dict
form and a list-based Taylor-shift form (pure additions) -- are cross-
checked coefficient-exactly in R0d; the fast form runs the heavy boxes.
Self-tests R0 pin the map on knowable examples (including the lo != 0
affine map -- the historical h**b bug class -- and exact-zero
accounting).  Semantics identical to the committed engine.

HONEST NOTES / SCOPE:
 *  The claim is the OPEN angle interval (0, pi): s = pi realness is
    the separate committed lemma (the t=2 face of core genuinely
    vanishes at the excluded pinch P1 and along the excluded edge E;
    both are OUTSIDE Omega).
 *  u = 1 (D = thetatilde4 * Dbar) IS included; u = 0 (D = 0), t = 0
    (s = 0), sg in {0,1} (p in {3/4, 0}) are excluded limits where the
    quotient construction itself degenerates (eta(w=1) = 0 forces
    rank 1: the t^7 factor; D = 0 likewise: u^14; p = 3/4 rank 1:
    sg^48; p = 0: the channel is the identity: (1-sg)^34).
 *  NO full-strip claim: R8b exhibits DISC < 0 above the curve at an
    interior angle (t = 39/20).  The A=3 full-strip statement simply
    fails to replicate at A=4; the curve-restricted statement is the
    correct A=4 analogue and still covers the entire Gray region.
 *  Simplicity (distinctness) is part of the certificate: disc > 0.
 *  R7/R7b are floating-point cross-checks of an exact object (belt
    and suspenders on top of the exact R6); every load-bearing
    certificate is exact integer arithmetic.
 *  Runtime ~14 min single-threaded (recorded run: 828 s; the largest
    single certificate is the HU chart, 7130955 Bernstein coefficients),
    peak RSS ~2.5 GB.

Deps: sympy, mpmath; imports the committed pipeline
lemma_7_34_a4_chi_certificates.py (same directory, read-only).
Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools
import os
import pickle
import random
import sys
import time
from fractions import Fraction as Fr
from math import comb, gcd

import sympy as sp
import mpmath as mp

sys.path.insert(0, "/Users/para/work/rnr/scripts/verify")
from lemma_7_34_a4_chi_certificates import (  # noqa: E402
    AVAL, D, p, w, t, sg, th, X, STATES, REPS, _pat,
    build_Q, sym_to_t_mid, PSUB, DSUB)

mp.mp.dps = 60
T0 = time.time()
PASS = True
CACHE = os.environ.get("A4RI_CACHE", "")   # optional dev cache dir
MSCALE = 9                                 # root scale: 2*MSCALE = 18


def rep(name, ok):
    global PASS
    PASS = PASS and bool(ok)
    print(f"  [{time.time()-T0:7.1f}s] {name:<64} {'PASS' if ok else 'FAIL'}",
          flush=True)
    return bool(ok)


def log(msg):
    print(f"  [{time.time()-T0:7.1f}s] {msg}", flush=True)


# ===========================================================================
# packed-key integer polynomial engine
# ===========================================================================
SH3 = (20, 10, 0)
MSK = 1023


def pmul(Aa, Bb):
    out = {}
    g = out.get
    if len(Aa) > len(Bb):
        Aa, Bb = Bb, Aa
    for kA, cA in Aa.items():
        for kB, cB in Bb.items():
            k = kA + kB
            out[k] = g(k, 0) + cA * cB
    return {k: v for k, v in out.items() if v}


def padd_into(acc, Bb, mulB):
    g = acc.get
    for k, v in Bb.items():
        nv = g(k, 0) + mulB * v
        if nv:
            acc[k] = nv
        elif k in acc:
            del acc[k]
    return acc


def ppow(Aa, e, cache):
    if e in cache:
        return cache[e]
    if e == 1:
        cache[1] = Aa
        return Aa
    half = ppow(Aa, e // 2, cache)
    r = pmul(half, half)
    if e % 2:
        r = pmul(r, Aa)
    cache[e] = r
    return r


def affine_axis_int(T, ax, lo, hi, dax):
    """exact integer affine map v -> lo + (hi-lo)v on axis ax; output is
    scaled by q^dax > 0 (q = common denominator) -- sign-preserving.
    v^a -> sum_b C(a,b) (q lo)^(a-b) (q(hi-lo))^b q^(dax-a) v^b."""
    lo = Fr(lo)
    hi = Fr(hi)
    q = lo.denominator * hi.denominator // gcd(lo.denominator,
                                               hi.denominator)
    nl = int(lo * q)
    nh = int((hi - lo) * q)
    PL = [1] * (dax + 1)
    PH = [1] * (dax + 1)
    PQ = [1] * (dax + 1)
    for i in range(1, dax + 1):
        PL[i] = PL[i - 1] * nl
        PH[i] = PH[i - 1] * nh
        PQ[i] = PQ[i - 1] * q
    sh = SH3[ax]
    out = {}
    g = out.get
    for k, c in T.items():
        a = (k >> sh) & MSK
        base = k - (a << sh)
        cq = c * PQ[dax - a]
        for b in range(a + 1):
            f = comb(a, b) * PL[a - b] * PH[b]
            if f == 0:
                continue
            key = base + (b << sh)
            nv = g(key, 0) + f * cq
            if nv:
                out[key] = nv
            elif key in out:
                del out[key]
    return out


def bern_axis(T, ax, dax):
    """scaled Bernstein transform along one axis (box [0,1]) -- the
    committed A=3 dict form."""
    sh = SH3[ax]
    out = {}
    g = out.get
    for k, cval in T.items():
        a2 = (k >> sh) & MSK
        base = k - (a2 << sh)
        for b2 in range(a2, dax + 1):
            f = comb(dax - a2, b2 - a2)
            key = base + (b2 << sh)
            out[key] = g(key, 0) + f * cval
    return {k: v for k, v in out.items() if v}


def bern_axis_fast(T, ax, dax):
    """same map via the reversal/Taylor-shift identity (pure additions on
    dense per-slice lists); cross-checked against bern_axis in R0d."""
    sh = SH3[ax]
    slices = {}
    for k, cval in T.items():
        a2 = (k >> sh) & MSK
        base = k - (a2 << sh)
        s = slices.get(base)
        if s is None:
            s = slices[base] = [0] * (dax + 1)
        s[a2] = cval
    out = {}
    for base, s in slices.items():
        c = s[::-1]
        for j in range(dax):
            for i in range(dax - 1, j - 1, -1):
                c[i] += c[i + 1]
        for b2 in range(dax + 1):
            v = c[dax - b2]
            if v:
                out[base + (b2 << sh)] = v
    return out


def bern_stats(T, degs, box, axes=(2, 1, 0)):
    cur = dict(T)
    for ax, (lo, hi) in enumerate(box):
        if Fr(lo) == 0 and Fr(hi) == 1:
            continue
        cur = affine_axis_int(cur, ax, lo, hi, degs[ax])
    for ax in axes:
        cur = bern_axis_fast(cur, ax, degs[ax])
    total = 1
    for ax in axes:
        total *= (degs[ax] + 1)
    negs = sum(1 for v in cur.values() if v < 0)
    nzeros = total - len(cur)
    return negs, nzeros, len(cur)


def no_negatives(T, degs, box, name, axes=(2, 1, 0)):
    negs, nzeros, nnz = bern_stats(T, degs, box, axes)
    log(f"{name}: nonzero {nnz}, exact zeros {nzeros}, negatives {negs}")
    return negs == 0 and nnz > 0


def closed_strict(T, degs, box, name, depth=0, maxdepth=8, budget=None):
    """closed-strict certificate with sound bisection fallback."""
    if budget is None:
        budget = [80]
    negs, nzeros, nnz = bern_stats(T, degs, box)
    if negs == 0 and nzeros == 0 and nnz > 0:
        if depth == 0:
            log(f"{name}: closed-strict at root box ({nnz} coeffs)")
        return True
    if depth >= maxdepth or budget[0] <= 0:
        log(f"{name}: closed-strict FAILED (negs {negs}, zeros {nzeros}, "
            f"depth {depth})")
        return False
    budget[0] -= 1
    wid = [(Fr(hi) - Fr(lo)) * (degs[i] + 1) for i, (lo, hi) in
           enumerate(box)]
    ax = max(range(3), key=lambda i: wid[i])
    lo, hi = box[ax]
    mid = (Fr(lo) + Fr(hi)) / 2
    b1 = list(box)
    b1[ax] = (lo, mid)
    b2 = list(box)
    b2[ax] = (mid, hi)
    log(f"{name}: bisect axis {ax} at {mid} (negs {negs}, zeros {nzeros})")
    return (closed_strict(T, degs, b1, name + "L", depth + 1, maxdepth,
                          budget)
            and closed_strict(T, degs, b2, name + "R", depth + 1, maxdepth,
                              budget))


def face_u1(T):
    out = {}
    for k, v in T.items():
        key = ((k >> 20) << 20) | (k & 1023)
        out[key] = out.get(key, 0) + v
    return {k: v for k, v in out.items() if v}


# ===========================================================================
# R0: engine self-tests
# ===========================================================================
def engine_selftest():
    # (a) Bernstein coefficients of x^2 on [1/2,1] must be {1/4, 1/2, 1}
    T = {(2 << 20): 1}
    cur = affine_axis_int(T, 0, Fr(1, 2), Fr(1), 2)
    cur = bern_axis(cur, 0, 2)
    got = {(k >> 20): v for k, v in cur.items()}
    okA = (Fr(got.get(0, 0), comb(2, 0) * 4) == Fr(1, 4)
           and Fr(got.get(1, 0), comb(2, 1) * 4) == Fr(1, 2)
           and Fr(got.get(2, 0), comb(2, 2) * 4) == 1)
    rep("R0a engine: Bernstein coeffs of x^2 on [1/2,1] == {1/4,1/2,1}", okA)
    # (b) f = x on [0,1]: exactly one exact-zero coefficient (vertex zero)
    cur = bern_axis({(1 << 20): 1}, 0, 1)
    okB = all(v > 0 for v in cur.values()) and (2 - len(cur)) == 1
    rep("R0b engine: f = x on [0,1] has exactly one exact-zero coeff", okB)
    # (c) 2xy - 1 must show a negative; x + y + 1 must be closed-strict
    cur = {(1 << 20) | (1 << 10): 2, 0: -1}
    for ax, d_ in ((1, 1), (0, 1)):
        cur = bern_axis(cur, ax, d_)
    okC1 = any(v < 0 for v in cur.values())
    cur = {(1 << 20): 1, (1 << 10): 1, 0: 1}
    for ax, d_ in ((1, 1), (0, 1)):
        cur = bern_axis(cur, ax, d_)
    okC2 = all(v > 0 for v in cur.values()) and len(cur) == 4
    rep("R0c engine: 2xy-1 shows negatives; x+y+1 closed-strict",
        okC1 and okC2)
    # (d) fast Taylor-shift axis map == dict axis map on a random tensor
    random.seed(12345)
    T = {}
    for _ in range(200):
        k = (random.randint(0, 9) << 20) | (random.randint(0, 7) << 10) \
            | random.randint(0, 5)
        T[k] = T.get(k, 0) + random.randint(-99, 99)
    T = {k: v for k, v in T.items() if v}
    okD = True
    for ax, dax in ((0, 9), (1, 7), (2, 5)):
        okD &= (bern_axis(T, ax, dax) == bern_axis_fast(T, ax, dax))
    rep("R0d engine: fast Taylor-shift == dict Bernstein axis map", okD)


# ===========================================================================
# R8 helper: w = -1 split and the t=2 face identity (A=4)
# ===========================================================================
def t2_face_gate(DISC):
    Xs = sp.Symbol('X')
    Aval = 4

    def _pat4(tr):
        x0, u_, up = tr
        if x0 == u_ == up:
            return 0
        if x0 == u_ and u_ != up:
            return 1
        if x0 == up and u_ != up:
            return 2
        if u_ == up and x0 != u_:
            return 3
        return 4

    Ew = -D / ((Aval - 1) * (1 - D))          # w = -1
    etaw = sp.cancel(((Aval - 1) * D * Ew + D - (Aval - 1) * Ew)
                     / ((Aval - 1) * D * Ew + D - (Aval - 1)))
    Tm = [[(1 - p) if i == j else p / (Aval - 1) for j in range(Aval)]
          for i in range(Aval)]
    states = list(itertools.product(range(Aval), repeat=3))
    reps_ = [None] * 5
    for k_, tr_ in enumerate(states):
        if reps_[_pat4(tr_)] is None:
            reps_[_pat4(tr_)] = k_
    Qm = sp.zeros(5, 5)
    for a_ in range(5):
        x0, xp_, xq = states[reps_[a_]]
        for (y, yp, yq) in states:
            term = Tm[xp_][yp] * Tm[xq][yq] / Tm[x0][y]
            if yp != y:
                term *= etaw
            if yq != y:
                term *= etaw
            Qm[a_, _pat4((y, yp, yq))] += term
    # c|_{w=-1} = p(p-1)(2D^2-4D+3)^2  (G Gt|_{w=-1} = -(2D^2-4D+3)^2)
    c_ = p * (p - 1) * (2 * D ** 2 - 4 * D + 3) ** 2
    Mm = sp.zeros(5, 5)
    for i in range(5):
        for j in range(5):
            n2, d2 = sp.fraction(sp.cancel(c_ * Qm[i, j]))
            Mm[i, j] = sp.expand(n2 / d2)
    chiM = sp.expand(Mm.charpoly(Xs).as_expr())
    lam1M = sp.expand(sp.Rational(2, 3) * D * (D - 1) * (4 * p - 3) ** 2
                      * (2 * D ** 2 - 4 * D + 3))
    f1p = Xs - lam1M
    q4, r4 = sp.div(chiM, f1p, Xs)
    if sp.expand(r4) != 0:
        return False
    f4p = sp.expand(q4)
    disc4 = sp.discriminant(sp.Poly(f4p, Xs))
    disc4E = sp.expand(disc4.as_expr() if isinstance(disc4, sp.Poly)
                       else disc4)
    Rres = sp.expand(f4p.subs(Xs, lam1M))
    face2 = {}
    for k, v in DISC.items():
        key = (k >> 16, (k >> 8) & 255)
        face2[key] = face2.get(key, 0) + v * (2 ** (k & 255))
    face2E = sum(sp.Integer(v) * p ** a * D ** b
                 for (a, b), v in face2.items() if v)
    return sp.expand(face2E - (2 * MSCALE) ** 20 * disc4E * Rres ** 2) == 0


# ===========================================================================
# R9 helpers: curve substitution and strips
# ===========================================================================
def substitute_curve(DISC, Ap, Bd):
    """target = 4^Ap 8^Bd (1+sg^2)^Bd DISC(p->(3/4)(1-sg^2),
    D->(3/8)u(2-sg^2)(1-sg)^2/(1+sg^2), t) as a packed (sg,u,t) dict."""
    pow1m = [[1]]                       # (1-sg)^m
    for m in range(1, 2 * Bd + 1):
        prev = pow1m[-1]
        cur = [0] * (len(prev) + 1)
        for i, v in enumerate(prev):
            cur[i] += v
            cur[i + 1] -= v
        pow1m.append(cur)
    pow1p2 = [[1]]                      # (1+sg^2)^m
    for m in range(1, Bd + 1):
        prev = pow1p2[-1]
        cur = [0] * (len(prev) + 2)
        for i, v in enumerate(prev):
            cur[i] += v
            cur[i + 2] += v
        pow1p2.append(cur)
    pow2m = [[1]]                       # (2-sg^2)^m
    for m in range(1, Bd + 1):
        prev = pow2m[-1]
        cur = [0] * (len(prev) + 2)
        for i, v in enumerate(prev):
            cur[i] += 2 * v
            cur[i + 2] -= v
        pow2m.append(cur)

    def conv(Aa, Bb):
        out = [0] * (len(Aa) + len(Bb) - 1)
        for i, x in enumerate(Aa):
            if x:
                for j, y in enumerate(Bb):
                    if y:
                        out[i + j] += x * y
        return out

    Vb = [conv(pow2m[b], conv(pow1m[2 * b], pow1p2[Bd - b]))
          for b in range(Bd + 1)]
    groups = {}
    for k, cf in DISC.items():
        a, b, cc = k >> 16, (k >> 8) & 255, k & 255
        groups.setdefault((a, b), []).append((cc, cf))
    P3 = [3 ** m for m in range(Ap + Bd + 1)]
    P2 = [2 ** m for m in range(2 * Ap + 3 * Bd + 1)]
    target = {}
    g = target.get
    for (a, b), lst in groups.items():
        V = Vb[b]
        W = [0] * (2 * a + len(V))
        for i in range(a + 1):          # (1-sg^2)^a
            u_ = comb(a, i) * (-1 if i & 1 else 1)
            base = 2 * i
            for j, y in enumerate(V):
                if y:
                    W[base + j] += u_ * y
        mulc = P3[a + b] * P2[2 * (Ap - a) + 3 * (Bd - b)]
        for (cc, cf) in lst:
            f = cf * mulc
            kb = (b << 10) | cc
            for e, y in enumerate(W):
                if y:
                    k2 = (e << 20) | kb
                    nv = g(k2, 0) + f * y
                    if nv:
                        target[k2] = nv
                    elif k2 in target:
                        del target[k2]
                    g = target.get
    return target


def strip_target(target):
    def to_slices(T):
        S = {}
        for k, v in T.items():
            a, b, cc = k >> 20, (k >> 10) & 1023, k & 1023
            S.setdefault((b, cc), {})[a] = v
        return S

    def from_slices(S):
        T = {}
        for (b, cc), d in S.items():
            for a, v in d.items():
                if v:
                    T[(a << 20) | (b << 10) | cc] = v
        return T

    def div_poly_sg(T, dcoef):
        """exact division of every sg-slice by dcoef (list, low->high);
        integer long division; None if any remainder != 0.  Leading
        coefficient of dcoef must be +-1."""
        dd = len(dcoef) - 1
        lead = dcoef[-1]
        assert lead in (1, -1)
        S = to_slices(T)
        out = {}
        for key, dsl in S.items():
            deg = max(dsl)
            if deg < dd:
                return None
            c = [dsl.get(i, 0) for i in range(deg + 1)]
            q = [0] * (deg - dd + 1)
            for i in range(deg, dd - 1, -1):
                qi = c[i] // lead
                if qi * lead != c[i]:
                    return None
                q[i - dd] = qi
                if qi:
                    for j2, dv in enumerate(dcoef):
                        c[i - dd + j2] -= qi * dv
            if any(c[:dd]):
                return None
            out[key] = {i: v for i, v in enumerate(q) if v}
        return from_slices(out)

    def div_mono(T, which):
        out = {}
        for k, v in T.items():
            a, b, cc = k >> 20, (k >> 10) & 1023, k & 1023
            if which == 'sg':
                if a == 0:
                    return None
                a -= 1
            elif which == 'u':
                if b == 0:
                    return None
                b -= 1
            else:
                if cc == 0:
                    return None
                cc -= 1
            out[(a << 20) | (b << 10) | cc] = v
        return out

    core = dict(target)
    strips = {}
    for name, fn in (("u", lambda T: div_mono(T, 'u')),
                     ("sg", lambda T: div_mono(T, 'sg')),
                     ("t", lambda T: div_mono(T, 't')),
                     ("1-sg", lambda T: div_poly_sg(T, [1, -1])),
                     ("1+sg", lambda T: div_poly_sg(T, [1, 1])),
                     ("1+sg^2", lambda T: div_poly_sg(T, [1, 0, 1])),
                     ("2-sg^2", lambda T: div_poly_sg(T, [2, 0, -1]))):
        while True:
            r = fn(core)
            if r is None:
                break
            core = r
            strips[name] = strips.get(name, 0) + 1
    g = 0
    for v in core.values():
        g = gcd(g, abs(v))
        if g == 1:
            break
    if g > 1:
        core = {k: v // g for k, v in core.items()}
    return core, strips, g


# ===========================================================================
# R13 helpers: the P1 blow-up
# ===========================================================================
def build_local(core, CU):
    """core_loc = 3^CU * core(s, 2/3 + m, 2 - v) as {(i,j,k): int}; also
    returns the stripped integer content g (core_loc*g == 3^CU core(...))."""
    loc = {}
    P3 = [3 ** e for e in range(CU + 1)]
    for k0, v in core.items():
        a, b, cc = k0 >> 20, (k0 >> 10) & 1023, k0 & 1023
        for j in range(b + 1):
            cb = comb(b, j) * (2 ** (b - j)) * P3[CU - (b - j)]
            vb = v * cb
            for kk in range(cc + 1):
                cu = comb(cc, kk) * (2 ** (cc - kk)) * ((-1) ** kk)
                key = (a, j, kk)
                val = loc.get(key, 0) + vb * cu
                if val:
                    loc[key] = val
                elif key in loc:
                    del loc[key]
    g = 0
    for v in loc.values():
        g = gcd(g, abs(v))
        if g == 1:
            break
    if g > 1:
        loc = {k: v // g for k, v in loc.items()}
    return loc, max(g, 1)


def build_chart(loc, kind):
    out = {}
    for (i, j, kk), v in loc.items():
        wd = i + j + 2 * kk - 4
        if kind == 'S':
            key = (wd << 20) | (j << 10) | kk
            sgn = 1
        elif kind == 'M+':
            key = (wd << 20) | (i << 10) | kk
            sgn = 1
        elif kind == 'M-':
            key = (wd << 20) | (i << 10) | kk
            sgn = (-1) ** j
        else:                    # 'U'
            key = (wd << 20) | (i << 10) | j
            sgn = 1
        nv = out.get(key, 0) + sgn * v
        if nv:
            out[key] = nv
        elif key in out:
            del out[key]
    return out


def chart_reconstruct(H, kind):
    """exact inverse of the chart monomial map -> local dict."""
    out = {}
    for k, v in H.items():
        wd, x, y = k >> 20, (k >> 10) & 1023, k & 1023
        if kind == 'S':
            j, kk = x, y
            i = wd + 4 - j - 2 * kk
            sgn = 1
        elif kind in ('M+', 'M-'):
            i, kk = x, y
            j = wd + 4 - i - 2 * kk
            sgn = (-1) ** j if kind == 'M-' else 1
        else:
            i, j = x, y
            r = wd + 4 - i - j
            if r % 2 != 0:
                return None
            kk = r // 2
            sgn = 1
        if i < 0 or j < 0 or kk < 0:
            return None
        out[(i, j, kk)] = sgn * v
    return out


def chart_degs(T):
    return [max(k >> 20 for k in T), max((k >> 10) & 1023 for k in T),
            max(k & 1023 for k in T)]


# ===========================================================================
# main
# ===========================================================================
def main():
    print("=" * 78)
    print("Interior-angle realness, A = 4: disc(psi) > 0 on the rationalized")
    print("CURVE domain for all t in (0,2)  --  exact Bernstein certificates")
    print("=" * 78)
    engine_selftest()

    # ---------------- R1 ----------------
    log("R1 building Q, M = cQ, charpoly, Etld_k ...")
    Q = build_Q()
    G = D ** 2 * w + (D - 3) * (1 - D)
    Gt = D ** 2 + (D - 3) * (1 - D) * w
    c = sp.expand(p * (1 - p) * G * Gt)
    M = sp.zeros(5, 5)
    ok_poly = True
    for i in range(5):
        for j in range(5):
            n_, d_ = sp.fraction(sp.cancel(c * Q[i, j]))
            ok_poly &= d_.is_number
            M[i, j] = sp.expand(n_ / d_)
    rep("R1a M = p(1-p) G Gt Q is a polynomial matrix", ok_poly)
    cp = M.charpoly(X)
    E = [sp.expand(co * (-1) ** k) for k, co in enumerate(cp.all_coeffs())]
    Etld, mks = [sp.Integer(1)], [0]
    for k in range(1, 6):
        ek, mk = sym_to_t_mid(E[k])
        Etld.append(ek)
        mks.append(mk)
    rep("R1b grading: E_k w-symmetric at central power (-w)^k",
        mks == [0, 1, 2, 3, 4, 5] and all(e is not None for e in Etld))
    ok_lau = True
    tsub = 1 - (w + 1 / w) / 2
    for k in range(1, 6):
        diff = sp.expand((E[k] - (-w) ** k * Etld[k].subs(t, tsub)) * w ** k)
        ok_lau &= (sp.expand(diff) == 0)
    rep("R1c exact Laurent identities E_k == (-w)^k Etld_k(t(w))", ok_lau)

    # ---------------- R2 ----------------
    Bpoly = (4 * D - 3) ** 2 + 2 * t * D ** 2 * (1 - D) * (3 - D)
    chat, mc = sym_to_t_mid(c)
    ok2 = (mc == 1 and sp.expand(chat + p * (1 - p) * Bpoly) == 0)
    ok2 &= (sp.expand(Gt - w * G.subs(w, 1 / w)) == 0)
    rep("R2a c = w p(1-p) B(D,t), B = (4D-3)^2 + 2tD^2(1-D)(3-D); "
        "Gt = w G(1/w)", ok2)
    # B > 0 on the curve domain: 3 - 4D = (3/2) N/(1+sg^2), N linear in u
    uu = sp.Symbol('u')
    DSUBC = sp.cancel(DSUB.subs(th, uu * (1 - sg ** 2 / 2)))
    ok2b = (sp.cancel(DSUBC - sp.Rational(3, 8) * uu * (2 - sg ** 2)
                      * (1 - sg) ** 2 / (1 + sg ** 2)) == 0)
    N = 2 * (1 + sg ** 2) - uu * (2 - sg ** 2) * (1 - sg) ** 2
    ok2b &= (sp.cancel((3 - 4 * DSUBC) * (1 + sg ** 2)
                       - sp.Rational(3, 2) * N) == 0)
    ok2b &= (sp.expand(N - ((1 - uu) * 2 * (1 + sg ** 2)
                            + uu * sg * (sg * (1 - sg) ** 2 + 4))) == 0)
    rep("R2b curve domain: D = (3/8)u(2-sg^2)(1-sg)^2/(1+sg^2); "
        "3-4D > 0 for sg>0", ok2b)

    # ---------------- R3 ----------------
    Ehat = [None]
    ok3 = True
    for k in range(1, 6):
        Pk = sp.Poly(sp.expand((2 * MSCALE) ** k * Etld[k]), p, D, t)
        Tk = {}
        for mon, cc in zip(Pk.monoms(), Pk.coeffs()):
            ok3 &= bool(cc.is_Integer)
            Tk[(mon[0] << 16) | (mon[1] << 8) | mon[2]] = int(cc)
        Ehat.append(Tk)
    rep("R3 Ehat_k = 18^k Etld_k have integer coefficients", ok3)

    # ---------------- R4 ----------------
    a1, a2_, a3, a4, a5, Z = sp.symbols('a1 a2 a3 a4 a5 Z')
    disc5 = sp.discriminant(sp.Poly(
        Z ** 5 + a1 * Z ** 4 + a2_ * Z ** 3 + a3 * Z ** 2 + a4 * Z + a5, Z))
    P5 = sp.Poly(disc5, a1, a2_, a3, a4, a5)
    terms5 = [(mon, int(cc)) for mon, cc in zip(P5.monoms(), P5.coeffs())]
    wts = set(sum((i + 1) * e for i, e in enumerate(mon))
              for mon, _ in terms5)
    rep(f"R4 universal quintic disc: {len(terms5)} terms, weight {wts}",
        len(terms5) == 59 and wts == {20})

    # ---------------- R5 ----------------
    cachef = os.path.join(CACHE, "disc_a4.pkl") if CACHE else None
    if cachef and os.path.exists(cachef):
        with open(cachef, "rb") as f:
            DISC = pickle.load(f)
            if isinstance(DISC, tuple):
                DISC = DISC[0]
        log(f"R5 loaded cached DISC ({len(DISC)} terms)")
    else:
        log("R5 assembling DISC = disc_Y(psi) (universal formula) ...")
        caches = [None] + [dict() for _ in range(5)]
        DISC = {}
        for mon, cc in terms5:
            factors = []
            for i, e in enumerate(mon):
                if e:
                    factors.append(ppow(Ehat[i + 1], e, caches[i + 1]))
            factors.sort(key=len)
            prod = {0: 1}
            for f in factors:
                prod = pmul(prod, f)
            padd_into(DISC, prod, cc)
        del caches
        if cachef:
            with open(cachef, "wb") as f:
                pickle.dump(DISC, f)
    dP = max(k >> 16 for k in DISC)
    dD = max((k >> 8) & 255 for k in DISC)
    dT = max(k & 255 for k in DISC)
    rep(f"R5 DISC: {len(DISC)} terms, degs ({dP},{dD},{dT})",
        (dP, dD, dT) == (54, 78, 20) and len(DISC) == 34173)

    # ---------------- R6 ----------------
    log("R6 exact univariate-disc validation, 16 random rational points")
    random.seed(20260704)
    ok6 = True
    Y = sp.Symbol('Y')
    for _ in range(16):
        pv = sp.Rational(random.randint(1, 63), 96)   # p in (0, 3/4)
        Dv = sp.Rational(random.randint(1, 63), 100)
        tv = sp.Rational(random.randint(0, 32), 16)
        psi = Y ** 5
        for k in range(1, 6):
            psi += (sp.Integer(2 * MSCALE) ** k * Etld[k].subs(
                {p: pv, D: Dv, t: tv})) * Y ** (5 - k)
        dref = sp.discriminant(sp.Poly(psi, Y))
        dval = sum(Fr(cc) * Fr(int(pv.p), int(pv.q)) ** (k >> 16)
                   * Fr(int(Dv.p), int(Dv.q)) ** ((k >> 8) & 255)
                   * Fr(int(tv.p), int(tv.q)) ** (k & 255)
                   for k, cc in DISC.items())
        ok6 &= (Fr(int(sp.fraction(dref)[0]), int(sp.fraction(dref)[1]))
                == dval)
    rep("R6 DISC == disc_Y(psi) exactly at 16 random rational points", ok6)

    # ---------------- R7 ----------------
    log("R7 mpmath eigen-disc cross-check (sign + ratio), dps 160")

    def Q_hp(pv, Dv, sv):
        """5x5 pattern-quotient Q at full mp precision (float inputs are
        promoted BEFORE any arithmetic -- a double-precision leak here once
        capped the achievable ratio agreement at ~1e-11)."""
        mpp = mp.mpf(pv)
        mpD = mp.mpf(Dv)
        Ee = (mpD / ((AVAL - 1) * (1 - mpD))) * mp.e ** (1j * sv)
        et = ((AVAL - 1) * mpD * Ee + mpD - (AVAL - 1) * Ee) / \
            ((AVAL - 1) * mpD * Ee + mpD - (AVAL - 1))
        etb = mp.conj(et)
        Tn = [[(1 - mpp) if i == j else mpp / (AVAL - 1)
               for j in range(AVAL)] for i in range(AVAL)]
        Qn = mp.zeros(5, 5)
        for a2 in range(5):
            x_, xp, xq = STATES[REPS[a2]]
            for (y, yp, yq) in STATES:
                Qn[a2, _pat((y, yp, yq))] += (Tn[xp][yp] * Tn[xq][yq]
                                              / Tn[x_][y]
                                              * (et if yp != y else 1)
                                              * (etb if yq != y else 1))
        return Qn

    random.seed(7)
    ok7 = True
    npts = 0
    worst = mp.mpf(0)
    dps_save = mp.mp.dps
    mp.mp.dps = 160
    THT4 = 1 - sg ** 2 / 2
    for _ in range(60):
        sgv = random.uniform(0.02, 0.98)
        uv = random.uniform(0.05, 1.0)
        tv = random.uniform(0.001, 1.999)
        pv = float(PSUB.subs(sg, sp.Float(sgv)))
        thv = uv * float(THT4.subs(sg, sp.Float(sgv)))
        Dv = float(DSUB.subs({sg: sp.Float(sgv), th: sp.Float(thv)}))
        if Dv < 1e-6:
            continue
        sv = mp.acos(1 - mp.mpf(tv))   # t <-> s link at working precision
        Qn = Q_hp(pv, Dv, sv)
        ev = mp.eig(Qn)[0]
        Gv = mp.mpf(Dv) ** 2 * mp.e ** (1j * sv) \
            + (mp.mpf(Dv) - (AVAL - 1)) * (1 - mp.mpf(Dv))
        rho = mp.mpf(pv) * (1 - mp.mpf(pv)) * abs(Gv) ** 2
        rts = [2 * MSCALE * rho * e for e in ev]
        dnum = mp.mpf(1)
        min_relgap = mp.mpf(10)
        for i in range(5):
            for j in range(i + 1, 5):
                dnum *= (rts[i] - rts[j]) ** 2
                rg = abs(rts[i] - rts[j]) / (abs(rts[i]) + abs(rts[j])
                                             + mp.mpf('1e-40'))
                min_relgap = min(min_relgap, rg)
        dnum = mp.re(dnum)
        if min_relgap < mp.mpf('1e-3'):   # conditioning filter
            continue
        pw_ = [mp.mpf(1)] * (dP + 1)
        Dw_ = [mp.mpf(1)] * (dD + 1)
        tw_ = [mp.mpf(1)] * (dT + 1)
        for i in range(1, dP + 1):
            pw_[i] = pw_[i - 1] * pv
        for i in range(1, dD + 1):
            Dw_[i] = Dw_[i - 1] * Dv
        for i in range(1, dT + 1):
            tw_[i] = tw_[i - 1] * tv
        dval = mp.mpf(0)
        for k, cc in DISC.items():
            dval += mp.mpf(cc) * pw_[k >> 16] * Dw_[(k >> 8) & 255] \
                * tw_[k & 255]
        npts += 1
        rerr = abs(dval / dnum - 1)
        worst = max(worst, rerr)
        ok7 &= (rerr < mp.mpf('1e-30')) and ((dval > 0) == (dnum > 0))
    mp.mp.dps = dps_save
    rep(f"R7 sign+ratio agree at {npts} pts (worst |ratio-1| ~ "
        f"{mp.nstr(worst, 3)})", ok7 and npts >= 25)

    # ---------------- R8 ----------------
    log("R8 t=2 face identity (w = -1 split rederived in-script) ...")
    rep("R8 DISC(p,D,2) == 18^20 disc(f4) f4(lam1M)^2", t2_face_gate(DISC))

    # ---------------- R8b: the interior-angle island ----------------
    sgv = Fr(13, 20)
    thv = Fr(39, 40)
    tv = Fr(39, 20)
    pv = Fr(3, 4) * (1 - sgv ** 2)
    Dbv = Fr(3, 4) * (1 - sgv) ** 2 / (1 + sgv ** 2)
    Dv = thv * Dbv
    curve = 1 - sgv ** 2 / 2
    dI = sum(Fr(cc) * pv ** (k >> 16) * Dv ** ((k >> 8) & 255)
             * tv ** (k & 255) for k, cc in DISC.items())
    rep("R8b island at INTERIOR angle: DISC < 0 at (13/20, 39/40, 39/20), "
        "th > curve", dI < 0 and thv > curve and thv < 1 and tv < 2)

    # ---------------- R9 ----------------
    log("R9 curve-domain substitution -> target(sg,u,t)")
    target = substitute_curve(DISC, dP, dD)
    random.seed(31)
    ok9a = True
    for _ in range(6):
        sgv = Fr(random.randint(1, 15), 16)
        uv = Fr(random.randint(1, 16), 16)
        tv = Fr(random.randint(0, 32), 16)
        pv = Fr(3, 4) * (1 - sgv ** 2)
        Dv = Fr(3, 8) * uv * (2 - sgv ** 2) * (1 - sgv) ** 2 \
            / (1 + sgv ** 2)
        dv = sum(Fr(cf) * pv ** (k >> 16) * Dv ** ((k >> 8) & 255)
                 * tv ** (k & 255) for k, cf in DISC.items())
        tv2 = sum(Fr(cf) * sgv ** (k >> 20) * uv ** ((k >> 10) & 1023)
                  * tv ** (k & 1023) for k, cf in target.items())
        ok9a &= (tv2 == Fr(4) ** dP * Fr(8) ** dD * (1 + sgv ** 2) ** dD
                 * dv)
    rep("R9a target == 4^54 8^78 (1+sg^2)^78 DISC(subst): 6 exact points",
        ok9a)
    core, strips, content = strip_target(target)
    log(f"R9 strips: {strips}; content bits {content.bit_length()}")
    ok9b = (strips.get('u', 0) == 14 and strips.get('sg', 0) == 48
            and strips.get('t', 0) == 7 and strips.get('1-sg', 0) == 34
            and strips.get('1+sg^2', 0) == 10
            and strips.get('2-sg^2', 0) == 14
            and set(strips) <= {'u', 'sg', 't', '1-sg', '1+sg', '1+sg^2',
                                '2-sg^2'}
            and content > 0)
    rep("R9b target == content u^14 sg^48 t^7 (1-sg)^34 (1+sg^2)^10 "
        "(2-sg^2)^14 core", ok9b)
    CS = max(k >> 20 for k in core)
    CU = max((k >> 10) & 1023 for k in core)
    CT = max(k & 1023 for k in core)
    degs = [CS, CU, CT]
    log(f"R9 core: {len(core)} terms, degs ({CS},{CU},{CT})")
    random.seed(41)
    ok9c = True
    for _ in range(4):
        sgv = Fr(random.randint(1, 15), 16)
        uv = Fr(random.randint(1, 16), 16)
        tv = Fr(random.randint(1, 31), 16)
        cv = sum(Fr(cf) * sgv ** (k >> 20) * uv ** ((k >> 10) & 1023)
                 * tv ** (k & 1023) for k, cf in core.items())
        tv2 = sum(Fr(cf) * sgv ** (k >> 20) * uv ** ((k >> 10) & 1023)
                  * tv ** (k & 1023) for k, cf in target.items())
        recon = (content * uv ** 14 * sgv ** 48 * tv ** 7
                 * (1 - sgv) ** 34 * (1 + sgv ** 2) ** 10
                 * (2 - sgv ** 2) ** 14 * cv)
        ok9c &= (recon == tv2)
    rep("R9c strip reconstruction identity at 4 exact points", ok9c)
    del target

    # ---------------- R9d: zero geography gates ----------------
    # (i) the edge (sg=0, u=1): core identically zero along it
    edgeE = {}
    for k, v in core.items():
        if (k >> 20) == 0:
            edgeE[k & 1023] = edgeE.get(k & 1023, 0) + v
    ok9d1 = all(v == 0 for v in edgeE.values())
    rep("R9d-i edge E = (sg=0, u=1, t): core identically zero", ok9d1)
    # (ii) pinch location: core(0,u,2) == (1-u)^4 (3u-2)^4 E54(u), E54 > 0
    uY = sp.Symbol('u')
    E2 = {}
    for k, v in core.items():
        if (k >> 20) == 0:
            a = (k >> 10) & 1023
            E2[a] = E2.get(a, 0) + v * (2 ** (k & 1023))
    Ep = sp.Poly(sum(sp.Integer(v) * uY ** a for a, v in E2.items()), uY)
    q54, r54 = sp.div(Ep, sp.Poly((1 - uY) ** 4 * (3 * uY - 2) ** 4, uY))
    ok9d2 = (r54 == 0)
    if ok9d2:
        Q54 = sp.Poly(q54, uY)
        ok9d2 &= (Q54.count_roots(0, 1) == 0
                  and Q54.eval(sp.Rational(1, 2)) > 0)
    rep("R9d-ii core(0,u,2) == (1-u)^4 (3u-2)^4 E54(u), E54 > 0 on [0,1]",
        ok9d2)

    # ---------------- R7b: adversarial numeric scan (cross-check) -------
    log("R7b dps-40 adversarial scan: sg in [0.55,0.85], u in [0.9,1], "
        "t in [1.7,2]")
    dps_save = mp.mp.dps
    mp.mp.dps = 40
    SL = {}
    for k, v in core.items():
        a, b, cc = k >> 20, (k >> 10) & 1023, k & 1023
        SL.setdefault((b, cc), []).append((a, v))
    neg7b = 0
    cnt7b = 0
    for iu in range(6):
        uv = mp.mpf('0.9') + mp.mpf('0.02') * iu
        for it in range(7):
            tv = mp.mpf('1.7') + mp.mpf('0.05') * it
            pu = [mp.mpf(1)] * (CU + 1)
            pt = [mp.mpf(1)] * (CT + 1)
            for i in range(1, CU + 1):
                pu[i] = pu[i - 1] * uv
            for i in range(1, CT + 1):
                pt[i] = pt[i - 1] * tv
            cs = [mp.mpf(0)] * (CS + 1)
            for (b, cc), lst in SL.items():
                wgt = pu[b] * pt[cc]
                for a, v in lst:
                    cs[a] += v * wgt
            for isg in range(31):
                sgv = mp.mpf('0.55') + mp.mpf('0.01') * isg
                acc = mp.mpf(0)
                for c2 in reversed(cs):
                    acc = acc * sgv + c2
                cnt7b += 1
                if acc <= 0:
                    neg7b += 1
    del SL
    mp.mp.dps = dps_save
    rep(f"R7b core > 0 at all {cnt7b} tight-margin points (numeric "
        f"cross-check)", neg7b == 0 and cnt7b == 1302)

    # ---------------- R10: C1 + the u=1 face ----------------
    ok10 = no_negatives(core, degs, [(0, 1), (0, 1), (0, Fr(15, 8))],
                        "R10 C1 [0,1]^2x[0,15/8]")
    rep("R10 C1 no negatives (core > 0 on the open C1 box)", ok10)
    fu1 = face_u1(core)
    ok10f = no_negatives(fu1, [CS, 0, CT], [(0, 1), (0, 1), (0, 2)],
                         "R10f u=1 face, t in [0,2]", axes=(2, 0))
    rep("R10f u=1 face, ALL t in [0,2]: no negatives (open-face > 0)",
        ok10f)
    del fu1

    # ---------------- R11: C2/C3/C4 ----------------
    ok11a = no_negatives(core, degs,
                         [(0, 1), (0, Fr(5, 8)), (Fr(7, 4), 2)],
                         "R11a C2 [0,1]x[0,5/8]x[7/4,2]")
    rep("R11a C2 no negatives", ok11a)
    ok11b = no_negatives(core, degs,
                         [(0, 1), (Fr(17, 24), 1), (Fr(7, 4), 2)],
                         "R11b C3 [0,1]x[17/24,1]x[7/4,2]")
    rep("R11b C3 no negatives", ok11b)
    ok11c = closed_strict(core, degs,
                          [(Fr(1, 64), 1), (Fr(5, 8), Fr(17, 24)),
                           (Fr(7, 4), 2)],
                          "R11c C4 [1/64,1]x[5/8,17/24]x[7/4,2]")
    rep("R11c C4 CLOSED-STRICT", ok11c)

    # ---------------- R12: wedge boxes ----------------
    ok12a = closed_strict(core, degs,
                          [(0, Fr(1, 64)), (Fr(5, 8), Fr(17, 24)),
                           (Fr(7, 4), Fr(127, 64))], "R12a W3a")
    rep("R12a W3a [0,1/64]x[5/8,17/24]x[7/4,127/64] CLOSED-STRICT", ok12a)
    ok12b = closed_strict(core, degs,
                          [(0, Fr(1, 64)), (Fr(5, 8), Fr(31, 48)),
                           (Fr(127, 64), 2)], "R12b W3b")
    rep("R12b W3b [0,1/64]x[5/8,31/48]x[127/64,2] CLOSED-STRICT", ok12b)
    ok12c = closed_strict(core, degs,
                          [(0, Fr(1, 64)), (Fr(11, 16), Fr(17, 24)),
                           (Fr(127, 64), 2)], "R12c W3c")
    rep("R12c W3c [0,1/64]x[11/16,17/24]x[127/64,2] CLOSED-STRICT", ok12c)

    # ---------------- R13: the P1 blow-up ----------------
    log("R13 local dict at P1 = (0, 2/3, 2) ...")
    loc, gloc = build_local(core, CU)
    # spot-gate the local dict against core at 2 exact rational points
    random.seed(97)
    okloc = True
    for _ in range(2):
        sv_ = Fr(random.randint(1, 9), 64)
        mv_ = Fr(random.randint(-9, 9), 100)
        vv_ = Fr(random.randint(1, 9), 64)
        lv = sum(Fr(cf) * sv_ ** i * mv_ ** j * vv_ ** kk
                 for (i, j, kk), cf in loc.items())
        cv = sum(Fr(cf) * sv_ ** (k >> 20)
                 * (Fr(2, 3) + mv_) ** ((k >> 10) & 1023)
                 * (2 - vv_) ** (k & 1023) for k, cf in core.items())
        okloc &= (Fr(gloc) * lv == Fr(3) ** CU * cv)
    rep("R13a g*core_loc == 3^64 core(s, 2/3+m, 2-v): 2 exact points",
        okloc)
    del core            # memory: R14/R15 use DISC/Etld only
    wmin = min(i + j + 2 * kk for (i, j, kk) in loc)
    rep(f"R13b minimal (1,1,2)-weighted degree of core_loc == 4 "
        f"(got {wmin})", wmin == 4)
    c002 = loc.get((0, 0, 2), 0)
    c021 = loc.get((0, 2, 1), 0)
    c111 = loc.get((1, 1, 1), 0)
    c201 = loc.get((2, 0, 1), 0)
    c040 = loc.get((0, 4, 0), 0)
    c130 = loc.get((1, 3, 0), 0)
    c220 = loc.get((2, 2, 0), 0)
    c310 = loc.get((3, 1, 0), 0)
    c400 = loc.get((4, 0, 0), 0)
    # v=0 quartic must be lam4*(81 m^4 - 432 m^3 s + 864 m^2 s^2
    #                           - 768 m s^3 + 3328 s^4)
    okq = (c040 % 81 == 0)
    lam4 = c040 // 81 if okq else 0
    okq = (okq and lam4 > 0
           and c130 == -432 * lam4 and c220 == 864 * lam4
           and c310 == -768 * lam4 and c400 == 3328 * lam4)
    # integer identity: 81m^4-432m^3s+864m^2s^2-768ms^3+3328s^4
    #                   == (3m-4s)^4 + 3072 s^4
    sQ, mQ = sp.symbols('sQ mQ')
    okq &= (sp.expand(81 * mQ ** 4 - 432 * mQ ** 3 * sQ
                      + 864 * mQ ** 2 * sQ ** 2 - 768 * mQ * sQ ** 3
                      + 3328 * sQ ** 4
                      - ((3 * mQ - 4 * sQ) ** 4 + 3072 * sQ ** 4)) == 0)
    okL4 = (c002 > 0 and c021 > 0
            and c111 ** 2 == 4 * c021 * c201
            and 3 * c111 == -8 * c021
            and okq)
    rep("R13c L4 = c002 v^2 + c021 (m-(4/3)s)^2 v + lam4[(3m-4s)^4 + "
        "3072 s^4]", okL4)
    charts = {}
    ok_coll = True
    for kind in ('S', 'M+', 'M-', 'U'):
        H = build_chart(loc, kind)
        charts[kind] = H
        ok_coll &= (len(H) == len(loc))
        back = chart_reconstruct(H, kind)
        ok_coll &= (back == loc)
    rep("R13d chart monomial maps bijective + exact inverse "
        "reconstruction", ok_coll)
    random.seed(5)
    qP = (1 << 61) - 1

    def loc_eval_q(s_, m_, u_):
        ps = {}
        pm = {}
        pu = {}
        acc = 0
        for (i, j, kk), v in loc.items():
            if i not in ps:
                ps[i] = pow(s_, i, qP)
            if j not in pm:
                pm[j] = pow(m_, j, qP)
            if kk not in pu:
                pu[kk] = pow(u_, kk, qP)
            acc = (acc + v * ps[i] % qP * pm[j] % qP * pu[kk]) % qP
        return acc

    def dict_eval_q(T, x, y, z):
        px = {}
        py = {}
        pz = {}
        acc = 0
        for k, v in T.items():
            i, j, kk = k >> 20, (k >> 10) & 1023, k & 1023
            if i not in px:
                px[i] = pow(x, i, qP)
            if j not in py:
                py[j] = pow(y, j, qP)
            if kk not in pz:
                pz[kk] = pow(z, kk, qP)
            acc = (acc + v * px[i] % qP * py[j] % qP * pz[kk]) % qP
        return acc

    # Schwartz-Zippel identity gates mod 2^61 - 1: the difference
    # polynomials have total degree <= ~800, so a nonzero difference
    # survives 4 uniform points with probability <= (800/q)^4 ~ 1e-70.
    ok_id = True
    for _ in range(4):
        s_ = random.randrange(qP)
        mh = random.randrange(qP)
        uh = random.randrange(qP)
        sh_ = random.randrange(qP)
        mu_ = random.randrange(qP)
        nu_ = random.randrange(qP)
        ok_id &= (loc_eval_q(s_, s_ * mh % qP, s_ * s_ % qP * uh % qP)
                  == pow(s_, 4, qP) * dict_eval_q(charts['S'], s_, mh, uh)
                  % qP)
        ok_id &= (loc_eval_q(mu_ * sh_ % qP, mu_,
                             mu_ * mu_ % qP * uh % qP)
                  == pow(mu_, 4, qP) * dict_eval_q(charts['M+'], mu_, sh_,
                                                   uh) % qP)
        ok_id &= (loc_eval_q(mu_ * sh_ % qP, (qP - mu_) % qP,
                             mu_ * mu_ % qP * uh % qP)
                  == pow(mu_, 4, qP) * dict_eval_q(charts['M-'], mu_, sh_,
                                                   uh) % qP)
        ok_id &= (loc_eval_q(nu_ * sh_ % qP, nu_ * mh % qP,
                             nu_ * nu_ % qP)
                  == pow(nu_, 4, qP) * dict_eval_q(charts['U'], nu_, sh_,
                                                   mh) % qP)
    rep("R13e chart blow-up identities (Schwartz-Zippel mod 2^61-1, "
        "4 pts)", ok_id)
    del loc                 # memory: charts carry everything from here
    jobs = [
        ("R13f HS  [0,1/64]x[-1,1]x[0,1]", 'S',
         [(0, Fr(1, 64)), (-1, 1), (0, 1)]),
        ("R13f HM+ [0,1/48]x[0,1]x[0,1]", 'M+',
         [(0, Fr(1, 48)), (0, 1), (0, 1)]),
        ("R13f HM- [0,1/48]x[0,1]x[0,1]", 'M-',
         [(0, Fr(1, 48)), (0, 1), (0, 1)]),
        ("R13f HU  [0,1/8]x[0,1]x[-1,1]", 'U',
         [(0, Fr(1, 8)), (0, 1), (-1, 1)]),
    ]
    for name, kind, box in jobs:
        H = charts[kind]
        dd = chart_degs(H)
        okj = closed_strict(H, dd, box, name)
        rep(f"{name} CLOSED-STRICT", okj)
        charts[kind] = None     # memory: free after certification

    # ---------------- R14: coverage ----------------
    okcov = (Fr(7, 4) < Fr(15, 8) < Fr(127, 64) < 2
             and Fr(5, 8) < Fr(31, 48) < Fr(2, 3) < Fr(11, 16)
             < Fr(17, 24) < 1
             and Fr(31, 48) == Fr(2, 3) - Fr(1, 48)
             and Fr(11, 16) == Fr(2, 3) + Fr(1, 48)
             and Fr(127, 64) == 2 - Fr(1, 64)
             and Fr(1, 8) ** 2 == Fr(1, 64)
             and Fr(1, 64) < Fr(1, 48))
    rep("R14 coverage arithmetic (t and u unions, window bounds, "
        "sqrt(1/64) = 1/8)", okcov)
    log("R14 union: u=1: F1 face;  t<15/8: C1;  t>=15/8: u<=5/8: C2;")
    log("    u in [5/8,17/24]: C4 (sg>=1/64) | W3a (t<=127/64) |")
    log("    W3b/W3c (|u-2/3|>=1/48) | W3d charts (window, sg>0);")
    log("    u>17/24: C3  ==> core > 0 on all of Omega")

    # ---------------- R15: anchor ----------------
    pv = sp.Rational(9, 16)          # (sg,u,t) = (1/2, 3/4, 1)
    Dv = sp.Rational(63, 640)
    okan = (sp.cancel(PSUB.subs(sg, sp.Rational(1, 2)) - pv) == 0)
    okan &= (sp.cancel(DSUB.subs(
        th, sp.Rational(3, 4) * (1 - sg ** 2 / 2)).subs(
        sg, sp.Rational(1, 2)) - Dv) == 0)
    Yv = sp.Symbol('Y')
    psi = Yv ** 5
    for k in range(1, 6):
        psi += (sp.Integer(2 * MSCALE) ** k
                * Etld[k].subs({p: pv, D: Dv, t: 1})) * Yv ** (5 - k)
    nreal = sp.Poly(psi, Yv).count_roots()
    danch = sum(Fr(cc) * Fr(9, 16) ** (k >> 16)
                * Fr(63, 640) ** ((k >> 8) & 255)
                for k, cc in DISC.items())    # t = 1: t-powers all 1
    rep("R15 anchor (sg,u,t)=(1/2,3/4,1) <=> (p,D)=(9/16,63/640): 5 "
        "distinct real roots, DISC > 0", okan and nreal == 5 and danch > 0)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("disc(psi) > 0 on Omega = (0,1)x(0,1]x(0,2) in (sg,u,t); anchor")
    print("has 5 distinct real roots; Omega connected => the A=4 replica")
    print("spectrum is REAL and SIMPLE for all 0<p<3/4,")
    print("0<D<=thetatilde4*Dbar, 0<s<pi.  With the committed s=pi lemma")
    print("(lemma_7_34_a4_realness_pi.py): all angles in (0,pi] on the")
    print("curve-restricted domain, which contains the whole A=4 Gray")
    print("region (cert5-on-gray certificate (B)).")
    return PASS


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
