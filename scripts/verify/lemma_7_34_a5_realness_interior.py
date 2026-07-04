#!/usr/bin/env python3
r"""
lemma_7_34_a5_realness_interior.py
============================================================================
BUG-009-D / Route B (A=5): REALNESS OF THE REPLICA SPECTRUM AT ALL
INTERIOR ANGLES 0 < s < pi, proven on the TWO-PIECE CURVE-RESTRICTED
Gray domain.

THEOREM (certified here, exact integer/rational arithmetic end to end).
  For A = 5, all sg in (0,1), all D with 0 < D <= curve*(sg) Dbar(sg),
  and ALL spectral angles s in (0, pi), the 5x5 pattern-quotient replica
  matrix Q(p, D, w = e^{is}) has 5 REAL, DISTINCT eigenvalues, where
      p = (4/5)(1 - sg^2)                  (all p in (0, 4/5)),
      Dbar(sg) = (4/5)(1-sg)^2/(1+sg^2),
      curve*(sg) = curveM(sg) := 1 - (3/4) sg^2      on (0, 7/10],
                   curveP(sg) := (9 - 5 sg^2)/10     on [7/10, 1)
  (curveP >= curveM on the overlap: curveP - curveM = (5sg^2-2)/20 > 0
  for sg >= 7/10 -- gated).  The claim is proven PER PIECE, each piece
  rationalized onto its own unit cube:
      piece M:  sg = (3/4) x,        x in (0,1)  ==> sg in (0, 3/4),
      piece P:  sg = (7 + 3x)/10,    x in [0,1)  ==> sg in [7/10, 1),
  and the union of the two pieces covers all sg in (0,1).  Together
  with the s = pi endpoint lemma (lemma_7_34_a5_realness_pi.py, same
  two-piece domain), realness holds for ALL angles s in (0, pi].

WHY TWO PIECES AND A CURVE (vs the A=3 full-strip interior lemma):
full-strip realness is FALSE at A=5 already at interior angles, and the
A=5 complexification island cuts DEEPER than at A=4 (bottom th ~ 0.727
at sg ~ 0.64, t = 2; at interior t = 15/8 the bottom is ~ 0.76).  This
script EXHIBITS the island at an interior angle exactly (gate R8b): at
(sg, th, t) = (16/25, 87/100, 15/8) -- strictly inside the open
superdomain strip (th < 1), strictly ABOVE curveM (87/100 > 433/625 =
curveM(16/25)), at an interior angle (t = 15/8 < 2) -- the discriminant
DISC below is < 0 (exact rational): two eigenvalues of Q are genuinely
complex there.  Also, no SINGLE quadratic curve 1 - a sg^2 clears the
island while threading the (A)/(B) corridor (the corridor pinches to
a in (0.599, 0.604) at sg -> 1 but the island needs a >= ~0.72): the
two-piece splice is necessary, not an artifact.

REDUCTION TO A POLYNOMIAL DISCRIMINANT (each step gated in-script):
  (1) M := c Q, c = p(1-p) G Gt, G = D^2 w + (D-4)(1-D), Gt = D^2 +
      (D-4)(1-D) w, is a POLYNOMIAL matrix (R1a); chiM = charpoly(M) =
      sum_k (-1)^k E_k X^{5-k} with E_k w-symmetric at the central power
      (-w)^k:  E_k = (-w)^k Etld_k(p,D,t)|_{t=1-(w+1/w)/2} -- an EXACT
      Laurent identity per k (R1c).  On |w| = 1 (t = 1 - cos s):
        psi(Y) := Y^5 + sum_k 8^k Etld_k(p,D,t) Y^{5-k}
      has INTEGER polynomial coefficients (R3; the A=5 psi-scale is
      8^k = (2*4)^k, replacing A=4's 18^k = (2*9)^k) and satisfies
      psi(Y) = (8/w)^5 chiM(wY/8): its roots are exactly
        8 p(1-p) B(D,t) * lambda_i(Q),
      B(D,t) = (5D-4)^2 + 2tD^2(1-D)(4-D) = G Gt / w  (R2).  B > 0 on
      the curve domain: 4 - 5D > 0 there via the exact per-piece
      identities (R2b)
        piece M:  64(16+9x^2) - u(64-27x^2)(4-3x)^2
                  = (1-u) 64(16+9x^2) + u 3x (512+144x-216x^2+81x^3),
        piece P:  200(149+42x+9x^2) - 9u(131-42x-9x^2)(1-x)^2
                  = (1-u) 200(149+42x+9x^2)
                    + u (28621+11136x-54x^2+216x^3+81x^4),
      with the u=1 cofactors root-free on [0,1] (Sturm) -- so the
      multiplier 8 p(1-p) B is strictly positive on the whole claimed
      domain, and realness/simplicity of Q's spectrum is EQUIVALENT to
      realness/simplicity of psi's roots.
  (2) DISC(p,D,t) := disc_Y(psi): ONE integer polynomial (34173 terms,
      deg (54,78,20) -- the same shape as A=4), assembled through the
      universal quintic discriminant (59 terms, isobaric weight 20,
      R4/R5) and validated twice BEFORE any box run: exact rational
      equality with sympy's univariate discriminant at 16 random points
      (R6), and an mpmath eigenvalue-based sign+ratio cross-check at
      >= 25 points, dps 160 (R7); plus a targeted dps-40 adversarial
      scan of the island-adjacent tight zone of piece M (R7b; numeric
      cross-check, not load-bearing).
  (3) t=2 face structure (R8): DISC(p,D,2) == 8^20 disc_X(f4) *
      f4(lam1M)^2, with chiM|_{w=-1} = f1 f4 the committed w=-1 split,
      f1 = X - lam1M, lam1M = (1/2) D(D-1)(5p-4)^2 (2D^2-5D+4) --
      rederived in-script (the s=pi lemma proves realness ON this face;
      here the identity is a structural gate only).
  (4) Rationalized CURVE domain PER PIECE (x,u,t) in [0,1]x[0,1]x[0,2]:
      piece M:  p = (4/5)(1 - (9/16)x^2),
                D = u (64-27x^2)(4-3x)^2 / (80 (16+9x^2)),
        target_M := 5^132 4^210 (16+9x^2)^78 DISC(subst) is an integer
        polynomial (R9a, exact Fraction gates), and by EXACT division
        (R9b/R9c)
          target_M = content u^14 x^48 t^7 (4-3x)^34 (16+9x^2)^2
                     (64-27x^2)^14 * coreM,
        every stripped factor > 0 on Omega_M := (0,1)x(0,1]x(0,2)
        (x=0 <=> p=4/5, u=0 <=> D=0, t=0 <=> s=0 are excluded limits);
        coreM: 153638 terms, degs (306, 64, 13).
      piece P:  p = (4/5)(1-sg^2), sg = (7+3x)/10,
                D = (9/250) u (131-42x-9x^2)(1-x)^2 / (149+42x+9x^2),
        target_P := 5^396 2^156 (149+42x+9x^2)^78 DISC(subst), and
          target_P = content u^14 t^7 (1-x)^34 (7+3x)^48
                     (149+42x+9x^2)^2 (131-42x-9x^2)^14 * coreP,
        every stripped factor > 0 on Omega_P := [0,1)x(0,1]x(0,2)
        (x=1 <=> sg=1 <=> p=0 is the excluded limit; x=0 <=> sg=7/10 is
        INTERIOR and included); coreP: 153982 terms, degs (306,64,13).
      Hence sign(disc psi) = sign(core) per piece and the theorem
      reduces to core > 0 on Omega_M resp. Omega_P.
  (5) CLOSED-BOX ZERO GEOGRAPHY (gated R9d; piece M has exactly the
      A=4 shape, piece P is zero-free):
      piece M:
        * the identically-zero edge  E := (x=0, u=1, t in [0,2])
          (R9d-i; sg=0 <=> p = pmax, the same degeneration as A=4), and
        * the boundary pinch  P1 := (x,u,t) = (0, 5/8, 2)
          (R9d-ii, the pinch-location identity: coreM(0,u,2) ==
          (1-u)^4 (8u-5)^4 E54(u) with E54 > 0 on [0,1] -- exact
          division + Sturm; u0 = 5/8 is the A=5 analogue of A=4's
          u0 = 2/3: both map to D = 1/2 at sg = 0).
      Both are OUTSIDE Omega_M, but P1 forces a blow-up exactly as at
      A=4 (and A=3's Z3).
      piece P: coreP has NO zeros on the closed box: the x=0 and x=1
      t=2 edge slices are root-free resp. a positive constant (R9d-P),
      and the ONE closed-strict root-box certificate (R10P) proves
      coreP > 0 on ALL of [0,1]x[0,1]x[0,2] -- no decomposition, no
      blow-ups.
  (6) coreM > 0 on Omega_M -- Bernstein decomposition (all certificates
      exact integer arithmetic; "no-negatives" = all scaled Bernstein
      coefficients >= 0, not all zero => f >= 0 on the closed box AND
      f > 0 on the OPEN box; "closed-strict" = min > 0 AND no exact-
      zero coefficients => f > 0 on the CLOSED box; closed-strict
      pieces may be bisected, that gluing is sound):
        C1   [0,1] x [0,1]       x [0,15/8]        no-negatives
        F1   the u=1 FACE, (x,t) in [0,1]x[0,2]    no-negatives
        C2   [0,1] x [0,7/12]    x [7/4,2]         no-negatives
        C3a  [0,5/8] x [2/3,1]   x [7/4,2]         no-negatives
        C3b  [5/8,1] x [2/3,1]   x [7/4,2]         CLOSED-STRICT
             (the island-adjacent split of A=4's C3: the A=5 island
             clears the u=1 face by only ~3.3% in u near x ~ 0.83, so
             the x >= 5/8 part is certified closed-strict away from
             the zero edge E, and the x <= 5/8 part -- which contains
             E -- gets the open-box certificate)
        C4   [1/64,1] x [7/12,2/3] x [7/4,2]       CLOSED-STRICT
        W3a  [0,1/64] x [7/12,2/3] x [7/4,127/64]  CLOSED-STRICT
        W3b  [0,1/64] x [7/12,29/48] x [127/64,2]  CLOSED-STRICT
        W3c  [0,1/64] x [31/48,2/3] x [127/64,2]   CLOSED-STRICT
        W3d  the blow-up window [0,1/64]x[29/48,31/48]x[127/64,2]
             around the boundary pinch P1 = (0, 5/8, 2).
      (In the recorded run every certificate passes AT ITS ROOT BOX --
      zero bisections; the difficulty was the domain geography, exactly
      as at A=3/A=4.)
  (7) W3d, quasi-homogeneous order-4 blow-up (R13).  Local coordinates
      s = x, m = u - 5/8, v = 2 - t (window |m| <= 1/48, s <= 1/64,
      v <= 1/64); core_loc := 8^64 * coreM(s, 5/8+m, 2-v) / g in
      Z[s,m,v] (g = integer content; spot-gated against coreM).
      VERIFIED: minimal (1,1,2)-weighted degree of core_loc is EXACTLY
      4 (R13b), and the weighted-leading form is
        L4 = c002 v^2 + c021 (m - (15/16)s)^2 v
             + (lam4/2) [ (16m - 15s)^4 + (25s)^4 ],
      c002, c021, lam4 > 0 -- EXACT identities (R13c; the middle
      coefficient is a perfect square c111^2 = 4 c021 c201 with
      c111/c021 = -15/8, and the v=0 quartic satisfies the integer
      identity 2(32768 m^4 - 122880 m^3 s + 172800 m^2 s^2
      - 108000 m s^3 + 220625 s^4) == (16m-15s)^4 + (25s)^4).
      L4 > 0 for v >= 0, (s,m,v) != 0.  Charts (monomial maps,
      bijective, exact inverse reconstruction gates R13d,
      Schwartz-Zippel identity gates mod 2^61-1 R13e):
        HS :  core_loc(s, s mh, s^2 vh)      = s^4  HS(s, mh, vh)
        HM+/-: core_loc(mu sh, +-mu, mu^2 vh) = mu^4 HM+-(mu, sh, vh)
        HU :  core_loc(nu sh, nu mh, nu^2)   = nu^4 HU(nu, sh, mh)
      certified CLOSED-STRICT on  [0,1/64]x[-1,1]x[0,1],
      [0,1/48]x[0,1]x[0,1] (both signs), [0,1/8]x[0,1]x[-1,1] (R13f).
      SOUNDNESS: for (s,m,v) != 0 in the window let rho = max(s, |m|,
      sqrt(v)); the corresponding chart point lies in its CLOSED box
      and the identity gives core = rho^4 * H(...) > 0.  Only P1
      itself (t = 2, outside Omega_M) is excluded.
  (8) COVERAGE (R14, exact rational comparisons).
      piece M:  u = 1: face F1 (open-face > 0 for x in (0,1), t in
      (0,2)).  u in (0,1), t in (0,15/8): C1 (open box).  u in (0,1),
      t in [15/8,2): u-union (0,7/12) (C2, open, t-range (7/4,2)
      contains [15/8,2)) u [7/12,2/3] u (2/3,1), where [7/12,2/3] is
      covered for x >= 1/64 by C4 (closed) and for x <= 1/64 by W3a
      (t <= 127/64), W3b/W3c (u outside the [29/48,31/48] window), and
      W3d-charts (the window, (s,m,v) != 0); and (2/3,1) is covered by
      C3a (open, x < 5/8) u C3b (closed, x >= 5/8).  t-overlap: 7/4 <
      15/8 < 127/64; window arithmetic 29/48 = 5/8 - 1/48, 31/48 =
      5/8 + 1/48, 7/12 = 5/8 - 1/24, 2/3 = 5/8 + 1/24, 127/64 =
      2 - 1/64, sqrt(1/64) = 1/8.   =>  coreM > 0 on all of Omega_M.
      piece P: the single R10P closed-box certificate covers Omega_P.
      UNION over sg: piece M covers sg in (0, 3/4), piece P covers
      sg in [7/10, 1); (0,3/4) u [7/10,1) = (0,1), and on the overlap
      curveP >= curveM (gated), so curve* is covered everywhere.
  (9) ANCHOR + CONNECTEDNESS (R15, per piece): piece M at (x,u,t) =
      (2/3, 3/4, 1) <=> (p,D) = (3/5, 39/400), s = pi/2: psi has 5
      DISTINCT REAL roots (exact Sturm) and DISC > 0 (exact); piece P
      at (x,u,t) = (1/3, 1/2, 1) <=> (p,D) = (36/125, 29/5125).  Each
      Omega is path-connected; disc psi > 0 on it; along any path the
      real-root count of a monic real quintic can change only through
      a collision (disc = 0), which never happens on Omega.  Hence 5
      distinct real roots everywhere on both Omegas, i.e. the replica
      spectrum is REAL and SIMPLE.                                 QED

COMBINED SPECTRAL STATE for A=5 after this certificate: realness holds
on the two-piece curve-restricted domain for all angles s in (0, pi]
(this script + the s=pi lemma lemma_7_34_a5_realness_pi.py).  Whether
the two-piece domain contains the whole A=5 Gray region (D_c < curve*
Dbar) is the (B)-side of the separate cert5-on-gray deliverable (the
scout's piecewise Sturm inclusion passed for both pieces; that lemma is
NOT this script).

ENGINE (packed-key integer tensors, key = a_x<<20 | b_u<<10 | c_t):
scaled Bernstein coefficients b^_beta = sum_alpha C(d-alpha, beta-alpha)
a_alpha = C(d,beta) b_beta (positive rescale); affine box maps
v -> lo + (hi-lo)v are exact-integer (global positive scale q^deg).
Two implementations of the Bernstein axis map -- the committed A=3/A=4
dict form and a list-based Taylor-shift form (pure additions) -- are
cross-checked coefficient-exactly in R0d; the fast form runs the heavy
boxes.  Self-tests R0 pin the map on knowable examples (including the
lo != 0 affine map -- the historical h**b bug class -- and exact-zero
accounting).  Semantics identical to the committed engine.

HONEST NOTES / SCOPE:
 *  The claim is the OPEN angle interval (0, pi): s = pi realness is
    the separate s=pi lemma (the t=2 face of coreM genuinely vanishes
    at the excluded pinch P1 and along the excluded edge E; both are
    OUTSIDE Omega_M).
 *  u = 1 (D = curve* Dbar) IS included; u = 0 (D = 0), t = 0 (s = 0),
    sg in {0,1} (p in {4/5, 0}) are excluded limits where the quotient
    construction itself degenerates (eta(w=1) = 0 forces rank 1: the
    t^7 factor; D = 0 likewise: u^14; p = 4/5 rank 1: x^48 in piece M;
    p = 0: the channel is the identity: (1-x)^34 in piece P).
 *  NO full-strip claim: R8b exhibits DISC < 0 above the curve at an
    interior angle (t = 15/8).  The A=3 full-strip statement fails to
    replicate at A=5 (worse than at A=4: the island clears curveM by
    only ~0.024 in th); the two-piece curve-restricted statement is
    the correct A=5 analogue.
 *  Simplicity (distinctness) is part of the certificate: disc > 0.
 *  R7/R7b are floating-point cross-checks of an exact object (belt
    and suspenders on top of the exact R6); every load-bearing
    certificate is exact integer arithmetic.
 *  The two pieces are certified SEQUENTIALLY (piece M fully freed
    before piece P is substituted): peak RSS ~2 GB.
 *  Runtime ~20 min single-threaded (recorded run: 1215 s, 53/53 gates
    PASS; EVERY Bernstein certificate -- including all four blow-up
    charts and the piece P whole-box run -- passed AT ITS ROOT BOX,
    zero bisections; the largest single certificate is the HU chart of
    the P1 blow-up, 7842315 Bernstein coefficients).

Deps: sympy, mpmath only; fully self-contained (the A=5 replica
machinery is derived in-script; nothing is imported from other lemma
scripts).  Optional dev cache: env A5RI_CACHE (directory for the DISC
pickle).  Python: /Users/para/.venvs/rnr/bin/python.
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

mp.mp.dps = 60
T0 = time.time()
PASS = True
CACHE = os.environ.get("A5RI_CACHE", "")   # optional dev cache dir
MSCALE = 4                                 # root scale: 2*MSCALE = 8
AVAL = 5

p, D, w, t = sp.symbols('p D w t')
X = sp.Symbol('X')
U = sp.Symbol('u')


def rep(name, ok):
    global PASS
    PASS = PASS and bool(ok)
    print(f"  [{time.time()-T0:7.1f}s] {name:<64} {'PASS' if ok else 'FAIL'}",
          flush=True)
    return bool(ok)


def log(msg):
    print(f"  [{time.time()-T0:7.1f}s] {msg}", flush=True)


# ===========================================================================
# A=5 replica machinery (self-contained; the a5 analogue of the committed
# A=4 pipeline imports)
# ===========================================================================
def _pat(tr):
    x, u_, up = tr
    if x == u_ == up: return 0
    if x == u_ and u_ != up: return 1
    if x == up and u_ != up: return 2
    if u_ == up and x != u_: return 3
    return 4


STATES = list(itertools.product(range(AVAL), repeat=3))
REPS = [None] * 5
for _k, _tr in enumerate(STATES):
    if REPS[_pat(_tr)] is None:
        REPS[_pat(_tr)] = _k

ETH = D / ((AVAL - 1) * (1 - D))


def _etaf(W):
    E = ETH * W
    return ((AVAL - 1) * D * E + D - (AVAL - 1) * E) \
        / ((AVAL - 1) * D * E + D - (AVAL - 1))


def build_Q():
    eta, etb = _etaf(w), _etaf(1 / w)
    T = [[(1 - p) if i == j else p / (AVAL - 1) for j in range(AVAL)]
         for i in range(AVAL)]
    Q = sp.zeros(5, 5)
    for a_ in range(5):
        x, xp, xq = STATES[REPS[a_]]
        for (y, yp, yq) in STATES:
            term = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
            if yp != y: term *= eta
            if yq != y: term *= etb
            Q[a_, _pat((y, yp, yq))] += term
    return Q


def laurent_sym_to_t(expr, K=60):
    """w-symmetric Laurent polynomial -> polynomial in t via Chebyshev in
    u = w + 1/w = 2 - 2t on |w| = 1.  Returns None if not symmetric."""
    ex = sp.expand(expr)
    npoly = sp.Poly(sp.expand(ex * w ** K), w)
    co = {}
    for (md,), c in npoly.terms():
        co[md - K] = co.get(md - K, 0) + c
    ks = sorted(co.keys())
    maxk = max(abs(k) for k in ks) if ks else 0
    C = {0: sp.Integer(2), 1: U}
    for k2 in range(2, maxk + 1):
        C[k2] = sp.expand(U * C[k2 - 1] - C[k2 - 2])
    out = sp.Integer(0)
    done = set()
    for k2 in ks:
        if k2 in done: continue
        if k2 == 0:
            out += co[0]; done.add(0)
        else:
            kk = abs(k2)
            if sp.simplify(co.get(kk, 0) - co.get(-kk, 0)) != 0:
                return None
            out += co.get(kk, 0) * C[kk]
            done.add(kk); done.add(-kk)
    return sp.expand(out.subs(U, 2 - 2 * t))


def sym_to_t_mid(expr, K=60):
    """Laurent in w -> (t-polynomial, central power m of (-w)); None if the
    central power is half-integer or the Laurent part is asymmetric."""
    ex = sp.expand(expr)
    poly = sp.Poly(sp.expand(ex * w ** K), w)
    degs = [md - K for (md,), cc in poly.terms()]
    mid = sp.Rational(min(degs) + max(degs), 2)
    if int(mid) != mid:
        return None, None
    r = laurent_sym_to_t(sp.cancel(ex / (-w) ** int(mid)), K)
    return r, (int(mid) if r is not None else None)


# ===========================================================================
# packed-key integer polynomial engine (committed A=4 form)
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
    scaled by q^dax > 0 (q = common denominator) -- sign-preserving."""
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
    committed A=3/A=4 dict form."""
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
# R8 helper: w = -1 split and the t=2 face identity (A=5)
# ===========================================================================
def t2_face_gate(DISC):
    Xs = sp.Symbol('X')
    K5 = 2 * D ** 2 - 5 * D + 4
    Ew = -D / ((AVAL - 1) * (1 - D))          # w = -1
    etaw = sp.cancel(((AVAL - 1) * D * Ew + D - (AVAL - 1) * Ew)
                     / ((AVAL - 1) * D * Ew + D - (AVAL - 1)))
    Tm = [[(1 - p) if i == j else p / (AVAL - 1) for j in range(AVAL)]
          for i in range(AVAL)]
    Qm = sp.zeros(5, 5)
    for a_ in range(5):
        x0, xp_, xq = STATES[REPS[a_]]
        for (y, yp, yq) in STATES:
            term = Tm[xp_][yp] * Tm[xq][yq] / Tm[x0][y]
            if yp != y:
                term *= etaw
            if yq != y:
                term *= etaw
            Qm[a_, _pat((y, yp, yq))] += term
    # c|_{w=-1} = p(p-1)(2D^2-5D+4)^2  (G Gt|_{w=-1} = -(2D^2-5D+4)^2)
    c_ = p * (p - 1) * K5 ** 2
    Mm = sp.zeros(5, 5)
    for i in range(5):
        for j in range(5):
            n2, d2 = sp.fraction(sp.cancel(c_ * Qm[i, j]))
            Mm[i, j] = sp.expand(n2 / d2)
    chiM = sp.expand(Mm.charpoly(Xs).as_expr())
    lam1M = sp.expand(sp.Rational(1, 2) * D * (D - 1) * (5 * p - 4) ** 2
                      * K5)
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
# R9 helpers: per-piece curve substitution (x-coordinates) and strips
# ===========================================================================
def polypow_list(base, n):
    """list of dense coefficient lists base^m, m = 0..n (integer conv)."""
    out = [[1]]
    for m in range(1, n + 1):
        prev = out[-1]
        cur = [0] * (len(prev) + len(base) - 1)
        for i, v in enumerate(prev):
            if v:
                for j, e in enumerate(base):
                    if e:
                        cur[i + j] += v * e
        out.append(cur)
    return out


def conv(Aa, Bb):
    out = [0] * (len(Aa) + len(Bb) - 1)
    for i, x in enumerate(Aa):
        if x:
            for j, y in enumerate(Bb):
                if y:
                    out[i + j] += x * y
    return out


def substitute_piece(DISC, pnum, Dnum_lists, mulc_fn):
    """target(x,u,t) = clear * DISC(p(x), D(x,u), t) as a packed dict.
    pnum[a]: coefficient list of the p-numerator^a (in x); Dnum_lists[b]:
    coefficient list for the D^b numerator INCLUDING the denominator-
    clearing complement; mulc_fn(a,b): the integer scale."""
    groups = {}
    for k, cf in DISC.items():
        a, b, cc = k >> 16, (k >> 8) & 255, k & 255
        groups.setdefault((a, b), []).append((cc, cf))
    target = {}
    g = target.get
    for (a, b), lst in groups.items():
        Wl = conv(pnum[a], Dnum_lists[b])
        mulc = mulc_fn(a, b)
        for (cc, cf) in lst:
            fmul = cf * mulc
            kb = (b << 10) | cc
            for e, y in enumerate(Wl):
                if y:
                    k2 = (e << 20) | kb
                    nv = g(k2, 0) + fmul * y
                    if nv:
                        target[k2] = nv
                    elif k2 in target:
                        del target[k2]
                    g = target.get
    return target


def strip_target(target, factors):
    """iterated exact division by monomials/x-polynomials; the polynomial
    divisor may have any nonzero integer leading coefficient (each long-
    division step checks exact integer divisibility)."""
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

    def div_poly_x(T, dcoef):
        dd = len(dcoef) - 1
        lead = dcoef[-1]
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
            if which == 'x':
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
    for name, spec in factors:
        while True:
            r = div_mono(core, spec) if isinstance(spec, str) \
                else div_poly_x(core, spec)
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
    return core, strips, max(g, 1)


# ===========================================================================
# R13 helpers: the P1 blow-up (piece M; u0 = 5/8, denominator 8)
# ===========================================================================
def build_local(core, CU):
    """core_loc = 8^CU * coreM(s, 5/8 + m, 2 - v) as {(i,j,k): int}; also
    returns the stripped integer content g (core_loc*g == 8^CU core(...))."""
    loc = {}
    P8 = [8 ** e for e in range(CU + 1)]
    for k0, v in core.items():
        a, b, cc = k0 >> 20, (k0 >> 10) & 1023, k0 & 1023
        for j in range(b + 1):
            cb = comb(b, j) * (5 ** (b - j)) * P8[CU - (b - j)]
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


def dict_eval_fr(T, xv, uv, tv):
    return sum(Fr(cf) * xv ** (k >> 20) * uv ** ((k >> 10) & 1023)
               * tv ** (k & 1023) for k, cf in T.items())


def disc_eval_fr(DISC, pv, Dv, tv):
    return sum(Fr(cf) * pv ** (k >> 16) * Dv ** ((k >> 8) & 255)
               * tv ** (k & 255) for k, cf in DISC.items())


# ===========================================================================
# main
# ===========================================================================
def main():
    print("=" * 78)
    print("Interior-angle realness, A = 5: disc(psi) > 0 on the TWO-PIECE")
    print("rationalized curve domain for all t in (0,2) -- exact Bernstein")
    print("=" * 78)
    engine_selftest()

    # ---------------- R1 ----------------
    log("R1 building Q, M = cQ, charpoly, Etld_k ...")
    Q = build_Q()
    G = D ** 2 * w + (D - 4) * (1 - D)
    Gt = D ** 2 + (D - 4) * (1 - D) * w
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
    Bpoly = (5 * D - 4) ** 2 + 2 * t * D ** 2 * (1 - D) * (4 - D)
    chat, mc = sym_to_t_mid(c)
    ok2 = (mc == 1 and sp.expand(chat + p * (1 - p) * Bpoly) == 0)
    ok2 &= (sp.expand(Gt - w * G.subs(w, 1 / w)) == 0)
    rep("R2a c = w p(1-p) B(D,t), B = (5D-4)^2 + 2tD^2(1-D)(4-D); "
        "Gt = w G(1/w)", ok2)
    # B > 0 on both pieces: 4 - 5D > 0 via exact positive-combination
    # identities in the piece coordinates
    xx, uu = sp.symbols('x u')
    DSUB_M = uu * (64 - 27 * xx ** 2) * (4 - 3 * xx) ** 2 \
        / (80 * (16 + 9 * xx ** 2))
    DSUB_P = sp.Rational(9, 250) * uu * (131 - 42 * xx - 9 * xx ** 2) \
        * (1 - xx) ** 2 / (149 + 42 * xx + 9 * xx ** 2)
    # piece coordinates match the global curve parametrization
    sgM = sp.Rational(3, 4) * xx
    sgP = (7 + 3 * xx) / sp.Integer(10)
    DBARg = sp.Rational(4, 5) * (1 - sp.Symbol('sigma')) ** 2 \
        / (1 + sp.Symbol('sigma') ** 2)
    sgv = sp.Symbol('sigma')
    ok2b = (sp.cancel(DSUB_M - (uu * (1 - sp.Rational(3, 4) * sgv ** 2)
                                * DBARg).subs(sgv, sgM)) == 0)
    ok2b &= (sp.cancel(DSUB_P - (uu * (9 - 5 * sgv ** 2) / 10
                                 * DBARg).subs(sgv, sgP)) == 0)
    rep("R2b piece D-parametrizations == u * curve * Dbar (exact)", ok2b)
    # piece M: (4-5D) 16(16+9x^2) == (1-u) 64(16+9x^2) + u 3x cubM
    cubM = 512 + 144 * xx - 216 * xx ** 2 + 81 * xx ** 3
    okM = (sp.cancel((4 - 5 * DSUB_M) * 16 * (16 + 9 * xx ** 2)
                     - ((1 - uu) * 64 * (16 + 9 * xx ** 2)
                        + uu * 3 * xx * cubM)) == 0)
    okM &= (sp.Poly(cubM, xx).count_roots(0, 1) == 0
            and cubM.subs(xx, 0) == 512)
    quaP = 28621 + 11136 * xx - 54 * xx ** 2 + 216 * xx ** 3 + 81 * xx ** 4
    okP = (sp.cancel((4 - 5 * DSUB_P) * 50 * (149 + 42 * xx + 9 * xx ** 2)
                     - ((1 - uu) * 200 * (149 + 42 * xx + 9 * xx ** 2)
                        + uu * quaP)) == 0)
    okP &= (sp.Poly(quaP, xx).count_roots(0, 1) == 0
            and quaP.subs(xx, 0) == 28621)
    rep("R2c 4-5D > 0 on both pieces (positive-combination identities + "
        "Sturm)", okM and okP)

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
    rep("R3 Ehat_k = 8^k Etld_k have integer coefficients (A=5 psi-scale)",
        ok3)

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
    cachef = os.path.join(CACHE, "disc_a5_scale8.pkl") if CACHE else None
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
        pv = sp.Rational(random.randint(1, 63), 80)   # p in (0, 4/5)
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
        promoted BEFORE any arithmetic -- the committed A=4 repair)."""
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
    for _ in range(60):
        sgv2 = random.uniform(0.02, 0.98)
        uv = random.uniform(0.05, 1.0)
        tv = random.uniform(0.001, 1.999)
        pv = 0.8 * (1 - sgv2 ** 2)
        if sgv2 <= 0.75:
            cv = 1 - 0.75 * sgv2 ** 2
        else:
            cv = (9 - 5 * sgv2 ** 2) / 10
        Dv = uv * cv * 0.8 * (1 - sgv2) ** 2 / (1 + sgv2 ** 2)
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
    rep("R8 DISC(p,D,2) == 8^20 disc(f4) f4(lam1M)^2", t2_face_gate(DISC))

    # ---------------- R8b: the interior-angle island ----------------
    sgvI = Fr(16, 25)
    thvI = Fr(87, 100)
    tvI = Fr(15, 8)
    pvI = Fr(4, 5) * (1 - sgvI ** 2)
    DbvI = Fr(4, 5) * (1 - sgvI) ** 2 / (1 + sgvI ** 2)
    DvI = thvI * DbvI
    curveI = 1 - Fr(3, 4) * sgvI ** 2
    dI = disc_eval_fr(DISC, pvI, DvI, tvI)
    rep("R8b island at INTERIOR angle: DISC < 0 at (16/25, 87/100, 15/8), "
        "th > curveM", dI < 0 and thvI > curveI and thvI < 1 and tvI < 2)

    # =================================================================
    # PIECE M
    # =================================================================
    log("R9M piece M curve substitution -> target(x,u,t), sg = (3/4)x")
    pnumM = polypow_list([16, 0, -9], dP)          # (16-9x^2)^a
    p64 = polypow_list([64, 0, -27], dD)
    p43 = polypow_list([4, -3], 2 * dD)
    p169 = polypow_list([16, 0, 9], dD)
    DnumM = [conv(p64[b], conv(p43[2 * b], p169[dD - b]))
             for b in range(dD + 1)]
    del p64, p43, p169
    target = substitute_piece(
        DISC, pnumM, DnumM,
        lambda a, b: 5 ** (dP - a) * 4 ** (dP - a) * 5 ** (dD - b)
        * 4 ** (2 * (dD - b)))
    del pnumM, DnumM
    random.seed(31)
    ok9a = True
    for _ in range(6):
        xv = Fr(random.randint(1, 15), 16)
        uv = Fr(random.randint(1, 16), 16)
        tv = Fr(random.randint(0, 32), 16)
        sgv2 = Fr(3, 4) * xv
        pv = Fr(4, 5) * (1 - sgv2 ** 2)
        Dv = Fr(1, 5) * uv * (4 - 3 * sgv2 ** 2) * (1 - sgv2) ** 2 \
            / (1 + sgv2 ** 2)
        dv = disc_eval_fr(DISC, pv, Dv, tv)
        tv2 = dict_eval_fr(target, xv, uv, tv)
        ok9a &= (tv2 == Fr(5) ** (dP + dD) * Fr(4) ** (dP + 2 * dD)
                 * (16 + 9 * xv ** 2) ** dD * dv)
    rep("R9Ma target_M == 5^132 4^210 (16+9x^2)^78 DISC(subst): 6 exact "
        "points", ok9a)
    core, strips, content = strip_target(
        target, [("u", 'u'), ("x", 'x'), ("t", 't'),
                 ("4-3x", [4, -3]), ("4+3x", [4, 3]),
                 ("16+9x^2", [16, 0, 9]), ("64-27x^2", [64, 0, -27])])
    log(f"R9M strips: {strips}; content bits {content.bit_length()}")
    ok9b = (strips.get('u', 0) == 14 and strips.get('x', 0) == 48
            and strips.get('t', 0) == 7 and strips.get('4-3x', 0) == 34
            and strips.get('16+9x^2', 0) == 2
            and strips.get('64-27x^2', 0) == 14
            and strips.get('4+3x', 0) == 0 and content > 0)
    rep("R9Mb target_M == content u^14 x^48 t^7 (4-3x)^34 (16+9x^2)^2 "
        "(64-27x^2)^14 coreM", ok9b)
    CS = max(k >> 20 for k in core)
    CU = max((k >> 10) & 1023 for k in core)
    CT = max(k & 1023 for k in core)
    degs = [CS, CU, CT]
    log(f"R9M coreM: {len(core)} terms, degs ({CS},{CU},{CT})")
    rep("R9Mc coreM: 153638 terms, degs (306,64,13)",
        len(core) == 153638 and (CS, CU, CT) == (306, 64, 13))
    random.seed(41)
    ok9c = True
    for _ in range(4):
        xv = Fr(random.randint(1, 15), 16)
        uv = Fr(random.randint(1, 16), 16)
        tv = Fr(random.randint(1, 31), 16)
        cv = dict_eval_fr(core, xv, uv, tv)
        tv2 = dict_eval_fr(target, xv, uv, tv)
        recon = (content * uv ** 14 * xv ** 48 * tv ** 7
                 * (4 - 3 * xv) ** 34 * (16 + 9 * xv ** 2) ** 2
                 * (64 - 27 * xv ** 2) ** 14 * cv)
        ok9c &= (recon == tv2)
    rep("R9Md strip reconstruction identity at 4 exact points", ok9c)
    del target

    # ---------------- R9d: piece M zero geography ----------------
    # (i) the edge (x=0, u=1): coreM identically zero along it
    edgeE = {}
    for k, v in core.items():
        if (k >> 20) == 0:
            edgeE[k & 1023] = edgeE.get(k & 1023, 0) + v
    ok9d1 = all(v == 0 for v in edgeE.values())
    rep("R9d-i edge E = (x=0, u=1, t): coreM identically zero", ok9d1)
    # (ii) pinch location: coreM(0,u,2) == (1-u)^4 (8u-5)^4 E54(u), E54 > 0
    uY = sp.Symbol('u')
    E2 = {}
    for k, v in core.items():
        if (k >> 20) == 0:
            a = (k >> 10) & 1023
            E2[a] = E2.get(a, 0) + v * (2 ** (k & 1023))
    Ep = sp.Poly(sum(sp.Integer(v) * uY ** a for a, v in E2.items()), uY)
    q54, r54 = sp.div(Ep, sp.Poly((1 - uY) ** 4 * (8 * uY - 5) ** 4, uY))
    ok9d2 = (r54 == 0)
    if ok9d2:
        Q54 = sp.Poly(q54, uY)
        ok9d2 &= (Q54.count_roots(0, 1) == 0
                  and Q54.eval(sp.Rational(1, 2)) > 0)
    rep("R9d-ii coreM(0,u,2) == (1-u)^4 (8u-5)^4 E54(u), E54 > 0 on [0,1]",
        ok9d2)

    # ---------------- R7b: adversarial numeric scan (cross-check) -------
    # the A=5 island-adjacent tight zone: x in [0.70,0.95] (sg in
    # [0.525,0.7125]), u in [0.9,1], t in [1.7,2] -- the island clears the
    # u=1 face by only ~3.3% in u near x ~ 0.83
    log("R7b dps-40 adversarial scan: x in [0.70,0.95], u in [0.9,1], "
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
            for ix in range(26):
                xv = mp.mpf('0.70') + mp.mpf('0.01') * ix
                acc = mp.mpf(0)
                for c2 in reversed(cs):
                    acc = acc * xv + c2
                cnt7b += 1
                if acc <= 0:
                    neg7b += 1
    del SL
    mp.mp.dps = dps_save
    rep(f"R7b coreM > 0 at all {cnt7b} tight-margin points (numeric "
        f"cross-check)", neg7b == 0 and cnt7b == 1092)

    # ---------------- R10: C1 + the u=1 face ----------------
    ok10 = no_negatives(core, degs, [(0, 1), (0, 1), (0, Fr(15, 8))],
                        "R10 C1 [0,1]^2x[0,15/8]")
    rep("R10 C1 no negatives (coreM > 0 on the open C1 box)", ok10)
    fu1 = face_u1(core)
    ok10f = no_negatives(fu1, [CS, 0, CT], [(0, 1), (0, 1), (0, 2)],
                         "R10f u=1 face, t in [0,2]", axes=(2, 0))
    rep("R10f u=1 face, ALL t in [0,2]: no negatives (open-face > 0)",
        ok10f)
    del fu1

    # ---------------- R11: C2/C3a/C3b/C4 ----------------
    ok11a = no_negatives(core, degs,
                         [(0, 1), (0, Fr(7, 12)), (Fr(7, 4), 2)],
                         "R11a C2 [0,1]x[0,7/12]x[7/4,2]")
    rep("R11a C2 no negatives", ok11a)
    ok11b = no_negatives(core, degs,
                         [(0, Fr(5, 8)), (Fr(2, 3), 1), (Fr(7, 4), 2)],
                         "R11b C3a [0,5/8]x[2/3,1]x[7/4,2]")
    rep("R11b C3a no negatives (contains the zero edge E)", ok11b)
    ok11c = closed_strict(core, degs,
                          [(Fr(5, 8), 1), (Fr(2, 3), 1), (Fr(7, 4), 2)],
                          "R11c C3b [5/8,1]x[2/3,1]x[7/4,2]")
    rep("R11c C3b CLOSED-STRICT (island-adjacent)", ok11c)
    ok11d = closed_strict(core, degs,
                          [(Fr(1, 64), 1), (Fr(7, 12), Fr(2, 3)),
                           (Fr(7, 4), 2)],
                          "R11d C4 [1/64,1]x[7/12,2/3]x[7/4,2]")
    rep("R11d C4 CLOSED-STRICT", ok11d)

    # ---------------- R12: wedge boxes ----------------
    ok12a = closed_strict(core, degs,
                          [(0, Fr(1, 64)), (Fr(7, 12), Fr(2, 3)),
                           (Fr(7, 4), Fr(127, 64))], "R12a W3a")
    rep("R12a W3a [0,1/64]x[7/12,2/3]x[7/4,127/64] CLOSED-STRICT", ok12a)
    ok12b = closed_strict(core, degs,
                          [(0, Fr(1, 64)), (Fr(7, 12), Fr(29, 48)),
                           (Fr(127, 64), 2)], "R12b W3b")
    rep("R12b W3b [0,1/64]x[7/12,29/48]x[127/64,2] CLOSED-STRICT", ok12b)
    ok12c = closed_strict(core, degs,
                          [(0, Fr(1, 64)), (Fr(31, 48), Fr(2, 3)),
                           (Fr(127, 64), 2)], "R12c W3c")
    rep("R12c W3c [0,1/64]x[31/48,2/3]x[127/64,2] CLOSED-STRICT", ok12c)

    # ---------------- R13: the P1 blow-up ----------------
    log("R13 local dict at P1 = (0, 5/8, 2) ...")
    loc, gloc = build_local(core, CU)
    # spot-gate the local dict against coreM at 2 exact rational points
    random.seed(97)
    okloc = True
    for _ in range(2):
        sv_ = Fr(random.randint(1, 9), 64)
        mv_ = Fr(random.randint(-9, 9), 100)
        vv_ = Fr(random.randint(1, 9), 64)
        lv = sum(Fr(cf) * sv_ ** i * mv_ ** j * vv_ ** kk
                 for (i, j, kk), cf in loc.items())
        cv = sum(Fr(cf) * sv_ ** (k >> 20)
                 * (Fr(5, 8) + mv_) ** ((k >> 10) & 1023)
                 * (2 - vv_) ** (k & 1023) for k, cf in core.items())
        okloc &= (Fr(gloc) * lv == Fr(8) ** CU * cv)
    rep("R13a g*core_loc == 8^64 coreM(s, 5/8+m, 2-v): 2 exact points",
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
    # v=0 quartic must be lam4*(32768 m^4 - 122880 m^3 s + 172800 m^2 s^2
    #                           - 108000 m s^3 + 220625 s^4)
    okq = (c040 % 32768 == 0)
    lam4 = c040 // 32768 if okq else 0
    okq = (okq and lam4 > 0
           and c130 == -122880 * lam4 and c220 == 172800 * lam4
           and c310 == -108000 * lam4 and c400 == 220625 * lam4)
    # integer identity: 2(32768m^4 - 122880m^3s + 172800m^2s^2 - 108000ms^3
    #                     + 220625s^4) == (16m-15s)^4 + (25s)^4
    sQ, mQ = sp.symbols('sQ mQ')
    okq &= (sp.expand(2 * (32768 * mQ ** 4 - 122880 * mQ ** 3 * sQ
                           + 172800 * mQ ** 2 * sQ ** 2
                           - 108000 * mQ * sQ ** 3 + 220625 * sQ ** 4)
                      - ((16 * mQ - 15 * sQ) ** 4 + (25 * sQ) ** 4)) == 0)
    okL4 = (c002 > 0 and c021 > 0
            and c111 ** 2 == 4 * c021 * c201
            and 8 * c111 == -15 * c021
            and okq)
    rep("R13c L4 = c002 v^2 + c021 (m-(15/16)s)^2 v + (lam4/2)[(16m-15s)^4"
        " + (25s)^4]", okL4)
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
    # polynomials have total degree <= ~900, so a nonzero difference
    # survives 4 uniform points with probability <= (900/q)^4 ~ 1e-70.
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

    # =================================================================
    # PIECE P
    # =================================================================
    log("R9P piece P curve substitution -> target(x,u,t), sg = (7+3x)/10")
    p1x17 = polypow_list(conv([1, -1], [17, 3]), dP)   # ((1-x)(17+3x))^a
    p131 = polypow_list([131, -42, -9], dD)
    p1x = polypow_list([1, -1], 2 * dD)
    p149 = polypow_list([149, 42, 9], dD)
    DnumP = [conv(p131[b], conv(p1x[2 * b], p149[dD - b]))
             for b in range(dD + 1)]
    del p131, p1x, p149
    targetP = substitute_piece(
        DISC, p1x17, DnumP,
        lambda a, b: 3 ** a * 5 ** (3 * (dP - a)) * 18 ** b
        * 2 ** (2 * (dD - b)) * 5 ** (3 * (dD - b)))
    del p1x17, DnumP
    random.seed(53)
    ok9aP = True
    for _ in range(6):
        xv = Fr(random.randint(1, 15), 16)
        uv = Fr(random.randint(1, 16), 16)
        tv = Fr(random.randint(0, 32), 16)
        sgv2 = Fr(7 + 3 * xv, 10)
        pv = Fr(4, 5) * (1 - sgv2 ** 2)
        Dv = Fr(2, 25) * uv * (9 - 5 * sgv2 ** 2) * (1 - sgv2) ** 2 \
            / (1 + sgv2 ** 2)
        dv = disc_eval_fr(DISC, pv, Dv, tv)
        tv2 = dict_eval_fr(targetP, xv, uv, tv)
        ok9aP &= (tv2 == Fr(5) ** (3 * (dP + dD)) * Fr(2) ** (2 * dD)
                  * (149 + 42 * xv + 9 * xv ** 2) ** dD * dv)
    rep("R9Pa target_P == 5^396 2^156 (149+42x+9x^2)^78 DISC(subst): 6 "
        "exact points", ok9aP)
    coreP, stripsP, contentP = strip_target(
        targetP, [("u", 'u'), ("t", 't'), ("1-x", [1, -1]),
                  ("7+3x", [7, 3]), ("17+3x", [17, 3]),
                  ("149+42x+9x^2", [149, 42, 9]),
                  ("131-42x-9x^2", [131, -42, -9])])
    log(f"R9P strips: {stripsP}; content bits {contentP.bit_length()}")
    ok9bP = (stripsP.get('u', 0) == 14 and stripsP.get('t', 0) == 7
             and stripsP.get('1-x', 0) == 34
             and stripsP.get('7+3x', 0) == 48
             and stripsP.get('149+42x+9x^2', 0) == 2
             and stripsP.get('131-42x-9x^2', 0) == 14
             and stripsP.get('17+3x', 0) == 0 and contentP > 0)
    rep("R9Pb target_P == content u^14 t^7 (1-x)^34 (7+3x)^48 "
        "(149+42x+9x^2)^2 (131-42x-9x^2)^14 coreP", ok9bP)
    CSP = max(k >> 20 for k in coreP)
    CUP = max((k >> 10) & 1023 for k in coreP)
    CTP = max(k & 1023 for k in coreP)
    degsP = [CSP, CUP, CTP]
    log(f"R9P coreP: {len(coreP)} terms, degs ({CSP},{CUP},{CTP})")
    rep("R9Pc coreP: 153982 terms, degs (306,64,13)",
        len(coreP) == 153982 and (CSP, CUP, CTP) == (306, 64, 13))
    random.seed(61)
    ok9cP = True
    for _ in range(4):
        xv = Fr(random.randint(1, 15), 16)
        uv = Fr(random.randint(1, 16), 16)
        tv = Fr(random.randint(1, 31), 16)
        cv = dict_eval_fr(coreP, xv, uv, tv)
        tv2 = dict_eval_fr(targetP, xv, uv, tv)
        recon = (contentP * uv ** 14 * tv ** 7 * (1 - xv) ** 34
                 * (7 + 3 * xv) ** 48 * (149 + 42 * xv + 9 * xv ** 2) ** 2
                 * (131 - 42 * xv - 9 * xv ** 2) ** 14 * cv)
        ok9cP &= (recon == tv2)
    rep("R9Pd strip reconstruction identity at 4 exact points", ok9cP)
    del targetP

    # ---------------- R9d-P: piece P zero geography ----------------
    uY = sp.Symbol('u')
    okg = True
    for xv, lbl, want_const in ((Fr(0), "x=0 (sg=7/10)", False),
                                (Fr(1), "x=1 (sg=1)", True)):
        c1 = {}
        for k, v in coreP.items():
            a, b, cc = k >> 20, (k >> 10) & 1023, k & 1023
            c1[b] = c1.get(b, 0) + v * xv ** a * Fr(2) ** cc
        den = 1
        for vv in c1.values():
            den = den * vv.denominator // gcd(den, vv.denominator)
        pv2 = sp.Poly(sum(sp.Integer(int(vv * den)) * uY ** b
                          for b, vv in c1.items()), uY)
        if want_const:
            okg &= (pv2.degree() <= 0 and pv2.eval(0) > 0)
        else:
            okg &= (pv2.degree() == 62 and pv2.count_roots(0, 1) == 0
                    and pv2.eval(sp.Rational(1, 2)) > 0)
    rep("R9d-P t=2 edge slices: x=0 deg-62 root-free; x=1 positive "
        "constant", okg)

    # ---------------- R10P: the single piece P certificate ----------
    fu1P = face_u1(coreP)
    ok10P = no_negatives(fu1P, [CSP, 0, CTP], [(0, 1), (0, 1), (0, 2)],
                         "R10Pf u=1 face, t in [0,2]", axes=(2, 0))
    rep("R10Pf piece P u=1 face: no negatives (open-face > 0; sanity -- "
        "the root box below subsumes it)", ok10P)
    del fu1P
    okP = closed_strict(coreP, degsP, [(0, 1), (0, 1), (0, 2)],
                        "R10P root box [0,1]x[0,1]x[0,2]")
    rep("R10P coreP CLOSED-STRICT on the WHOLE box [0,1]x[0,1]x[0,2]", okP)
    del coreP

    # ---------------- R14: coverage ----------------
    okcov = (Fr(7, 4) < Fr(15, 8) < Fr(127, 64) < 2
             and Fr(7, 12) < Fr(29, 48) < Fr(5, 8) < Fr(31, 48)
             < Fr(2, 3) < 1
             and Fr(29, 48) == Fr(5, 8) - Fr(1, 48)
             and Fr(31, 48) == Fr(5, 8) + Fr(1, 48)
             and Fr(7, 12) == Fr(5, 8) - Fr(1, 24)
             and Fr(2, 3) == Fr(5, 8) + Fr(1, 24)
             and Fr(127, 64) == 2 - Fr(1, 64)
             and Fr(1, 8) ** 2 == Fr(1, 64)
             and Fr(1, 64) < Fr(1, 48))
    rep("R14 coverage arithmetic (t and u unions, window bounds, "
        "sqrt(1/64) = 1/8)", okcov)
    # two-piece union over sg
    sgS = sp.Symbol('sigma')
    okuni = (sp.expand((9 - 5 * sgS ** 2) / sp.Integer(10)
                       - (1 - sp.Rational(3, 4) * sgS ** 2)
                       - (5 * sgS ** 2 - 2) / sp.Integer(20)) == 0
             and 5 * Fr(7, 10) ** 2 - 2 == Fr(9, 20)
             and Fr(7, 10) < Fr(3, 4))
    rep("R14 union: curveP - curveM == (5sg^2-2)/20 > 0 on overlap; "
        "(0,3/4) u [7/10,1) = (0,1)", okuni)
    log("R14 piece M union: u=1: F1 face;  t<15/8: C1;  t>=15/8:")
    log("    u<7/12: C2 | u in [7/12,2/3]: C4 (x>=1/64) | W3a (t<=127/64)")
    log("    | W3b/W3c (|u-5/8|>=1/48) | W3d charts (window, rho>0) |")
    log("    u>2/3: C3a (x<5/8) + C3b (x>=5/8)  ==> coreM > 0 on Omega_M")
    log("R14 piece P: R10P covers the whole closed box  ==> coreP > 0 on "
        "Omega_P")

    # ---------------- R15: anchors ----------------
    # piece M: (x,u,t) = (2/3, 3/4, 1) -> sg = 1/2, p = 3/5, D = 39/400
    pvM = sp.Rational(3, 5)
    DvM = sp.Rational(39, 400)
    okanM = (sp.cancel(DSUB_M.subs({xx: sp.Rational(2, 3),
                                    uu: sp.Rational(3, 4)}) - DvM) == 0
             and sp.Rational(4, 5) * (1 - (sp.Rational(3, 4)
                                           * sp.Rational(2, 3)) ** 2) == pvM)
    # piece P: (x,u,t) = (1/3, 1/2, 1) -> sg = 4/5, p = 36/125, D = 29/5125
    pvP = sp.Rational(36, 125)
    DvP = sp.Rational(29, 5125)
    okanP = (sp.cancel(DSUB_P.subs({xx: sp.Rational(1, 3),
                                    uu: sp.Rational(1, 2)}) - DvP) == 0
             and sp.Rational(4, 5) * (1 - sp.Rational(4, 5) ** 2) == pvP)
    Yv = sp.Symbol('Y')
    okan = okanM and okanP
    for pv, Dv in ((pvM, DvM), (pvP, DvP)):
        psi = Yv ** 5
        for k in range(1, 6):
            psi += (sp.Integer(2 * MSCALE) ** k
                    * Etld[k].subs({p: pv, D: Dv, t: 1})) * Yv ** (5 - k)
        nreal = sp.Poly(psi, Yv).count_roots()
        danch = sum(Fr(cc) * Fr(int(pv.p), int(pv.q)) ** (k >> 16)
                    * Fr(int(Dv.p), int(Dv.q)) ** ((k >> 8) & 255)
                    for k, cc in DISC.items())    # t = 1: t-powers all 1
        okan &= (nreal == 5 and danch > 0)
    rep("R15 anchors M (2/3,3/4,1)->(3/5,39/400), P (1/3,1/2,1)->"
        "(36/125,29/5125): 5 distinct real roots, DISC > 0", okan)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("disc(psi) > 0 on Omega_M and Omega_P; anchors have 5 distinct")
    print("real roots; each Omega connected => the A=5 replica spectrum is")
    print("REAL and SIMPLE for all 0<p<4/5, 0<D<=curve* Dbar, 0<s<pi,")
    print("curve* = curveM = 1-(3/4)sg^2 on (0,7/10] / curveP = (9-5sg^2)/10")
    print("on [7/10,1).  With the s=pi lemma (lemma_7_34_a5_realness_pi.py):")
    print("all angles in (0,pi] on the two-piece curve-restricted domain.")
    return PASS


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
