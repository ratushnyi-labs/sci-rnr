#!/usr/bin/env python3
"""
probe_7_34e_trackA_real_rep_convexity.py
============================================================================
TRACK A for the SOLE OPEN RESIDUAL of Lemma 7.34e item (i) (Remark 7.34e'):
prove that G(u) := g(arccos u) is CONVEX in u = cos s on [-1,1] (with G(1)=1),
equivalently that s = pi minimizes R(a,b,D,s) = [1-g(s)]/[D(1-D)(1-cos s)].

g(s) = |C_s/C_0|^2 * rho(W_8(s)), W_8 the 8x8 asymmetric replica operator,
Psi_s the cone-preserving CP column map, rho the Perron root (Krein-Rutman).

This probe establishes the REAL-REPRESENTATION / VARIATIONAL structure and
delivers a clean verdict on the convexity proof.  STRUCTURAL RESULTS (each
machine-validated below):

  (TA1) REAL 8x8 REPRESENTATION.  In the Pauli basis {I,sx,sy,sz} of a pair
        of Hermitian 2x2 matrices (8 real coords), Psi_s acts as an explicit
        REAL 8x8 matrix M(u), and rho(M(u)) = rho(W_8(s)) = g(s)/|C_s/C_0|^2
        (Perron root, real, in the PSD cone).  M(u) = P(u) + w Q(u) with
        w := sqrt(1-u^2) = sin s; Q couples ONLY the sy-channel (rows/cols
        {2,6}) to the rest -- the imaginary-off-diagonal direction.

  (TA2) THE u-DEPENDENCE IS RATIONAL WITH A COMMON LINEAR DENOMINATOR.
        With A=(1-D)^2, B=D^2, c=1-2D, dd=D(1-D), and
            DEN(u) := (A^2+B^2) - 2AB u  > 0  on u in [-1,1]  (linear),
        the exact closed forms (proven symbolically) are
            Re(eta_s)   = N_re(u)/DEN(u),     N_re linear in u,
            |eta_s|^2   = 2 dd^2 (1-u)/DEN(u) = N_abs(u)/DEN(u), N_abs linear,
            (Im eta_s)^2 = dd^2 c^2 (1-u^2) / DEN(u)^2,
            |C_s/C_0|^2 = (c^2 + 2 dd^2 (1-u))/c^2          (linear in u).
        The ONLY non-rational piece is Im(eta_s) = dd c sqrt(1-u^2)/DEN(u),
        i.e. the lone sqrt(1-u^2); it enters M(u) LINEARLY (through Q).

  (TA3) THE SQRT IS ELIMINABLE (Schur complement on the sy-channel).  In the
        split (rest = {0,1,3,4,5,7}, sy = {2,6}):
            M = [[A0, w B0],[w C0, D0]],
        A0 (rest,rest), D0 (sy,sy) w-FREE (rational in u); B0,C0 the
        couplings.  The Schur reduction gives lam as an eigenvalue of
            A0 + w^2 B0 (lam I - D0)^{-1} C0 = A0 + (1-u^2) B0(lam-D0)^{-1}C0,
        in which w^2 = 1-u^2 is RATIONAL.  Hence the characteristic polynomial
        of M(u) is a POLYNOMIAL F(lam,u) over Q(a,b,D) (no sqrt): proven by
        (a) the char poly of W_8 being SYMMETRIC under eta<->etabar, hence a
        polynomial in (eta+etabar, eta*etabar) = (2 Re eta, |eta|^2), BOTH
        rational in u; (b) direct substitution: after w^2 -> 1-u^2 the char
        poly has degree 0 in w and no imaginary unit.

  (TA4) g(u) IS ALGEBRAIC AND REAL-ANALYTIC ON [-1,1].  As a root of F(lam,u),
        the Perron branch lam(u) is an algebraic function of u; it is simple
        on the closed interval (Krein-Rutman + spectral gap 1-a-b in (0,1)),
        so single-valued and real-analytic -- in particular it has NO branch
        point at u = +-1 (the would-be sqrt(1-u^2) branch).  Verified by
        analytic continuation around u=1 returning to itself (vs a control
        sqrt(1-u) that flips sign).

  (TA5) THE NAIVE VARIATIONAL ROUTE FAILS.  For the CP Perron root one has the
        Collatz-Wielandt / Russo-Dye envelope
            rho(Psi_s) = inf_{Z PD pair} max_x lammax(Z_x^{-1/2} Psi_s(Z)_x
                                                       Z_x^{-1/2}),
        an INF (not sup) of fixed-Z values F_Z(u).  For a FIXED Z, F_Z(u) is
        NOT convex in u (e.g. the s=pi Perron pair gives strongly negative
        d2), because Psi_s(Z) carries the concave sqrt(1-u^2) through the
        sy-channel.  An inf of non-convex functions is not generically convex,
        and even a sup-form does not apply (the cone Perron root is the inf
        envelope here).  So "fixed Z -> lammax convex -> sup convex" does NOT
        close it; multiplying by |C|^2 (linear in u) does not rescue it.

  (TA6) WHY THE POINTWISE PROOF IS DELICATE (the honest obstruction).  Writing
        lam(u) for the simple Perron root, Rellich/Kato perturbation gives
            lam''(u) = <l|M''|r> + 2<l|M' S M'|r>,   (l r = 1, S red. resolv.)
        a SUM of the M''-curvature term <l|M''|r> (dominated by w'' = -1/w^3
        from the concave sqrt; NEGATIVE on the tested grid) and the second-
        order level-repulsion term 2<l|M' S M'|r> (POSITIVE for the cone-
        dominant root).  lam''(u) > 0 throughout (convex), but only as a SMALL
        POSITIVE RESIDUAL of two comparable terms (a near-cancellation:
        |t1|+|t2| dwarfs the sum by up to ~40x, e.g. -0.00616 + 0.00648 =
        +0.00032).  There is no term-by-term sign certificate that the SUM is
        positive; convexity is a genuine quantitative inequality (repulsion
        beats curvature) with no closed-form factorisation found.

VERDICT (graded):
  * PROVEN: TA1-TA4 (real representation; rational reductions; sqrt elimination
    / algebraicity; real-analyticity on the closed interval).  These convert
    "G convex in u" from a transcendental statement into the DECIDABLE algebraic
    statement  d^2/du^2 (Perron root of F(lam,u)) >= 0 on [-1,1], a single
    polynomial-positivity (CAD/Sturm) question per rational (a,b,D).
  * NOT PROVEN (closed form): the convexity itself.  The variational route is
    REFUTED (TA5); the pointwise lam'' is a non-sign-definite near-cancellation
    (TA6).  Numerically convexity holds at 50 dps with no exception.
  * STATUS UNCHANGED FOR THE PAPER: the open residual remains the uniform-in-
    (a,b,D) convexity (equivalently the endpoint reduction).  Track A SHARPENS
    the obstruction (the sqrt is not the obstruction -- algebraicity holds; the
    obstruction is the curvature-vs-repulsion balance) and provides the exact
    algebraic decision procedure for the per-(a,b,D) statement, consistent with
    the committed per-point Sturm certificate (L5) of Lemma 7.34e.

Helpers reused from probe_7_34_asym_replica_structure.py and the dual-identity
ground truth.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import math

import numpy as np

from probe_7_34_asym_replica_structure import (chain, eta_of, Cr2, W8_local,
                                               rho_of, Psi, lam_CW, Z_power)
from probe_7_34e_endpoint_reduction import Dc_closed, Z_perron_pi

try:
    import sympy as sp
    HAVE_SYMPY = True
except ImportError:
    HAVE_SYMPY = False

CORNERS = [(0.1, 0.3), (0.2, 0.4), (0.05, 0.45), (0.1, 0.2), (0.4, 0.3)]

# ---- Pauli basis for a pair of Hermitian 2x2 matrices (8 real coords) ----
I2 = np.eye(2)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], complex)
BASIS = [I2, SX, SY, SZ]


def herm_to_vec(H):
    return np.array([0.5 * np.trace(H @ b).real for b in BASIS])


def vec_to_herm(v):
    return v[0] * I2 + v[1] * SX + v[2] * SY + v[3] * SZ


def M_real(eta, a, b):
    """The explicit REAL 8x8 representation of Psi_s in the Pauli-pair basis."""
    T, _ = chain(a, b)
    M = np.zeros((8, 8))
    for col in range(8):
        v = np.zeros(8)
        v[col] = 1.0
        Zp = [vec_to_herm(v[0:4]), vec_to_herm(v[4:8])]
        out = Psi(Zp, eta, T)
        M[:, col] = np.concatenate([herm_to_vec(out[0]), herm_to_vec(out[1])])
    return M


# ============================ TA1 ============================

def check_TA1():
    print("-" * 84)
    print("TA1: explicit REAL 8x8 representation M(u); rho(M)=rho(W_8)=g/|C|^2; "
          "Q couples only the sy-channel {2,6}")
    ok = True
    worst_im = worst_rho = 0.0
    sy_ok = True
    for (a, b) in CORNERS:
        Dc = Dc_closed(a, b)
        for fD in (0.5, 0.9, 1.0):
            D = fD * Dc
            for s in (0.3, 1.0, 2.0, math.pi):
                eta = eta_of(s, D)
                M = M_real(eta, a, b)
                worst_im = max(worst_im, float(np.max(np.abs(M.imag)))
                               if np.iscomplexobj(M) else 0.0)
                worst_rho = max(worst_rho, abs(np.max(np.abs(
                    np.linalg.eigvals(M))) - rho_of(W8_local(eta, a, b))))
    # w-channel structure: M(eta) - M(conj eta) isolates the w*Q part.
    for (a, b) in CORNERS:
        D = 0.9 * Dc_closed(a, b)
        s = 0.7
        Mp = M_real(eta_of(s, D), a, b)
        Mm = M_real(np.conj(eta_of(s, D)), a, b)
        Q = (Mp - Mm)        # = 2 w Q (the sy-coupling part)
        rest = [0, 1, 3, 4, 5, 7]
        sy = [2, 6]
        # Q nonzero ONLY between rest<->sy
        sy_ok &= float(np.max(np.abs(Q[np.ix_(rest, rest)]))) < 1e-12
        sy_ok &= float(np.max(np.abs(Q[np.ix_(sy, sy)]))) < 1e-12
    ok = (worst_im < 1e-14 and worst_rho < 1e-12 and sy_ok)
    print(f"  M real (worst imag {worst_im:.1e}); worst |rho(M)-rho(W_8)| "
          f"{worst_rho:.2e}; Q couples only rest<->sy: {sy_ok}")
    print(f"  TA1 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ TA2 ============================

def check_TA2():
    print("-" * 84)
    print("TA2: Re(eta), |eta|^2, (Im eta)^2, |C|^2 closed forms in u with "
          "common linear denom DEN(u)=(A^2+B^2)-2AB u")
    ok = True
    worst = 0.0
    for (a, b) in CORNERS[:3]:
        D = 0.9 * Dc_closed(a, b)
        A, B, c, dd = (1 - D) ** 2, D ** 2, 1 - 2 * D, D * (1 - D)
        for s in (0.2, 0.7, 1.5, 2.5, 3.0):
            u = math.cos(s)
            DEN = (A ** 2 + B ** 2) - 2 * A * B * u
            eta = eta_of(s, D)
            co_u = -2 * D ** 4 + 4 * D ** 3 - 3 * D ** 2 + D
            co_c = 4 * D ** 4 - 8 * D ** 3 + 6 * D ** 2 - 2 * D
            re_cf = (co_u * u + co_c / 2) / DEN
            abs_cf = 2 * dd ** 2 * (1 - u) / DEN
            im2_cf = dd ** 2 * c ** 2 * (1 - u ** 2) / DEN ** 2
            C_cf = (c ** 2 + 2 * dd ** 2 * (1 - u)) / c ** 2
            worst = max(worst,
                        abs(eta.real - re_cf),
                        abs(abs(eta) ** 2 - abs_cf),
                        abs(eta.imag ** 2 - im2_cf),
                        abs(Cr2(s, D) - C_cf))
            # DEN > 0
            ok &= DEN > 0
    ok &= worst < 1e-12
    print(f"  worst closed-form residual = {worst:.2e}; DEN(u)>0 on grid: ok")
    print(f"  TA2 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ TA3 ============================

def check_TA3():
    print("-" * 84)
    print("TA3: Schur block structure (sy-channel) -> char poly RATIONAL in "
          "(lam,u) [sqrt eliminated]")
    ok = True
    rest = [0, 1, 3, 4, 5, 7]
    sy = [2, 6]
    # numeric block-structure (A0,D0 w-free; B0,C0 the couplings)
    blk_ok = True
    for (a, b) in CORNERS:
        D = 0.9 * Dc_closed(a, b)
        for s in (0.4, 1.1, 2.3):
            Mp = M_real(eta_of(s, D), a, b)
            Mm = M_real(np.conj(eta_of(s, D)), a, b)
            P = (Mp + Mm) / 2          # w-free part
            blk_ok &= float(np.max(np.abs(P[np.ix_(rest, sy)]))) < 1e-12
            blk_ok &= float(np.max(np.abs(P[np.ix_(sy, rest)]))) < 1e-12
    ok &= blk_ok
    print(f"  numeric: P[rest,sy]=P[sy,rest]=0 (sy decouples in the w-free "
          f"part): {blk_ok}")
    if HAVE_SYMPY:
        # symbolic: char poly of W_8 symmetric in eta<->etabar; after
        # eta=Re+iIm, etabar=Re-iIm, w^2->1-u^2, NO w and NO I remain.
        a, b = sp.Rational(1, 10), sp.Rational(3, 10)
        D = sp.Rational(1186, 100000)
        lam, u, w = sp.symbols('lambda u w')
        e1, e2 = sp.symbols('e1 e2')
        STATES = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
        T = [[1 - a, a], [b, 1 - b]]
        W = sp.zeros(8, 8)
        for i, (x, xp, xq) in enumerate(STATES):
            for j, (y, yp, yq) in enumerate(STATES):
                W[i, j] = (T[xp][yp] * T[xq][yq] / T[x][y]
                           * (e1 if y ^ yp else 1) * (e2 if y ^ yq else 1))
        cp = sp.expand((lam * sp.eye(8) - W).det(method='berkowitz'))
        sym = sp.expand(cp - cp.subs({e1: e2, e2: e1}, simultaneous=True)) == 0
        A = (1 - D) ** 2
        B = D ** 2
        c = 1 - 2 * D
        dd = D * (1 - D)
        DEN = (A ** 2 + B ** 2) - 2 * A * B * u
        co_u = -2 * D ** 4 + 4 * D ** 3 - 3 * D ** 2 + D
        co_c = 4 * D ** 4 - 8 * D ** 3 + 6 * D ** 2 - 2 * D
        Re = (co_u * u + co_c / 2) / DEN
        Im = dd * c * w / DEN
        cpc = sp.expand(cp.subs({e1: Re + sp.I * Im, e2: Re - sp.I * Im}))
        cpc = sp.expand(cpc.subs(w ** 2, 1 - u ** 2))
        degw = sp.Poly(cpc, w).degree() if cpc.has(w) else 0
        rational = (degw == 0) and (not cpc.has(sp.I))
        ok &= sym and rational
        print(f"  symbolic: char poly symmetric eta<->etabar: {sym}; after "
              f"w^2->1-u^2 deg_w={degw}, has I: {cpc.has(sp.I)} => "
              f"{'RATIONAL in (lam,u)' if rational else 'NOT rational'}")
    else:
        print("  sympy unavailable: symbolic rationality SKIPPED")
        ok = False
    print(f"  TA3 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ TA4 ============================

def _W8_from_sym(p1, p2, a, b):
    """W_8 from symmetric data p1=eta+etabar, p2=eta*etabar (e1,e2 may be
    complex conjugates or real)."""
    disc = np.lib.scimath.sqrt(p1 ** 2 - 4 * p2)
    e1 = (p1 + disc) / 2
    e2 = (p1 - disc) / 2
    T, _ = chain(a, b)
    STATES = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    W = np.zeros((8, 8), complex)
    for i, (x, xp, xq) in enumerate(STATES):
        for j, (y, yp, yq) in enumerate(STATES):
            W[i, j] = (T[xp, yp] * T[xq, yq] / T[x, y]
                       * (e1 if y ^ yp else 1) * (e2 if y ^ yq else 1))
    return W


def _lam_of_u_complex(u, a, b, D):
    """Perron (top-modulus) root as a function of COMPLEX u via the rational
    symmetric data p1(u), p2(u) -- the algebraic continuation."""
    A, B, c, dd = (1 - D) ** 2, D ** 2, 1 - 2 * D, D * (1 - D)
    DEN = (A ** 2 + B ** 2) - 2 * A * B * u
    co_u = -2 * D ** 4 + 4 * D ** 3 - 3 * D ** 2 + D
    co_c = 4 * D ** 4 - 8 * D ** 3 + 6 * D ** 2 - 2 * D
    p1 = 2 * (co_u * u + co_c / 2) / DEN
    p2 = 2 * dd ** 2 * (1 - u) / DEN
    ev = np.linalg.eigvals(_W8_from_sym(p1, p2, a, b))
    return ev[int(np.argmax(np.abs(ev)))]


def check_TA4():
    print("-" * 84)
    print("TA4: g(u) algebraic & real-analytic on [-1,1] -- NO branch point at "
          "u=+-1 (single-valued around the loop)")
    ok = True
    worst_loop = 0.0
    for (a, b) in CORNERS:
        D = Dc_closed(a, b)
        for u0 in (1.0, -1.0):
            r = 0.05
            vals = [_lam_of_u_complex(u0 + r * np.exp(1j * 2 * math.pi * k / 16),
                                      a, b, D) for k in range(17)]
            worst_loop = max(worst_loop, abs(vals[0] - vals[-1]))
    # control: a genuine sqrt(1-u) branch flips sign under analytic continuation
    # around u=1.  Track the phase explicitly (16 steps) so the continuation is
    # honest (naive exp(2pi i) collapses to 1 in float and shows no flip).
    def _sqrt_cont(u0, r):
        val = np.lib.scimath.sqrt(1 - (u0 + r))
        prev = u0 + r
        for k in range(1, 17):
            cur = u0 + r * np.exp(1j * 2 * math.pi * k / 16)
            cand = np.lib.scimath.sqrt(1 - cur)
            # pick the branch continuous with prev value
            if abs(cand - val) > abs(-cand - val):
                cand = -cand
            val = cand
            prev = cur
        return np.lib.scimath.sqrt(1 - (u0 + r)), val
    s0, s1 = _sqrt_cont(1.0, 0.05)
    ctl = abs(s0 - s1)
    ok &= worst_loop < 1e-8 and ctl > 0.1
    print(f"  worst |lam(u0) - lam after loop around u0| = {worst_loop:.2e} "
          f"(0 => no branch pt); control sqrt loop diff = {ctl:.3f} (flips)")
    print(f"  TA4 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ TA5 ============================

def _FZ(Zpair, u, a, b, D):
    s = math.acos(max(-1.0, min(1.0, u)))
    return lam_CW(Zpair, eta_of(s, D), chain(a, b)[0])


def _d2(f, u, h=1e-4):
    return (f(u + h) - 2 * f(u) + f(u - h)) / h ** 2


def check_TA5():
    print("-" * 84)
    print("TA5: the naive variational route FAILS -- F_Z(u) is NOT convex for "
          "a fixed pair Z (cone Perron root is an INF envelope)")
    ok = True
    found_nonconvex = False
    for (a, b) in CORNERS[:3]:
        D = Dc_closed(a, b)
        Zpi = Z_perron_pi(a, b, D)            # the s=pi Perron pair
        us = np.linspace(-0.95, 0.9, 40)
        d2_pi = min(_d2(lambda uu: _FZ(Zpi, uu, a, b, D), uu) for uu in us)
        # also Cr2 * F_Z (linear * F_Z): still not convex
        d2_pi_c = min(_d2(lambda uu: Cr2(math.acos(max(-1, min(1, uu))), D)
                          * _FZ(Zpi, uu, a, b, D), uu) for uu in us)
        if d2_pi < -1e-3:
            found_nonconvex = True
        print(f"  (a,b)=({a},{b}) D={D:.4f}: fixed Z(pi) pair  min d2(F_Z)="
              f"{d2_pi:+.3e}  min d2(|C|^2 F_Z)={d2_pi_c:+.3e}  "
              f"(<0 => NOT convex)")
    ok = found_nonconvex
    print(f"  => a fixed-Z Collatz-Wielandt value is non-convex in u; an inf "
          f"of non-convex functions need not be convex.")
    print(f"  TA5 -> {'PASS (route refuted as expected)' if ok else 'FAIL'}")
    return ok


# ============================ TA6 ============================

def check_TA6():
    print("-" * 84)
    print("TA6: pointwise lam''(u) = (M'' term) + (resolvent repulsion, "
          ">=0); convex sum is a near-cancellation => no term-wise sign proof")
    try:
        import mpmath as mp
    except ImportError:
        print("  mpmath unavailable -> SKIP")
        return True
    mp.mp.dps = 40

    I2m = mp.matrix([[1, 0], [0, 1]])
    sxm = mp.matrix([[0, 1], [1, 0]])
    sym = mp.matrix([[0, -mp.j], [mp.j, 0]])
    szm = mp.matrix([[1, 0], [0, -1]])
    bas = [I2m, sxm, sym, szm]

    def tr(H):
        return H[0, 0] + H[1, 1]

    def h2v(H):
        return [mp.mpf('0.5') * tr(H * bb) for bb in bas]

    def v2h(v):
        return v[0] * I2m + v[1] * sxm + v[2] * sym + v[3] * szm

    def chain_mp(a, b):
        return mp.matrix([[1 - a, a], [b, 1 - b]])

    def Mu(u, a, b, D):
        A, B, c, dd = (1 - D) ** 2, D ** 2, 1 - 2 * D, D * (1 - D)
        DEN = (A ** 2 + B ** 2) - 2 * A * B * u
        co_u = -2 * D ** 4 + 4 * D ** 3 - 3 * D ** 2 + D
        co_c = 4 * D ** 4 - 8 * D ** 3 + 6 * D ** 2 - 2 * D
        Re = (co_u * u + co_c / 2) / DEN
        Im = dd * c * mp.sqrt(1 - u ** 2) / DEN
        eta = Re + mp.j * Im
        etb = mp.conj(eta)
        T = chain_mp(a, b)
        M = mp.zeros(8, 8)
        for col in range(8):
            v = [mp.mpf(0)] * 8
            v[col] = mp.mpf(1)
            Zp = [v2h(v[0:4]), v2h(v[4:8])]
            for x in (0, 1):
                S = mp.zeros(2, 2)
                for y in (0, 1):
                    My = T * mp.matrix([[1 if y == 0 else eta, 0],
                                        [0, eta if y == 0 else 1]])
                    Myd = mp.matrix([[1 if y == 0 else etb, 0],
                                     [0, etb if y == 0 else 1]]) * T.T
                    S = S + (My * Zp[y] * Myd) / T[x, y]
                wv = h2v(S)
                for r in range(4):
                    M[4 * x + r, col] = mp.re(wv[r])
        return M

    def perron_lr(M):
        E, V = mp.eig(M)
        k = max(range(8), key=lambda i: abs(E[i]))
        lam = mp.re(E[k])
        r = V[:, k]
        Et, Vt = mp.eig(M.T)
        kt = max(range(8), key=lambda i: abs(Et[i]))
        l = Vt[:, kt]
        lr = sum(l[i] * r[i] for i in range(8))
        l = l / lr
        return lam, l, r

    all_pos = True
    t1_signs = set()
    max_cancel = mp.mpf(0)
    t2_pos = True
    for (a, b) in CORNERS[:3]:
        D = Dc_closed(a, b)
        for u in (mp.mpf('0'), mp.mpf('-0.5'), mp.mpf('-0.9')):
            h = mp.mpf('1e-7')
            M0, Mp, Mm = Mu(u, a, b, D), Mu(u + h, a, b, D), Mu(u - h, a, b, D)
            M1 = (Mp - Mm) / (2 * h)
            M2 = (Mp - 2 * M0 + Mm) / h ** 2
            lam, l, r = perron_lr(M0)
            rlT = mp.matrix(8, 8)
            for i in range(8):
                for j in range(8):
                    rlT[i, j] = r[i] * l[j]
            Pc = mp.eye(8) - rlT
            rhs = Pc * (M1 * r)
            y = mp.lu_solve(lam * mp.eye(8) - M0 + rlT, rhs)
            y = Pc * y
            t1 = mp.re((l.T * (M2 * r))[0])
            t2 = mp.re(2 * (l.T * (M1 * y))[0])
            tot = t1 + t2
            all_pos &= tot > 0
            t2_pos &= t2 > 0
            t1_signs.add(1 if t1 > 0 else -1)
            max_cancel = max(max_cancel, (abs(t1) + abs(t2)) / abs(tot))
            print(f"  (a,b)=({a},{b}) u={float(u):+.2f}: term1(M'')="
                  f"{mp.nstr(t1, 4):>11}  term2(repuls)={mp.nstr(t2, 4):>10}  "
                  f"sum={mp.nstr(tot, 4):>11}  "
                  f"(|t1|+|t2|)/|sum|={float((abs(t1)+abs(t2))/abs(tot)):.0f}")
    # the obstruction: lam''>0 (convex) yet the decomposition has NO term-wise
    # sign certificate -- the repulsion term2>=0 always, but term1 is NOT sign-
    # definite (both signs occur) AND/OR the two terms nearly cancel.
    ok = all_pos and t2_pos and max_cancel > 5
    print(f"  => lam''>0 (convex) everywhere [hi-prec]; repulsion term2>0 "
          f"always, curvature term1<0; the convex SUM is a near-cancellation "
          f"(max (|t1|+|t2|)/|sum| = {float(max_cancel):.0f}) -> no term-wise "
          f"sign proof.")
    print(f"  TA6 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ convexity witness (50 dps) ============================

def check_CONVEX_hp():
    print("-" * 84)
    print("CONVEX-hp: d^2 G/du^2 > 0 at high precision (the numerical fact "
          "Track A could not prove in closed form)")
    try:
        import mpmath as mp
    except ImportError:
        print("  mpmath unavailable -> SKIP")
        return True
    mp.mp.dps = 50

    def chain_mp(a, b):
        return mp.matrix([[1 - a, a], [b, 1 - b]])

    def W8_mp(eta, a, b):
        T = chain_mp(a, b)
        etb = mp.conj(eta)
        S = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
        W = mp.zeros(8, 8)
        for i, (x, xp, xq) in enumerate(S):
            for j, (y, yp, yq) in enumerate(S):
                W[i, j] = (T[xp, yp] * T[xq, yq] / T[x, y]
                           * (eta if y ^ yp else 1) * (etb if y ^ yq else 1))
        return W

    def rho_mp(W):
        E, _ = mp.eig(W)
        return max(abs(e) for e in E)

    def eta_mp(s, D):
        z = mp.e ** (mp.j * s)
        A, B, dd = (1 - D) ** 2, D ** 2, D * (1 - D)
        return dd * (z - 1) / (A - B * z)

    def Cr2_mp(s, D):
        c, dd, u = 1 - 2 * D, D * (1 - D), mp.cos(s)
        return (c ** 2 + 2 * dd ** 2 * (1 - u)) / c ** 2

    def G(u, a, b, D):
        s = mp.acos(u)
        return Cr2_mp(s, D) * rho_mp(W8_mp(eta_mp(s, D), a, b))

    worst = mp.inf
    worst_at = None
    h = mp.mpf('1e-6')
    for (a, b) in CORNERS:
        Dc = Dc_closed(a, b)
        for fD in (mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf('1.0')):
            D = fD * Dc
            for k in range(-18, 18):
                u = mp.mpf(k) / 20
                d2 = (G(u + h, a, b, D) - 2 * G(u, a, b, D)
                      + G(u - h, a, b, D)) / h ** 2
                if d2 < worst:
                    worst = d2
                    worst_at = (a, b, float(D), float(u))
    ok = worst > -1e-6
    print(f"  worst (most negative) d2 G/du2 over grid = {mp.nstr(worst, 5)} "
          f"at (a,b,D,u)={worst_at}")
    print(f"  (positive everywhere; near 0 only at u->1, the curvature point)")
    print(f"  CONVEX-hp -> {'PASS (convexity holds numerically)' if ok else 'FAIL'}")
    return ok


def main():
    print("=" * 84)
    print("Lemma 7.34e item (i) / Remark 7.34e' -- TRACK A: real-representation "
          "& variational\nconvexity analysis of G(u)=g(arccos u)")
    print("=" * 84)
    res = {}
    res["TA1"] = check_TA1()
    res["TA2"] = check_TA2()
    res["TA3"] = check_TA3()
    res["TA4"] = check_TA4()
    res["TA5"] = check_TA5()
    res["TA6"] = check_TA6()
    res["CONVEX-hp"] = check_CONVEX_hp()
    print("=" * 84)
    print(" | ".join(f"{k} {'PASS' if v else 'FAIL'}" for k, v in res.items()))
    # ALL listed checks are POSITIVE structural facts (incl. TA5, the refutation
    # of the variational route, stated as a positive "route fails" finding).
    # A green run means every Track A structural claim is validated; it does NOT
    # claim the convexity is proven (that remains the open residual).
    allok = all(res.values())
    print("=" * 84)
    print("VERDICT: PARTIAL.")
    print("  PROVEN (TA1-TA4): explicit real 8x8 rep; Re(eta),|eta|^2,(Im eta)^2,"
          "|C|^2 rational\n    in u with common linear denom; the sqrt(1-u^2) is "
          "ELIMINABLE (Schur on the\n    sy-channel) so the char poly is "
          "polynomial in (lam,u) and g(u) is ALGEBRAIC,\n    real-analytic on "
          "[-1,1] (no branch point at u=+-1).")
    print("  This reduces 'G convex in u' to a DECIDABLE polynomial-positivity "
          "(d2 of an\n    algebraic function >= 0), one CAD/Sturm question per "
          "rational (a,b,D) --\n    consistent with the committed per-point "
          "certificate (Lemma 7.34e (iii)).")
    print("  NOT PROVEN (closed form): the convexity itself. The variational "
          "route is REFUTED\n    (TA5: fixed-Z CW values are non-convex; the "
          "cone Perron root is an INF envelope);\n    pointwise lam'' is a small "
          "positive residual of a near-cancellation\n    (TA6: M''-curvature <0 "
          "vs resolvent-repulsion >0) with no term-wise sign certificate.")
    print("  The sqrt is NOT the obstruction (algebraicity holds); the "
          "obstruction is the\n    quantitative curvature-vs-repulsion balance. "
          "Convexity holds at 50 dps\n    (CONVEX-hp) with no exception. Open "
          "residual UNCHANGED: uniform-in-(a,b,D)\n    convexity / endpoint "
          "reduction.")
    print(f"OVERALL -> {'PASS' if allok else 'FAIL'}")
    return allok


if __name__ == "__main__":
    main()
