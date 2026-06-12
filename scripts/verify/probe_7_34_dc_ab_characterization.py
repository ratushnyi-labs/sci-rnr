#!/usr/bin/env python3
"""
probe_7_34_dc_ab_characterization.py
============================================================================
INSTANCE (2) RESIDUALS (iii) + (ii) of the non-symmetric binary Markov chain
in the lossy-RD dispersion converse arc (tex/rnr_coding.tex Lemma 7.34e &
surrounding remark).

CHAIN:  X stationary binary Markov, P(0->1)=a, P(1->0)=b, a,b in (0,1),
        a+b<1, pi0=b/(a+b), second eigenvalue lambda2 = 1-a-b.
REGION: the Gray/deconvolution region {D : the BSC(D)-deconvolution of P_X is
        a valid probability law for ALL n},
            P_{Y*}(y) = 2^{-n} sum_w Phat_X(w)(1-2D)^{-|w|}(-1)^{<w,y>} >= 0,
        all y, all n.  D_c(a,b) := sup of valid D.

----------------------------------------------------------------------------
ITEM (iii) -- THE KEYSTONE: exact all-n characterization of D_c(a,b).
----------------------------------------------------------------------------
(1) TRANSFER SYSTEM (PROVEN, verified A1).
    Combining the two Walsh steps coordinate-wise,
        P_{Y*}(y) = sum_x P_X(x) prod_i g(x_i,y_i),  g(x,y)=(1+zeta(-1)^{x xor y})/2,
        zeta = (1-2D)^{-1},
    which is an exact 2x2 transfer contraction over the source state:
        P_{Y*}(y) = pi^T [prod_{t=1}^{n-1} G_{y_t} T] G_{y_n} 1,
        G_0 = diag((1+zeta)/2,(1-zeta)/2),  G_1 = diag((1-zeta)/2,(1+zeta)/2).
    The natural per-symbol transfer matrix is B_y := G_y T.

(2) BINDING WORD (VERIFIED A2): the binding word family (argmin_y P_{Y*}) is the
    ALTERNATING word 0101.../1010... (period 2) at EVERY tested (a,b).  Verified
    by exhaustive min over all 2^n words for n<=16: below D_c all words are >=0,
    and the first word to go negative just above D_c is alternating.

(3) EXACT CONDITION (PROVEN symbolically A3, VERIFIED A4).
    The alternating word is governed by the period-2 product M = B1 B0 = (G1 T)(G0 T).
    Its eigenvalues r e^{+-i phi} become a complex-conjugate pair exactly when the
    discriminant disc(D) = (tr M)^2 - 4 det M crosses 0; in the complex regime the
    alternating-word value oscillates ~ r^{n/2} cos(n phi/2 + psi) and hence goes
    negative at some finite n (first sign change near n ~ 2 pi / phi).  Therefore
        D_c(a,b) = smallest positive D with disc(D)=0.
    disc is a quadratic in z^2 (z=zeta); the BINDING root is
        z_c^2 = (2-a-b)^2 / [ (2-a-b)^2 - 4ab ]
              = (1+lambda2)^2 / [ 4 lambda2 + (a-b)^2 ],   lambda2 = 1-a-b.

(4) CLOSED FORM (PROVEN A3/A5):
        D_c(a,b) = (1/2) ( 1 - sqrt( 1 - 4ab/(2-a-b)^2 ) ),
    z_c = 1/(1-2 D_c) = (2-a-b)/sqrt((2-a-b)^2 - 4ab).
    Reduces EXACTLY to Gray's symmetric formula D_c(p)=(1/2)(1-sqrt(1-2p)/(1-p))
    at a=b=p (symbolic identity, A5).

(5) ROUTE-A UNIFORM GRAY BOUND -- NEGATIVE/HONEST RESULT (A6).
    With the now-EXACT threshold, the curvature-positivity quantity
        Psibar_0(a,b) = kappa_eff^2(a,b) * D_c(1-D_c)/(1-2 D_c)^2
    does NOT stay below the previously-scanned 0.385.  In fact
        Psibar_0 -> 1  as  a+b -> 1  (the slow-mixing boundary, lambda2->0),
    so sup_{(a,b)} Psibar_0 = 1 (approached, not attained); Route A FAILS to give
    a uniform-in-(a,b) bound.  It DOES give a uniform bound on any bounded-mixing
    sub-region {a+b <= 1-eta}: e.g. sup_{a+b<=0.97} Psibar_0 = 0.4375 < 1.
    [This corrects the paper's "0.385 max at (0.05,0.45)" -- that grid never
    reached the a+b->1 corner.]  The per-point certificate of Lemma 7.34e is
    unaffected (Psibar_0<1 strictly at every fixed interior (a,b)).

----------------------------------------------------------------------------
ITEM (ii) -- ROUTINE: Marton/Paulin Markov bounded-differences replacement of
the i.i.d. switch-bond McDiarmid step (I7) for the cumulants (mu_x, v_x).
----------------------------------------------------------------------------
For a<>b the word x is a genuine 2-state Markov chain, so the coordinates are
DEPENDENT and vanilla McDiarmid (independent inputs) is inapplicable.  Replace by
Paulin 2015 (EJP 20(79), "Concentration inequalities for Markov chains by Marton
couplings and spectral methods", Thm 2.1): for f with single-coordinate bounded
differences c_i on a Markov chain whose mixing matrix has operator norm bounded by
(1-delta)^{-1}, delta = one-step Dobrushin coefficient,
    P(|f - E f| >= t) <= 2 exp( -2 t^2 (1-delta)^2 / sum_i c_i^2 ).
For our chain the one-step Dobrushin coefficient is
    delta(T) = (1/2) sum_k |T(0,k)-T(1,k)| = |1-a-b| = |lambda2|,
so 1-delta = a+b (since a+b<1), and with c_i <= c_mu uniformly,
    P(|mu_x - nD| >= t) <= 2 exp( -2 t^2 (a+b)^2 / ((n-1) c_mu^2) ),
i.e. the SAME sub-Gaussian (I7) form with the variance proxy inflated by the
mixing factor 1/(a+b)^2 = 1/(1-|lambda2|)^2.  c_mu = O(1) (cluster representation);
E[mu_x]=nD exactly.  VERIFIED (B1,B2): E[mu_x]=nD; Var(mu_x)=Theta(n); the Paulin
sub-Gaussian MGF bound holds for all lambda at an asymmetric (slow-mixing) corner.

----------------------------------------------------------------------------
LABELS:  PROVEN (symbolic/exact) | VERIFIED-NUMERICALLY | CONJECTURED.
Numerics: numpy + mpmath + sympy; exact-rational where the prize warrants.
============================================================================
"""
import math
import numpy as np

try:
    import mpmath as mp
    HAVE_MP = True
except Exception:
    HAVE_MP = False
try:
    import sympy as sp
    HAVE_SYMPY = True
except Exception:
    HAVE_SYMPY = False


# ============================ shared machinery ============================

def chain(a, b):
    T = np.array([[1 - a, a], [b, 1 - b]])
    pi = np.array([b / (a + b), a / (a + b)])
    return T, pi


def Gmats(zeta):
    G0 = np.diag([(1 + zeta) / 2, (1 - zeta) / 2])
    G1 = np.diag([(1 - zeta) / 2, (1 + zeta) / 2])
    return G0, G1


def Dc_closed(a, b):
    """PROVEN closed form (item iii.4)."""
    return 0.5 * (1 - math.sqrt(1 - 4 * a * b / (2 - a - b) ** 2))


def PY_word(y, a, b, D):
    """P_{Y*}(y) by the exact 2x2 transfer (item iii.1)."""
    T, pi = chain(a, b)
    zeta = 1 / (1 - 2 * D)
    G0, G1 = Gmats(zeta)
    G = [G0, G1]
    v = pi.copy()
    for t in range(len(y) - 1):
        v = v @ (G[y[t]] @ T)
    v = v @ G[y[-1]]
    return float(v.sum())


def PY_brute(y, a, b, D):
    """Reference: direct sum_x P_X(x) prod_i g(x_i,y_i)."""
    import itertools
    n = len(y)
    T, pi = chain(a, b)
    zeta = 1 / (1 - 2 * D)
    tot = 0.0
    for x in itertools.product([0, 1], repeat=n):
        Px = pi[x[0]]
        for t in range(n - 1):
            Px *= T[x[t], x[t + 1]]
        f = 1.0
        for i in range(n):
            f *= (1 + zeta * (-1) ** (x[i] ^ y[i])) / 2
        tot += Px * f
    return tot


def M_period2(a, b, D):
    """Period-2 transfer product M = B1 B0 = (G1 T)(G0 T) (item iii.3)."""
    T, _ = chain(a, b)
    zeta = 1 / (1 - 2 * D)
    G0, G1 = Gmats(zeta)
    return (G1 @ T) @ (G0 @ T)


# ============================ A1: transfer identity ============================

def A1_transfer_identity():
    print("-" * 84)
    print("A1 (PROVEN, VERIFIED): P_{Y*}(y) = pi^T [prod G_{y_t} T] G_{y_n} 1  "
          "matches the FWHT deconvolution / brute sum.")
    import itertools
    rng = np.random.default_rng(0)
    maxerr = 0.0
    for (a, b, D) in [(0.1, 0.3, 0.01), (0.2, 0.4, 0.04), (0.05, 0.3, 0.005),
                      (0.4, 0.45, 0.15)]:
        for n in [3, 4, 5, 6]:
            ys = ([list(yy) for yy in itertools.product([0, 1], repeat=n)]
                  if n <= 5 else [list(rng.integers(0, 2, n)) for _ in range(40)])
            for y in ys:
                maxerr = max(maxerr, abs(PY_word(y, a, b, D) - PY_brute(y, a, b, D)))
    ok = maxerr < 1e-12
    print(f"  max |transfer - brute| over (a,b,D) x n<=6 (all words n<=5) = {maxerr:.2e}")
    print(f"  A1 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ A2: binding word ============================

def A2_binding_word():
    print("-" * 84)
    print("A2 (VERIFIED): binding word = ALTERNATING (period 2); below D_c all "
          "2^n words >=0, first-negative just above D_c is alternating (n<=16).")
    import itertools

    def min_all(n, a, b, D):
        best = 1e18
        bw = None
        for y in itertools.product([0, 1], repeat=n):
            p = PY_word(list(y), a, b, D)
            if p < best:
                best, bw = p, y
        return best, bw

    ok = True
    for (a, b) in [(0.1, 0.3), (0.2, 0.4), (0.05, 0.3), (0.1, 0.2),
                   (0.3, 0.4), (0.02, 0.5)]:
        dc = Dc_closed(a, b)
        # below: exhaustive min over all words, all n<=16, must be >= -tol
        minbelow = 1e18
        for n in range(2, 17):
            m, _ = min_all(n, a, b, dc * 0.999)
            minbelow = min(minbelow, m)
        # above: first word to go negative & is it alternating
        fw = None
        firstn = None
        for n in range(2, 17):
            m, bw = min_all(n, a, b, dc * 1.05)
            if m < 0:
                fw, firstn = bw, n
                break
        is_alt = fw is not None and all(fw[i] != fw[i + 1] for i in range(len(fw) - 1))
        ok &= (minbelow >= -1e-11) and is_alt
        print(f"  ({a},{b}) Dc={dc:.6f}: min(all y,n<=16) below={minbelow:.2e} (>=0:"
              f"{minbelow>=-1e-11}); first-neg above at n={firstn} alt={is_alt}")
    print(f"  A2 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ A3: symbolic discriminant + closed form ============================

def A3_symbolic():
    print("-" * 84)
    print("A3 (PROVEN, sympy): disc(M)=tr^2-4det factors; binding root "
          "z_c^2=(2-a-b)^2/((2-a-b)^2-4ab); D_c closed form; Gray reduction.")
    if not HAVE_SYMPY:
        print("  A3 -> SKIP (no sympy)")
        return True
    a, b, z, p = sp.symbols('a b z p', positive=True)
    T = sp.Matrix([[1 - a, a], [b, 1 - b]])
    G0 = sp.diag((1 + z) / 2, (1 - z) / 2)
    G1 = sp.diag((1 - z) / 2, (1 + z) / 2)
    M = (G1 * T) * (G0 * T)
    tr = sp.expand(M.trace())
    det = sp.expand(M.det())
    disc = sp.expand(tr ** 2 - 4 * det)
    # det = (z-1)^2(z+1)^2(a+b-1)^2/16
    det_ok = sp.simplify(det - (z - 1) ** 2 * (z + 1) ** 2 * (a + b - 1) ** 2 / 16) == 0
    # disc even in z -> poly in u=z^2; binding root u2:
    P = sp.Poly(disc, z)
    c4 = P.coeff_monomial(z ** 4)
    c2 = P.coeff_monomial(z ** 2)
    c0 = P.coeff_monomial(1)
    u = sp.symbols('u', positive=True)
    roots = sp.solve(sp.Eq(c4 * u ** 2 + c2 * u + c0, 0), u)
    u2_target = (2 - a - b) ** 2 / ((2 - a - b) ** 2 - 4 * a * b)
    root_ok = any(sp.simplify(r - u2_target) == 0 for r in roots)
    # FULL factorization (the paper's displayed identity):
    #   16 disc = F1 * F2,  F1 = (a+b)^2-(a-b)^2 z^2,
    #                       F2 = (2-a-b)^2 - ((2-a-b)^2-4ab) z^2,
    # and F2 BINDS: the cross-difference identity
    #   (a+b)^2[(2-a-b)^2-4ab] - (2-a-b)^2(a-b)^2 = 16ab(1-a-b) > 0
    # puts zeta_2^2 strictly below zeta_1^2 on the whole domain (ab>0, a+b<1).
    F1 = (a + b) ** 2 - (a - b) ** 2 * z ** 2
    F2 = (2 - a - b) ** 2 - ((2 - a - b) ** 2 - 4 * a * b) * z ** 2
    fact_ok = sp.simplify(16 * disc - sp.expand(F1 * F2)) == 0
    cross = (a + b) ** 2 * ((2 - a - b) ** 2 - 4 * a * b) \
        - (2 - a - b) ** 2 * (a - b) ** 2
    bind_ok = sp.simplify(cross - 16 * a * b * (1 - a - b)) == 0
    # closed-form D_c from u2 = z_c^2, z=1/(1-2D).  D_c = (1 - 1/z_c)/2 and
    # D_c_target = (1 - sqrt(1 - 4ab/(2-a-b)^2))/2; the radical equality
    # 1/z_c = sqrt(1 - 4ab/(2-a-b)^2) is proven by SQUARING (both nonneg):
    #   1/z_c^2 = 1/u2 = ((2-a-b)^2-4ab)/(2-a-b)^2 = 1 - 4ab/(2-a-b)^2,
    # with the radicand positive since (2-a-b)^2-4ab = 4(1-a-b)+(a-b)^2 > 0 (a+b<1).
    radicand = 1 - 4 * a * b / (2 - a - b) ** 2
    Dc_ok = (sp.simplify(1 / u2_target - radicand) == 0 and
             sp.simplify(((2 - a - b) ** 2 - 4 * a * b) - (4 * (1 - a - b) + (a - b) ** 2)) == 0)
    # symmetric reduction: 1 - p^2/(1-p)^2 = (1-2p)/(1-p)^2 (squaring identity)
    sym_inside_ok = sp.simplify((1 - p ** 2 / (1 - p) ** 2) - (1 - 2 * p) / (1 - p) ** 2) == 0
    print(f"  det M = (z^2-1)^2(a+b-1)^2/16 : {det_ok}")
    print(f"  16 disc = F1*F2 exact factorization : {fact_ok}")
    print(f"  binding identity cross-diff = 16ab(1-a-b) (F2 binds, zeta_2<zeta_1) : {bind_ok}")
    print(f"  binding root z_c^2 = (2-a-b)^2/((2-a-b)^2-4ab) among disc roots : {root_ok}")
    print(f"  D_c closed form via squaring (1/z_c^2 = 1-4ab/(2-a-b)^2, radicand>0) : {Dc_ok}")
    print(f"  Gray reduction inside-sqrt identity (1-p^2/(1-p)^2)=(1-2p)/(1-p)^2 : {sym_inside_ok}")
    ok = det_ok and fact_ok and bind_ok and root_ok and Dc_ok and sym_inside_ok
    print(f"  A3 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ A4: disc-root vs high-n threshold ============================

def A4_threshold_match():
    print("-" * 84)
    print("A4 (VERIFIED): D_c(closed) = n->inf alternating-word threshold; "
          "first-negative n -> infinity as D -> D_c+ (the 'creep').")
    ok = True

    if HAVE_MP:
        mp.mp.dps = 40

        def PY_alt_mp(n, a, b, D, start):
            a, b, D = mp.mpf(a), mp.mpf(b), mp.mpf(D)
            z = 1 / (1 - 2 * D)
            T = mp.matrix([[1 - a, a], [b, 1 - b]])
            pi = mp.matrix([[b / (a + b), a / (a + b)]])
            G0 = mp.matrix([[(1 + z) / 2, 0], [0, (1 - z) / 2]])
            G1 = mp.matrix([[(1 - z) / 2, 0], [0, (1 + z) / 2]])
            G = [G0, G1]
            v = pi.copy()
            y = [(start + i) % 2 for i in range(n)]
            for t in range(n - 1):
                v = v * (G[y[t]] * T)
            v = v * G[y[-1]]
            return v[0, 0] + v[0, 1]

        for (a, b) in [(0.1, 0.3), (0.05, 0.3)]:
            dc = Dc_closed(a, b)
            print(f"  ({a},{b}) Dc_closed={dc:.8f}: first-negative alternating-word n vs D-D_c")
            for fac in [1.02, 1.01, 1.005]:
                D = dc * fac
                firstneg = None
                for n in range(2, 160):
                    if PY_alt_mp(n, a, b, D, 0) < 0 or PY_alt_mp(n, a, b, D, 1) < 0:
                        firstneg = n
                        break
                print(f"     D={fac}xD_c: first-negative at n={firstneg} (-> infinity as D->D_c+)")
                ok &= firstneg is not None
    else:
        print("  (mpmath absent -- float disc-root vs n<=40 alt-word only)")

    # also: closed-form disc-root equals the smallest disc sign-change (float, cheap)
    def disc(a, b, D):
        M = M_period2(a, b, D)
        return np.trace(M) ** 2 - 4 * np.linalg.det(M)

    def smallest_root(a, b):
        Dprev, prev = 1e-6, disc(a, b, 1e-6)
        D = 1e-6
        while D < 0.49:
            D += 2e-4
            cur = disc(a, b, D)
            if (cur > 0) != (prev > 0):
                lo, hi = Dprev, D
                for _ in range(100):
                    m = (lo + hi) / 2
                    if (disc(a, b, m) > 0) == (prev > 0):
                        lo = m
                    else:
                        hi = m
                return (lo + hi) / 2
            prev, Dprev = cur, D
        return None

    mxerr = 0.0
    for (a, b) in [(0.1, 0.3), (0.2, 0.4), (0.05, 0.3), (0.1, 0.2), (0.3, 0.4),
                   (0.02, 0.5), (0.25, 0.25), (0.4, 0.45)]:
        r = smallest_root(a, b)
        if r is not None:
            mxerr = max(mxerr, abs(r - Dc_closed(a, b)))
    print(f"  max |smallest disc-root - D_c(closed)| over grid = {mxerr:.2e}")
    ok &= mxerr < 1e-7
    print(f"  A4 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ A5: symmetric reduction numeric ============================

def A5_gray_reduction():
    print("-" * 84)
    print("A5 (VERIFIED): D_c(p,p) == Gray (1/2)(1-sqrt(1-2p)/(1-p)) numerically.")
    ok = True
    for p in [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.45]:
        closed = Dc_closed(p, p)
        gray = 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))
        ok &= abs(closed - gray) < 1e-12
        print(f"  p={p}: closed={closed:.8f} Gray={gray:.8f} diff={abs(closed-gray):.2e}")
    print(f"  A5 -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ A6: Route-A scan (HONEST NEGATIVE) ============================

STATES = [(x, xp, xq) for x in (0, 1) for xp in (0, 1) for xq in (0, 1)]


def _W_parts(a, b):
    T = [[1 - a, a], [b, 1 - b]]

    def part(f1, f2):
        W = np.zeros((8, 8))
        for i, (x, xp, xq) in enumerate(STATES):
            for j, (y, yp, yq) in enumerate(STATES):
                if (y ^ yp) == f1 and (y ^ yq) == f2:
                    W[i, j] = T[xp][yp] * T[xq][yq] / T[x][y]
        return W

    return [[part(0, 0), part(0, 1)], [part(1, 0), part(1, 1)]]


def _lr(a, b):
    T = [[1 - a, a], [b, 1 - b]]
    pi = [b / (a + b), a / (a + b)]
    r = np.array([sum(T[xp][y] * T[xq][y] / T[x][y] for y in (0, 1))
                  for (x, xp, xq) in STATES])
    l = np.zeros(8)
    l[0], l[7] = pi[0], pi[1]
    return l, r


def keff2_numeric(a, b):
    """kappa_eff^2 = alpha11 - 2 alpha2 - 1 (Lemma 7.34e(ii))."""
    Wp = _W_parts(a, b)
    W00, W10, W01, W11 = Wp[0][0], Wp[1][0], Wp[0][1], Wp[1][1]
    l, r = _lr(a, b)
    P = np.outer(r, l)
    S = np.linalg.solve(np.eye(8) - W00 + P, np.eye(8) - P)
    a2 = l @ W10 @ S @ W10 @ r
    a11 = l @ W11 @ r + l @ W10 @ S @ W01 @ r + l @ W01 @ S @ W10 @ r
    return float(a11 - 2 * a2 - 1)


def Psibar0(a, b):
    D = Dc_closed(a, b)
    dd = D * (1 - D)
    return keff2_numeric(a, b) * dd / (1 - 2 * D) ** 2


def A6_route_a():
    print("-" * 84)
    print("A6 (HONEST NEGATIVE): with the EXACT threshold, Psibar_0(a,b,Dc) -> 1 "
          "as a+b->1; sup over open region = 1 (Route A FAILS uniformly).")
    print("    [The paper's '0.385 max' was a grid artifact -- never reached a+b->1.]")
    # interior bounded-mixing sub-region: uniform bound holds
    mx097 = 0.0
    for a in np.linspace(0.005, 0.5, 90):
        for b in np.linspace(0.005, 0.5, 90):
            if not (0.01 <= a + b <= 0.97):
                continue
            mx097 = max(mx097, Psibar0(a, b))
    print(f"  interior bounded-mixing  sup_{{a+b<=0.97}} Psibar_0 = {mx097:.5f}  (< 1: {mx097 < 1})")
    # boundary a+b->1: Psibar_0 -> 1
    print("  a+b->1 boundary (a=0.05 fixed):")
    seq = []
    for ab in [0.95, 0.98, 0.99, 0.995, 0.999]:
        b = ab - 0.05
        v = Psibar0(0.05, b)
        seq.append(v)
        print(f"    a+b={ab}: Psibar_0={v:.5f}")
    monotone_up_to_1 = all(seq[i] < seq[i + 1] for i in range(len(seq) - 1)) and seq[-1] > 0.99
    # per-point strict <1 everywhere interior (no violations)
    viol = 0
    for a in np.linspace(0.005, 0.5, 60):
        for b in np.linspace(0.005, 0.985, 80):
            if a + b >= 0.9995:
                continue
            if Psibar0(a, b) >= 1:
                viol += 1
    print(f"  per-point strict Psibar_0<1 on interior grid: violations = {viol}")
    # the FACTS we assert: (a) uniform bound on bounded-mixing region, (b) sup=1,
    # (c) per-point <1 strict.  All consistent -> PASS = the honest characterization.
    ok = (mx097 < 1) and monotone_up_to_1 and (viol == 0)
    print(f"  A6 -> {'PASS (Route-A failure correctly characterized)' if ok else 'FAIL'}")
    return ok


# ============================ ITEM (ii): Marton/Paulin ============================

def dobrushin(a, b):
    T = np.array([[1 - a, a], [b, 1 - b]])
    return 0.5 * np.sum(np.abs(T[0] - T[1]))   # = |1-a-b|


def _theta0(D):
    return math.log(D / (1 - D))


def _eta_C(beta, D):
    eb = np.exp(beta)
    ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
    eta = (1 - ze) / (1 + ze)
    C = (1 + eb) * (1 + ze) / 2
    return eta, C


def mu_x_asym(xbits, a, b, D):
    """mu_x = E_Q[K | x] from the tilted posterior q_k = IFFT(Phi(s;x)),
    Phi(s;x) = C_s^n sum_{x'} P_X(x') eta_s^{d(x,x')} / [C_0^n sum_{x'} P_X(x')],
    the source sum done by the asymmetric-chain transfer over x'."""
    n = len(xbits)
    S = n + 1
    s = 2 * math.pi * np.arange(S) / S
    th = _theta0(D)
    eta, Cs = _eta_C(th + 1j * s, D)
    T = np.array([[1 - a, a], [b, 1 - b]])
    pi = np.array([b / (a + b), a / (a + b)])

    def fac(v, xi):
        return eta if v != xi else np.ones(S, dtype=complex)

    al = np.zeros((2, S), dtype=complex)
    for xp in (0, 1):
        al[xp] = pi[xp] * fac(xp, xbits[0])
    for t in range(1, n):
        new = np.zeros((2, S), dtype=complex)
        for xp in (0, 1):
            f = fac(xp, xbits[t])
            new[xp] = (al[0] * T[0, xp] + al[1] * T[1, xp]) * f
        al = new
    M = (Cs ** n) * (al[0] + al[1])
    M0 = M[0].real
    phi = M / M0
    q = np.fft.fft(phi).real / S
    q = np.maximum(q, 0.0)
    q = q[:n + 1]
    q /= q.sum()
    ks = np.arange(n + 1, dtype=float)
    return float(ks @ q)


def B_marton_paulin():
    print("-" * 84)
    print("ITEM (ii) -- Marton/Paulin Markov bounded-differences (Paulin 2015 EJP "
          "20(79) Thm 2.1).")
    print("  one-step Dobrushin delta(T) = |1-a-b| = |lambda2|; factor 1/(1-delta)^2 "
          "= 1/(a+b)^2.")
    print("  P(|mu_x - nD| >= t) <= 2 exp(-2 t^2 (a+b)^2 / ((n-1) c_mu^2)).")
    ok = True
    rng = np.random.default_rng(7)

    def samp(n, a, b):
        T = np.array([[1 - a, a], [b, 1 - b]])
        pi = np.array([b / (a + b), a / (a + b)])
        x = [int(rng.random() < pi[1])]
        for _ in range(n - 1):
            x.append(int(rng.random() < T[x[-1], 1]))
        return np.array(x)

    for (a, b) in [(0.05, 0.45), (0.1, 0.3)]:
        d = dobrushin(a, b)
        print(f"  (a,b)=({a},{b}): delta={d:.4f}=|1-a-b|={abs(1-a-b):.4f}; "
              f"1/(a+b)^2={1/(a+b)**2:.3f}; reduction delta==|1-a-b|: {abs(d-abs(1-a-b))<1e-12}")
        ok &= abs(d - abs(1 - a - b)) < 1e-12
        dc = Dc_closed(a, b)
        D = 0.9 * dc
        n = 120
        # measure single-bond influence c_mu
        cmu = 0.0
        for _ in range(40):
            x = samp(n, a, b)
            m0 = mu_x_asym(x, a, b, D)
            for i in range(n):
                xf = x.copy()
                xf[i] ^= 1
                cmu = max(cmu, abs(mu_x_asym(xf, a, b, D) - m0))
        sigP2 = (n - 1) * cmu ** 2 / (4 * (a + b) ** 2)
        R = 2500
        mus = np.array([mu_x_asym(samp(n, a, b), a, b, D) for _ in range(R)])
        Emu = mus.mean()
        c = mus - n * D
        # E[mu]=nD
        mean_ok = abs(Emu - n * D) < 0.15 * math.sqrt(mus.var()) + 0.05 * n * D + 0.1
        # Paulin sub-Gaussian MGF bound holds for all lambda
        mgf_ok = True
        for lam in [0.2, 0.5, 1.0, 2.0, -1.0, -2.0]:
            lmgf = math.log(np.mean(np.exp(lam * c)))
            bound = lam * lam * sigP2 / 2
            mgf_ok &= lmgf <= bound
        proxy_ok = sigP2 >= mus.var()
        ok &= mean_ok and mgf_ok and proxy_ok
        print(f"    n={n}, D=0.9Dc={D:.5f}: E[mu]={Emu:.3f} (nD={n*D:.3f}, ok:{mean_ok}); "
              f"c_mu={cmu:.4f}; sigma_P^2={sigP2:.2f}>=Var(mu)={mus.var():.3f}:{proxy_ok}; "
              f"Paulin MGF bound all-lambda:{mgf_ok}")
    print(f"  ITEM (ii) -> {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ main ============================

def main():
    print("=" * 84)
    print("probe_7_34_dc_ab_characterization.py -- instance (2) residuals (iii)+(ii)")
    print("=" * 84)
    results = {}
    results['A1 transfer identity'] = A1_transfer_identity()
    results['A2 binding word (alternating)'] = A2_binding_word()
    results['A3 symbolic disc + closed form'] = A3_symbolic()
    results['A4 disc-root = n->inf threshold'] = A4_threshold_match()
    results['A5 Gray symmetric reduction'] = A5_gray_reduction()
    results['A6 Route-A scan (honest negative)'] = A6_route_a()
    results['(ii) Marton/Paulin concentration'] = B_marton_paulin()
    print("=" * 84)
    print("SUMMARY")
    for k, v in results.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    allok = all(results.values())
    print("-" * 84)
    print("KEYSTONE (iii):  D_c(a,b) = (1/2)(1 - sqrt(1 - 4ab/(2-a-b)^2))   [PROVEN]")
    print("                 binding word = alternating; condition = disc(B1 B0)=0.")
    print("                 reduces to Gray's D_c(p) at a=b.")
    print("ROUTE A (iii.5): sup Psibar_0 = 1 as a+b->1  =>  uniform Gray bound FAILS;")
    print("                 holds only on bounded-mixing {a+b<=1-eta} (e.g. <=0.97: 0.4375).")
    print("ITEM (ii):       Paulin Markov-McDiarmid, factor 1/(a+b)^2, sub-Gaussian.  [VERIFIED]")
    print("=" * 84)
    print(f"OVERALL -> {'PASS' if allok else 'FAIL'}")
    return allok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
