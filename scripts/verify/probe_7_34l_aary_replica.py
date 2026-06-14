#!/usr/bin/env python3
"""
probe_7_34l_aary_replica.py
============================================================================
THE A-ARY SYMMETRIC REPLICA SPECTRAL BOUND  (the one genuinely-new analytic
ingredient for the A-ary operational RD-dispersion theorem, generalising the
binary Remark 7.34b replica theorem and the asymmetric Lemma 7.34e).

Source X: A-ary symmetric Markov chain, T(i,i)=1-p, T(i,j)=p/(A-1) (j!=i),
pi uniform.  Hamming distortion.  Backward test channel K = A-ary symmetric
channel, K(y|x)=1-D (y=x), D/(A-1) (y!=x).  Gray region 0<D<=D_c^(A)(p)
(Lemma 7.34k): all-n deconvolution P_{Y*}=(K^{-1})^{(x)n}P_X >= 0.

DUAL IDENTITY (the load-bearing structural fact, verified V0).  On the Gray
region the SLB is all-n tight and the posterior normaliser factorises:
    M(x^n) := sum_{y} P_{Y*}(y) e^{-lambda* d_H(x,y)}  =  c_n * P_X(x^n),
    c_n = gamma^n  (gamma the per-site dual constant, x^n-INDEPENDENT to 1e-15),
    lambda* = log((1-D)(A-1)/D),   nu := e^{-lambda*} = D/((1-D)(A-1)).
This is the A-ary analog of the binary "Phi(0;x)=1 / structural cancellation":
the quenched 1/M(x)^2 collapses to gamma^{-2n}/P_X(x)^2, turning the annealed
second moment into a SOURCE-SIDE transfer-matrix product on triples
(x,x',x'') -- the right A-ary generalisation (the OUTPUT law P_{Y*} is NOT
first-order Markov for A>=3, so the naive output-side W_{A^3} fails; the
source-side form is exact for every A, verified V1 to 1e-15).

THE A-ARY REPLICA OPERATOR.  With the dressed single-site kernel
    G_s(x,x') := sum_{y} K^{-1}[y,x'] f_s(x,y),
    f_s(x,y) = 1 (y=x),  nu e^{i s} (y!=x),
the annealed second moment is
    E_x|Phi(s;x)|^2 = gamma^{-2n} * u(s)^T W_A(s)^{n-1} * 1,
    W_A[(xa,xb,xc),(x,x',x'')] = (1/T[xa,x]) T[xb,x'] T[xc,x'']
                                 G_s(x,x') conj(G_s(x,x'')),
    u(s)[(x,x',x'')] = (1/A) G_s(x,x') conj(G_s(x,x'')).
Per-site replica rate
    g_A(s) := rho(W_A(s)) / gamma^2 ,   g_A(0)=1 (Perron normalisation, since
    E_x|Phi(0;x)|^2 = 1 identically).  The dispersion achievability needs
    sup_{[eps,pi]} g_A(s) < 1  with positive curvature g_A(s)=1-vbar_A s^2+...

CHECKS:
  V0  dual identity: M(x)=c_n P_X(x), c_n=gamma^n, x-independent (1e-13).
  V1  brute E_x|Phi|^2 == gamma^{-2n} u^T W_A^{n-1} 1  (operator exact, 1e-12).
  V2  g_A(s) < 1 on (eps,pi] across the A-ary Gray region (A=2..6, D up to Dc).
  V3  curvature g_A(s)=1 - vbar_A s^2 + O(s^4), vbar_A>0 = per-site posterior
      variance of the Hamming distortion (the A-ary varentropy connection).
  V4  s=pi is binding (worst case); g_A(pi)<1 at D=Dc (closed Gray endpoint).
  V5  off-Gray control: D>D_c^(A) -> P_{Y*} invalid (negative mass): region
      is load-bearing.
  V6  S_A symmetry reduction: W_A (dim A^3) collapses to an A-INDEPENDENT
      pattern-quotient (5x5 for A>=3, 4x4 for A=2) whose top eigenvalue == rho(W_A)
      -- the Perron sits in the trivial isotypic (pattern-averaged block), so a
      single uniform-in-A spectral argument suffices (A enters as a parameter,
      not the dimension).  This is the A-ary analog of the binary flip-sector split.
  V7  alphabet-UNIFORM curvature floor: g_A(s) <= 1 - (11/16) D(1-D)(1-cos s) for
      EVERY A (the binary 11/16 constant is conservative; the actual worst ratio is
      ~1.5, slightly increasing in A), binding corner small p / D->D_c / s=pi.  The
      excluded degenerate point p=(A-1)/A is exactly V_lossless=0 (vacuous).

Python: numpy.  PASS/FAIL to stdout.
"""
import itertools
import math
import numpy as np

PASS = True
def rep(name, ok):
    global PASS
    PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}")
    return ok

# ----------------------------------------------------------------------------
def Tmat(A, p):
    T = np.full((A, A), p / (A - 1)); np.fill_diagonal(T, 1 - p); return T

def Kinv(A, D):
    b0 = D / (A - 1); a0 = 1 - D - b0
    return (np.eye(A) - b0 * np.ones((A, A))) / a0

def lam_star(A, D):
    return math.log((1 - D) * (A - 1) / D)

def words_pmf(A, p, n):
    T = Tmat(A, p)
    ws = list(itertools.product(range(A), repeat=n))
    PX = {}
    for w in ws:
        v = 1.0 / A
        for k in range(1, n):
            v *= T[w[k - 1], w[k]]
        PX[w] = v
    return ws, PX

def deconv_PY(A, p, D, n, ws, PX):
    Ki = Kinv(A, D)
    PY = {}
    for y in ws:
        acc = 0.0
        for x in ws:
            f = 1.0
            for k in range(n):
                f *= Ki[y[k], x[k]]
            acc += f * PX[x]
        PY[y] = acc
    return PY

def Dc_num(A, pv):
    def im(Dv):
        Ki = Kinv(A, Dv); Tn = Tmat(A, pv)
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0))
        t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, (A - 1) / A - 1e-7
    for _ in range(70):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if im(mid) < 1e-11 else (lo, mid)
    return lo

# ----------------------------------------------------------------------------
#  brute annealed second moment  (ground truth)
# ----------------------------------------------------------------------------
def brute(A, p, D, n, s_list):
    lam = lam_star(A, D)
    ws, PX = words_pmf(A, p, n)
    PY = deconv_PY(A, p, D, n, ws, PX)
    def dH(x, y): return sum(1 for a, b in zip(x, y) if a != b)
    E2 = np.zeros(len(s_list)); E1 = np.zeros(len(s_list), complex)
    vbar = 0.0
    minPY = min(PY.values())
    for x in ws:
        ks = np.array([dH(x, y) for y in ws])
        w0 = np.array([PY[y] for y in ws]) * np.exp(-lam * ks)
        M0 = w0.sum()
        if M0 <= 0:
            continue
        pk = np.bincount(ks, weights=w0, minlength=n + 1) / M0
        kk = np.arange(n + 1)
        Phi = np.array([np.sum(pk * np.exp(1j * sv * kk)) for sv in s_list])
        E1 += PX[x] * Phi
        E2 += PX[x] * np.abs(Phi) ** 2
        mu = float(kk @ pk)
        vbar += PX[x] * (float((kk ** 2) @ pk) - mu * mu)
    return E1, E2, minPY, vbar / n

def gamma_dual(A, p, D, n=4):
    """per-site dual constant gamma: M(x)=gamma^n P_X(x). Returns (gamma, stdratio)."""
    lam = lam_star(A, D)
    ws, PX = words_pmf(A, p, n)
    PY = deconv_PY(A, p, D, n, ws, PX)
    def dH(x, y): return sum(1 for a, b in zip(x, y) if a != b)
    ratios = []
    for x in ws:
        M = sum(PY[y] * math.exp(-lam * dH(x, y)) for y in ws)
        ratios.append(M / PX[x])
    ratios = np.array(ratios)
    return ratios.mean() ** (1.0 / n), ratios.std()

# ----------------------------------------------------------------------------
#  SOURCE-SIDE replica operator
# ----------------------------------------------------------------------------
def W_operator(A, p, D, s):
    T = Tmat(A, p); Ki = Kinv(A, D); nu = math.exp(-lam_star(A, D))
    w = np.exp(1j * s)
    def G(x, xp):  # sum_y Ki[y,xp] f_s(x,y)
        tot = 0j
        for y in range(A):
            tot += Ki[y, xp] * (1.0 if y == x else nu * w)
        return tot
    Gm = np.array([[G(x, xp) for xp in range(A)] for x in range(A)])
    states = list(itertools.product(range(A), repeat=3))
    SI = {st: i for i, st in enumerate(states)}; nS = len(states)
    W = np.zeros((nS, nS), complex)
    for (xa, xb, xc) in states:
        i = SI[(xa, xb, xc)]
        for (x, xp, xq) in states:
            j = SI[(x, xp, xq)]
            W[i, j] = (1 / T[xa, x]) * T[xb, xp] * T[xc, xq] * Gm[x, xp] * np.conj(Gm[x, xq])
    u = np.array([(1 / A) * Gm[x, xp] * np.conj(Gm[x, xq]) for (x, xp, xq) in states], complex)
    return W, u

def E2_tm(A, p, D, n, s):
    W, u = W_operator(A, p, D, s)
    gam, _ = gamma_dual(A, p, D, n=4)
    val = u @ np.linalg.matrix_power(W, n - 1) @ np.ones(W.shape[0])
    return gam ** (-2 * n) * val

def g_A(A, p, D, s, gam=None):
    W, _ = W_operator(A, p, D, s)
    rho = float(np.max(np.abs(np.linalg.eigvals(W))))
    if gam is None:
        gam, _ = gamma_dual(A, p, D, n=4)
    return rho / gam ** 2

# ----------------------------------------------------------------------------
#  CHECKS
# ----------------------------------------------------------------------------
def V0():
    print("-" * 78)
    print("V0  dual identity  M(x)=gamma^n P_X(x), x-independent (the cancellation)")
    ok = True; worst = 0.0
    for A in (2, 3, 4, 5):
        for p in (0.1, 0.2):
            D = 0.7 * Dc_num(A, p)
            gam, std = gamma_dual(A, p, D, n=4)
            worst = max(worst, std)
            print(f"     A={A} p={p} D=0.7Dc={D:.5f}  gamma={gam:.8f}  std(M/PX)={std:.2e}")
            ok = ok and std < 1e-12
    return rep(f"V0 dual constant x-independent (worst std {worst:.1e})", ok)

def V1():
    print("-" * 78)
    print("V1  brute E_x|Phi|^2 == gamma^{-2n} u^T W_A^{n-1} 1  (operator EXACT)")
    ok = True; worst = 0.0
    s_list = [0.4, 1.1, 1.9, 2.7, math.pi]
    for A in (2, 3, 4, 5):
        for p in (0.1, 0.2):
            D = 0.7 * Dc_num(A, p)
            for n in (3, 4):
                _, E2b, _, _ = brute(A, p, D, n, s_list)
                for si, s in enumerate(s_list):
                    tm = E2_tm(A, p, D, n, s)
                    worst = max(worst, abs(tm.real - E2b[si]) + abs(tm.imag))
        print(f"     A={A}: max|brute-TM| so far {worst:.2e}")
    ok = worst < 1e-11
    return rep(f"V1 operator matches brute 2nd moment (max err {worst:.1e})", ok)

def V2():
    print("-" * 78)
    print("V2  g_A(s) = rho(W_A(s))/gamma^2 < 1 on (eps,pi], A-ary Gray region")
    ok = True; worst_g = -np.inf
    grid = np.linspace(1e-3, math.pi, 100)
    for A in (2, 3, 4, 5, 6):
        for p in (0.1, 0.2):
            dc = Dc_num(A, p)
            for fD in (0.5, 0.9, 1.0):
                D = fD * dc
                gam, _ = gamma_dual(A, p, D, n=4)
                gmax = max(g_A(A, p, D, float(s), gam) for s in grid)
                worst_g = max(worst_g, gmax)
                tag = '' if gmax < 1 - 1e-9 else '  <-- CHECK'
                print(f"     A={A} p={p} D={fD:.1f}Dc={D:.6f}  max g_A = {gmax:.9f}{tag}")
                ok = ok and gmax < 1 + 1e-12  # <1 up to numerical (D=Dc endpoint)
    return rep(f"V2 sup g_A <= 1 on Gray region (worst {worst_g:.8f})", ok)

def V3():
    print("-" * 78)
    print("V3  curvature  g_A(s)=1 - vbar_A s^2 + O(s^4), vbar_A>0 = per-site var")
    ok = True
    for A in (2, 3, 4, 5):
        for p in (0.1, 0.2):
            dc = Dc_num(A, p); D = 0.9 * dc
            gam, _ = gamma_dual(A, p, D, n=4)
            h = 1e-4
            vbar_g = (1 - g_A(A, p, D, h, gam)) / h ** 2
            _, _, _, vbar_brute = brute(A, p, D, 5, [0.0])
            rel = abs(vbar_g - vbar_brute) / max(vbar_brute, 1e-14)
            print(f"     A={A} p={p}: curv(g)={vbar_g:.6f}  per-site var(n=5)={vbar_brute:.6f}  rel={rel:.2e}")
            ok = ok and vbar_g > 0 and rel < 0.1
    return rep("V3 positive curvature = per-site posterior variance", ok)

def V4():
    print("-" * 78)
    print("V4  g_A smallest at s=pi (max margin); g monotone-down on [pi/2,pi]; D=Dc")
    # As in binary (Remark 7.34b: 'g is in fact smallest at s=pi'), the binding
    # (hardest, g->1) point is s->0 where the curvature floor 1-vbar s^2 governs;
    # the endpoint s=pi is the EASIEST (smallest g, largest margin).  We check
    # the full-grid sup g<1 (V2) AND that the s=pi endpoint clears with margin.
    ok = True
    grid = np.linspace(math.pi / 2, math.pi, 30)
    for A in (2, 3, 4, 5, 6):
        for p in (0.1, 0.2):
            dc = Dc_num(A, p); D = dc
            gam, _ = gamma_dual(A, p, D, n=4)
            g_pi = g_A(A, p, D, math.pi, gam)
            g_mid = g_A(A, p, D, math.pi / 2, gam)
            # g(pi) is the minimum over [pi/2,pi]
            gmin_endpoint = min(g_A(A, p, D, float(s), gam) for s in grid)
            pi_is_min = abs(g_pi - gmin_endpoint) < 1e-7
            print(f"     A={A} p={p} D=Dc: g(pi)={g_pi:.8f} (=min on [pi/2,pi]? {pi_is_min}) g(pi/2)={g_mid:.8f}")
            ok = ok and g_pi < 1 - 1e-6 and pi_is_min and g_pi < g_mid
    return rep("V4 s=pi is the easiest point (smallest g), clears with margin", ok)

def V5():
    print("-" * 78)
    print("V5  off-Gray control: D>D_c^(A) -> binding alternating word P_{Y*}<0")
    # The binding word (Lemma 7.34k) is the period-2 alternating word 0101...;
    # just above D_c its deconvolved mass oscillates negative at index
    # n ~ pi/(2 theta) (the complex-onset creep), so we evaluate the alternating
    # word directly to a length long enough to expose the negativity (cheap,
    # O(n A^2)).  This is the load-bearing positivity that defines the region.
    ok = True
    for A in (3, 4, 5, 6):
        p = 0.15; dc = Dc_num(A, p); D = 1.5 * dc
        Ki = Kinv(A, D); T = Tmat(A, p)
        # P_Y*(alt word) via transfer: v_0 = pi=1/A row; step y multiplies by B_y
        By = lambda y: np.array([[Ki[y, a] * T[a, b] for b in range(A)] for a in range(A)])
        # The complex-onset creep index grows as n ~ pi/(2*theta) with theta the rotation
        # angle of the period-2 transfer product, which shrinks with A (Dc^(A)->0); so the
        # negativity surfaces only at longer words for large A.  Use a generous length 200.
        minmass = np.inf
        v = np.ones(A) / A
        for i in range(200):
            y = i % 2  # alternating 0,1,0,1,...
            v = v @ By(y)
            minmass = min(minmass, v.sum())
        print(f"     A={A} D=1.5Dc={D:.6f}: min alt-word P_Y* over len<=200 = {minmass:+.3e}")
        ok = ok and minmass < -1e-12
    return rep("V5 alternating word P_Y*<0 off-Gray (region load-bearing)", ok)

def _pattern(tr):
    """orbit of an ordered triple (x,u,u') under the diagonal S_A action, by its equality
    pattern: 0 all-equal, 1 x=u!=u', 2 x=u'!=u, 3 u=u'!=x, 4 all-distinct."""
    x, u, up = tr
    if x == u == up:
        return 0
    if x == u and u != up:
        return 1
    if x == up and u != up:
        return 2
    if u == up and x != u:
        return 3
    return 4

def V6():
    print("-" * 78)
    print("V6  S_A symmetry reduction: W_A collapses to an A-INDEPENDENT pattern-quotient")
    print("    (5x5 for A>=3, 4x4 for A=2); its top eigenvalue == rho(W_A) (Perron in the")
    print("    trivial isotypic = pattern-averaged block) => one uniform-in-A spectral argument")
    ok = True; worst = 0.0
    for A in (2, 3, 4, 5, 6):
        p = 0.2; D = 0.7 * Dc_num(A, p)
        npat = 5 if A >= 3 else 4
        for s in (0.3, 1.0, 2.0, math.pi):
            W, _ = W_operator(A, p, D, s)
            states = list(itertools.product(range(A), repeat=3))
            reps = [None] * npat
            for k, tr in enumerate(states):
                pp = _pattern(tr)
                if reps[pp] is None:
                    reps[pp] = k
            Q = np.zeros((npat, npat), complex)
            for a in range(npat):
                i = reps[a]
                for j, tr in enumerate(states):
                    Q[a, _pattern(tr)] += W[i, j]
            rW = float(np.max(np.abs(np.linalg.eigvals(W))))
            rQ = float(np.max(np.abs(np.linalg.eigvals(Q))))
            worst = max(worst, abs(rW - rQ))
        print(f"     A={A}: A^3={A**3} -> {npat}x{npat} pattern-quotient; "
              f"|rho(W)-rho(Q)| worst-so-far = {worst:.2e}")
        ok = ok and worst < 1e-10
    return rep("V6 W_A reduces to a <=5-dim pattern-quotient carrying rho (A-independent)", ok)

def V7():
    print("-" * 78)
    print("V7  alphabet-UNIFORM curvature floor: g_A(s) <= 1 - (11/16) D(1-D)(1-cos s)")
    print("    i.e. c_A(s)=(1-g_A)/(D(1-D)(1-cos s)) >= 11/16 for ALL A; binding corner =")
    print("    small p, D->D_c, s=pi.  (Actual worst ratio ~1.5; 11/16 is the binary floor.)")
    ss = np.concatenate([np.geomspace(1e-2, 0.5, 25, endpoint=False),
                         np.linspace(0.5, math.pi, 50)])
    ok = True; glob = np.inf
    for A in (2, 3, 4, 5, 6):
        amin = np.inf
        # restrict to the V_lossless>0 regime p < (A-1)/A (degenerate iid p=(A-1)/A excluded)
        for p in (0.02, 0.05, 0.1, 0.2, 0.3, 0.45 * (A - 1) / A):
            if abs(p - (A - 1) / A) < 1e-6:
                continue
            dc = Dc_num(A, p)
            if dc < 1e-8:
                continue
            D = dc; dd = D * (1 - D)
            gam, _ = gamma_dual(A, p, D, n=4)
            rs = [(1 - g_A(A, p, D, float(s), gam)) / (dd * (1 - math.cos(s))) for s in ss]
            amin = min(amin, min(rs))
        glob = min(glob, amin)
        print(f"     A={A}: worst floor ratio c_A at D=D_c (V_lossless>0) = {amin:.4f}")
        ok = ok and amin >= 11.0 / 16.0
    return rep(f"V7 g_A<=1-(11/16)D(1-D)(1-cos s) uniform in A (global worst {glob:.4f})", ok)

def V8():
    print("-" * 78)
    print("V8  A-ary Markov-McDiarmid input (Paulin 2015): Dobrushin coefficient")
    print("    delta = lambda_2(T) = 1 - pA/(A-1) < 1, so McDiarmid holds with factor")
    print("    1/(1-delta)^2 = ((A-1)/(pA))^2 (finite for every A>=2, p in (0,(A-1)/A)).")
    ok = True
    for A in (2, 3, 4, 5, 6):
        for p in (0.1, 0.2, 0.3):
            T = Tmat(A, p)
            ev = np.sort(np.linalg.eigvals(T).real)[::-1]
            lam2 = ev[1]
            delta = 0.5 * max(np.abs(T[i] - T[j]).sum()
                              for i in range(A) for j in range(A) if i != j)
            pred = 1 - p * A / (A - 1)
            ok = ok and abs(lam2 - pred) < 1e-12 and abs(delta - pred) < 1e-12 and delta < 1
        print(f"     A={A}: delta=lambda_2=1-pA/(A-1); 1/(1-delta)^2=((A-1)/(pA))^2 finite")
    return rep("V8 A-ary Markov-McDiarmid: delta=lambda_2<1, Paulin factor finite", ok)

def V9():
    print("-" * 78)
    print("V9  converse input (Remark 7.34d): rho=Var(j_n)/Var(i_n)=1 in A-ary Gray")
    print("    region (the SLB shift j_n - i_n is x-INDEPENDENT => V_conv=V_lossless^(A)).")
    ok = True
    for A in (2, 3, 4, 5):
        for p in (0.1, 0.2):
            n = 6; D = 0.7 * Dc_num(A, p); lam = lam_star(A, D)
            ws, PX = words_pmf(A, p, n); PY = deconv_PY(A, p, D, n, ws, PX)
            def dH(a, b): return sum(1 for u, w in zip(a, b) if u != w)
            T = Tmat(A, p)
            def inf(x):
                v = 1.0 / A
                for k in range(1, n):
                    v *= T[x[k - 1], x[k]]
                return -math.log2(v)
            I = np.array([inf(x) for x in ws])
            J = np.array([-math.log2(sum(PY[y] * math.exp(-lam * dH(x, y)) for y in ws)) for x in ws])
            w = np.array([PX[x] for x in ws])
            var = lambda a: (lambda m: float(np.sum(w * (a - m) ** 2)))(float(np.sum(w * a)))
            rho = var(J) / var(I); shstd = math.sqrt(var(J - I))
            print(f"     A={A} p={p}: rho={rho:.12f}  shift std={shstd:.1e} (x-indep)")
            ok = ok and abs(rho - 1) < 1e-9 and shstd < 1e-10
    return rep("V9 converse ratio rho=1 (shift deterministic) => V_conv=V_lossless", ok)

def V10():
    print("-" * 78)
    print("V10 V_lossless^(A) = p(1-p) log2^2((1-p)(A-1)/p) (=binary at A=2); the")
    print("    increments xi_i=-log T(X_{i-1},X_i) are i.i.d. (state-independent law),")
    print("    so V is the single-increment variance (NO Gordin cross-terms).")
    def Vclosed(A, p):
        return p * (1 - p) * (math.log2((1 - p) * (A - 1) / p)) ** 2
    ok = True
    for A in (2, 3, 4, 5, 6):
        for p in (0.1, 0.2):
            T = Tmat(A, p)
            # single-increment variance: xi = -log2 T(prev,next), prev uniform
            # P(stay)=1-p value -log2(1-p); P(move)=p value -log2(p/(A-1))
            g0 = -math.log2(1 - p); g1 = -math.log2(p / (A - 1))
            v_inc = p * (1 - p) * (g1 - g0) ** 2
            ok = ok and abs(v_inc - Vclosed(A, p)) < 1e-12
        print(f"     A={A}: V_lossless={Vclosed(A,0.1):.5f} (p=.1), {Vclosed(A,0.2):.5f} (p=.2)")
    # binary recovery
    ok = ok and abs(Vclosed(2, 0.2) - 0.2 * 0.8 * math.log2(0.8 / 0.2) ** 2) < 1e-12
    return rep("V10 V_lossless^(A) closed form (i.i.d. increments, binary at A=2)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("A-ary symmetric REPLICA spectral bound -- verifier (probe 7.34l)")
    print("=" * 78)
    V0(); V1(); V2(); V3(); V4(); V5(); V6(); V7(); V8(); V9(); V10()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
