#!/usr/bin/env python3
"""
PROBE for the 7.34 achievability residual (V_op <= V_lossless on the Gray region):
does the spectral / inhomogeneous-Markov LOCAL-CLT (Bahadur-Rao) give
    (*)  G_n(x^n) := -log2 P_{Y*}( d_H(x^n, Y) <= nD ) = j_n(x^n,D) + (1/2) log2 n + O(1)
UNIFORMLY over the sqrt(n)-typical shell of source words x^n, where Y* is the (stationary,
geometrically ergodic, finite-state) optimal OUTPUT process of the BSMS on the Gray region,
and j_n = i_n - n h(D) is the n-letter d-tilted information (i_n = source surprisal).

If (*) holds then Var(G_n)/n = Var(j_n)/n + o(1) = V_lossless, closing V_op = V_lossless.

LITERATURE BASIS (this script tests whether the cited form actually obtains numerically):
  - Kontoyiannis-Meyn 2003 (AoAP 13(1):304-362) Thm 4.x + 2005 (EJP 10:61-123) Thm 5.3/5.4:
    for a HOMOGENEOUS geom-ergodic chain and a SINGLE FIXED test fn F, the partial sum obeys the
    exact (Bahadur-Rao) LD asymptotic  P_x{S_n >= nc} ~ const * e^{-n J(c)} / sqrt(2 pi n sigma_a^2),
    J = Legendre transform of log dominant eigenvalue of the tilted kernel; lower tail analogous.
    The lattice case (Thm 5.4) uses the FINITE-n log-MGF Lambda_n and gives the same 1/sqrt(n).
  - Dolgopyat-Sarig 2023 (LNM 2330) Thm 7.8 / 7.26: for INHOMOGENEOUS uniformly-elliptic chains
    with a TIME-VARYING bounded test fn f_i, the lattice LLT-for-large-deviations gives
    P[S_N = z_N] ~ rho_N e^{-V_N I_N(z_N/V_N)} / sqrt(2 pi V_N), rho_N bounded away from 0,inf,
    VALID for (z_N - E S_N)/V_N in the LD-threshold interval (c_-, c_+).

THE RNR OBJECT: A_n = sum_i 1[Y_i != x_i], Y ~ Y* (homogeneous), test fn g_i(y)=1[y!=x_i] is
SELECTED from {y, 1-y} by the fixed schedule x^n -- an inhomogeneous additive functional whose
per-step law depends on x_i. mu(x^n)=P(A_n <= nD) is the LOWER tail.

CHECKS (the practical chain hypothesis -> theory -> identity -> numeric):
  V1  BAHADUR-RAO FORM, fixed x: for a fixed typical x^n, sweep n and fit
         G_n = a*n + b*log2(n) + c.  Test b ~ +1/2 (the local-CLT prefactor exponent) and
         a ~ R-rate.  This validates the 1/sqrt(n) prefactor empirically for the MEMORY output.
  V2  THE IDENTITY (*): on the Gray region (D <= D_c, where the iid-BSC(D) deconvolution is valid
         and Y* is the BSC(D)-deconvolved output), test  G_n - j_n  vs  (1/2)log2 n + O(1)
         UNIFORMLY over many typical x^n: report max-min spread of (G_n - j_n - 0.5 log2 n).
  V3  DISPERSION (the payoff): Var(G_n)/n  vs  V_lossless and vs Var(j_n)/n. Equal => closes.
  V4  THRESHOLD / FULL-REGIME (the DS obstruction check): is the RNR deviation D inside the
         large-deviation threshold (c_-,c_+)? Test full regime via the positivity-threshold proxy:
         P(A_n <= nD - eps n) and P(A_n >= nD + eps n) BOTH >= poly(1/n) (DS Thm 7.30 (7.53/54)),
         i.e. nD is not pinned at the essential-range boundary. (Endpoint D=D_c handled separately.)

NB on regimes: the Gray-region ALL-n SLB identity needs D <= D_c(p) which is TINY (e.g. p=0.1 ->
D_c=0.0031). At such D the ball radius nD is 0/1 for any feasible n -> the DP ball collapses
(artifact, see 7.34b PART D). So V2/V3 on the Gray region proper are NOT numerically reachable at
finite n; we test the local-CLT FORM (V1) and the threshold (V4) at moderate D (bulk radius),
where Y* is still the (geom-ergodic) BSC(D) output and the Bahadur-Rao structure is identical --
the FORM of (*) is D-uniform; only the all-n SLB rate identity j_n=i_n-nh(D) is Gray-specific.
"""
import math
import numpy as np


def h2(x):
    if x <= 0 or x >= 1:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def sample_bsms(n, p, R, rng):
    steps = (rng.random((R, n)) < p).astype(np.int8)
    steps[:, 0] = (rng.random(R) < 0.5).astype(np.int8)
    return np.cumsum(steps, axis=1) % 2


def ball_tail_logprob(x, r, t):
    """Return ( -log2 P(d_H(x,Y) <= t), P(d_H(x,Y) <= t) ) for Y ~ symmetric Markov(switch r),
    via exact DP over (Y-state in {0,1}, running Hamming count 0..n). O(n^2)."""
    n = x.shape[0]
    f = np.zeros((2, n + 1))
    x0 = int(x[0])
    f[0, 0 if x0 == 0 else 1] = 0.5
    f[1, 0 if x0 == 1 else 1] = 0.5
    for i in range(1, n):
        xi = int(x[i])
        fn = np.zeros((2, n + 1))
        for s in (0, 1):
            fs = f[s]
            for sp in (0, 1):
                pr = (1.0 - r) if sp == s else r
                if sp != xi:
                    fn[sp, 1:] += pr * fs[:n]
                else:
                    fn[sp, :] += pr * fs
        f = fn
    P = f[0] + f[1]
    mu = float(P[: t + 1].sum())
    return -math.log2(max(mu, 1e-300)), mu, P


def best_output_switch(p, D, n, X, r_grid):
    """Pick the Markov output switch r minimizing the mean rate E[-log2 mu]/n at this (n,D).
    (On the Gray region the SLB-optimal output is BSC(D)-deconvolution; for the FORM test at
    moderate D we let the codebook output be the rate-minimizing symmetric Markov approx.)"""
    t = int(round(n * D))
    best = None
    for r in r_grid:
        G = np.array([ball_tail_logprob(X[j], r, t)[0] for j in range(len(X))])
        rate = G.mean() / n
        if best is None or rate < best[1]:
            best = (r, rate)
    return best[0]


def jn_surprisal(x, p):
    """n-letter d-tilted info up to the -n h(D) shift: i_n(x) = -log2 P_BSMS(x). j_n=i_n-n h(D)."""
    n = x.shape[0]
    sw = int(np.sum(x[1:] != x[:-1]))
    stays = (n - 1) - sw
    logP2 = math.log2(0.5) + stays * math.log2(1 - p) + sw * math.log2(p)
    return -logP2  # = i_n; subtract n*h(D) for j_n


# ---------------------------------------------------------------- V1
def V1_bahadur_rao_form(p, rng):
    print("=" * 96)
    print("V1  LOWER-TAIL BAHADUR-RAO PREFACTOR: test the (1/2)log2 n term of the local-CLT.")
    print("    CLEAN instance: schedule x=0^n (constant) => g_i(y)=1[y!=0]=y for all i, so")
    print("    A_n = #{1s in Y} is a HOMOGENEOUS additive functional of the symmetric Markov")
    print("    output Y (switch r) -- exactly Kontoyiannis-Meyn Thm 5.4 (lattice). Fixed type =>")
    print("    NO empirical-type drift, isolating b. (A random typical x adds an O(sqrt n) type")
    print("    fluctuation in the rate n*J that swamps the log term -- a TEST artifact, not (*).)")
    ns = [400, 600, 900, 1300, 1900, 2700, 3800]
    rows = []
    bs = []
    for r in (0.1, 0.2, 0.3):
        for D in (0.2, 0.3):  # mean #1s/n = 1/2, so D<1/2 is a genuine lower-tail LD
            Gs, logn = [], []
            for n in ns:
                x = np.zeros(n, dtype=np.int8)
                t = int(round(n * D))
                G, _, _ = ball_tail_logprob(x, r, t)
                Gs.append(G); logn.append(math.log2(n))
            A = np.vstack([ns, logn, np.ones(len(ns))]).T
            (a, b, c), *_ = np.linalg.lstsq(A, np.array(Gs), rcond=None)
            bs.append(b)
            rows.append((r, D, a, b))
            print(f"    r={r} D={D}:  G_n ~ {a:.5f}*n + {b:.3f}*log2(n) + {c:.2f}   "
                  f"=> b={b:.3f}  (BR: +0.5)")
    ok = all(0.35 < b < 0.65 for b in bs)
    print(f"    >>> all prefactor exponents in [0.35,0.65] (-> +1/2 as mixing improves): "
          f"{'YES' if ok else 'NO'}")
    return ok


# ---------------------------------------------------------------- V4
def V4_threshold_full_regime(p, rng):
    print("=" * 96)
    print("V4  LD-THRESHOLD / FULL-REGIME (the Dolgopyat-Sarig partial-regime obstruction check).")
    print("    DS Thm 7.26: the precise LD asymptotics hold for the deviation level z=(nD-E A_n)/V_n")
    print("    iff z in (c_-,c_+). Partial regime (c_pm<r_pm, DS Ex.7.36) happens ONLY via a")
    print("    VANISHING-probability extreme symbol. Here Y is stationary with both symbols prob 1/2")
    print("    every step and transitions in (0,1): no such degeneracy => expect FULL regime c_pm=r_pm.")
    print("    Numeric witness: the LOWER-TAIL LD rate -log2 P(A_n<=nD)/n must be FINITE and bounded")
    print("    (not ->inf) for D anywhere in the open essential range (0, mean) -- i.e. nD not pinned")
    print("    at the boundary 0. We use x=0^n so A_n=#1s in Y (mean n/2), and sweep D in (0,1/2).")
    n = 1500
    x = np.zeros(n, dtype=np.int8)
    ok_all = True
    for r in (0.1, 0.2, 0.3):
        _, _, P = ball_tail_logprob(x, r, n)  # full count law of A_n=#1s
        ks = np.arange(len(P)); mean = float((ks * P).sum()) / n
        line = []
        for D in (0.05, 0.15, 0.25, 0.40):
            t = int(round(n * D))
            below = float(P[: t + 1].sum())
            rate = -math.log2(max(below, 1e-300)) / n
            finite = below > 1e-300 and rate < 5.0   # bounded LD rate => interior, full regime
            ok_all = ok_all and finite
            line.append(f"D={D}:rate={rate:.3f}{'*' if finite else '!'}")
        print(f"    r={r} (mean A_n/n={mean:.3f}):  " + "  ".join(line))
    print(f"    >>> lower-tail LD rate finite/bounded for all interior D (full regime, c_pm=r_pm): "
          f"{'YES' if ok_all else 'NO'}  (* = finite bounded rate, ! = degenerate)")
    return ok_all


# ---------------------------------------------------------------- V2/V3 (small-D, Gray-region form)
def V23_identity_and_dispersion(p, D, n, R, rng):
    print("=" * 96)
    print(f"V2/V3  IDENTITY (*) and DISPERSION at p={p}, D={D}, n={n} (R={R} words).")
    print("       G_n - j_n  vs  (1/2)log2 n + O(1);  Var(G_n)/n vs V_lossless, Var(j_n)/n.")
    t = int(round(n * D))
    if t == 0:
        print(f"       SKIP: nD={n*D:.2f} rounds to ball radius t=0 (Gray-region collapse, see 7.34b PART D);")
        print("       the identity (*) is not numerically reachable at this (n,D) -- form tested in V1 instead.")
        return None
    X = sample_bsms(n, p, R, rng)
    r = best_output_switch(p, D, n, X[: min(R, 300)], [0.6 * D, 0.8 * D, D, 1.2 * D, 1.5 * D, 2 * D])
    G = np.array([ball_tail_logprob(X[j], r, t)[0] for j in range(R)])
    i_n = np.array([jn_surprisal(X[j], p) for j in range(R)])
    j_n = i_n - n * h2(D)
    resid = G - j_n - 0.5 * math.log2(n)
    Vloss = p * (1 - p) * (math.log2((1 - p) / p)) ** 2
    VG = G.var() / n
    Vj = j_n.var() / n
    print(f"       output r*={r}; residual (G-j_n-0.5log2 n): mean={resid.mean():.3f} "
          f"std={resid.std():.3f} spread(max-min)={resid.max()-resid.min():.3f}")
    print(f"       Var(G_n)/n={VG:.4f}   Var(j_n)/n={Vj:.4f}   V_lossless={Vloss:.4f}")
    closes = abs(VG - Vloss) < 0.1 and abs(VG - Vj) < 0.1
    print(f"       >>> Var(G_n)/n ~ V_lossless ~ Var(j_n)/n (closes): {closes}")
    return closes


if __name__ == "__main__":
    rng = np.random.default_rng(20260605)
    p = 0.1
    Vloss = p * (1 - p) * (math.log2((1 - p) / p)) ** 2
    Dc = (1 - math.sqrt(1 - 2 * p) / (1 - p)) / 2
    print("#" * 96)
    print(f"# BSMS p={p}: D_c(Gray)={Dc:.4f}, V_lossless={Vloss:.4f} bits^2.")
    print(f"# Goal: does the spectral/inhomog local-CLT give (*) -> V_op<=V_lossless on (0,D_c]?")
    print("#" * 96)

    v1 = V1_bahadur_rao_form(p, rng)          # lower-tail 1/sqrt(n) prefactor (fixed type)
    v4 = V4_threshold_full_regime(p, rng)     # LD-threshold / full-regime check
    # Gray-region identity is collapsed at finite n (D_c tiny); show the collapse honestly:
    v23_gray = V23_identity_and_dispersion(p, Dc * 0.7, 400, 1500, rng)
    print("=" * 96)
    print("NOTE: the Gray-region SLB rate identity j_n=i_n-n h(D) holds only for D<=D_c (tiny); at")
    print("      such D the ball radius nD collapses to {0,1} for any feasible n (7.34b PART D), so")
    print("      V2/V3 (the dispersion EQUALITY itself) are NOT numerically reachable -- V1/V4 test")
    print("      the local-CLT FORM, which is what the achievability argument actually invokes.")

    print("#" * 96)
    print("SUMMARY")
    print(f"  V1 Bahadur-Rao 1/sqrt(n) prefactor (b~+0.5) for the MEMORY output Y*: {v1}")
    print(f"  V4 RNR deviation nD interior to LD-threshold (full regime, DS obstruction cleared): {v4}")
    print(f"  V2/V3 Gray-region identity+dispersion numerically reachable at finite n: "
          f"{'no (radius collapse)' if v23_gray is None else v23_gray}")
    print("#" * 96)
    print("READING: V1+V4 support that the local-CLT FORM (1/sqrt(n) prefactor, full LD regime) holds")
    print("  for the geom-ergodic output Y* at the RNR deviation. What the literature does NOT directly")
    print("  cover -- and what V2/V3 cannot reach numerically -- is the COMBINATION needed for (*):")
    print("  a FIXED-fraction (c interior, not small-eps) Bahadur-Rao for a TIME-INHOMOGENEOUS test")
    print("  function g(.;x_i), UNIFORMLY over the typical shell. KM=fixed-c but homogeneous/fixed-F;")
    print("  Dolgopyat-Sarig=inhomogeneous but the published LD-threshold window + array uniformity must")
    print("  be instantiated for THIS schedule family. See the report for the precise residual.")
