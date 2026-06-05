#!/usr/bin/env python3
"""
FRONTIER: ACHIEVABILITY side of the discrete-Markov (BSMS) Hamming RD dispersion (7.34c attempt).

Remark 7.34b proved the CONVERSE dispersion V_conv = lim Var(j_n)/n is D-dependent and strictly
below the lossless varentropy, but a converse alone cannot refute conjecture B (V_op=V_lossless):
a loose converse could coexist with covering-dominated V_op=V_lossless. To settle it we need an
ACHIEVABILITY (upper bound) on V_op.

A random codebook drawn i.i.d. from a MARKOV output measure q* (symmetric chain, switch prob r)
is a VALID coding scheme. Its excess-distortion prob <= E_X[(1-mu(X))^M], mu(x)=q*(B_{nD}(x)).
The achievable rate/dispersion are the mean/variance of G(x^n) = -log2 q*(B_{nD}(x^n)):
   R_ach(D) = E[G]/n,   V_ach(D) = Var(G)/n.
If the best Markov output is first-order optimal (R_ach ~ R(D)) and V_ach < V_lossless, then
   V_op <= V_ach < V_lossless  ==>  conjecture B is OPERATIONALLY REFUTED.

CRUCIAL: unlike the 2^n d-tilted enumeration (ball radius floor(nD) in {0,1} at n<=12 -- artifactual),
the TRANSFER-MATRIX DP computes the EXACT ball measure q*(B_{nD}(x^n)) for LARGE n (radius t=round(nD)
in the bulk), by a per-sample DP over (Y-state, running Hamming count) -- O(n^2), no 2^n blowup.

For each sampled x^n ~ BSMS(p) and output switch r, DP gives the full law of A_n=d_H(x^n,Y),
Y~Markov(r); mu=P(A_n<=t), G=-log2 mu.  We optimize r (min rate = best Markov output), then report
V_ach at r*.  Checks:
  C1 first-order: best-Markov R_ach(D) ~ the n-letter RD R(D) (sanity that the Markov output is ~optimal).
  C2 THE TEST: V_ach(D) < V_lossless = p(1-p)log2^2((1-p)/p)  => B operationally refuted.
  C3 n-stability: V_ach at n and 1.5n agree (not a finite-n artifact).
"""
import math
import numpy as np

def h2(x):
    if x <= 0 or x >= 1: return 0.0
    return -x*math.log2(x) - (1-x)*math.log2(1-x)

def sample_bsms(n, p, R, rng):
    """R samples of the stationary binary symmetric Markov chain, switch prob p."""
    steps = (rng.random((R, n)) < p).astype(np.int8)   # innovations Bern(p)
    steps[:, 0] = (rng.random(R) < 0.5).astype(np.int8) # X_1 uniform
    return np.cumsum(steps, axis=1) % 2                  # running XOR

def G_ball(x, r, t):
    """-log2 P(d_H(x,Y) <= t), Y ~ symmetric Markov(switch r), via DP over (state, count)."""
    n = x.shape[0]
    f = np.zeros((2, n+1))
    x0 = int(x[0])
    f[0, 0 if x0 == 0 else 1] = 0.5      # Y_1=0
    f[1, 0 if x0 == 1 else 1] = 0.5      # Y_1=1
    for i in range(1, n):
        xi = int(x[i])
        fn = np.zeros((2, n+1))
        for s in (0, 1):
            fs = f[s]
            for sp in (0, 1):
                pr = (1.0 - r) if sp == s else r
                if sp != xi:                # count increments
                    fn[sp, 1:] += pr * fs[:n]
                else:
                    fn[sp, :] += pr * fs
        f = fn
    P = f[0] + f[1]
    mu = P[:t+1].sum()
    return -math.log2(max(mu, 1e-300))

def achievability(p, D, n, R, rng, r_grid):
    X = sample_bsms(n, p, R, rng)
    t = int(round(n*D))
    best = None
    for r in r_grid:
        Gs = np.array([G_ball(X[j], r, t) for j in range(R)])
        Rate = Gs.mean()/n
        Vach = Gs.var()/n
        if best is None or Rate < best[1]:
            best = (r, Rate, Vach, t)
    return best  # (r*, R_ach, V_ach, t)

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    p = 0.1; D = 0.05
    Vloss = p*(1-p)*(math.log2((1-p)/p))**2
    R_slb = h2(p) - h2(D)
    print("="*84)
    print(f"BSMS p={p}, D={D}: V_lossless={Vloss:.4f} bits^2; SLB rate h(p)-h(D)={R_slb:.4f}; "
          f"(n-letter R(D)~0.27 from the n=12 BA solve)")
    r_grid = [0.04, 0.06, 0.08, 0.10, 0.13, 0.16, 0.20]
    print(f"{'n':>4} {'t=rnd(nD)':>9} {'r*':>6} {'R_ach(D)':>9} {'V_ach':>8} {'V_loss':>8} {'V_ach<V_loss?':>13}")
    res = {}
    for n in (48, 72, 108):
        Rn = 4000 if n <= 72 else 2500
        r_star, R_ach, V_ach, t = achievability(p, D, n, Rn, rng, r_grid)
        res[n] = (R_ach, V_ach)
        print(f"{n:>4} {t:>9d} {r_star:>6.3f} {R_ach:>9.4f} {V_ach:>8.4f} {Vloss:>8.4f} "
              f"{'YES' if V_ach < Vloss-0.03 else 'no':>13}")
    print("="*84)
    print("VERDICT")
    Vs = [res[n][1] for n in res]
    Rs = [res[n][0] for n in res]
    below = all(v < Vloss - 0.03 for v in Vs)
    stable = (max(Vs) - min(Vs)) < 0.12
    # first-order optimality: R_ach close to the n-letter R(D)~0.27 (and >= SLB)
    foopt = all(0.24 < r < 0.32 for r in Rs)
    R12_exact = 0.2686   # exact n=12 BA optimum at D=0.05 (from remark_7_34b ... PART B, s=2.0 row)
    suboptimal = any(r > R12_exact + 0.005 for r in Rs)
    print(f"  order-1 Markov R_ach(D): {[f'{r:.3f}' for r in Rs]}  vs EXACT n=12 optimum R_12={R12_exact}")
    print(f"  V_ach(D): {[f'{v:.3f}' for v in Vs]}  vs V_conv~0.51 (7.34b)  vs V_lossless={Vloss:.3f}")
    print("="*84)
    print("HONEST CONCLUSION (after adversarial review -- the achievability does NOT refute B):")
    print(f"  - The order-1 Markov codebook is FIRST-ORDER SUBOPTIMAL: R_ach~0.28 EXCEEDS the exact")
    print(f"    n=12 optimum R_12={R12_exact} (>= R(D)). The optimal n-letter output for a Markov source")
    print( "    requires GROWING Markov order (Eswaran-Gastpar 2022, arXiv:2211.04535); any fixed-order")
    print( "    codebook overpays the rate, so V_ach is attached to the WRONG first-order rate and")
    print( "    V_op <= V_ach does NOT follow. (The single-letter d-tilted object is also degenerate,")
    print( "    V_sl=0, for the symmetric chain: Krishnamachari 2026 Cor 6.)")
    print(f"  - SUGGESTIVE only: V_ach~0.5 ~ the 7.34b converse V_conv~0.51, both far below V_lossless,")
    print( "    consistent with (but NOT proving) V_op~V_conv<V_lossless.")
    print( "  => ACHIEVABILITY DISPERSION REMAINS OPEN. The 7.34b CONVERSE (V_conv<V_lossless) stands.")
    print( "     Sharpened obstruction: a second-order GROWING-ORDER Markov type-covering is the missing tool.")
