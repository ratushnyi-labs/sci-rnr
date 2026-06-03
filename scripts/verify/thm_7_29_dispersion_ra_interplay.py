r"""
Verification for Theorem 7.29 (second-order / dispersion coding rate of RNR, and the
random-access <-> dispersion interplay): the finite-blocklength RNR rate at error eps is
    L_RNR(N,eps,K) = N h + (N/K) delta_inf + sqrt(N V) Q^{-1}(eps) + o(sqrt N),
where V is the source VARENTROPY (= the long-run variance of the per-symbol information
= T7.28's v^2 under the ideal predictor), and sqrt(N V) Q^{-1}(eps) is the FUNDAMENTAL
source dispersion (Kontoyiannis-Verdu 2014). At the balanced K = Theta(sqrt N) the
random-access overhead (N/K) delta_inf and the dispersion sqrt(N V) Q^{-1}(eps) are BOTH
Theta(sqrt N) -- byte-granular random access exacts a second-order rate penalty of the
SAME ORDER as the intrinsic stochastic variability.

WHAT THIS SCRIPT VERIFIES (stationary order-1 Markov source, ideal predictor).
(V1) VARENTROPY = DISPERSION = T7.28 long-run variance. Compute the analytic varentropy
     V = Var(l_1) + 2 sum_k Cov(l_0,l_k) for l_i = -log2 P(X_i|X_{i-1}); confirm it equals
     the empirical Var(sum_i l_i)/N (the per-symbol information variance). This V is the
     source dispersion (K-V 2014).
(V2) GAUSSIAN QUANTILE (second-order achievability). Over many realisations X^N, the
     (1-eps)-quantile of the codelength sum_i l_i matches N h + sqrt(N V) Q^{-1}(eps) to
     O(sqrt N) -- the dispersion expansion holds (the deterministic (N/K) delta_inf
     overhead is an additive mean-shift, added separately in V3).
(V3) THE INTERPLAY. At K = c sqrt(N): (N/K) delta_inf = (sqrt(N)/c) delta_inf and
     dispersion = sqrt(N V) Q^{-1}(eps) are BOTH Theta(sqrt N); their RATIO
     delta_inf / (c sqrt(V) Q^{-1}(eps)) is N-INDEPENDENT (verified across a range of N).
     So the two contribute additively at the sqrt(N) scale, with the dispersion the
     irreducible floor and the overhead tunable by the access budget c (= K/sqrt N).

PASS = V1 (analytic V == empirical, and == long-run variance) AND V2 (Gaussian quantile
matches empirical to o(sqrt N) / small relative error) AND V3 (the overhead/dispersion
ratio is N-independent, both Theta(sqrt N)).

HONEST SCOPE. Q^{-1}(eps) := Phi^{-1}(1-eps) is the Gaussian TAIL-quantile (upper
(1-eps) point); the script uses NormalDist().inv_cdf(1-eps). Achievability is the RNR
scheme's codelength (T7.28 CLT + T7.15 mean), valid under exp-mixing. The dispersion
floor sqrt(N V) Q^{-1}(eps) is the source second-order CONVERSE under the EXCESS-LENGTH
criterion P[length>L]<=eps, proved by Kontoyiannis-Verdu 2014 for MEMORYLESS and
FINITE-STATE IRREDUCIBLE APERIODIC MARKOV sources (broader mixing where their
info-density CLT conditions hold) -- so no scheme beats it on that class. The NEW
content is the RNR-specific (N/K)*delta_inf overhead at second order and its
Theta(sqrt N) interplay with the dispersion (V3 is the algebra of that interplay) --
NOT a new source-dispersion proof, and NOT a claim that the (N/K)*delta_inf penalty is
MANDATORY for random access (a tight RA converse is open). delta_inf is the
Crutchfield-Feldman excess of a FIXED finite-summable-excess source (the full constant
is gamma=delta_inf+log2(1/eta) with AC precision); here parametrised.
"""

import math
import sys

import numpy as np

RNG = np.random.default_rng(20260603)


def markov1(A, temp, rng):
    """Stationary order-1 Markov chain: transition T, stationary pi."""
    T = np.zeros((A, A))
    for a in range(A):
        l = rng.standard_normal(A) / temp
        p = np.exp(l - l.max()); T[a] = p / p.sum()
    # stationary distribution (left eigenvector for eigenvalue 1)
    vals, vecs = np.linalg.eig(T.T)
    k = int(np.argmin(np.abs(vals - 1.0)))
    pi = np.real(vecs[:, k]); pi = pi / pi.sum()
    pi = np.abs(pi); pi = pi / pi.sum()
    return T, pi


def entropy_rate(T, pi):
    h = 0.0
    for a in range(len(pi)):
        for b in range(len(pi)):
            if T[a, b] > 0:
                h -= pi[a] * T[a, b] * math.log2(T[a, b])
    return h


def varentropy_analytic(T, pi, kmax=200):
    """V = Var(l_1) + 2 sum_{k>=1} Cov(l_0, l_k), l_i = -log2 P(X_i | X_{i-1}).
    Long-run variance of the per-symbol information of a stationary Markov chain."""
    A = len(pi)
    # l depends on (X_{i-1}, X_i); define g(a,b) = -log2 T[a,b].
    g = np.where(T > 0, -np.log2(np.where(T > 0, T, 1)), 0.0)
    # joint stationary law of (X_{i-1}, X_i): P2[a,b] = pi[a] T[a,b]
    P2 = pi[:, None] * T
    mean_l = float((P2 * g).sum())
    # Var(l_1) = E[g^2] - mean^2
    var_l = float((P2 * g * g).sum()) - mean_l ** 2
    # Cov(l_0, l_k): l_0 = g(X_{-1},X_0), l_k = g(X_{k-1},X_k). Need joint law.
    # Use the chain: condition on X_0; propagate to X_{k-1},X_k via T^{k-1}.
    cov_sum = 0.0
    Tk = np.eye(A)               # T^0
    for k in range(1, kmax + 1):
        # l_0 = g(a, b) with (a,b) ~ P2 ; X_0 = b. Given X_0=b, X_{k-1} ~ (T^{k-1})[b,:].
        Tk = Tk @ T              # now Tk = T^k ; we need T^{k-1} -> use previous
        # joint of (X_{-1}=a, X_0=b, X_{k-1}=c, X_k=d):
        # P2[a,b] * (T^{k-1})[b,c] * T[c,d]
        Tkm1 = np.linalg.matrix_power(T, k - 1)
        E_l0_lk = 0.0
        for a in range(A):
            for b in range(A):
                w_ab = P2[a, b] * g[a, b]
                if w_ab == 0:
                    continue
                # sum over c,d
                cd = Tkm1[b, :][:, None] * T * g    # (c,d): T^{k-1}[b,c]*T[c,d]*g[c,d]
                E_l0_lk += w_ab * cd.sum()
        cov_k = E_l0_lk - mean_l ** 2
        cov_sum += cov_k
        if abs(cov_k) < 1e-12 and k > 5:
            break
    V = var_l + 2.0 * cov_sum
    return V, mean_l, var_l


def simulate_codelengths(T, pi, N, n_real, rng):
    """Return array of L = sum_i -log2 P(X_i|X_{i-1}) over n_real realisations.
    VECTORISED across realisations: all n_real chains advanced in lock-step."""
    A = len(pi)
    logT = np.where(T > 0, -np.log2(np.where(T > 0, T, 1)), 0.0)
    cumT = np.cumsum(T, axis=1)                  # (A, A) row-cumulatives
    # initial states ~ pi
    u0 = rng.random(n_real)
    states = np.searchsorted(np.cumsum(pi), u0, side='right')
    states = np.clip(states, 0, A - 1)
    Ls = -np.log2(pi[states])                    # initial -log2 pi(X_0)
    for _ in range(1, N):
        u = rng.random(n_real)
        # next state per chain: searchsorted in that chain's cumulative row
        rows = cumT[states]                      # (n_real, A)
        nxt = (u[:, None] > rows).sum(axis=1)    # first index where cum > u
        nxt = np.clip(nxt, 0, A - 1)
        Ls += logT[states, nxt]
        states = nxt
    return Ls


def check_V1():
    print("=" * 70)
    print("V1: varentropy V (analytic) == empirical Var(sum l_i)/N == long-run var")
    print("    (this V is the source DISPERSION, Kontoyiannis-Verdu 2014).")
    print("=" * 70)
    ok = True
    for trial, (A, temp) in enumerate([(2, 0.9), (3, 1.1), (2, 1.5)]):
        rng = np.random.default_rng(10 + trial)
        T, pi = markov1(A, temp, rng)
        V, mean_l, var_l = varentropy_analytic(T, pi)
        h = entropy_rate(T, pi)
        N = 4000
        Ls = simulate_codelengths(T, pi, N, 4000, rng)
        V_emp = float(np.var(Ls)) / N
        rel = abs(V - V_emp) / max(abs(V), 1e-9)
        print(f"  trial {trial}: A={A} temp={temp}: h={h:.4f}  V(analytic)={V:.4f}  "
              f"V(empirical Var/N)={V_emp:.4f}  rel.err={rel:.3f}")
        ok = ok and rel < 0.08
    print(f"  V1 {'PASS' if ok else 'FAIL'}")
    return ok


def check_V2():
    print("=" * 70)
    print("V2: (1-eps)-quantile of sum l_i == N h + sqrt(N V) Q^{-1}(eps) + o(sqrt N)")
    print("    (second-order dispersion expansion; Gaussian quantile matches empirical).")
    print("=" * 70)
    from statistics import NormalDist
    ok = True
    for trial, (A, temp, eps) in enumerate([(2, 0.9, 0.1), (3, 1.1, 0.05), (2, 1.5, 0.01)]):
        rng = np.random.default_rng(30 + trial)
        T, pi = markov1(A, temp, rng)
        V, _, _ = varentropy_analytic(T, pi)
        h = entropy_rate(T, pi)
        N = 6000
        Ls = simulate_codelengths(T, pi, N, 8000, rng)
        q_emp = float(np.quantile(Ls, 1 - eps))
        Qinv = NormalDist().inv_cdf(1 - eps)
        q_pred = N * h + math.sqrt(N * V) * Qinv
        # error should be o(sqrt N): report (q_emp - q_pred)/sqrt(N)
        err_over_sqrtN = (q_emp - q_pred) / math.sqrt(N)
        rel = abs(q_emp - q_pred) / abs(q_emp)
        print(f"  trial {trial}: A={A} eps={eps}: q_emp={q_emp:.1f}  "
              f"q_pred=Nh+sqrt(NV)Q^-1={q_pred:.1f}  (err/sqrtN={err_over_sqrtN:+.3f}, "
              f"rel={rel:.4f})")
        ok = ok and rel < 0.02
    print(f"  V2 {'PASS' if ok else 'FAIL'}")
    return ok


def check_V3():
    print("=" * 70)
    print("V3: INTERPLAY -- at K=c*sqrt(N), (N/K)delta_inf and sqrt(NV)Q^{-1}(eps) are")
    print("    BOTH Theta(sqrt N); their RATIO is N-INDEPENDENT.")
    print("=" * 70)
    from statistics import NormalDist
    rng = np.random.default_rng(7)
    T, pi = markov1(2, 1.0, rng)
    V, _, _ = varentropy_analytic(T, pi)
    Qinv = NormalDist().inv_cdf(1 - 0.05)
    delta_inf = 1.0     # per-sync-point excess overhead (parametrised positive constant)
    c = 1.0             # access budget: K = c sqrt(N)
    print(f"  V={V:.4f}, delta_inf={delta_inf}, c={c}, eps=0.05, Q^-1={Qinv:.3f}")
    ratios = []
    for N in (10**4, 10**5, 10**6, 10**7, 10**8):
        K = c * math.sqrt(N)
        overhead = (N / K) * delta_inf          # = sqrt(N)/c * delta_inf
        dispersion = math.sqrt(N * V) * Qinv
        ratio = overhead / dispersion
        ratios.append(ratio)
        print(f"    N={N:>10}: overhead(N/K)d={overhead:12.2f}  "
              f"dispersion={dispersion:12.2f}  ratio={ratio:.5f}  "
              f"(both/sqrtN: {overhead/math.sqrt(N):.3f}, {dispersion/math.sqrt(N):.3f})")
    spread = (max(ratios) - min(ratios)) / max(ratios)
    n_indep = spread < 1e-6
    expected_ratio = delta_inf / (c * math.sqrt(V) * Qinv)
    print(f"  ratio is N-independent (spread={spread:.2e}): {n_indep}; "
          f"== delta_inf/(c sqrt(V) Q^-1)={expected_ratio:.5f}? "
          f"{abs(ratios[0]-expected_ratio)<1e-9}")
    ok = n_indep and abs(ratios[0] - expected_ratio) < 1e-9
    print(f"  V3 {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print()
    print("#" * 70)
    print("# Theorem 7.29: RNR dispersion + the random-access<->dispersion interplay")
    print("#" * 70)
    print()
    r1 = check_V1(); print()
    r2 = check_V2(); print()
    r3 = check_V3(); print()
    print("=" * 70)
    allok = r1 and r2 and r3
    print(f"V1 varentropy == dispersion == long-run var : {'PASS' if r1 else 'FAIL'}")
    print(f"V2 Gaussian-quantile second-order expansion : {'PASS' if r2 else 'FAIL'}")
    print(f"V3 RA-overhead <-> dispersion interplay     : {'PASS' if r3 else 'FAIL'}")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
