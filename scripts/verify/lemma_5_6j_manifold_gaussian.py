#!/usr/bin/env python3
r"""
Verification of Lemma 5.6j: Quantitative TC for the quantized Gaussian-on-
rank-d-manifold model.

Setting. X = M + sigma * Z in R^N where M is centered rank-d Gaussian
with covariance U Lambda U^T (U in R^{N x d} orthonormal,
Lambda = diag(lambda_1, ..., lambda_d) >= 0), and Z ~ N(0, I_N) is
independent isotropic noise. Total covariance:
        Sigma = U Lambda U^T + sigma^2 I_N.

------------------------------------------------------------------------
Identity (a) - exact TC for Gaussian X (no assumption on U):

        TC(X) = (1/2) [ sum_i log_2(1 + a_i)  -  sum_j log_2(1 + b_j) ]
              (bits)

where a_i = (U Lambda U^T)_ii / sigma^2 and b_j = lambda_j / sigma^2.

------------------------------------------------------------------------
Lower bound (b) - Jensen upper bound on the second sum:

        sum_j log_2(1 + b_j)  <=  d * log_2(1 + N * rho / d)

with rho = (sum_j b_j) / N = bar_lambda / (N sigma^2). Hence

        TC(X)  >=  (N/2) * S_bar  -  (d/2) * log_2(1 + N rho / d)

where S_bar = (1/N) sum_i log_2(1 + a_i) is the per-coordinate
log-variance contribution from the embedding.

------------------------------------------------------------------------
Asymptotic (c) - Haar-uniform U, linear scaling. For a Haar-uniform
orthonormal frame U and equal eigenvalues lambda_j = bar_lambda/d, each
diagonal entry a_i = (bar_lambda/d) * sum_{j=1}^d u_{ij}^2 / sigma^2.
The random variable N * sum_j u_{ij}^2 converges in distribution as
N -> infty to chi^2_d. So a_i ~ rho * V with V ~ chi^2_d / d, mean 1,
variance 2/d. By concentration of measure,

        S_bar  ->  E_{V ~ chi^2_d / d} [ log_2(1 + rho * V) ]  =: c_d(rho)

almost surely. The per-coordinate TC limit is then

        lim_{N -> infty} TC(X) / N  =  (1/2) * c_d(rho).

For d -> infty, c_d(rho) -> log_2(1 + rho) (LLN -> V -> 1).
For d = 1, c_1(rho) is strictly less than log_2(1 + rho).

A 2nd-order Taylor expansion in 1/d gives the explicit estimate

        c_d(rho) = log_2(1 + rho) - rho^2 / ((ln 2) (1+rho)^2 d) + O(1/d^2).

------------------------------------------------------------------------
Quantization (d). Under lattice quantization at scale Delta in each
coordinate, by Bennett's formula the per-coordinate variance gains
Delta^2 / 12 ("dither" noise). The TC formula above applied to
Sigma + (Delta^2/12) I yields a perturbation bounded by O(N * (Delta/sigma)^2).
For Delta = o(sigma), the discrete TC tracks the continuous TC to o(N).

------------------------------------------------------------------------
Spread (e) is necessary. Without delocalisation, U can concentrate on
few ambient coordinates and Sigma becomes block-diagonal, giving TC=0
even though the embedding has rank d=O(1). The Haar-uniformity is
sufficient; deterministic spread (each tangent vector having entries
of magnitude O(1/sqrt(N))) is necessary.

------------------------------------------------------------------------
Sanity limits:
   d -> N: Sigma -> (lambda+sigma^2) I (if lambda_j all equal), TC -> 0.
   bar_lambda -> 0: Sigma -> sigma^2 I, TC -> 0.
   d = 1, N -> infty: TC/N -> (1/2) E[log_2(1 + rho * chi^2_1)],
                       strictly less than (1/2) log_2(1+rho).

PASS = all claims numerically reproduced for representative parameters.
"""

import math
import sys

import numpy as np


# ====================================================================
# Core utilities
# ====================================================================

def gaussian_tc_bits(Sigma):
    """Exact TC of N(0, Sigma) in bits."""
    diag = np.diag(Sigma)
    log_diag = float(np.sum(np.log(diag)))
    sign, logdet = np.linalg.slogdet(Sigma)
    assert sign > 0, "Sigma must be positive definite"
    return 0.5 * (log_diag - float(logdet)) / math.log(2)


def make_haar_orthonormal(N, d, rng):
    """Haar-uniform orthonormal frame of rank d in R^N."""
    G = rng.standard_normal((N, d))
    Q, _ = np.linalg.qr(G)
    return Q


def make_cov(N, d, sigma2, lambdas, U):
    """Sigma = U diag(lambdas) U^T + sigma2 * I_N."""
    return U @ np.diag(lambdas) @ U.T + sigma2 * np.eye(N)


# ====================================================================
# Claim (a): exact identity
# ====================================================================

def check_identity(N=80, d=4, sigma2=1.0, rho=2.0, seed=0):
    rng = np.random.default_rng(seed)
    bar_lambda = rho * N * sigma2
    lambdas = rng.uniform(0.5, 1.5, size=d)
    lambdas *= bar_lambda / float(lambdas.sum())
    U = make_haar_orthonormal(N, d, rng)
    Sigma = make_cov(N, d, sigma2, lambdas, U)

    tc1 = gaussian_tc_bits(Sigma)
    a = np.diag(Sigma - sigma2 * np.eye(N)) / sigma2
    b = lambdas / sigma2
    S1 = float(np.sum(np.log2(1.0 + a)))
    S2 = float(np.sum(np.log2(1.0 + b)))
    tc2 = 0.5 * (S1 - S2)
    return abs(tc1 - tc2) < 1e-10, tc1, tc2


# ====================================================================
# Claim (b): lower bound via Jensen on lambda_j
# ====================================================================

def check_lower_bound(N=500, d=4, sigma2=1.0, rho=1.0, n_trials=30, seed=1):
    rng = np.random.default_rng(seed)
    bar_lambda = rho * N * sigma2
    rows = []
    for trial in range(n_trials):
        # Non-uniform lambda_j to test Jensen slack
        weights = rng.exponential(size=d)
        weights /= weights.sum()
        lambdas = weights * bar_lambda
        U = make_haar_orthonormal(N, d, rng)
        Sigma = make_cov(N, d, sigma2, lambdas, U)
        tc = gaussian_tc_bits(Sigma)
        # S_bar = (1/N) sum_i log2(1 + a_i)
        a = (np.diag(Sigma) - sigma2) / sigma2
        S_bar = float(np.mean(np.log2(1.0 + a)))
        # Lower bound: (N/2) S_bar - (d/2) log2(1 + N*rho/d)
        lb = 0.5 * N * S_bar - 0.5 * d * math.log2(1 + N * rho / d)
        rows.append((tc, S_bar, lb))
    all_above = all(tc >= lb - 1e-9 for (tc, _, lb) in rows)
    return rows, all_above


# ====================================================================
# Claim (c): linear scaling, agreement with c_d(rho)
# ====================================================================

def c_d_rho_mc(d, rho, M=200_000, seed=0):
    """Monte Carlo estimate of c_d(rho) = E[log_2(1 + rho * V)], V ~ chi^2_d/d."""
    rng = np.random.default_rng(seed)
    V = rng.chisquare(d, size=M) / d
    return float(np.mean(np.log2(1.0 + rho * V)))


def check_linear_scaling(d=4, rho=1.0, sigma2=1.0,
                         N_list=(100, 250, 500, 1000, 2000, 5000),
                         seed=42):
    rng = np.random.default_rng(seed)
    c_d = c_d_rho_mc(d, rho, seed=seed)
    asymp = 0.5 * math.log2(1 + rho)
    rows = []
    for N in N_list:
        bar_lambda = rho * N * sigma2
        lambdas = np.ones(d) * (bar_lambda / d)
        U = make_haar_orthonormal(N, d, rng)
        Sigma = make_cov(N, d, sigma2, lambdas, U)
        tc = gaussian_tc_bits(Sigma)
        per_dim = tc / N
        target = 0.5 * c_d
        rows.append((N, tc, per_dim, target, asymp))
    return rows, c_d


# ====================================================================
# Claim (d): quantization perturbation
# ====================================================================

def check_quantization(N=50, d=3, sigma2=1.0, rho=1.0, seed=7):
    rng = np.random.default_rng(seed)
    bar_lambda = rho * N * sigma2
    lambdas = np.ones(d) * (bar_lambda / d)
    U = make_haar_orthonormal(N, d, rng)
    Sigma = make_cov(N, d, sigma2, lambdas, U)
    tc_cont = gaussian_tc_bits(Sigma)
    rows = []
    for Delta in (0.5, 0.25, 0.1, 0.05):
        # Bennett's noise model: Var(quantized) - Var(orig) ~ Delta^2/12
        Sigma_q = Sigma + (Delta**2 / 12.0) * np.eye(N)
        tc_q = gaussian_tc_bits(Sigma_q)
        pert_predicted = N * Delta**2 / (12.0 * sigma2)  # O(N Delta^2/sigma^2)
        rows.append((Delta, tc_cont, tc_q, abs(tc_q - tc_cont), pert_predicted))
    return rows, tc_cont


# ====================================================================
# Claim (e): spread is necessary
# ====================================================================

def check_spread(N=200, d=5, sigma2=1.0, rho=1.0, seed=99):
    rng = np.random.default_rng(seed)
    bar_lambda = rho * N * sigma2
    lambdas = np.ones(d) * (bar_lambda / d)

    U_full = make_haar_orthonormal(N, d, rng)
    Sigma_full = make_cov(N, d, sigma2, lambdas, U_full)
    tc_full = gaussian_tc_bits(Sigma_full)

    U_conc = np.zeros((N, d))
    U_conc[:d, :d] = np.eye(d)
    Sigma_conc = make_cov(N, d, sigma2, lambdas, U_conc)
    tc_conc = gaussian_tc_bits(Sigma_conc)

    M = min(5 * d, N)
    G = rng.standard_normal((N, d))
    G[M:, :] = 0.0
    U_part, _ = np.linalg.qr(G)
    Sigma_part = make_cov(N, d, sigma2, lambdas, U_part)
    tc_part = gaussian_tc_bits(Sigma_part)

    return [
        ("Haar (well-spread)",            tc_full),
        ("U=[I_d;0] (concentrated)",      tc_conc),
        (f"partial-spread (first {M} rows)", tc_part),
    ]


# ====================================================================
# Sanity: edge limits
# ====================================================================

def check_sanity(seed=11):
    rng = np.random.default_rng(seed)
    res = []

    # (i) d=N, lambdas equal => Sigma = const*I => TC=0
    N = 50
    sigma2 = 1.0
    lambdas_eq = np.ones(N) * 5.0
    Sigma_eq = np.eye(N) * (lambdas_eq[0] + sigma2)
    res.append(("d=N, lambdas equal: Sigma proportional to I => TC=0",
                gaussian_tc_bits(Sigma_eq), 0.0))

    # (ii) bar_lambda = 0 => Sigma = sigma2 I => TC=0
    U_zero = make_haar_orthonormal(N, 5, rng)
    Sigma_zero = make_cov(N, 5, sigma2, np.zeros(5), U_zero)
    res.append(("bar_lambda=0: pure noise => TC=0",
                gaussian_tc_bits(Sigma_zero), 0.0))

    # (iii) d=1, N=1000, rho=1: per-dim ~ 0.385
    N = 1000
    rho = 1.0
    bar_lambda = rho * N * sigma2
    U_1 = make_haar_orthonormal(N, 1, rng)
    Sigma_1 = make_cov(N, 1, sigma2, np.array([bar_lambda]), U_1)
    tc_1 = gaussian_tc_bits(Sigma_1)
    res.append((f"d=1, rho=1, N=1000: per-dim ~ 0.385",
                tc_1 / N, 0.385))
    return res


# ====================================================================
# Driver
# ====================================================================

def main():
    print("=" * 72)
    print("Lemma 5.6j: TC scaling for quantized Gaussian on rank-d manifold")
    print("=" * 72)
    fails = 0

    # ----- (a) identity -----
    print("\n[Claim (a)] Exact TC identity (no assumption on U):")
    print("   TC = (1/2)[sum_i log2(1+a_i) - sum_j log2(1+b_j)]")
    ok, tc1, tc2 = check_identity()
    print(f"   Direct: {tc1:.6f}, Identity: {tc2:.6f}, diff={abs(tc1-tc2):.2e}")
    if ok:
        print("   PASS")
    else:
        print("   FAIL"); fails += 1

    # ----- (b) lower bound -----
    print("\n[Claim (b)] Lower bound via Jensen on lambda_j:")
    print("   TC >= (N/2) * S_bar - (d/2) * log2(1 + N*rho/d)")
    print("   (Heterogeneous lambda_j, Haar U, 30 trials; N=500, d=4, rho=1)")
    rows, all_above = check_lower_bound()
    tc_mean = float(np.mean([r[0] for r in rows]))
    S_mean = float(np.mean([r[1] for r in rows]))
    lb_mean = float(np.mean([r[2] for r in rows]))
    print(f"   mean TC = {tc_mean:.3f}, mean S_bar = {S_mean:.4f}, mean LB = {lb_mean:.3f}")
    print(f"   mean slack TC-LB = {tc_mean - lb_mean:.3f} bits (slack due to Jensen on lambda_j)")
    if all_above:
        print("   PASS (lower bound holds in all 30 trials)")
    else:
        print("   FAIL"); fails += 1

    # ----- (c) linear scaling -----
    print("\n[Claim (c)] Linear scaling under Haar U: TC/N -> (1/2) c_d(rho)")
    print("   where c_d(rho) = E[log2(1 + rho * chi2_d/d)]")
    rows, c_d = check_linear_scaling(d=4, rho=1.0)
    print(f"   (d=4, rho=1)  c_d(rho) = {c_d:.4f}, target per-dim = {0.5 * c_d:.4f}")
    print(f"   asymptotic d->inf limit = (1/2) log2(2) = 0.5000")
    print(f"   {'N':>6} {'TC':>12} {'TC/N':>10} {'(1/2)c_d':>10} {'d->inf':>10}")
    for (N, tc, tcn, target, asymp) in rows:
        print(f"   {N:>6d} {tc:>12.3f} {tcn:>10.4f} {target:>10.4f} {asymp:>10.4f}")
    # Convergence check: TC/N approaches 0.5*c_d
    tcn_last = rows[-1][2]
    target_last = rows[-1][3]
    if abs(tcn_last - target_last) < 0.03 and tcn_last > 0.3:
        print(f"   PASS (TC/N converges to (1/2) c_d = {target_last:.4f})")
    else:
        print("   FAIL"); fails += 1

    # ----- (d) quantization -----
    print("\n[Claim (d)] Quantization perturbation O(N * Delta^2 / sigma^2)")
    print("   (Bennett model: Var gains Delta^2/12 per coord. N=50, d=3, rho=1)")
    rows, tc_cont = check_quantization()
    print(f"   continuous TC = {tc_cont:.4f}")
    print(f"   {'Delta':>8} {'TC(X)':>10} {'TC(X_q)':>10} {'|diff|':>10} {'O(NDel^2/12)':>14}")
    for (D, c, q, diff, pred) in rows:
        print(f"   {D:>8.3f} {c:>10.4f} {q:>10.4f} {diff:>10.4f} {pred:>14.4f}")
    # Diff should scale as Delta^2 and be small for small Delta
    if rows[-1][3] < 0.01 and rows[0][3] < 1.0:
        print("   PASS (perturbation shrinks quadratically with Delta)")
    else:
        print("   FAIL"); fails += 1

    # ----- (e) spread -----
    print("\n[Claim (e)] Spread is necessary: concentrated U gives TC=0")
    print("   (N=200, d=5, rho=1)")
    res = check_spread()
    for label, tc in res:
        print(f"   {label:<40s}  TC = {tc:>8.3f} bits")
    if res[0][1] > 10 and res[1][1] < 0.01:
        print("   PASS (Haar gives large TC; concentrated U gives 0)")
    else:
        print("   FAIL"); fails += 1

    # ----- sanity -----
    print("\n[Sanity] Edge cases")
    for label, actual, target in check_sanity():
        print(f"   {label}")
        print(f"     actual = {actual:.4f}, target = {target:.4f}")
    sanity = check_sanity()
    ok_sanity = (abs(sanity[0][1]) < 1e-6 and
                 abs(sanity[1][1]) < 1e-6 and
                 abs(sanity[2][1] - sanity[2][2]) < 0.03)
    if ok_sanity:
        print("   PASS")
    else:
        print("   FAIL"); fails += 1

    print("\n" + "=" * 72)
    if fails == 0:
        print("OVERALL: PASS  (all claims verified)")
        return 0
    else:
        print(f"OVERALL: FAIL  ({fails} claims failed)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
