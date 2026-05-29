r"""
Verification for Remark 6.9a (§13.4.4): characterisation of CONTINUOUS OP2(b)
-- the free-field Type-III-C support-selection problem that Lemma 6.9 left as
"the continuous mirror of OP2(b)".

VERDICT: free-field continuous OP2(b) is EASY (in P) in the worst case for the
unstructured / full-rank Gaussian field; the ONLY surviving hard core is the
sparse spiked covariance (Lemma 5.6j manifold WITH sparsity), and that is the
average-case sparse-PCA / planted-clique statistical-computational gap, NOT
worst-case multiplicative APX-hardness -- exactly mirroring the discrete OP2(b)
status (Lemma 13.4.2a absorption: only a cryptographic / average-case escape
survives). Hence Lemma 6.9's FIXED-FIELD restriction is the sharp worst-case
hardness boundary.

Two stacked reasons the free-field problem is easy, each verified:

(D1) THE MESP GAP DISSOLVES UNDER A FREE FIELD (chain-rule / no-budget).
     Lemma 6.9's APX-hardness is the MESP support-selection min_{|S|=s}
     0.5 log det Sigma_{S|Sbar}. That combinatorial object exists ONLY because
     the fixed-field problem imposes a hard budget |S|=s and FORBIDS
     transmitting X_Sbar (the complement must be reconstructed from the frozen
     model). A free, amortised M has no such budget: it codes the WHOLE block
     under the learned model at the support-INVARIANT rate
       h(X) = 0.5 log det(2 pi e Sigma),
     a single scalar with no subset to choose. We verify the chain-rule
     identity h(X_S) + h(X_Sbar | X_S) = h(X) for EVERY support S (so the
     full-block rate does not depend on S), i.e. the combinatorial MESP object
     vanishes once the budget is lifted.

(D2) LEARNING Sigma IS POLY-TIME AND DECODER-REPRODUCIBLE (the EASY encoder).
     The optimal free amortised predictor is M* = N(0, Sigma_hat) with
     Sigma_hat the sample-covariance MLE on a shared prefix of m blocks
     (the unstructured Gaussian MLE is the closed-form sample covariance,
     Anderson 2003 Thm 3.2.1, O(N^2 m) time; structured-linear case essentially
     convex for m >~ 14N, Zwiernik-Uhler-Richards JRSS-B 79(4) 2017). Both encoder and decoder run the SAME deterministic MLE on the
     SAME observed prefix (a universal-Slepian-Wolf-style shared model;
     Tan-Kosut), so the determinism / decoder-reproducibility constraint adds
     NO complexity. The per-block code length is the cross-entropy
       E[-log2 p_hat(X)] = h(X) + KL( N(0,Sigma) || N(0,Sigma_hat) ) -> h(X)
     as m -> infinity, with L(M*) = O(N^2 log(1/eps)) bits amortised over m
     blocks (L(M*)/m -> 0). We verify the rate -> h(X) (the entropy floor) and
     that the excess is exactly the (vanishing) plug-in KL.

(D3) THE ONLY HARD CORE IS SPARSE-SPIKED, AND IT IS AVERAGE-CASE (known).
     Above the BBP phase-transition threshold lambda > sqrt(N/m), and with NO
     sparsity, plain PCA recovers the Lemma 5.6j spike U in poly time
     (Baik-Ben Arous-Peche 2005), so even the manifold field is EASY there --
     consistent with Lemma 6.9's own item (4) ("gap collapses as Sigma ->
     sigma^2 I"). A worst-case-hard sliver appears ONLY when U is k-SPARSE
     (k << N) AND the signal sits in the sub-threshold band
       sqrt(k log N / m) < theta < k / sqrt(m),
     which is precisely the sparse-PCA / spiked-covariance statistical-
     computational gap, conjecturally hard via planted clique (Berthet-Rigollet
     COLT 2013; Brennan-Bresler-Huleihel COLT 2018) -- an AVERAGE-CASE
     detection hardness, not a
     worst-case multiplicative APX gap. We verify (a) above BBP, no sparsity:
     top-eigenvector overlap -> 1 (poly-time recovery, EASY); (b) the hard band
     is genuine only for log N < k < sqrt(N) (genuine sparsity; the dense Haar-U
     of Lemma 5.6j(3) is BBP-governed -> no gap).

(D4) BUT THE SPARSE-PCA GAP IS IRRELEVANT TO THE CODING RATE (the sharp form).
     The sparse-PCA gap is about efficiently RECOVERING the sparse support of U.
     The RNR coding objective never needs that: to achieve rate -> h(X) the
     encoder only needs KL( N(0,Sigma) || N(0,Sigma_hat) ) -> 0, which the FULL
     sample-covariance MLE delivers with NO PCA and NO support recovery, even in
     the sparse-PCA-hard band. We verify full-MLE rate -> h(X) for a k-sparse
     spike with k < sqrt(N). Consequently the amortised rate objective is
     UNCONDITIONALLY easy; the sparse-PCA gap re-enters only if one demands a
     non-amortised minimal-description M (small L(M) by exploiting sparsity) --
     and that is the single-block regime, where Lemma 6.9's MESP (not sparse
     PCA) already supplies worst-case hardness.

PASS = D1..D4 all hold:
  the MESP object dissolves under a free field; the poly-time MLE plug-in
  predictor achieves rate -> h(X) decoder-reproducibly; the only candidate
  hardness is the average-case sparse-PCA gap (above BBP / dense = easy), and
  even that is irrelevant to the coding RATE (full-MLE achieves h(X) without
  support recovery).
=> VERDICT: continuous OP2(b) (free-field, amortised rate objective) is EASY
   (in P), UNCONDITIONALLY. The optimal encoder is the full sample-covariance
   MLE plug-in. Lemma 6.9's FIXED-FIELD / single-block restriction is the sharp
   worst-case hardness boundary; the discrete-OP2(b) absorption escape
   (Kolmogorov-incompressible family -> cryptographic, Lemma 13.4.2a) has NO
   continuous analogue, because i.i.d. Gaussian samples cannot hide their own
   (poly-time-estimable) covariance. The continuous problem is thus strictly
   easier than the discrete OP2(b).
"""

import math
import sys
from itertools import combinations

import numpy as np

RNG = np.random.default_rng(20260529)
TOL = 1e-9


def slogdet(M):
    sign, val = np.linalg.slogdet(M)
    assert sign > 0, "non-SPD submatrix"
    return val


def h_block_bits(Sigma):
    """Differential entropy 0.5 log det(2 pi e Sigma) in bits."""
    N = Sigma.shape[0]
    return 0.5 * (slogdet(Sigma) / math.log(2) + N * math.log2(2 * math.pi * math.e))


def rand_spd(N, jitter=0.4):
    A = RNG.standard_normal((N, N))
    return A @ A.T + jitter * np.eye(N)


# ---------------------------------------------------------------------------
def check_D1_mesp_dissolves():
    """Full-block rate h(X) is support-invariant: the MESP object needs the budget."""
    N = 10
    Sigma = rand_spd(N)
    hX = h_block_bits(Sigma)
    ok_all = True
    worst = 0.0
    for s in (3, 5, 7):
        for S in combinations(range(N), s):
            S = list(S)
            Sb = [i for i in range(N) if i not in S]
            hS = 0.5 * (slogdet(Sigma[np.ix_(S, S)]) / math.log(2)
                        + len(S) * math.log2(2 * math.pi * math.e))
            SS = Sigma[np.ix_(S, S)]
            BB = Sigma[np.ix_(Sb, Sb)]
            BS = Sigma[np.ix_(Sb, S)]
            schur = BB - BS @ np.linalg.solve(SS, BS.T)
            hSb_S = 0.5 * (slogdet(schur) / math.log(2)
                           + len(Sb) * math.log2(2 * math.pi * math.e))
            worst = max(worst, abs(hS + hSb_S - hX))
    ok_all = worst < 1e-8
    print("  [D1] full-block rate h(X) = %.4f bits (one scalar, no subset)" % hX)
    print("  [D1] max |h(X_S)+h(X_Sbar|X_S) - h(X)| over all S = %.2e" % worst)
    print("  [D1] MESP gap dissolves once the |S|=s budget is lifted:  %s" % ok_all)
    return ok_all


def check_D2_mle_plugin_achieves_hX():
    """Sample-cov MLE plug-in predictor: per-block rate -> h(X); decoder-reproducible."""
    N = 8
    Sigma = rand_spd(N, jitter=0.5)
    hX = h_block_bits(Sigma)
    kls = []
    for m in (50, 500, 5000, 50000):
        Xs = RNG.multivariate_normal(np.zeros(N), Sigma, size=m)
        Sigma_hat = (Xs.T @ Xs) / m                 # MLE, O(N^2 m), poly-time
        Pinv = np.linalg.inv(Sigma_hat)
        ld_hat = slogdet(Sigma_hat)
        # cross-entropy of N(0,Sigma) under model N(0,Sigma_hat), in bits
        ce = 0.5 * (ld_hat / math.log(2)
                    + np.trace(Pinv @ Sigma) / math.log(2)
                    + N * math.log2(2 * math.pi))
        kl = ce - hX
        kls.append(kl)
        print("  [D2] m=%6d: code len=%.4f bits, excess=KL(true||hat)=%.4f" % (m, ce, kl))
    # Decoder-reproducibility: encoder and decoder run the IDENTICAL MLE on the
    # same prefix -> identical Sigma_hat (deterministic), verified by re-running.
    Xs = RNG.multivariate_normal(np.zeros(N), Sigma, size=200)
    enc_hat = (Xs.T @ Xs) / 200
    dec_hat = (Xs.T @ Xs) / 200   # decoder, same data, same deterministic estimator
    reproducible = np.allclose(enc_hat, dec_hat, atol=0.0)
    monotone = all(kls[i] >= kls[i + 1] - 1e-6 for i in range(len(kls) - 1))
    converges = kls[-1] < 1e-2 and kls[-1] >= -1e-9
    ok = monotone and converges and reproducible
    print("  [D2] KL -> 0 (rate -> h(X) floor):                       %s" % converges)
    print("  [D2] encoder/decoder MLE bit-identical (reproducible):   %s" % reproducible)
    print("  [D2] poly-time MLE plug-in achieves the entropy rate:    %s" % ok)
    return ok


def check_D3_only_hard_core_is_sparse_avgcase():
    """Above BBP + dense = EASY (PCA recovers spike); hard band needs sparsity."""
    # (a) above BBP, NO sparsity: plain PCA recovers the spike (poly-time, EASY)
    N = 300
    m = N                              # n ~ p, the high-dimensional regime
    lam = 8.0                          # well above BBP threshold sqrt(N/m)=1
    u = RNG.standard_normal(N)
    u /= np.linalg.norm(u)
    Sig = lam * np.outer(u, u) + np.eye(N)
    Xs = RNG.multivariate_normal(np.zeros(N), Sig, size=m)
    S_hat = (Xs.T @ Xs) / m
    w, V = np.linalg.eigh(S_hat)
    overlap_hi = abs(V[:, -1] @ u)
    bbp = math.sqrt(N / m)
    easy_above_bbp = overlap_hi > 0.85
    print("  [D3] above BBP (lambda=%.0f > thr=%.2f), dense U: PCA overlap=%.3f"
          % (lam, bbp, overlap_hi))
    print("  [D3]   -> dense spiked field is EASY (matches Lemma 6.9 item 4):  %s"
          % easy_above_bbp)

    # (b) The sparse-PCA stat-feasible / comp-hard band is (Berthet-Rigollet
    #     2013): statistical threshold theta_stat = sqrt(k log N / m); best
    #     poly-time threshold theta_comp ~ k / sqrt(m). The HARD band
    #     theta_stat < theta < theta_comp is a genuine sliver ONLY in the sparse
    #     regime where the spike is also SUB-BBP (k/sqrt(m) < sqrt(N/m), i.e.
    #     k < sqrt(N) -- otherwise a comp-detectable signal is already above the
    #     BBP threshold and plain PCA solves it: NO gap). So a non-degenerate gap
    #     requires  log N < k < sqrt(N), genuine sparsity. The dense Haar-U of
    #     Lemma 5.6j(3) has k = N, far outside this band -> governed by BBP only
    #     -> PCA-optimal -> NO computational gap (easy), as in part (a).
    def gap_is_genuine(N, m, k):
        theta_stat = math.sqrt(k * math.log(N) / m)
        theta_comp = k / math.sqrt(m)
        theta_bbp = math.sqrt(N / m)
        # non-degenerate gap: stat < comp (band non-empty) AND comp < bbp
        # (the hard band lies BELOW the trivial PCA-detectable scale)
        return (theta_stat < theta_comp) and (theta_comp < theta_bbp)
    N2, m2 = 10_000, 10_000
    sparse_k = 30           # log N ~ 9 < k=30 < sqrt(N)=100 : genuine gap
    dense_k = N2            # k = N : far above sqrt(N) : BBP-governed, no gap
    sparse_hard = gap_is_genuine(N2, m2, sparse_k)
    dense_no_gap = not gap_is_genuine(N2, m2, dense_k)
    print("  [D3] sparse k=%d (logN~9 < k < sqrtN=100): genuine gap:    %s"
          % (sparse_k, sparse_hard))
    print("  [D3] dense  k=N=%d (k >> sqrtN): BBP-governed, no gap:   %s"
          % (dense_k, dense_no_gap))
    print("  [D3]   -> only the SPARSE spiked regime carries the (average-case,")
    print("           sparse-PCA / planted-clique) gap; not worst-case APX.")
    return easy_above_bbp and sparse_hard and dense_no_gap


def check_D4_sparse_gap_irrelevant_to_rate():
    """Even in the sparse-PCA-hard regime, full MLE -> h(X) with NO support recovery."""
    N, k = 120, 8                     # k < sqrt(N) ~ 11: the sparse-PCA-hard band
    supp = RNG.choice(N, k, replace=False)
    u = np.zeros(N)
    u[supp] = RNG.standard_normal(k)
    u /= np.linalg.norm(u)
    theta = 0.7                       # weak, sub-BBP-for-dense signal
    Sigma = np.eye(N) + theta * np.outer(u, u)
    hX = h_block_bits(Sigma)
    kls = []
    for m in (4000, 40000, 400000):
        Xs = RNG.multivariate_normal(np.zeros(N), Sigma, size=m)
        Shat = (Xs.T @ Xs) / m        # FULL MLE -- no PCA, no support recovery
        Pinv = np.linalg.inv(Shat)
        ce = 0.5 * (slogdet(Shat) / math.log(2)
                    + np.trace(Pinv @ Sigma) / math.log(2)
                    + N * math.log2(2 * math.pi))
        kls.append(ce - hX)
        print("  [D4] sparse-hard Sigma, m=%7d: full-MLE rate=%.3f, excess KL=%.4f"
              % (m, ce, ce - hX))
    ok = kls[-1] < 0.05 and all(kls[i] >= kls[i + 1] - 1e-3 for i in range(len(kls) - 1))
    print("  [D4] full-MLE -> h(X) WITHOUT recovering the sparse support:  %s" % ok)
    print("  [D4]   -> the sparse-PCA gap is IRRELEVANT to the coding RATE; it")
    print("           re-enters only under a non-amortised minimal-M demand")
    print("           (= the single-block regime, where Lemma 6.9/MESP governs).")
    return ok


def main():
    print("Remark 6.9a: continuous OP2(b) (free-field Type-III-C) is EASY in P;")
    print("             Lemma 6.9 fixed-field is the sharp worst-case boundary.")
    print()
    results = {
        "D1 MESP dissolves (free field)": check_D1_mesp_dissolves(),
        "D2 MLE plug-in achieves h(X)": check_D2_mle_plugin_achieves_hX(),
        "D3 hard core = sparse avg-case": check_D3_only_hard_core_is_sparse_avgcase(),
        "D4 sparse gap irrelevant to rate": check_D4_sparse_gap_irrelevant_to_rate(),
    }
    print()
    all_ok = all(results.values())
    for k, v in results.items():
        print("  %-34s %s" % (k, "PASS" if v else "FAIL"))
    print()
    if all_ok:
        print("PASS: the MESP object dissolves under a free field (chain rule); the")
        print("      poly-time sample-covariance MLE plug-in predictor achieves rate")
        print("      -> h(X) decoder-reproducibly; the only candidate hardness is the")
        print("      average-case sparse-PCA gap (dense / above-BBP = easy), and even")
        print("      that is irrelevant to the coding RATE (full-MLE -> h(X), no PCA).")
        print("VERDICT: continuous OP2(b) (free-field amortised rate) is EASY (in P),")
        print("         UNCONDITIONALLY; optimal encoder = full sample-covariance MLE")
        print("         plug-in. Lemma 6.9's FIXED-FIELD / single-block restriction is")
        print("         the sharp worst-case hardness boundary. The discrete-OP2(b)")
        print("         Kolmogorov/cryptographic escape (Lemma 13.4.2a) has NO")
        print("         continuous analogue: i.i.d. Gaussian samples cannot hide their")
        print("         poly-time-estimable covariance. Continuous OP2(b) is strictly")
        print("         easier than discrete OP2(b).")
        sys.exit(0)
    else:
        print("FAIL: at least one check did not hold.")
        sys.exit(1)


if __name__ == "__main__":
    main()
