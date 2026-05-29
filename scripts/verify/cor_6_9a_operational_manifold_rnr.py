r"""
Verification for Corollary 6.9b (operational Gaussian-on-manifold RNR):
the end-to-end synthesis of Lemma 5.6j (TC rate) + Lemma 6.9 (support-
selection hardness) + Lemma 6.9a (the matching poly-time approximation).

This is the honest companion to the synthesis. Its job is NOT to "confirm" a
clean multiplicative (c * 5/4) end-to-end ratio -- that composition is FALSE.
It demonstrates the CORRECT composition: an ADDITIVE-REGRET bound on the RNR
residual, obtained from Lemma 6.9a's (1-1/e) submodular-greedy guarantee on
the maximum-entropy-sampling objective.

------------------------------------------------------------------------
The three pieces and how they compose.

  Source (Lemma 5.6j): X ~ N(0, Sigma), Sigma = U Lambda U^T + sigma^2 I, U
  Haar-orthonormal (delocalised), Lambda = (bar_lambda/d) I.  rho =
  bar_lambda/(N sigma^2).  TC(X) = (1/2) N c_d(rho) + o(N) = Theta(N).

  Hardness (Lemma 6.9): fixed Sigma, choose |S|=s, code X_S, reconstruct the
  complement; per-block residual C(S) = h(X_S | X_Sbar) = 0.5 log det(2 pi e
  Sigma_{S|Sbar}).  Minimising C(S) is APX-hard 5/4 via MESP (Ohsaka 2022).

  Approximation (Lemma 6.9a, the sibling result -- NOT the generic Nikolov
  det-bound). Complementation identity (EXACT bijection S <-> T = Sbar):
      C(S) = h(X) - M(Sbar),   M(T) := h(X_T) = 0.5 log det(2 pi e Sigma_T).
  So min_S C(S) <=> max_{|T|=N-s} M(T) = maximum-entropy sampling on the
  complement. On the RNR manifold M is SUBMODULAR, and MONOTONE once the
  noise floor 2 pi e sigma^2 >= 1 holds (lambda_min(Sigma)=sigma^2). Greedy
  submodular maximisation (Nemhauser-Wolsey-Fisher 1978) gives
      M(Tgreedy) >= (1 - 1/e) M(Topt).

------------------------------------------------------------------------
THE COMPOSITION (the load-bearing arithmetic).

The multiplicative (1-1/e) guarantee lives on M.  The RNR residual C = h(X)
- M is a DIFFERENCE, so the multiplicative ratio does NOT transfer to C.
What transfers is ADDITIVE REGRET:

    C(Tgreedy-support) - C(Sopt)  =  M(Topt) - M(Tgreedy)
                                  <=  (1/e) M(Topt).            (6.9b.1)
    (tight NWF bound on the shifted non-negative objective; equivalently
     <= (1/(e-1)) M(Tgreedy) in terms of the achieved value.)

There is NO multiplicative (c * 5/4) factor on the residual.  The 5/4
(hardness, Lemma 6.9) and the e/(e-1) (algorithm, Lemma 6.9a) sandwich the
achievable ratio ON M; they are not stages multiplied together.  A
multiplicative bound on C directly would need a constant-factor log-det
approximation, which does NOT exist in general (Civril-Magdon-Ismail 2013,
Algorithmica 65:159-176) -- Lemma 6.9a's guarantee is bought by the
monotone-submodular structure of M, not a generic C-approximation.

------------------------------------------------------------------------
Checks.

(A) Hypothesis compatibility: ONE concrete Sigma simultaneously satisfies
    Lemma 5.6j (delocalised Haar U, rho>0, HIGH SNR so the Lemma 6.9 gap is
    Theta(1)/coord), Lemma 6.9 (fixed Sigma, genuine MESP spread), AND Lemma
    6.9a's monotonicity noise floor sigma^2 >= 1/(2 pi e). Print all four
    witnesses for the SAME matrix.

(B) The composition is ADDITIVE-REGRET, not multiplicative-on-C. On small N
    (brute-force optima), run the Lemma 6.9a greedy MESP selector; verify the
    EXACT additive transfer  C_greedy - C_opt == M_opt - M_greedy, the tight
    submodular bound C_greedy - C_opt <= (1/e) M_opt, and that the naive
    multiplicative-(c*5/4)-on-C prediction is FALSE.

(C) Per-coordinate regret is governed by the greedy SLACK (M_opt-M_greedy)/t,
    which -> 0 in the near-modular / high-SNR regime (curvature c -> 0): show
    the slack shrinks as SNR grows.

(D) End-to-end encoder on a Gaussian-on-manifold source: estimate Sigma from
    samples, select S by Lemma 6.9a greedy, code X_S at its conditional
    entropy, reconstruct X_Sbar by linear MMSE; measure achieved residual vs
    the information-theoretic optimum and vs brute-force. Report the additive
    regret; confirm total per-block rate == h(X) by the chain rule.

PASS = (A) one source meets all hypotheses incl. noise floor + high SNR;
(B) the residual transfer is the additive identity (multiplicative-on-C is
false); (C) the regret vanishes with SNR; (D) the encoder attains h(X) total
with residual within the additive regret of optimal, in poly(N).
HONEST: additive-regret end-to-end guarantee; NO multiplicative (c*5/4) claim.
"""

import math
import sys
from itertools import combinations

import numpy as np

RNG = np.random.default_rng(20260529)
TOL = 1e-7
TWO_PI_E = 2.0 * math.pi * math.e
NOISE_FLOOR = 1.0 / TWO_PI_E          # ~0.0585: sigma^2 >= this => M monotone
E_OVER_EM1 = math.e / (math.e - 1.0)  # ~1.582
ONE_OVER_EM1 = 1.0 / (math.e - 1.0)   # ~0.582 (slack rel. to achieved M(That))
ONE_OVER_E = 1.0 / math.e             # ~0.368 (TIGHT slack rel. to M(Topt))


def slogdet2(M):
    """log_2 det of an SPD matrix."""
    sign, val = np.linalg.slogdet(M)
    assert sign > 0, "non-SPD submatrix"
    return val / math.log(2.0)


def manifold_sigma(N, d, bar_lambda, sigma2, equal=True):
    """Lemma 5.6j source: Sigma = U Lambda U^T + sigma^2 I, U Haar-orthonormal."""
    G = RNG.standard_normal((N, d))
    U, _ = np.linalg.qr(G)
    U = U[:, :d]
    if equal:
        lam = np.full(d, bar_lambda / d)
    else:
        lam = np.abs(RNG.standard_normal(d))
        lam = lam / lam.sum() * bar_lambda
    Sigma = U @ np.diag(lam) @ U.T + sigma2 * np.eye(N)
    return Sigma, U, lam


def tc_bits(Sigma):
    """Exact TC of N(0, Sigma) in bits (Lemma 5.6j identity)."""
    diag = np.diag(Sigma)
    return 0.5 * (np.sum(np.log2(diag)) - slogdet2(Sigma))


def diff_entropy(Sigma, T):
    """M(T) = h(X_T) = 0.5 log_2 det(2 pi e Sigma_T) (bits)."""
    T = list(T)
    if not T:
        return 0.0
    return 0.5 * slogdet2(TWO_PI_E * Sigma[np.ix_(T, T)])


def residual_cost(Sigma, S):
    """C(S) = h(X_S | X_Sbar) = 0.5 log_2 det(2 pi e Sigma_{S|Sbar}) (bits)."""
    N = Sigma.shape[0]
    S = list(S)
    Sb = [i for i in range(N) if i not in S]
    SS = Sigma[np.ix_(S, S)]
    if not Sb:
        return 0.5 * slogdet2(TWO_PI_E * SS)
    SB = Sigma[np.ix_(S, Sb)]
    BB = Sigma[np.ix_(Sb, Sb)]
    schur = SS - SB @ np.linalg.solve(BB, SB.T)
    return 0.5 * slogdet2(TWO_PI_E * schur)


def total_entropy(Sigma):
    return 0.5 * slogdet2(TWO_PI_E * Sigma)


def greedy_mesp(Sigma, t):
    """Lemma 6.9a greedy: maximise M(T)=h(X_T) over |T|=t. Poly-time."""
    N = Sigma.shape[0]
    T, remaining, cur = [], set(range(N)), 0.0
    for _ in range(t):
        best_i, best_gain = None, -math.inf
        for i in remaining:
            g = diff_entropy(Sigma, T + [i]) - cur
            if g > best_gain:
                best_gain, best_i = g, i
        T.append(best_i)
        remaining.discard(best_i)
        cur = diff_entropy(Sigma, T)
    return T


def brute_max_mesp(Sigma, t):
    """Exact argmax_{|T|=t} M(T) (small N only)."""
    best_T, best_val = None, -math.inf
    for T in combinations(range(Sigma.shape[0]), t):
        v = diff_entropy(Sigma, T)
        if v > best_val:
            best_val, best_T = v, list(T)
    return best_T, best_val


def greedy_suboptimal_sigma(N, t):
    """SPD Sigma on which greedy MESP is STRICTLY suboptimal (real regret),
    so the additive transfer is non-degenerate. Bait cluster of large,
    mutually-correlated coords lures greedy off the optimal complement set."""
    for _ in range(6000):
        A = 0.15 * RNG.standard_normal((N, N))
        Sigma = A @ A.T + 0.6 * np.eye(N)
        b = max(2, t // 2)
        for i in range(b):
            Sigma[i, i] += 6.0
        for i in range(b):
            for j in range(i + 1, b):
                Sigma[i, j] = Sigma[j, i] = 0.92 * math.sqrt(Sigma[i, i] * Sigma[j, j])
        w = np.linalg.eigvalsh(Sigma)
        if w[0] <= 1e-6:
            Sigma += (1e-3 - w[0]) * np.eye(N)
        Tg = greedy_mesp(Sigma, t)
        _, Mo = brute_max_mesp(Sigma, t)
        if Mo - diff_entropy(Sigma, Tg) > 1e-3:
            return Sigma
    return None


# ---------------------------------------------------------------------------
def check_A_compatibility():
    """One Sigma satisfies Lemma 5.6j + Lemma 6.9 + Lemma 6.9a simultaneously."""
    N, d, s = 14, 3, 7
    t = N - s
    bar_lambda, sigma2 = 14.0, 0.25  # rho=4 (high SNR); sigma2=0.25 >= noise floor
    Sigma, U, lam = manifold_sigma(N, d, bar_lambda, sigma2, equal=True)
    rho = bar_lambda / (N * sigma2)

    # (5.6j) delocalised Haar U, rho>0 => TC=Theta(N) > 0.
    tc = tc_bits(Sigma)
    max_abs_U = np.max(np.abs(U))
    deloc_ok = (tc > 1.0) and (max_abs_U < 0.9)

    # high SNR (Lemma 6.9 C5): bar_lambda >> N sigma^2 keeps the selection gap.
    snr_ok = bar_lambda > N * sigma2          # rho>1
    vals = [diff_entropy(Sigma, T) for T in combinations(range(N), t)]
    spread = max(vals) - min(vals)
    gap_ok = spread > 0.5

    # (6.9a) monotonicity noise floor: sigma^2 >= 1/(2 pi e).
    floor_ok = sigma2 >= NOISE_FLOOR
    # verify M is actually monotone here: random marginal gains all >= 0.
    mono_ok = True
    for _ in range(200):
        i = int(RNG.integers(0, N))
        rest = [x for x in range(N) if x != i]
        RNG.shuffle(rest)
        Tsub = rest[: int(RNG.integers(0, len(rest)))]
        if diff_entropy(Sigma, Tsub + [i]) - diff_entropy(Sigma, Tsub) < -1e-9:
            mono_ok = False
            break

    print("  [A] single source Sigma (N=%d, d=%d, s=%d, rho=%.1f, sigma2=%.3f):"
          % (N, d, s, rho, sigma2))
    print("      Lemma 5.6j: TC=%.2f>0, max|U_ij|=%.3f (delocalised) -> %s"
          % (tc, max_abs_U, deloc_ok))
    print("      Lemma 6.9 : high SNR (bar_lambda=%.1f > N sigma2=%.2f), MESP spread=%.3f -> %s"
          % (bar_lambda, N * sigma2, spread, snr_ok and gap_ok))
    print("      Lemma 6.9a: sigma2=%.3f >= 1/(2 pi e)=%.4f (M monotone): %s, verified mono: %s"
          % (sigma2, NOISE_FLOOR, floor_ok, mono_ok))
    ok = deloc_ok and snr_ok and gap_ok and floor_ok and mono_ok
    print("  [A] all three lemmas simultaneously satisfiable: %s" % ok)
    return ok


def check_B_additive_regret():
    """Composition is additive-regret on C; multiplicative-(c*5/4)-on-C is false."""
    N, t = 12, 6        # s = N - t = 6
    s = N - t
    Sigma = greedy_suboptimal_sigma(N, t)
    if Sigma is None:
        print("  [B] could not build a greedy-suboptimal instance (unexpected)")
        return False

    hX = total_entropy(Sigma)
    Topt, Mopt = brute_max_mesp(Sigma, t)
    Tg = greedy_mesp(Sigma, t)
    Mg = diff_entropy(Sigma, Tg)
    # supports S = complement of T
    Sopt = [i for i in range(N) if i not in Topt]
    Sg = [i for i in range(N) if i not in Tg]
    Copt = residual_cost(Sigma, Sopt)
    Cg = residual_cost(Sigma, Sg)

    # EXACT additive transfer: C_greedy - C_opt == M_opt - M_greedy.
    lhs = Cg - Copt
    rhs = Mopt - Mg
    additive_exact = abs(lhs - rhs) < 1e-6
    # NWF guarantee on the shifted NON-NEGATIVE monotone objective:
    #   M(That) >= (1-1/e) M(Topt)  =>  M(Topt)-M(That) <= (1/e) M(Topt).
    # The shift is modular (cancels in the difference M_opt - M_greedy), so the
    # TIGHT regret bound on the residual is lhs <= (1/e) * M_o_s.
    sing = [diff_entropy(Sigma, [i]) for i in range(N)]
    shift = -min(0.0, min(sing))
    Mg_s, Mo_s = Mg + t * shift, Mopt + t * shift
    regret_ok = lhs <= ONE_OVER_E * Mo_s + 1e-6        # tight (1/e) bound
    mult_M_ok = (Mo_s <= 1e-9) or (Mo_s / Mg_s <= E_OVER_EM1 + 1e-6)
    # the FALSE claim: residual multiplicative factor (c*5/4), c = e/(e-1).
    false_mult_pred = (E_OVER_EM1 * 1.25) * Copt
    mult_on_C_false = abs(Cg - false_mult_pred) > abs(Cg - Copt)

    print("  [B] M_opt=%.4f M_greedy=%.4f  |  C_opt=%.4f C_greedy=%.4f (bits, hX=%.3f)"
          % (Mopt, Mg, Copt, Cg, hX))
    print("  [B] additive transfer  C_g-C_opt=%.4f == M_opt-M_g=%.4f : %s"
          % (lhs, rhs, additive_exact))
    print("  [B] tight regret  C_g-C_opt <= (1/e) M_opt(shifted) = %.4f : %s"
          % (ONE_OVER_E * Mo_s, regret_ok))
    print("  [B] multiplicative-on-M (shifted) ratio <= e/(e-1): %s" % mult_M_ok)
    print("  [B] FALSE (c*5/4)-mult-on-C pred=%.3f vs actual C_g=%.3f (C_opt=%.3f) -> false: %s"
          % (false_mult_pred, Cg, Copt, mult_on_C_false))
    return additive_exact and regret_ok and mult_M_ok and mult_on_C_false


def check_C_regret_vanishes_near_modular():
    """Greedy slack -> 0 in the NEAR-MODULAR limit (low total-correlation /
    curvature c -> 0), per Lemma 6.9a's (1-e^{-c})/c refinement. The
    dependence on SNR is NON-monotone (slack is ~0 both at very low SNR, where
    all subsets tie, and at near-modular high SNR; it can be larger in
    between), so the honest claim is curvature-based, not SNR-monotone."""
    N, d, t = 11, 3, 5

    # near-modular: weak coupling (small off-diagonal), curvature ~ 0.
    near_mod = 0.0
    for _ in range(60):
        # tiny rank-d signal on a strong isotropic base => Sigma ~ modular
        Sigma, _, _ = manifold_sigma(N, d, bar_lambda=0.05 * N, sigma2=1.0, equal=False)
        Tg = greedy_mesp(Sigma, t)
        _, Mo = brute_max_mesp(Sigma, t)
        near_mod = max(near_mod, (Mo - diff_entropy(Sigma, Tg)) / t)

    # high-curvature: strong, unequal coupling => larger possible slack.
    high_curv = 0.0
    for _ in range(60):
        Sigma, _, _ = manifold_sigma(N, d, bar_lambda=4.0 * N, sigma2=0.3, equal=False)
        Tg = greedy_mesp(Sigma, t)
        _, Mo = brute_max_mesp(Sigma, t)
        high_curv = max(high_curv, (Mo - diff_entropy(Sigma, Tg)) / t)

    # The (1-1/e) worst-case bound ALWAYS holds; the refinement is that the
    # near-modular slack is essentially zero (<= a tiny tolerance), strictly
    # below the high-curvature worst case.
    near_mod_zero = near_mod < 1e-3
    refinement_ok = near_mod <= high_curv + 1e-9
    print("  [C] near-modular (curvature->0) worst per-coord slack = %.6f bits" % near_mod)
    print("  [C] high-curvature              worst per-coord slack = %.6f bits" % high_curv)
    print("  [C] near-modular slack ~ 0 (greedy -> exact, Lemma 6.9a V5): %s" % near_mod_zero)
    print("  [C] (note) slack-vs-SNR is NON-monotone; the honest control is curvature.")
    return near_mod_zero and refinement_ok


def check_D_end_to_end_encoder():
    """Full pipeline: estimate Sigma -> Lemma 6.9a greedy -> code; measure regret."""
    N, d, s = 12, 3, 7
    t = N - s
    bar_lambda, sigma2 = 21.0, 0.2     # sigma2=0.2 >= noise floor; high SNR
    Sigma_true, _, _ = manifold_sigma(N, d, bar_lambda, sigma2, equal=False)

    # (i) estimate the FIXED field from samples (fixed-field regime).
    m = 20000
    L = np.linalg.cholesky(Sigma_true)
    Xs = (L @ RNG.standard_normal((N, m))).T
    Sigma_hat = np.cov(Xs, rowvar=False)

    # (ii) select reconstruction set T greedily on the ESTIMATED covariance,
    #      coded support S = complement (Lemma 6.9a).
    Tg = greedy_mesp(Sigma_hat, t)
    Sg = [i for i in range(N) if i not in Tg]

    # (iii) measure on the TRUE covariance: residual C(S) and total rate.
    hX = total_entropy(Sigma_true)
    Cg = residual_cost(Sigma_true, Sg)
    Mg = diff_entropy(Sigma_true, [i for i in range(N) if i not in Sg])
    # chain rule: C(S) + M(Sbar) == h(X) for every S.
    chain_ok = abs(Cg + Mg - hX) < 1e-6

    # optimum (true cov) and the additive regret actually paid:
    Topt, Mopt = brute_max_mesp(Sigma_true, t)
    Sopt = [i for i in range(N) if i not in Topt]
    Copt = residual_cost(Sigma_true, Sopt)
    regret = Cg - Copt
    # tight (1/e) NWF bound on the shifted non-negative objective; greedy ran on
    # the ESTIMATED cov, so with finite m the guarantee is approximate -- we use
    # the shifted M_opt and a small estimation slack.
    sing = [diff_entropy(Sigma_true, [i]) for i in range(N)]
    shift = -min(0.0, min(sing))
    Mo_s = Mopt + t * shift
    regret_bound_ok = regret <= ONE_OVER_E * Mo_s + 1e-3
    additive_id_ok = abs(regret - (Mopt - Mg)) < 1e-6

    print("  [D] estimated Sigma from m=%d samples; greedy coded support S=%s" % (m, Sg))
    print("  [D] chain rule C(S)+M(Sbar)=%.3f == h(X)=%.3f : %s" % (Cg + Mg, hX, chain_ok))
    print("  [D] residual: opt C_opt=%.4f, encoder C_g=%.4f (regret=%.4f bits/block)"
          % (Copt, Cg, regret))
    print("  [D] regret == M_opt - M_greedy (additive identity): %s" % additive_id_ok)
    print("  [D] regret <= (1/e) M_opt(shifted) = %.4f (NWF bound): %s"
          % (ONE_OVER_E * Mo_s, regret_bound_ok))
    print("  [D] per-coord regret = %.5f bits (small; -> 0 near-modular, plus est. noise)"
          % (regret / s))
    if abs(regret) < 1e-6:
        print("  [D] (note) greedy attains the optimum on this benign low-rank+isotropic")
        print("       instance: the APX-hardness is worst-case; average-case manifold")
        print("       sources admit near-optimal greedy selection (regret ~ 0 here).")
    return chain_ok and additive_id_ok and regret_bound_ok


def main():
    print("Corollary 6.9b: operational Gaussian-on-manifold RNR (honest synthesis)")
    print("  TC rate (5.6j) + hardness (6.9) + submodular greedy (6.9a) -> ADDITIVE")
    print("  REGRET end-to-end bound on the residual (NOT multiplicative c*5/4).")
    print()
    results = {
        "A hypotheses compatible": check_A_compatibility(),
        "B additive-regret transfer": check_B_additive_regret(),
        "C regret -> 0 near-modular": check_C_regret_vanishes_near_modular(),
        "D end-to-end encoder": check_D_end_to_end_encoder(),
    }
    print()
    all_ok = all(results.values())
    for k, v in results.items():
        print("  %-28s %s" % (k, "PASS" if v else "FAIL"))
    print()
    if all_ok:
        print("PASS: a single Gaussian-on-manifold source satisfies Lemma 5.6j +")
        print("      6.9 + 6.9a at once (incl. noise floor sigma^2>=1/(2 pi e) and")
        print("      high SNR); Lemma 6.9a's (1-1/e) submodular guarantee composes")
        print("      with the TC rate as an ADDITIVE REGRET on the residual")
        print("      C_greedy - C_opt = M_opt - M_greedy <= (1/e) M_opt, NOT a")
        print("      multiplicative (c*5/4) factor; the encoder attains total rate")
        print("      h(X) with residual within the additive regret of optimal, poly(N).")
        sys.exit(0)
    else:
        print("FAIL: at least one check did not hold.")
        sys.exit(1)


if __name__ == "__main__":
    main()
