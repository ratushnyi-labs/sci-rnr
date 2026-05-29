r"""
Verification for Remark 6.9d: the Sigma >= I / noise-floor calibration of the
greedy e/(e-1) guarantee (Corollary 6.9b, Remark 6.9c) is NECESSARY, not
cosmetic. Below the entropy threshold (2 pi e lambda_min(Sigma) < 1) the RAW
log-det/entropy objective M(T) = 0.5 log2 det(2 pi e Sigma_T) can be negative,
and the Nemhauser-Wolsey-Fisher multiplicative bound M_greedy >= (1-1/e) M_opt
FAILS on raw M; it is restored only after the calibration shift to the
non-negative g.

This complements (does not duplicate) lemma_6_9c_typeiiic_threshold.py, which
calibrates Sigma >= I from the start and therefore never exhibits the failure
the calibration prevents. Here we exhibit precisely that failure.

CHECKS.
 (R1) The N=3 counterexample: raw bound FAILS, calibrated bound HOLDS.
 (R2) Across random sub-threshold Sigma, the raw (1-1/e) bound fails on a
      positive fraction of instances; the calibrated bound holds on all.
 (R3) The exact monotonicity condition is "all Var(X_i|X_T) >= 1/(2 pi e)";
      the proxy lambda_min(Sigma) >= 1/(2 pi e) is SUFFICIENT but NOT NECESSARY
      (a strongly-correlated Sigma below the eigenvalue threshold can still be
      monotone). Includes the explicit N=2 witness from the remark.
 (R4) The additive residual deficit C_greedy - C_opt = M_opt - M_greedy is
      shift-invariant (equal-size sets) and bounded by (1/e) g_opt on the
      CALIBRATED scale, but NOT by (1/e) M_opt on the raw scale.

PASS = R1..R4 all hold.
"""

import math
import sys
from itertools import combinations

import numpy as np

RNG = np.random.default_rng(20260529)
TWO_PI_E = 2.0 * math.pi * math.e
THR = 1.0 / TWO_PI_E  # ~0.05855


def slogdet2(M):
    sign, val = np.linalg.slogdet(M)
    assert sign > 0, "non-SPD submatrix"
    return val / math.log(2.0)


def Mraw(Sigma, T):
    """Raw differential entropy M(T) = 0.5 log2 det(2 pi e Sigma_T) (bits)."""
    T = list(T)
    if not T:
        return 0.0
    return 0.5 * slogdet2(TWO_PI_E * Sigma[np.ix_(T, T)])


def residual_cost(Sigma, S):
    """C(S) = h(X_S | X_Sbar) = 0.5 log2 det(2 pi e Sigma_{S|Sbar})."""
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


def greedy_max(Sigma, t, shift_c=0.0):
    """Greedy maximisation of g(T) = M(T) + |T|*shift_c over |T| = t.
    shift_c is the per-element calibration shift (0 = raw M). Returns T (list).
    The shift adds the same per-step constant to every candidate gain, so the
    greedy CHOICES are identical for any shift_c -- verified implicitly by R1."""
    N = Sigma.shape[0]
    T = []
    remaining = set(range(N))
    cur = 0.0
    for _ in range(t):
        best_i, best_gain = None, -math.inf
        for i in remaining:
            g = (Mraw(Sigma, T + [i]) + (len(T) + 1) * shift_c) - cur
            if g > best_gain:
                best_gain, best_i = g, i
        T.append(best_i)
        remaining.discard(best_i)
        cur = Mraw(Sigma, T) + len(T) * shift_c
    return T


def brute_max(Sigma, t):
    N = Sigma.shape[0]
    return max((Mraw(Sigma, T) for T in combinations(range(N), t)))


def manifold_sigma(N, d, lam_scale, sig2):
    U, _ = np.linalg.qr(RNG.standard_normal((N, d)))
    U = U[:, :d]
    lam = np.abs(RNG.standard_normal(d)) * lam_scale
    return U @ np.diag(lam) @ U.T + sig2 * np.eye(N)


def calib_shift(Sigma):
    """Per-element shift c so that lambda_min(alpha * 2 pi e Sigma) >= 1, i.e.
    g(T) = M(T) + |T| c is non-negative monotone. alpha = 1/lambda_min(2pe Sigma)."""
    lam_min_scaled = TWO_PI_E * float(np.min(np.linalg.eigvalsh(Sigma)))
    alpha = max(1.0, 1.0 / lam_min_scaled)
    return 0.5 * math.log2(alpha)


# ---------------------------------------------------------------------------
def check_R1_n3_counterexample():
    """The explicit N=3 sub-threshold counterexample from the remark."""
    A = np.array([[1.1, 0.0, 1.0 / math.sqrt(2)],
                  [0.0, 1.1, 1.0 / math.sqrt(2)],
                  [1.0 / math.sqrt(2), 1.0 / math.sqrt(2), 1.2]])
    Sigma = A / TWO_PI_E  # so that 2 pi e Sigma = A
    ev = np.linalg.eigvalsh(Sigma)
    sub_threshold = ev.min() < THR
    t = 2
    # raw greedy
    Tg = greedy_max(Sigma, t, shift_c=0.0)
    Mg = Mraw(Sigma, Tg)
    Mo = brute_max(Sigma, t)
    raw_fails = not (Mg >= (1.0 - 1.0 / math.e) * Mo - 1e-9)
    # calibrated greedy
    c = calib_shift(Sigma)
    Tgc = greedy_max(Sigma, t, shift_c=c)
    gg = Mraw(Sigma, Tgc) + t * c
    go = Mo + t * c   # opt of g over size-t == opt of M plus the constant shift
    calib_holds = gg >= (1.0 - 1.0 / math.e) * go - 1e-9
    # greedy choices identical under shift (R: shift preserves choices)
    same_choice = set(Tg) == set(Tgc)
    ok = sub_threshold and raw_fails and calib_holds and same_choice
    print("  [R1] N=3: eig(Sigma)=%s sub-threshold(lam_min<1/2pe)=%s"
          % (np.round(ev, 4).tolist(), sub_threshold))
    print("  [R1] RAW: M_greedy=%.4f  (1-1/e)M_opt=%.4f  -> raw bound FAILS=%s"
          % (Mg, (1 - 1 / math.e) * Mo, raw_fails))
    print("  [R1] CALIB: g_greedy=%.4f  (1-1/e)g_opt=%.4f -> calib bound HOLDS=%s"
          % (gg, (1 - 1 / math.e) * go, calib_holds))
    print("  [R1] shift preserves greedy choices (%s == %s): %s"
          % (sorted(Tg), sorted(Tgc), same_choice))
    return ok


def check_R2_random_subthreshold():
    """Raw bound fails on a positive fraction; calibrated holds on all."""
    raw_fail = 0
    calib_fail = 0
    n = 0
    for _ in range(400):
        N = int(RNG.integers(5, 9))
        t = int(RNG.integers(2, N - 1))
        d = int(RNG.integers(1, N))
        # sub-threshold: sig2 small so lambda_min < 1/(2 pi e); keep some signal
        Sigma = manifold_sigma(N, d, lam_scale=float(RNG.uniform(0.05, 0.6)),
                               sig2=float(RNG.uniform(0.005, 0.04)))
        if TWO_PI_E * np.min(np.linalg.eigvalsh(Sigma)) >= 1.0:
            continue  # only count genuine sub-threshold instances
        n += 1
        Tg = greedy_max(Sigma, t, shift_c=0.0)
        Mg, Mo = Mraw(Sigma, Tg), brute_max(Sigma, t)
        # raw bound only meaningful when M_opt > 0; failure = greedy below (1-1/e)M_opt
        if Mo > 1e-9 and Mg < (1.0 - 1.0 / math.e) * Mo - 1e-9:
            raw_fail += 1
        c = calib_shift(Sigma)
        Tgc = greedy_max(Sigma, t, shift_c=c)
        gg, go = Mraw(Sigma, Tgc) + t * c, Mo + t * c
        if go > 1e-9 and gg < (1.0 - 1.0 / math.e) * go - 1e-6:
            calib_fail += 1
    print("  [R2] sub-threshold instances tested: %d" % n)
    print("  [R2] RAW (1-1/e) bound failures: %d (>0 expected): %s" % (raw_fail, raw_fail > 0))
    print("  [R2] CALIBRATED bound failures:  %d (0 expected): %s" % (calib_fail, calib_fail == 0))
    return (raw_fail > 0) and (calib_fail == 0)


def check_R3_exact_condition_vs_proxy():
    """Conditional-variance condition is exact; lambda_min proxy sufficient not necessary."""
    # Explicit N=2 witness from the remark: lambda_min < 1/(2pe) but monotone.
    Sig = np.array([[100.0, 99.95], [99.95, 100.0]])
    lam_min = float(np.min(np.linalg.eigvalsh(Sig)))
    cond_var = Sig[1, 1] - Sig[1, 0] ** 2 / Sig[0, 0]  # Var(X_2|X_1)
    proxy_fails = lam_min < THR
    monotone = (0.5 * math.log2(TWO_PI_E * cond_var)) >= 0.0  # marginal gain >= 0
    witness_ok = proxy_fails and monotone and (cond_var >= THR)
    # Conversely: when ALL conditional variances >= 1/(2pe), M is monotone (random check)
    all_ok = True
    for _ in range(300):
        N = int(RNG.integers(4, 8))
        Sigma = manifold_sigma(N, int(RNG.integers(1, N)), 1.0, float(RNG.uniform(0.01, 1.0)))
        # check exact condition <=> monotone marginal gains
        prec = np.linalg.inv(Sigma)
        min_cv_full = min(1.0 / prec[i, i] for i in range(N))  # min over Var(X_i|rest)
        # monotone iff every marginal gain >= 0; test a sample of (i,T)
        worst_gain = math.inf
        for _ in range(8):
            i = int(RNG.integers(0, N))
            rest = [x for x in range(N) if x != i]
            RNG.shuffle(rest)
            T = rest[: int(RNG.integers(0, len(rest)))]
            worst_gain = min(worst_gain, Mraw(Sigma, T + [i]) - Mraw(Sigma, T))
        # if min conditional variance over the FULL conditioning >= THR, all gains >=0
        if min_cv_full >= THR and worst_gain < -1e-9:
            all_ok = False
    print("  [R3] N=2 witness: lam_min=%.4f (<1/2pe=%.4f: %s) but Var(X2|X1)=%.4f>=1/2pe -> monotone=%s"
          % (lam_min, THR, proxy_fails, cond_var, monotone))
    print("  [R3] => proxy lambda_min>=1/(2pe) is SUFFICIENT not NECESSARY: %s" % witness_ok)
    print("  [R3] when min Var(X_i|rest) >= 1/(2pe), sampled marginal gains all >=0: %s" % all_ok)
    return witness_ok and all_ok


def check_R4_additive_deficit_shift_invariant():
    """C_greedy - C_opt = M_opt - M_greedy (shift-invariant), bounded by g_opt/e."""
    ok_invariant = True
    ok_calib_bound = True
    # decisive raw additive-bound failure on the N=3 instance (from R1):
    A = np.array([[1.1, 0.0, 1.0 / math.sqrt(2)],
                  [0.0, 1.1, 1.0 / math.sqrt(2)],
                  [1.0 / math.sqrt(2), 1.0 / math.sqrt(2), 1.2]])
    Sg3 = A / TWO_PI_E
    Tg3 = greedy_max(Sg3, 2, shift_c=0.0)
    Mo3, Mg3 = brute_max(Sg3, 2), Mraw(Sg3, Tg3)
    # raw additive bound C_g-C_opt = M_opt-M_g <= M_opt/e ?  (should FAIL: 0.281 vs 0.050)
    ok_raw_bound_can_fail = (Mo3 > 0) and ((Mo3 - Mg3) > Mo3 / math.e + 1e-6)
    for _ in range(300):
        N = int(RNG.integers(5, 9))
        s = int(RNG.integers(1, N - 1))
        t = N - s
        d = int(RNG.integers(1, N))
        Sigma = manifold_sigma(N, d, lam_scale=float(RNG.uniform(0.1, 3.0)),
                               sig2=float(RNG.uniform(0.005, 0.5)))
        c = calib_shift(Sigma)
        Tgc = greedy_max(Sigma, t, shift_c=c)   # calibrated greedy choices
        Sg = [i for i in range(N) if i not in Tgc]
        Cg = residual_cost(Sigma, Sg)
        # brute optimum residual
        So = min(combinations(range(N), s), key=lambda S: residual_cost(Sigma, S))
        Co = residual_cost(Sigma, So)
        Mo = brute_max(Sigma, t)
        Mg = Mraw(Sigma, Tgc)
        # shift-invariance: residual deficit == entropy slack
        if abs((Cg - Co) - (Mo - Mg)) > 1e-6:
            ok_invariant = False
        # calibrated additive bound: slack <= g_opt/e
        go = Mo + t * c
        if (Mo - Mg) > go / math.e + 1e-6:
            ok_calib_bound = False
    print("  [R4] residual deficit == entropy slack (shift-invariant): %s" % ok_invariant)
    print("  [R4] calibrated additive bound slack <= g_opt/e holds:     %s" % ok_calib_bound)
    print("  [R4] raw additive bound slack <= M_opt/e VIOLATED on N=3 (%.3f > %.3f): %s"
          % (Mo3 - Mg3, Mo3 / math.e, ok_raw_bound_can_fail))
    return ok_invariant and ok_calib_bound and ok_raw_bound_can_fail


def main():
    print("Remark 6.9d: the Sigma>=I / noise-floor calibration of the greedy")
    print("             e/(e-1) guarantee is NECESSARY (raw bound fails below threshold)")
    print()
    results = {
        "R1 N=3 counterexample": check_R1_n3_counterexample(),
        "R2 random sub-threshold": check_R2_random_subthreshold(),
        "R3 exact cond vs proxy": check_R3_exact_condition_vs_proxy(),
        "R4 additive deficit": check_R4_additive_deficit_shift_invariant(),
    }
    print()
    all_ok = all(results.values())
    for k, v in results.items():
        print("  %-26s %s" % (k, "PASS" if v else "FAIL"))
    print()
    if all_ok:
        print("PASS: below the entropy threshold the RAW greedy e/(e-1) bound FAILS")
        print("      (explicit N=3, and a positive fraction of random instances); the")
        print("      calibration to the non-negative g restores it. The monotonicity")
        print("      proxy lambda_min>=1/(2 pi e) is sufficient not necessary (exact")
        print("      condition: all conditional variances >= 1/(2 pi e)). The additive")
        print("      residual deficit is shift-invariant and bounded by g_opt/e on the")
        print("      calibrated scale, NOT by M_opt/e on the raw scale. Hence the")
        print("      calibration hypothesis of Corollary 6.9b / Remark 6.9c is")
        print("      load-bearing and cannot be silently dropped.")
        sys.exit(0)
    else:
        print("FAIL: at least one check did not hold.")
        sys.exit(1)


if __name__ == "__main__":
    main()
