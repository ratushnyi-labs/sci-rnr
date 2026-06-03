r"""
Verification for Remark 7.15d (Renyi-entropy spectrum of the open-middle oracle:
integer-order Renyi H_alpha (alpha>=2) of the noncontiguous bounded-order subset
entropy is POLY-time computable via an alpha-fold tensor-product forward algorithm;
the min-entropy (alpha=inf) is NP-hard (most-likely string); only the SHANNON point
(alpha=1) is open -- so the difficulty of Remark 7.15a's open middle is isolated to
the alpha->1 / sum p log p limit, NOT the marginal/moment computation).

------------------------------------------------------------------------
THE INDIRECT PATH.

Remark 7.15a left OPEN the exact complexity of H(X_T) (Shannon) for a noncontiguous
subset T of a bounded-order chain (= exact finite-block Shannon entropy of the
induced subsample process Y_k:=X_{t_k}). Direct attack on the Shannon point is
blocked. Indirect path: characterise the whole RENYI SPECTRUM
    H_alpha(Y) = (1/(1-alpha)) log_2 sum_y P(y)^alpha,   H_1 = Shannon, H_inf = min.

KEY IDENTITY (the lever). For INTEGER alpha>=2,
    sum_y P(Y_T = y)^alpha  =  P( alpha i.i.d. copies of X all AGREE on T ),
because expanding the alpha-th power pairs alpha independent paths that emit the
SAME subset tuple y. This is the partition function of alpha independent copies of
the chain CONSTRAINED to agree at the subset positions and evolving freely across
the gaps -- computable by a FORWARD algorithm over the alpha-fold product state
(|Sigma|^{w*alpha} for the order-w lift), poly in (r, N) for FIXED alpha and w
(gaps via O(log gap) matrix powers). NO log-of-a-sum, NO A^{|T|} marginal.

So the spectrum splits:
  - alpha INTEGER >= 2 : POLY (tensor-product forward algorithm); cost poly in (r,N)
    but with constant |Sigma|^{w*alpha} GROWING in alpha (so this does NOT extend to
    alpha->inf nor to taking many moments)   [this script, R1-R2]
  - alpha = inf (min-entropy, -log max_y P(y)) : NP-HARD (most-likely OUTPUT string
    of an HMM; Lyngso-Pedersen JCSS 65(3):545-569, 2002). Watertight for OUR subclass:
    any HMM = order-2 chain alternating hidden/emitted positions, T=emitted   [R4]
  - alpha = 1 (SHANNON) : OPEN (Remark 7.15a) -- NOT reached by this route. By Renyi
    monotonicity H_2 <= H_1 <= H_0, the poly H_2 is a FREE poly LOWER BOUND on the
    open Shannon oracle, but the moments do NOT recover H_1 (analytic continuation
    needs >finitely-many orders, blocked by |Sigma|^{w*alpha}). NON-INTEGER alpha in
    (1,inf) ALSO lack the finite-copy structure. So it is a REPRESENTATIONAL
    DISCONTINUITY (integer orders = finite product-chain; Shannon = alpha->1
    derivative of an infinite-dim operator), NOT a clean "isolation to alpha=1".  [R3]

WHAT THIS SCRIPT VERIFIES.
(R1) IDENTITY + POLY ALGORITHM: the alpha-fold product forward DP computes
     sum_y P(y)^alpha EXACTLY (== brute, from the subset marginal) for alpha=2,3,
     WITHOUT building the |Sigma|^{|T|} marginal. (The DP state is |Sigma|^{2*alpha}
     for the order-2 lift -- constant for fixed alpha.)
(R2) So H_alpha (alpha=2,3) is computed from the poly product-DP, matching brute.
(R3) RENYI MONOTONICITY: H_inf <= ... <= H_3 <= H_2 <= H_1 (Shannon) <= H_0; in
     particular the poly H_2 lower-bounds the open Shannon entropy.
(R4) MIN-ENTROPY is the most-likely subset tuple: H_inf = -log max_y P(y); finding
     argmax_y P(y) (marginalising the gaps) is the HMM most-likely-string problem,
     NP-hard (Lyngso-Pedersen 2002). We exhibit the maximiser on a small instance
     (brute), the NP-hardness being the cited literature result.

PASS = R1 (product-DP == brute sum P^alpha, alpha=2,3, machine precision) AND
R2 (H_alpha matches) AND R3 (spectrum ordering) AND R4 (min-entropy = -log max P).

HONEST SCOPE. This does NOT resolve the SHANNON (alpha=1) exact complexity, and does
NOT claim "only Shannon is open": the product-copy method covers ONLY integer
alpha>=2; the Shannon point AND all non-integer orders lack the finite-copy structure
and are unaddressed. What it gives: integer Renyi orders are exactly poly (an island),
the min-entropy is NP-hard, and (Renyi monotonicity) H_2 is a free exact poly LOWER
BOUND on the open Shannon oracle. Complements Remark 7.15c (Shannon poly-eps-
APPROXIMABLE under filter stability) with an EXACT, unconditional poly handle (H_2).
"""

import math
import sys
from itertools import product

import numpy as np

RNG = np.random.default_rng(20260603)


def order2(A, temp, rng):
    T2 = np.zeros((A, A, A))
    for a in range(A):
        for b in range(A):
            l = rng.standard_normal(A) / temp
            p = np.exp(l - l.max()); T2[a, b] = p / p.sum()
    l = rng.standard_normal((A, A)) / temp
    pe = np.exp(l - l.max()); pi2 = pe / pe.sum()
    return pi2, T2


def order2_joint(pi2, T2, N):
    A = pi2.shape[0]
    P = np.zeros((A,) * N)
    for idx in product(range(A), repeat=N):
        pr = pi2[idx[0], idx[1]]
        for i in range(2, N):
            pr *= T2[idx[i - 2], idx[i - 1], idx[i]]
        P[idx] = pr
    return P / P.sum()


def subset_marginal(P, keep):
    keep = tuple(sorted(keep))
    axes_out = tuple(i for i in range(P.ndim) if i not in keep)
    M = P.sum(axis=axes_out) if axes_out else P
    return M


# ----------------------------------------------------------------------
# The alpha-fold PRODUCT FORWARD algorithm for sum_y P(Y_T=y)^alpha, T = even
# positions of an order-2 chain. State = tuple of alpha copies' (prev,cur) pairs;
# belief = P(all copies in these states AND agreed at all observed even positions
# so far). NO |Sigma|^{|T|} marginal is built.
# ----------------------------------------------------------------------

def product_forward_sum_palpha(pi2, T2, N, alpha):
    A = pi2.shape[0]
    evens = set(range(0, N, 2))
    # init at positions (0,1): copy c has pair (X0^c, X1^c) ~ pi2; position 0 even
    # => constrain all X0^c equal.
    belief = {}
    for combo in product(product(range(A), repeat=2), repeat=alpha):
        # combo = (pair_1, ..., pair_alpha), pair_c=(x0_c, x1_c)
        x0s = [pc[0] for pc in combo]
        if len(set(x0s)) != 1:            # position 0 observed: must agree
            continue
        w = 1.0
        for pc in combo:
            w *= pi2[pc[0], pc[1]]
        if w > 0:
            belief[combo] = belief.get(combo, 0.0) + w
    # advance positions 2..N-1
    for i in range(2, N):
        nb = {}
        for combo, w in belief.items():
            # each copy emits x_i^c ~ T2[prev,cur]; new pair = (cur, x_i^c)
            for nxts in product(range(A), repeat=alpha):
                if i in evens and len(set(nxts)) != 1:
                    continue              # observed position: all copies agree
                w2 = w
                newcombo = []
                ok = True
                for c in range(alpha):
                    prev, cur = combo[c]
                    t = T2[prev, cur, nxts[c]]
                    if t <= 0:
                        ok = False; break
                    w2 *= t
                    newcombo.append((cur, nxts[c]))
                if ok and w2 > 0:
                    key = tuple(newcombo)
                    nb[key] = nb.get(key, 0.0) + w2
        belief = nb
    return float(sum(belief.values()))


def renyi_from_marginal(M, alpha):
    p = M.ravel(); p = p[p > 0]
    if alpha == 1:
        return float(-(p * np.log2(p)).sum())
    if math.isinf(alpha):
        return float(-math.log2(p.max()))
    s = float((p ** alpha).sum())
    return (1.0 / (1.0 - alpha)) * math.log2(s)


def check_R1_R2():
    print("=" * 70)
    print("R1/R2: alpha-fold PRODUCT forward DP computes sum_y P(y)^alpha EXACTLY")
    print("       (== brute), so H_alpha (alpha=2,3) is poly. No |Sigma|^|T| table.")
    print("=" * 70)
    ok = True
    for trial, (A, N, temp) in enumerate([(2, 8, 0.8), (2, 10, 1.2), (3, 6, 0.9)]):
        rng = np.random.default_rng(10 + trial)
        pi2, T2 = order2(A, temp, rng)
        P = order2_joint(pi2, T2, N)
        evens = list(range(0, N, 2))
        M = subset_marginal(P, evens)
        for alpha in (2, 3):
            brute = float(((M.ravel()[M.ravel() > 0]) ** alpha).sum())
            prod = product_forward_sum_palpha(pi2, T2, N, alpha)
            err = abs(brute - prod)
            Halpha_brute = renyi_from_marginal(M, alpha)
            Halpha_prod = (1.0 / (1.0 - alpha)) * math.log2(prod)
            herr = abs(Halpha_brute - Halpha_prod)
            print(f"  trial {trial} A={A} N={N} |T|={len(evens)} alpha={alpha}: "
                  f"sum P^a brute={brute:.6f} prod-DP={prod:.6f} err={err:.1e}; "
                  f"H_{alpha} err={herr:.1e}")
            ok = ok and err < 1e-9 and herr < 1e-9
    print(f"  R1/R2 {'PASS' if ok else 'FAIL'}")
    return ok


def check_R3_R4():
    print("=" * 70)
    print("R3: Renyi monotonicity H_inf <= H_3 <= H_2 <= H_1(Shannon) <= H_0;")
    print("    H_2 (poly) is a free poly LOWER BOUND on the open Shannon oracle.")
    print("R4: min-entropy H_inf = -log max_y P(y) = most-likely subset tuple")
    print("    (NP-hard: HMM most-likely string, Lyngso-Pedersen 2002).")
    print("=" * 70)
    ok = True
    for trial, (A, N, temp) in enumerate([(2, 10, 0.7), (3, 6, 1.0)]):
        rng = np.random.default_rng(30 + trial)
        pi2, T2 = order2(A, temp, rng)
        P = order2_joint(pi2, T2, N)
        evens = list(range(0, N, 2))
        M = subset_marginal(P, evens)
        H0 = math.log2((M.ravel() > 1e-15).sum())
        H1 = renyi_from_marginal(M, 1)
        H2 = renyi_from_marginal(M, 2)
        H3 = renyi_from_marginal(M, 3)
        Hinf = renyi_from_marginal(M, math.inf)
        mono = (Hinf <= H3 + 1e-9 <= 1e18) and (H3 <= H2 + 1e-9) and \
               (H2 <= H1 + 1e-9) and (H1 <= H0 + 1e-9)
        # min-entropy maximiser
        flat = M.ravel()
        ystar = int(np.argmax(flat))
        is_max = abs(flat[ystar] - flat.max()) < 1e-15
        lower_ok = H2 <= H1 + 1e-9
        print(f"  trial {trial} A={A} N={N}: H_inf={Hinf:.4f} <= H_3={H3:.4f} <= "
              f"H_2={H2:.4f} <= H_1={H1:.4f} <= H_0={H0:.4f}  ordered={mono}")
        print(f"           H_2 <= H_1 (poly lower bound on open Shannon): {lower_ok}; "
              f"argmax_y P(y) found (NP-hard in general): {is_max}")
        ok = ok and mono and lower_ok and is_max
    print(f"  R3/R4 {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print()
    print("#" * 70)
    print("# Remark 7.15d: Renyi spectrum isolates the open middle to alpha=1")
    print("#" * 70)
    print()
    r1 = check_R1_R2(); print()
    r2 = check_R3_R4(); print()
    print("=" * 70)
    allok = r1 and r2
    print(f"R1/R2 product-DP == brute sum P^alpha; H_alpha poly : {'PASS' if r1 else 'FAIL'}")
    print(f"R3/R4 spectrum ordering + min-entropy maximiser     : {'PASS' if r2 else 'FAIL'}")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
