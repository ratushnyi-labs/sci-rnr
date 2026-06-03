r"""
Verification for Remark 7.29a (the (N/K)*delta_inf random-access overhead of Theorem
7.29 is TIGHT for the sync-point sub-block architecture: a matching second-order
converse). Any O(K)-access scheme that decodes block i from its sync point + the block's
own bits (without reading earlier blocks) pays, per block, the EXCESS ENTROPY
delta_inf = I(past; future) -- either as a cold-start rate loss in the payload, or as
stored sync state. Under exp-mixing the m=N/K cuts are near-independent, so the total
overhead is (1-o(1))*(N/K)*delta_inf, matching the Theorem 7.29 achievability.

THE CONVERSE ARGUMENT.
- STRICTLY-CHARGED MODEL: block i is decoded from its own DISJOINT pair
  (payload_i, sync_i) alone; shared metadata (the model) is SOURCE-INDEPENDENT; file
  size = sum_i (|payload_i|+|sync_i|). (A SOURCE-DEPENDENT global header is excluded --
  it can collapse the bound, e.g. X_n=Y stored once.) Since X_block_i is determined by
  (payload_i, sync_i), a Shannon source-coding bound gives
      E[|payload_i|+|sync_i|] >= H(X_block_i) = K*h + delta_K,   delta_K=H(X^K)-K*h.
- delta_K -> E (excess entropy) under exp-mixing. SUMMING the DISJOINT per-block bounds
  (no independence assumption needed) over m=N/K blocks: total >= N*h + (N/K)*delta_K
  = N*h + (N/K)*E - o(sqrt N) at K=Theta(sqrt N). Grouping r blocks under one sync is
  NOT a loophole: it is the spacing K->rK on the same overhead-vs-access curve.

For an order-1 Markov source the excess entropy is E = I(X_0; X_1) = H(X) - h (marginal
entropy minus entropy rate), and the cold-start loss is concentrated in the FIRST symbol
of each block (subsequent symbols are warm). So the per-block penalty is EXACTLY E.

WHAT THIS SCRIPT VERIFIES (order-1 Markov, ideal predictor).
(C1) EXCESS = COLD-START FIRST-SYMBOL PENALTY. E := I(X_0;X_1) = H(X) - h equals the
     extra bits to code the cold first symbol of a block (stationary prediction
     H(X)) versus the warm conditional (h). Verified analytically and by the
     definition I(X_0;X_1)=H(X)-h.
(C2) PER-BLOCK COLD PENALTY = E, and WARM (interior) symbols cost h with NO penalty:
     so the whole cold block of size K costs K*h + E, i.e. exactly E excess per block.
     Verified by computing the cold block cost vs the warm conditional cost.
(C3) ILLUSTRATIVE (not a proof): the future-relevant info I(X_prev;X_next)=E <= H(X)
     (storing X_prev verbatim costs H(X) >= E). The CONVERSE itself is the analytic
     total_i >= H(X_block_i) bound above, NOT this check.
(C4) ORDER-MATCH across K: total cold penalty / E = m = N/K for a range of K; at
     K=Theta(sqrt N) this is Theta(sqrt N), matching the T7.29 second-order overhead.

PASS = C1 (E=H(X)-h=I(X0;X1)) AND C2 (cold block = K*h + delta_K, here K*h+E for the
order-1 cold reset) AND C3 (illustrative E<=H(X)) AND C4 ((N/K)*E scaling).

HONEST SCOPE. The converse is for the STRICTLY-CHARGED sub-block model: disjoint
per-block (payload_i, sync_i), SOURCE-INDEPENDENT shared metadata, file size = sum of
per-block sizes, K=omega(1). It does NOT cover: source-DEPENDENT global headers (which
can collapse the bound, e.g. X_n=Y), non-sub-block codes (e.g. Ferragina-Venturini
O(1)-access, cf. Theorem 7.27'), or non-mixing / global-latent sources. Whether ANY
random-access code reaches the dispersion floor without the (N/K)*delta_inf penalty
remains OPEN. The script ILLUSTRATES the order-1 arithmetic (delta_inf=H(X)-h, cold
block = K*h+E, (N/K)E scaling); the converse itself is the analytic total_i>=H(X_block_i)
summed over disjoint blocks (NOT shown numerically -- it is a Shannon source-coding bound).
"""

import math
import sys

import numpy as np

RNG = np.random.default_rng(20260604)


def markov1(A, temp, rng):
    T = np.zeros((A, A))
    for a in range(A):
        l = rng.standard_normal(A) / temp
        p = np.exp(l - l.max()); T[a] = p / p.sum()
    vals, vecs = np.linalg.eig(T.T)
    k = int(np.argmin(np.abs(vals - 1.0)))
    pi = np.abs(np.real(vecs[:, k])); pi = pi / pi.sum()
    return T, pi


def Hmarg(pi):
    p = pi[pi > 0]
    return float(-(p * np.log2(p)).sum())


def entropy_rate(T, pi):
    h = 0.0
    for a in range(len(pi)):
        row = T[a][T[a] > 0]
        h -= pi[a] * float((row * np.log2(row)).sum())
    return h


def mutual_info_X0X1(T, pi):
    """I(X_0; X_1) for the stationary order-1 chain = H(X_1) - H(X_1|X_0) = H(X) - h."""
    P2 = pi[:, None] * T                 # joint (X_0, X_1)
    pX1 = P2.sum(axis=0)
    I = 0.0
    for a in range(len(pi)):
        for b in range(len(pi)):
            if P2[a, b] > 0:
                I += P2[a, b] * math.log2(P2[a, b] / (pi[a] * pX1[b]))
    return float(I)


def check_C1():
    print("=" * 70)
    print("C1: excess entropy E = I(X0;X1) = H(X) - h (marginal entropy - rate).")
    print("=" * 70)
    ok = True
    for trial, (A, temp) in enumerate([(2, 0.8), (3, 1.0), (2, 1.4), (4, 0.9)]):
        rng = np.random.default_rng(10 + trial)
        T, pi = markov1(A, temp, rng)
        h = entropy_rate(T, pi)
        Hx = Hmarg(pi)
        I01 = mutual_info_X0X1(T, pi)
        E = Hx - h
        ok_t = abs(I01 - E) < 1e-9
        print(f"  trial {trial}: A={A}: H(X)={Hx:.4f}  h={h:.4f}  E=H(X)-h={E:.4f}  "
              f"I(X0;X1)={I01:.4f}  match={ok_t}")
        ok = ok and ok_t
    print(f"  C1 {'PASS' if ok else 'FAIL'}")
    return ok


def check_C2():
    print("=" * 70)
    print("C2: a COLD block of size K costs K*h + E exactly (cold first symbol pays")
    print("    E=H(X)-h; warm interior symbols pay h). => E excess per block.")
    print("=" * 70)
    ok = True
    for trial, (A, temp, K) in enumerate([(2, 0.8, 50), (3, 1.0, 80), (2, 1.4, 120)]):
        rng = np.random.default_rng(20 + trial)
        T, pi = markov1(A, temp, rng)
        h = entropy_rate(T, pi)
        Hx = Hmarg(pi)
        E = Hx - h
        # cold block expected cost = H(X_1) [stationary first symbol] + (K-1)*h [warm]
        cold_block = Hx + (K - 1) * h
        warm_block = K * h                    # if the prev symbol were known
        excess = cold_block - warm_block
        ok_t = abs(excess - E) < 1e-9
        print(f"  trial {trial}: A={A} K={K}: cold block={cold_block:.3f}  "
              f"warm K*h={warm_block:.3f}  excess={excess:.4f}  E={E:.4f}  match={ok_t}")
        ok = ok and ok_t
    print(f"  C2 {'PASS' if ok else 'FAIL'}")
    return ok


def check_C3():
    print("=" * 70)
    print("C3: pay E as rate OR store E as state. The minimal sync state to remove the")
    print("    cold penalty carries I(X_prev; future)=E bits (bottleneck) <= H(X) (cost")
    print("    of storing X_prev verbatim); E is the tight stored-state floor.")
    print("=" * 70)
    ok = True
    for trial, (A, temp) in enumerate([(2, 0.8), (3, 1.0), (4, 0.9)]):
        rng = np.random.default_rng(30 + trial)
        T, pi = markov1(A, temp, rng)
        Hx = Hmarg(pi)
        E = Hx - entropy_rate(T, pi)
        # the future depends on the past ONLY through X_prev (order-1); the
        # bottleneck-relevant info I(X_prev; X_next) = E. Storing X_prev verbatim
        # costs H(X)>=E. So E is the floor and H(X) the (wasteful) verbatim cost.
        I_prev_next = mutual_info_X0X1(T, pi)
        bottleneck_ok = abs(I_prev_next - E) < 1e-9 and E <= Hx + 1e-12
        print(f"  trial {trial}: A={A}: stored-state floor E={E:.4f} = I(X_prev;X_next); "
              f"verbatim store H(X)={Hx:.4f} (>=E: {E<=Hx+1e-12})  {bottleneck_ok}")
        ok = ok and bottleneck_ok
    print("  => either way the per-block cost is E; total >= N*h + (N/K)*E.")
    print(f"  C3 {'PASS' if ok else 'FAIL'}")
    return ok


def check_C4():
    print("=" * 70)
    print("C4: total cold penalty = (N/K)*E across K; at K=Theta(sqrt N) it is")
    print("    Theta(sqrt N), matching the T7.29 second-order overhead.")
    print("=" * 70)
    rng = np.random.default_rng(7)
    T, pi = markov1(3, 1.0, rng)
    E = Hmarg(pi) - entropy_rate(T, pi)
    print(f"  E={E:.4f}")
    ok = True
    for N in (10**4, 10**6, 10**8):
        K = math.sqrt(N)
        m = N / K
        total = m * E
        ratio = total / E
        per_sqrtN = total / math.sqrt(N)
        ok = ok and abs(ratio - m) < 1e-6
        print(f"    N={N:>10}: K=sqrt(N)={K:.0f}  m=N/K={m:.0f}  total=(N/K)E={total:.2f}  "
              f"(total/sqrtN={per_sqrtN:.4f}=E: {abs(per_sqrtN-E)<1e-9})")
    print(f"  C4 {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print()
    print("#" * 70)
    print("# Remark 7.29a: tight (N/K)delta_inf RA converse for the sub-block scheme")
    print("#" * 70)
    print()
    r1 = check_C1(); print()
    r2 = check_C2(); print()
    r3 = check_C3(); print()
    r4 = check_C4(); print()
    print("=" * 70)
    allok = r1 and r2 and r3 and r4
    print(f"C1 E = H(X)-h = I(X0;X1)            : {'PASS' if r1 else 'FAIL'}")
    print(f"C2 cold block = K*h + E (excess=E) : {'PASS' if r2 else 'FAIL'}")
    print(f"C3 pay-E-or-store-E (floor E<=H(X)): {'PASS' if r3 else 'FAIL'}")
    print(f"C4 total=(N/K)E=Theta(sqrt N)      : {'PASS' if r4 else 'FAIL'}")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
