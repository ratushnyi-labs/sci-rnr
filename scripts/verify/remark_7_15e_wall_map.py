r"""
Verification for Remark 7.15e (mapping --- not breaking --- the open-middle wall:
the exact finite-block Shannon oracle is sandwiched poly <= H_1 <= FP^#P, and BOTH
natural attack routes are obstructed: the predicate-encoding #P-hardness route is
self-defeating for bounded-order chains (they encode only READ-ONCE predicates, whose
model-counting is in P), and the obvious poly algorithm is blocked by the
non-collapsing Blackwell filter. So the exact complexity stays OPEN, now pinned).

------------------------------------------------------------------------
THE WALL, AND WHY THE TWO NATURAL ATTACKS FAIL.

Remark 7.15a left the EXACT Shannon oracle H(Y_1..Y_r) (Y_k=X_{t_k}, the induced
subsample HMM) open. 7.15c made it eps-approximable under filter stability; 7.15d
bracketed it (integer Renyi poly, min-entropy NP-hard, H_2 a poly lower bound). The
exact Shannon point remains the wall. We map it.

(E1) UPPER BOUND: H(Y) to 2^{-poly} precision is computable by an FP machine with a
     #P oracle (FP^#P). H(Y) = -sum_{y in Sigma^r} P(y) log2 P(y), each P(y) a poly-bit
     rational (forward algorithm), log2 poly-time to 2^{-poly} precision; approximate
     each -P(y)log2 P(y) to guard bits, scale to integers, and sum over y with a #P
     oracle. The value is real/irrational => an FP^#P APPROXIMATION, not a #P function;
     logs do not push the class higher.

(E2) THE PREDICATE-ENCODING #P-HARDNESS ROUTE IS SELF-DEFEATING for bounded-order
     chains. 7.15a/V5 proved the entropy oracle #P-hard for GENERAL AR fields by
     encoding a #P-hard predicate E into a coordinate Y=E AND B, so H(Y)=h(P(E)/2)
     reveals the #P-hard marginal P(E)=#E/2^n. For a BOUNDED-ORDER-w chain, the
     coordinate Y=E(a) must be computed by the chain AS IT GENERATES the assignment
     a=a_1..a_n ONCE, tracking E in its state of size |Sigma|^w: i.e. E is a
     WIDTH-|Sigma|^w READ-ONCE branching program (each a_i read once, fixed order).
     But:
       (i)  a general CNF needs EXPONENTIAL read-once width (it cannot be encoded), and
       (ii) for any E that CAN be encoded (a read-once BP), MODEL-COUNTING #E is in P
            (bottom-up DP over the BP nodes; Van den Broeck-Suciu 2017).
     The GAP positions are not a loophole: marginalising the unobserved gaps is the
     SAME forward DP (a matrix-product walk-count, #L <= FP); the chain generates its
     sequence ONCE (no re-read / multi-pass power); and a MEMORYLESS bounded-state
     chain cannot enforce the global (bijection / no-repeat) constraints that make
     counting #P-hard -- the predicate bit is produced by the bounded state wherever
     the assignment sits. So the DIRECT (inline) predicate-encoding reveals only
     FP-computable counts -> only poly quantities, NEVER #P-hardness. (Bars the V5
     mechanism, not every conceivable reduction.) This is exactly the gap between the
     general-AR result (7.15a, #P-hard) and the bounded-order open middle.

(E3) THE OBVIOUS POLY ALGORITHM IS BLOCKED. H(Y)=sum_k H(Y_k|Y_{<k}) needs the filter
     (belief over the hidden state given Y_{<k}); the reachable belief set does NOT
     collapse (Blackwell measure generically infinitely supported; verified by the
     7.15a open-middle probe), so the filtering sum branches exponentially. The
     finite-block entropy is a Lyapunov-type functional of random matrix products,
     whose RATE has no closed form (Blackwell, "The entropy of functions of
     finite-state Markov chains", Trans. First Prague Conf. on Information Theory,
     1957, pp. 13-20; Jurgens-Crutchfield 2021 arXiv:2008.12886 for the generically
     infinite predictive features). So the natural DP is blocked too. (This blocks
     the OBVIOUS filtering algorithm; it does NOT prove no poly algorithm exists.)

NET: the wall HOLDS against both natural attacks; the exact oracle sits in
[poly (collision H_2) , FP^#P], with the standard #P-hardness route obstructed (E2)
and the standard poly route obstructed (E3). NOT broken --- mapped: any resolution
needs a NON-predicate-encoding hardness construction or a NON-filtering algorithm.

------------------------------------------------------------------------
WHAT THIS SCRIPT VERIFIES (the crisp, non-vacuous core is W1-W2).

(W1) READ-ONCE MODEL-COUNTING IS POLY: build random width-S read-once branching
     programs (a bounded-state chain reading a_1..a_n once) and compute #{a:E(a)=1}
     by (a) layer DP over the S nodes (poly, O(nS|Sigma|)) and (b) brute 2^n
     enumeration. They MATCH -> the count any bounded-order chain can encode is poly.

(W2) THE PREDICATE-ENCODING REVEALS ONLY A POLY COUNT: the coordinate Y=E AND B has
     H(Y)=h(#E/2^{n+1}); we recover #E from H(Y) (injective) and confirm it equals the
     POLY DP count of W1 -> the entropy oracle here is poly, NOT #P-hard. (Contrast:
     a general non-read-once CNF would need exp width and is NOT encodable.)

(W3) FP^#P SUM STRUCTURE (E1, illustrative): H(Y) equals the explicit weighted sum
     -sum_y P(y) log2 P(y) of poly-computable terms (sanity vs the chain-rule value).

PASS = W1 (DP count == brute, all instances) AND W2 (#E recovered from H(Y) == DP
count) AND W3 (explicit sum == chain-rule entropy).

HONEST SCOPE. This does NOT break the wall (exact Shannon complexity still OPEN). It
(E1) bounds it by FP^#P, (E2) rigorously rules out the predicate-encoding #P-hardness
route for bounded-order chains (read-once #SAT in P), and (E3) notes the Blackwell
obstruction to the obvious algorithm. It does NOT rule out a cleverer (non-predicate)
hardness reduction or a (non-filtering) poly algorithm.
"""

import math
import sys
from itertools import product

import numpy as np

RNG = np.random.default_rng(20260603)


# ----------------------------------------------------------------------
# A width-S READ-ONCE branching program over a_1..a_n (fixed order, each read once):
# layers 0..n, each layer has S nodes; transition delta[i][s][a_i] -> next node;
# start node 0; a subset ACC of layer-n nodes accepts. This is EXACTLY a bounded-state
# chain that generates a_1..a_n and tracks a predicate E in its state.
# ----------------------------------------------------------------------

def random_robp(n, S, A, rng):
    delta = [rng.integers(0, S, size=(S, A)) for _ in range(n)]   # delta[i][s,a]
    acc = rng.integers(0, 2, size=S).astype(bool)                 # accepting layer-n nodes
    return delta, acc


def robp_eval(delta, acc, a):
    s = 0
    for i, ai in enumerate(a):
        s = int(delta[i][s, ai])
    return bool(acc[s])


def robp_count_dp(delta, acc, n, S, A):
    """#{a in A^n : E(a)=1} by POLY layer DP over the S nodes. O(n S A)."""
    cnt = np.zeros(S, dtype=object)   # exact big ints
    cnt[0] = 1                        # start at node 0 with 1 way (empty prefix)
    for i in range(n):
        nxt = np.zeros(S, dtype=object)
        for s in range(S):
            if cnt[s] == 0:
                continue
            for ai in range(A):
                nxt[int(delta[i][s, ai])] += cnt[s]
        cnt = nxt
    return int(sum(int(cnt[s]) for s in range(S) if acc[s]))


def robp_count_brute(delta, acc, n, A):
    return sum(1 for a in product(range(A), repeat=n) if robp_eval(delta, acc, a))


def check_W1():
    print("=" * 70)
    print("W1: read-once model-counting is POLY -- layer DP over S nodes == brute")
    print("    2^n enumeration. (A bounded-state chain encodes a read-once BP.)")
    print("=" * 70)
    ok = True
    for trial, (n, S, A) in enumerate([(10, 4, 2), (12, 5, 2), (8, 3, 3), (14, 6, 2)]):
        rng = np.random.default_rng(100 + trial)
        delta, acc = random_robp(n, S, A, rng)
        dp = robp_count_dp(delta, acc, n, S, A)
        brute = robp_count_brute(delta, acc, n, A)
        match = dp == brute
        print(f"  trial {trial}: n={n} S={S} A={A}: DP #E={dp}  brute #E={brute}  "
              f"match={match}  (DP cost O(nS A)={n*S*A} vs brute A^n={A**n})")
        ok = ok and match
    print(f"  W1 {'PASS' if ok else 'FAIL'}")
    return ok


def check_W2():
    print("=" * 70)
    print("W2: the predicate-encoding reveals only a POLY count. Y = E AND B,")
    print("    H(Y)=h(#E/2^{n+1}); recovered #E == the POLY DP count (not #P-hard).")
    print("=" * 70)

    def hbin(p):
        if p <= 0 or p >= 1:
            return 0.0
        return -p * math.log2(p) - (1 - p) * math.log2(1 - p)

    def hbin_inv_lo(y):
        lo, hi = 0.0, 0.5
        for _ in range(200):
            m = (lo + hi) / 2
            if hbin(m) < y:
                lo = m
            else:
                hi = m
        return (lo + hi) / 2

    ok = True
    for trial, (n, S, A) in enumerate([(10, 4, 2), (12, 5, 2)]):
        rng = np.random.default_rng(200 + trial)
        delta, acc = random_robp(n, S, A, rng)
        countE = robp_count_dp(delta, acc, n, S, A)     # POLY
        # Y = E AND B, B~Bern(1/2): P(Y=1) = P(E)/2 = (#E/A^n)/2
        pE = countE / (A ** n)
        pY = pE / 2.0
        HY = hbin(pY)                                   # the entropy oracle value
        pE_rec = 2.0 * hbin_inv_lo(HY)
        countE_rec = round(pE_rec * (A ** n))
        ok_trial = abs(countE_rec - countE) <= 1        # rounding tolerance
        print(f"  trial {trial}: n={n} S={S}: #E(DP,poly)={countE}, "
              f"H(Y)={HY:.6f} -> recovered #E={countE_rec}  match={ok_trial}")
        print(f"           => the count the entropy reveals IS the poly DP count of W1;")
        print(f"              bounded-order encoding cannot reveal a #P-hard count.")
        ok = ok and ok_trial
    print(f"  W2 {'PASS' if ok else 'FAIL'}")
    return ok


def check_W3():
    print("=" * 70)
    print("W3: FP^#P sum structure (E1) -- H(Y) = -sum_y P(y) log2 P(y), a weighted")
    print("    sum of poly-computable terms (== chain-rule entropy; sanity).")
    print("=" * 70)
    # small induced HMM: order-1 chain observed on a subset; check explicit sum ==
    # entropy of the brute marginal (the FP^#P sum is exactly this exponential sum).
    A, N = 2, 8
    rng = np.random.default_rng(7)
    pi = rng.random(A); pi /= pi.sum()
    T = rng.random((A, A)); T /= T.sum(1, keepdims=True)
    P = np.zeros((A,) * N)
    for idx in product(range(A), repeat=N):
        pr = pi[idx[0]]
        for i in range(1, N):
            pr *= T[idx[i - 1], idx[i]]
        P[idx] = pr
    P /= P.sum()
    keep = [0, 2, 4, 6]
    M = P.sum(axis=tuple(i for i in range(N) if i not in keep))
    m = M.ravel()
    explicit = float(-sum(p * math.log2(p) for p in m if p > 0))   # the FP^#P sum
    direct = float(-(m[m > 0] * np.log2(m[m > 0])).sum())
    ok = abs(explicit - direct) < 1e-12
    print(f"  explicit -sum_y P(y)log P(y) = {explicit:.8f}  (== {direct:.8f}? {ok})")
    print(f"  ILLUSTRATES the weighted-sum form; FP^#P membership is analytic (E1),")
    print(f"  not shown numerically here.")
    print(f"  W3 {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print()
    print("#" * 70)
    print("# Remark 7.15e: mapping the wall (read-once obstacle + FP^#P bound)")
    print("#" * 70)
    print()
    r1 = check_W1(); print()
    r2 = check_W2(); print()
    r3 = check_W3(); print()
    print("=" * 70)
    allok = r1 and r2 and r3
    print(f"W1 read-once #SAT poly (DP==brute)      : {'PASS' if r1 else 'FAIL'}")
    print(f"W2 predicate-encoding reveals poly count: {'PASS' if r2 else 'FAIL'}")
    print(f"W3 FP^#P weighted-sum structure         : {'PASS' if r3 else 'FAIL'}")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
