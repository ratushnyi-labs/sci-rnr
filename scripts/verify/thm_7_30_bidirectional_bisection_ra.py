r"""
Verification for Remark 7.30 (fill-in-the-middle / bisection coding: a NEURAL realisation
of the Theorem 7.27' (Ferragina-Venturini) escape, and its limits). NOTE: an earlier draft
claimed a "Theorem 7.30: O(log N) random access at rate H(X)"; TWO adversarial reviews
(codex formal + Gemini structural) cut it down. What survives: (1) ideal rate = H(X) for
OBSERVABLE finite-order Markov; (2) it removes the sub-block delta_inf by SHARING anchors
(non-sub-block, like FV's 7.27' escape, realised with a learned bidirectional predictor).
What was BROKEN: "finite-state" scope (HMM does NOT screen -- V5); "O(log N) access" (a
single arithmetic-coded stream is sequential, far query = Theta(N) -- V3); compressed
random access at entropy is classical (GGV/FV) not new.

THE IDEA (user-proposed). A causal transformer predicts X_i | X_{<i} (forward only). A
DIRECTION-FLAGGED / fill-in-the-middle (FIM) predictor predicts X_i given context on
EITHER or BOTH sides -- the "bridge" P(X_mid | X_left, X_right). This enables BISECTION
order coding: store the endpoints, then code the midpoint of each interval conditioned on
its two endpoints, recursively. Two facts:
  (A) RATE = H(X) (observable order-m). The bisection order is a valid chain-rule
      reordering of the joint: for an OBSERVABLE order-m Markov chain P(X_mid | all-coded)
      = P(X_mid | nearest left, nearest right coded) (Markov screening with m-wide anchors)
      = the bridge, so the bridge product telescopes to P(X). No rate penalty.
  (B) delta_inf REMOVED by SHARING anchors (non-sub-block). Adjacent blocks share each
      boundary symbol (coded ONCE, as data), not two disjoint sync states -- so Remark
      7.29a's DISJOINT-charging converse does not bind. This is the SAME escape as Theorem
      7.27' (Ferragina-Venturini already attains Nh+o(N) with fast access for finite-order
      Markov); the FIM/bisection scheme is the NEURAL realisation, not a new asymptote.

WHAT WAS BROKEN by adversarial review (kept here as honest checks):
  - "finite-state" scope: an HMM does NOT screen with symbol anchors (V5).
  - "O(log N) access": logical tree DEPTH is not bitstream ACCESS; a single arithmetic
    stream is sequential, a far pre-order query needs Theta(N) predecessors (V3).
  - novelty: compressed random access at entropy is classical (GGV wavelet trees SODA
    2003; Ferragina-Venturini TCS 2007); FIM is prior art (OpenAI 2022; InCoder 2022).

WHAT THIS SCRIPT VERIFIES (order-1 Markov unless noted, ideal direction-flagged predictor).
(V1) BRIDGE CONDITIONAL. P(X_mid | X_left=a, X_right=b) at left/right distances dl,dr is
     proportional to (T^dl)[a,mid]*(T^dr)[mid,b]; verify it equals the exact conditional
     from the joint, and that summing over mid gives 1 (proper distribution).
(V2) RATE = H(X). The total bisection codelength sum over the tree of
     -log2 P(X_node | endpoints), plus the endpoints, equals -log2 P(X^N) for every
     realisation (telescoping). Hence E[codelength] = H(X^N): NO rate penalty.
(V3) DEPTH IS NOT ACCESS. Tree depth = ceil(log2(N-1)) is the LOGICAL bridge-depth; but in
     a single sequential arithmetic stream the bitstream rank of a far-right query's second
     ancestor is Theta(N), so naive random access is Theta(N), NOT O(log N). Real random
     access needs block-granular entry points (O(K) access). Both are tabulated.
(V4) ANCHOR WIDTH = MEMORY ORDER (the real scope boundary). The rate=H(X) telescoping
     needs the two anchors to d-SEPARATE the midpoint from the rest. For an order-1 source
     two single-symbol anchors screen (excess 0). For an ORDER-2 source single-symbol
     anchors do NOT screen: the bisection rate sum_nodes H(X_mid|X_L,X_R) EXCEEDS H(X^N)
     by a positive amount (verified by brute force) -- so an order-m source needs
     m-symbol anchors (=> O(m log N) access at rate H(X)), and an UNBOUNDED-memory source
     pays a genuine bisection EXCESS (the access-rate tradeoff resurfaces). NOTE: this is a
     RATE statement, not losslessness -- ANY normalised bridge codes losslessly (encoder
     and decoder agree); the excess is the price of using too-narrow anchors.
(V5) FINITE-STATE (HMM) DOES NOT SCREEN (codex counterexample). A hidden Markov source
     (hidden Y_t, P(Y_t=Y_{t-1})=0.9, noisy emission P(X_t=Y_t)=0.8) is finite-STATE but
     not observable-finite-ORDER: the sufficient statistic is the FILTER state, not any
     fixed window of symbols. So P(X_1|X_0,X_2) != P(X_1|X_0,X_2,X_4) and single-anchor
     bisection pays a positive excess (~0.028 bits at N=5). Hence the rate=H claim is
     scoped to OBSERVABLE finite-order sources; "finite-state" was an overclaim.
(V6) delta_inf REMOVED by SHARING anchors (refutes "pruned leaves become causal sub-blocks
     so delta_inf snaps back"). A PRUNED bisection -- coarse anchors every K, interiors
     decoded from BOTH endpoints -- telescopes to H(X^N)=Nh+O(1) (no per-block penalty);
     decoding the SAME leaves ONE-SIDED (causal sub-block) reintroduces a positive penalty
     = the shared-future info. So sharing the boundary (coded once as data) genuinely
     removes the sub-block (N/K)delta_inf. (Shown for order-1: one boundary symbol; an
     order-m source shares the full m-wide boundary context.)

PASS = V1 (bridge exact & normalised) AND V2 (order-1 bisection rate == -log2 P(X)) AND
V3 (depth = ceil(log2(N-1)) but far-query bitstream rank = Theta(N): depth != access) AND
V4 (order-1 single-anchor excess 0, order-2 > 0) AND V5 (HMM single-anchor excess > 0) AND
V6 (two-sided pruned bisection == H(X), one-sided causal leaves > H(X)).

HONEST SCOPE. RATE = H(X) (V2) holds for OBSERVABLE finite-order-m Markov with m-wide
anchors; finite-STATE/HMM (V5) and longer-memory/neural laws pay a bisection excess. The
construction REMOVES the sub-block (N/K)delta_inf by SHARING anchors (non-sub-block) -- the
SAME escape as Theorem 7.27' (Ferragina-Venturini), realised with a learned bidirectional
predictor; compressed random access at entropy is classical (GGV/FV), not new here. ACCESS
is O(K) at block-granular entry points, NOT O(log N) in a single sequential stream (V3).
Whether a GENUINELY NEURAL (exp-mixing, not fixed-order) source admits the escape with
w=Theta(log N) anchors at o(N) excess is the OPEN direction (toward Remark 7.27d).
"""

import math
import sys
from itertools import product

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


def bridge(T, a, b, dl, dr):
    """P(X_mid | X_left=a (dl steps left), X_right=b (dr steps right)), order-1 chain.
    prop to (T^dl)[a, mid] * (T^dr)[mid, b]."""
    A = T.shape[0]
    Tl = np.linalg.matrix_power(T, dl)
    Tr = np.linalg.matrix_power(T, dr)
    w = Tl[a, :] * Tr[:, b]
    s = w.sum()
    return w / s if s > 0 else np.ones(A) / A


def joint_logp(T, pi, x):
    lp = math.log2(pi[x[0]])
    for i in range(1, len(x)):
        lp += math.log2(T[x[i - 1], x[i]])
    return lp     # log2 P(x)


def bisection_order(lo, hi):
    """Yield (mid, lo, hi) for the bisection tree of the open interval (lo,hi);
    endpoints lo,hi are assumed already coded. Recurse on (lo,mid),(mid,hi)."""
    if hi - lo <= 1:
        return
    mid = (lo + hi) // 2
    yield (mid, lo, hi)
    yield from bisection_order(lo, mid)
    yield from bisection_order(mid, hi)


def check_V1():
    print("=" * 70)
    print("V1: bridge P(X_mid|X_left,X_right) prop (T^dl)[a,mid](T^dr)[mid,b];")
    print("    matches the exact conditional and is normalised.")
    print("=" * 70)
    ok = True
    for trial, (A, temp, N) in enumerate([(2, 0.9, 5), (3, 1.1, 5)]):
        rng = np.random.default_rng(10 + trial)
        T, pi = markov1(A, temp, rng)
        # exact joint over X^N
        P = np.zeros((A,) * N)
        for idx in product(range(A), repeat=N):
            pr = pi[idx[0]]
            for i in range(1, N):
                pr *= T[idx[i - 1], idx[i]]
            P[idx] = pr
        P /= P.sum()
        # check bridge at mid=2, left=0 (dl=2), right=4 (dr=2)
        mid, lo, hi = 2, 0, 4
        max_err = 0.0
        for a in range(A):
            for b in range(A):
                # exact P(X_2 | X_0=a, X_4=b) from the joint
                axes_out = tuple(i for i in range(N) if i not in (0, 2, 4))
                M = P.sum(axis=axes_out)            # P(X_0,X_2,X_4)
                denom = M[a, :, b].sum()
                if denom <= 0:
                    continue
                exact = M[a, :, b] / denom
                br = bridge(T, a, b, mid - lo, hi - mid)
                max_err = max(max_err, float(np.abs(exact - br).max()))
        norm_ok = abs(bridge(T, 0, 0, 2, 2).sum() - 1.0) < 1e-12
        print(f"  trial {trial}: A={A} N={N}: max|bridge-exact|={max_err:.2e}  "
              f"normalised={norm_ok}")
        ok = ok and max_err < 1e-9 and norm_ok
    print(f"  V1 {'PASS' if ok else 'FAIL'}")
    return ok


def check_V2():
    print("=" * 70)
    print("V2: bisection codelength == -log2 P(X^N) every realisation (rate = H(X)).")
    print("=" * 70)
    ok = True
    for trial, (A, temp, N) in enumerate([(2, 0.9, 9), (3, 1.1, 9), (2, 1.4, 17)]):
        rng = np.random.default_rng(20 + trial)
        T, pi = markov1(A, temp, rng)
        nodes = list(bisection_order(0, N - 1))   # interior bisection points
        max_err = 0.0
        for _ in range(2000):
            x = np.empty(N, dtype=int)
            x[0] = rng.choice(A, p=pi)
            for i in range(1, N):
                x[i] = rng.choice(A, p=T[x[i - 1]])
            # bisection codelength: endpoints X_0, X_{N-1}, then bridges
            L = -math.log2(pi[x[0]])
            L += -math.log2((np.linalg.matrix_power(T, N - 1))[x[0], x[N - 1]])  # X_{N-1}|X_0
            for (mid, lo, hi) in nodes:
                br = bridge(T, x[lo], x[hi], mid - lo, hi - mid)
                L += -math.log2(br[x[mid]])
            true_L = -joint_logp(T, pi, x)         # -log2 P(x)
            max_err = max(max_err, abs(L - true_L))
        print(f"  trial {trial}: A={A} N={N}: max |bisection L - (-log2 P(x))| = "
              f"{max_err:.2e}  ({len(nodes)+2} coded units)")
        ok = ok and max_err < 1e-7
    print("  => bisection coding is LOSSLESS at rate H(X): no first-order penalty.")
    print(f"  V2 {'PASS' if ok else 'FAIL'}")
    return ok


def _preorder_rank_of_far_query(N):
    """In the pre-order single-stream layout (root, then LEFT subtree fully, then RIGHT),
    decoding a far-right leaf needs the root's RIGHT child, whose pre-order rank is
    1 + (size of the whole left subtree) = Theta(N). Return that rank for N=2^q+1."""
    # nodes emitted before the right child = 1 (root) + all interior nodes of left half
    # left half is interval (0, mid) of size mid = (N-1)//2 -> (mid-1) interior nodes
    mid = (N - 1) // 2
    left_interior = max(mid - 1, 0)
    return 1 + left_interior  # rank (1-indexed) of the root's right child in the stream


def check_V3():
    print("=" * 70)
    print("V3: DEPTH != ACCESS. Logical bridge-depth is O(log N), but the BITSTREAM rank")
    print("    of a far-right query in a single sequential stream is Theta(N).")
    print("=" * 70)
    ok = True
    for q in (10, 20, 30, 40):
        N = 2**q + 1
        depth = math.ceil(math.log2(N - 1))           # logical tree depth
        rank = _preorder_rank_of_far_query(N)          # bitstream predecessors (sequential)
        print(f"    N={N-1:>14}: logical depth={depth:>3}   but far-query bitstream "
              f"rank={rank:>14}  (Theta(N): {rank/(N-1):.2f}*N)")
        # depth must be logarithmic AND the sequential bitstream rank must be ~N/2 (Theta N)
        ok = ok and depth <= q + 1 and rank >= (N - 1) // 4
    print("  => O(log N) is the LOGICAL bridge-depth, NOT bitstream access; a single")
    print("     arithmetic stream is sequential. Real random access needs block-granular")
    print("     entry points -> O(K) access (the O(sqrt N) regime, but WITHOUT delta_inf).")
    print(f"  V3 {'PASS' if ok else 'FAIL'}")
    return ok


def _entropy(p):
    p = np.asarray(p, dtype=float).ravel()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def _marg_entropy(P, keep):
    """Entropy of the marginal over coordinates in `keep`."""
    keep = tuple(sorted(keep))
    axes_out = tuple(i for i in range(P.ndim) if i not in keep)
    M = P.sum(axis=axes_out) if axes_out else P
    return _entropy(M)


def _cond_entropy(P, target, given):
    """H(target | given) = H(target u given) - H(given)."""
    tg = set(target) | set(given)
    Hg = _marg_entropy(P, given) if given else 0.0
    return _marg_entropy(P, tg) - Hg


def _bisection_rate(P, N):
    """Ideal bisection codelength with SINGLE-symbol anchors (two endpoints per node),
    from the true joint P: H(X0) + H(X_{N-1}|X0) + sum_nodes H(X_mid | X_lo, X_hi)."""
    r = _cond_entropy(P, {0}, set()) + _cond_entropy(P, {N - 1}, {0})
    for (mid, lo, hi) in bisection_order(0, N - 1):
        r += _cond_entropy(P, {mid}, {lo, hi})
    return r


def _joint_order1(A, rng, N):
    T, pi = markov1(A, 1.0, rng)
    P = np.zeros((A,) * N)
    for idx in product(range(A), repeat=N):
        pr = pi[idx[0]]
        for i in range(1, N):
            pr *= T[idx[i - 1], idx[i]]
        P[idx] = pr
    return P / P.sum()


def _joint_order2(A, rng, N):
    # random order-2 kernel K[a,b,c]=P(x_i=c|x_{i-2}=a,x_{i-1}=b) and pair init P2[a,b]
    K = rng.random((A, A, A)) + 0.05
    K /= K.sum(axis=2, keepdims=True)
    P2 = rng.random((A, A)) + 0.05
    P2 /= P2.sum()
    P = np.zeros((A,) * N)
    for idx in product(range(A), repeat=N):
        pr = P2[idx[0], idx[1]]
        for i in range(2, N):
            pr *= K[idx[i - 2], idx[i - 1], idx[i]]
        P[idx] = pr
    return P / P.sum()


def check_V4():
    print("=" * 70)
    print("V4: ANCHOR WIDTH = MEMORY ORDER. order-1 single-anchor screens (excess 0);")
    print("    order-2 single-anchor does NOT (bisection rate > H(X), excess > 0).")
    print("=" * 70)
    ok = True
    N = 5
    for trial in range(3):
        rng = np.random.default_rng(40 + trial)
        A = 2 + (trial % 2)
        P1 = _joint_order1(A, rng, N)
        H1 = _marg_entropy(P1, range(N))
        bis1 = _bisection_rate(P1, N)
        exc1 = bis1 - H1
        P2 = _joint_order2(A, rng, N)
        H2 = _marg_entropy(P2, range(N))
        bis2 = _bisection_rate(P2, N)
        exc2 = bis2 - H2
        print(f"  trial {trial}: A={A} N={N}")
        print(f"    order-1: H={H1:.4f}  bisection={bis1:.4f}  excess={exc1:+.2e}  "
              f"(single-symbol anchors SCREEN)")
        print(f"    order-2: H={H2:.4f}  bisection={bis2:.4f}  excess={exc2:+.4f}  "
              f"(single-symbol anchors do NOT screen)")
        ok = ok and abs(exc1) < 1e-9 and exc2 > 1e-3
    print("  => anchor width must match memory order m: order-m needs m-symbol anchors")
    print("     (access O(m log N) at rate H); unbounded memory pays a positive excess.")
    print("     (This is a RATE statement; any normalised bridge is still LOSSLESS.)")
    print(f"  V4 {'PASS' if ok else 'FAIL'}")
    return ok


def _joint_hmm(N, a=0.9, em=0.8):
    """Binary HMM joint over observations X^N: hidden Y_t with P(Y_t=Y_{t-1})=a, emission
    P(X_t=Y_t)=em. Finite-STATE but not observable-finite-ORDER. (codex counterexample)"""
    P = np.zeros((2,) * N)
    for x in product(range(2), repeat=N):
        total = 0.0
        for y in product(range(2), repeat=N):
            pr = 0.5
            for t in range(1, N):
                pr *= a if y[t] == y[t - 1] else 1 - a
            for yt, xt in zip(y, x):
                pr *= em if xt == yt else 1 - em
            total += pr
        P[x] = total
    return P / P.sum()


def check_V5():
    print("=" * 70)
    print("V5: FINITE-STATE (HMM) does NOT screen with symbol anchors (codex case).")
    print("    'finite-state' was an overclaim; scope is OBSERVABLE finite-order.")
    print("=" * 70)
    N = 5
    P = _joint_hmm(N)
    H = _marg_entropy(P, range(N))
    bis = _bisection_rate(P, N)
    exc = bis - H

    def cond(target, given_idx):
        # P(X_target=1 | X_i=1 for i in given_idx)
        num = den = 0.0
        for x in product(range(2), repeat=N):
            if all(x[i] == 1 for i in given_idx):
                den += P[x]
                if x[target] == 1:
                    num += P[x]
        return num / den

    p_local = cond(1, (0, 2))         # P(X1=1 | X0=1, X2=1)
    p_plusfar = cond(1, (0, 2, 4))    # P(X1=1 | X0=1, X2=1, X4=1)
    print(f"    P(X1=1 | X0=1, X2=1)       = {p_local:.10f}   (local anchors)")
    print(f"    P(X1=1 | X0=1, X2=1, X4=1) = {p_plusfar:.10f}   (+ far coded endpoint)")
    print(f"    => differ by {abs(p_plusfar - p_local):.4f}: anchors do NOT screen.")
    print(f"    H(X^5)={H:.6f}  single-anchor bisection={bis:.6f}  excess={exc:+.6f} bits")
    ok = abs(p_plusfar - p_local) > 1e-3 and exc > 1e-3
    print("  => the sufficient statistic is the hidden FILTER state, not a symbol window")
    print("     (the §7.15 'open middle'); rate=H holds only for OBSERVABLE order-m.")
    print(f"  V5 {'PASS' if ok else 'FAIL'}")
    return ok


def check_V6():
    print("=" * 70)
    print("V6: delta_inf REMOVED by SHARING anchors (refutes 'pruned leaves go causal').")
    print("    Pruned bisection = coarse anchors every K + interior blocks decoded TWO-SIDED")
    print("    telescopes to H(X) (no penalty); decoding those SAME leaves ONE-SIDED")
    print("    (causal sub-block) reintroduces a positive penalty = the shared-future info.")
    print("=" * 70)
    ok = True
    # order-1 Markov, N=9, anchors at 0,4,8 -> 2 interior blocks {1,2,3},{5,6,7}, K=4
    N, K = 9, 4
    anchors = [0, 4, 8]
    blocks = [((1, 2, 3), (0, 4)), ((5, 6, 7), (4, 8))]
    for trial in range(3):
        rng = np.random.default_rng(70 + trial)
        A = 2
        P = _joint_order1(A, rng, N)
        H_full = _marg_entropy(P, range(N))
        coarse = _marg_entropy(P, anchors)                          # H(X_0,X_4,X_8)
        two_sided = coarse + sum(_cond_entropy(P, set(intr), set(ep))
                                 for intr, ep in blocks)            # both endpoints
        one_sided = coarse + sum(_cond_entropy(P, set(intr), {ep[0]})
                                 for intr, ep in blocks)            # LEFT endpoint only
        gap = one_sided - H_full
        print(f"  trial {trial}: A={A} N={N} K={K}: H(X^N)={H_full:.5f}")
        print(f"    two-sided shared-anchor bisection = {two_sided:.5f}  "
              f"(penalty {two_sided-H_full:+.2e}: NONE -- = H(X))")
        print(f"    one-sided causal leaves           = {one_sided:.5f}  "
              f"(penalty {gap:+.5f}: the delta the SHARED right anchor removes)")
        ok = ok and abs(two_sided - H_full) < 1e-9 and gap > 1e-3
    print("  => sharing the boundary symbol (coded ONCE as data) and decoding interiors")
    print("     TWO-SIDED attains H(X); the leaves do NOT 'become causal sub-blocks'.")
    print("     So the sub-block (N/K)delta_inf is genuinely removed (non-sub-block).")
    print(f"  V6 {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print()
    print("#" * 70)
    print("# Remark 7.30: FIM/bisection RNR -- neural realisation of the 7.27' escape")
    print("#  (rate=H for OBSERVABLE order-m; depth != access; HMM does not screen)")
    print("#" * 70)
    print()
    r1 = check_V1(); print()
    r2 = check_V2(); print()
    r3 = check_V3(); print()
    r4 = check_V4(); print()
    r5 = check_V5(); print()
    r6 = check_V6(); print()
    print("=" * 70)
    allok = r1 and r2 and r3 and r4 and r5 and r6
    print(f"V1 bridge conditional exact+normalised : {'PASS' if r1 else 'FAIL'}")
    print(f"V2 bisection rate = H(X) (observable)  : {'PASS' if r2 else 'FAIL'}")
    print(f"V3 depth != bitstream access (Theta N) : {'PASS' if r3 else 'FAIL'}")
    print(f"V4 anchor width = memory order (scope) : {'PASS' if r4 else 'FAIL'}")
    print(f"V5 HMM does NOT screen (finite-state)  : {'PASS' if r5 else 'FAIL'}")
    print(f"V6 delta_inf removed by shared anchors : {'PASS' if r6 else 'FAIL'}")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
