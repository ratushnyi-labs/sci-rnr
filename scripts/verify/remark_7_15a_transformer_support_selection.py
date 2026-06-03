r"""
Verification for Remark 7.15a (transformer-predictor support selection:
the §6.9 Gaussian-MESP framework extends to the discrete autoregressive
neural-LM field of Theorem 7.15, with a sharper submodular structure and a
new representation-dependent oracle dichotomy).

------------------------------------------------------------------------
THE FRONTIER QUESTION.

§6.9 (Lemma 6.9 + Remarks 6.9a/c/d + Cor 6.9b) characterises Type-III-C
SUPPORT SELECTION for a GAUSSIAN entropy field X ~ N(0, Sigma):

    residual of coding S, reconstructing Sbar  =  C(S) = h(X_S | X_Sbar)
    complementation identity (chain rule):       C(S) = h(X) - M(Sbar),
        M(T) := h(X_T) = 0.5 log det(2 pi e Sigma_T)   (log-det)

  - M is SUBMODULAR (log-det; Kelmans-Kimelfeld / Krause-Guestrin 2008),
    MONOTONE and NON-NEGATIVE only after the calibration sigma^2 >= 1/(2 pi e)
    (Remark 6.9d: this floor is NECESSARY, not cosmetic);
  - greedy gives M(Tgreedy) >= (1-1/e) M(Topt) (Nemhauser-Wolsey-Fisher
    1978), hence the ADDITIVE-REGRET residual bound (Cor 6.9b, eq 6.9b.1);
  - fixed-field is APX-hard 5/4 via MESP (Lemma 6.9); free-field amortised
    is EASY in P (Remark 6.9a, sample-covariance MLE plug-in).

Theorem 7.15 is the END-TO-END NEURAL-LM RNR: a transformer predictor M'
produces a per-position field Q_i(. | x_{<i}); the field is DISCRETE,
NON-Gaussian, autoregressive. DOES the §6.9 picture transfer?

------------------------------------------------------------------------
THE ANSWER THIS SCRIPT VERIFIES (Remark 7.15a).

Formalisation. Ground set V = {1,...,N} of sub-block positions (T7.5:
inside one sub-block of size K = Theta(sqrt N)). Encoder picks a CODED set
S; codes X_S directly and RECONSTRUCTS X_Sbar from M'. The lossless
residual is
    C(S) = H(X_S) + Hrecon(Sbar | S, M'),
and with the IDEAL predictor M' = the true conditional law, the
information-theoretic reconstruction cost is Hrecon = H(X_Sbar | X_S), so
    C(S) = H(X_S) + H(X_Sbar | X_S) = H(X)         (chain rule, EXACT)
        = H(X) - M(Sbar) + M(Sbar)  ... i.e. the SELECTION lever is
    C_sel(S) := H(X_S | X_Sbar) = H(X) - M(Sbar),  M(T) := H(X_T).
This is the discrete mirror of Cor 6.9b's identity, with DIFFERENTIAL
entropy / log-det replaced by SHANNON entropy.

Three structural facts, each checked numerically:

(V1) M(T) = H(X_T) is monotone, non-negative, SUBMODULAR for ANY discrete
     joint law -- including the transformer/AR field -- by FUJISHIGE 1978
     (entropy is a polymatroid rank). UNLIKE the Gaussian case this needs
     NO calibration: discrete entropy is automatically >= 0 (probabilities
     <= 1) and automatically monotone (H(X_{T+i}) - H(X_T) = H(X_i|X_T) >=
     0). So Remark 6.9d's necessary noise floor sigma^2 >= 1/(2 pi e)
     DISAPPEARS -- the discrete result is STRICTLY CLEANER.

(V2) Greedy maximisation of M gives M(Tgreedy) >= (1-1/e) M(Topt)
     (Nemhauser-Wolsey-Fisher 1978), and the residual transfers ADDITIVELY
     exactly as in Cor 6.9b: C_sel(greedy) - C_sel(opt) = M(Topt) -
     M(Tgreedy) <= (1/e) M(Topt). Verified vs BRUTE-FORCE optimum on small
     discrete AR fields.

(V3) THE NEW DICHOTOMY (representation-dependent oracle, absent in the
     Gaussian setting where Sigma_T is read off in O(N^2)). Greedy needs an
     efficient VALUE ORACLE for H(X_T). For an AR field this is marginal
     inference P(x_T) = sum_{x_Tbar} prod_i Q_i(x_i|x_{<i}):
       (a) ORDER-1 MARKOV, ARBITRARY subset T: POLY-TIME. Marginalising the
           gaps between consecutive T-elements is a transfer-matrix product
           (Chapman-Kolmogorov), so H(X_T) for ANY T costs O(|gaps| * |Sig|^3).
           => greedy support selection is FULLY poly-time here.
       (b) BOUNDED-ORDER-W MARKOV, CONTIGUOUS block T={a..b}: POLY-TIME via
           the forward / transfer-matrix algorithm, O((b-a) |Sig|^{W+1}).
       (c) GENERAL AR (arbitrary Q_i) or HIGH-ORDER, ARBITRARY subset T:
           the entropy ORACLE is #P-hard. Marginal-inference #P-hardness
           (Roth 1996; exact inference already NP-hard, Cooper 1990) is about
           marginal PROBABILITIES, not entropy per se;
           the explicit reduction PROBABILITY -> ENTROPY-ORACLE is V5 below
           (Y = E AND fair-bit B makes H(Y)=h(P(E)/2) injective, so an
           entropy oracle computes the #P-hard marginal). So general AR adds
           a SEPARATE #P-hard value-oracle obstacle ON TOP of the (already
           NP-hard, V4) selection -- two obstacles, not a relocation. The
           |Sigma|^{W+1} poly-time in (a)/(b) is genuine only for CONSTANT
           alphabet; for W=O(log N) it is poly only with |Sigma|=O(1).
     We CHECK (a) that a GENUINELY poly-time O(r |Sig|^2) chain-rule entropy
     -- H(X_T)=H(X_{t_1})+sum_k H(X_{t_k}|X_{t_{k-1}}), never materialising the
     A^{|T|} table -- equals the brute-force marginal entropy on every subset
     (the induced subset process of an order-1 chain is again order-1: the
     observed previous kept symbol IS the full state, so no hidden-state belief
     branching arises -- which is exactly why ARBITRARY subsets are tractable
     at order 1; at order>=2 with gaps the induced process is hidden-state and
     the (a) method fails, so (b) uses contiguity as a SUFFICIENT condition --
     whether noncontiguous bounded-order is poly is left OPEN), and exhibit
     (c)'s exponential blow-up of brute marginalisation.

(V4) The selection problem itself, GIVEN an oracle, is NP-hard already in
     the discrete setting -- and over FREE position choice, faithful to
     Lemma 6.9's MESP (not a candidate-restricted variant). Gadget: U latent
     independent fair coins E_1..E_U (NOT selectable); one POSITION per
     candidate set C_i, with value X_i := (E_e)_{e in C_i}. For any subset S
     of positions H(X_S)=|union_{i in S} C_i| (distinct covered coins), so
     max_{|S|=s} H(X_S) IS max-coverage (NP-hard, Feige 1998). So hardness
     TRANSFERS from §6.9 (it is not softened by discreteness); what changes is
     (i) the cleaner unconditional submodularity (V1) and (ii) general AR ADDS
     a #P-hard oracle obstacle (V3/V5). Here the oracle H(X_S)=|union| is
     POLY-time, so the SELECTION is NP-hard even with an easy oracle -- the two
     obstacles are genuinely separate. The reduction needs a COMPACT field
     representation (the coin gadget is compact: O(N) coords). We verify on a
     small instance that the explicit independent-coin joint reproduces
     H(X_S)=|union| on every subset and that free max-entropy == max-coverage.

(V5) The entropy ORACLE is #P-hard for general AR fields -- via an explicit
     reduction from marginal probability (which is #P-hard, Roth 1996) to
     the entropy oracle, patching the "inference-hardness != oracle-hardness"
     gap. For any event E, set Y = E AND B with B ~ Bern(1/2) independent;
     then P(Y=1) = P(E)/2 in [0,1/2] and H(Y) = h(P(E)/2). Since binary
     entropy h is STRICTLY INCREASING on [0,1/2], H(Y) is INJECTIVE in P(E),
     so an entropy oracle recovers the #P-hard marginal P(E) = 2 h^{-1}(H(Y)).
     Verified: h injective on [0,1/2] and P(E) recovered to bisection
     precision from H(Y) alone.

PASS = V1 (submodular+monotone+nonneg, NO calibration) AND V2 (greedy
(1-1/e) + additive-regret identity vs brute force) AND V3 (transfer-matrix
oracle exact for order-1 arbitrary subsets + contiguous bounded-order; brute
blows up) AND V4 (max-coverage reduction: max-entropy subset == max-cover)
AND V5 (entropy-oracle #P-hardness reduction is injective/exact).

HONEST SCOPE. Like Cor 6.9b this is the FIXED-FIELD lever (the model M' is
given, as in T7.15 where M' is a pre-trained shared predictor). In the
free-field amortised regime the Remark 6.9a collapse applies mutatis
mutandis (a free M' codes at H(X) for every S, selection vacuous). The new
content over §6.9 is: (1) submodularity is unconditional/cleaner for the
discrete field; (2) a representation-dependent ORACLE dichotomy that has no
Gaussian analogue and locates the true practical bottleneck for the T7.15
encoder.
"""

import math
import sys
from itertools import combinations, product

import numpy as np

RNG = np.random.default_rng(20260602)
TOL = 1e-9


# ======================================================================
# Discrete autoregressive (transformer-like) field utilities.
# A "field" is given by per-position conditionals; we build the full joint
# table for small N so we can compute EXACT entropies / marginals and
# brute-force the optimum. This is the discrete analogue of a fixed Sigma.
# ======================================================================

def random_markov_field(N, A, order=1, rng=RNG, temp=1.0):
    """Return joint pmf table P[x_0,...,x_{N-1}] (shape A^N) of an
    order-`order` Markov chain over alphabet {0,..,A-1}. `temp` < 1 makes
    it peakier (lower entropy), > 1 flatter. This stands in for a
    transformer's autoregressive softmax field."""
    shape = (A,) * N
    P = np.zeros(shape)
    # conditional tables: Q[ctx -> distribution over next symbol]
    # ctx is a tuple of up to `order` previous symbols.
    cond = {}

    def get_cond(ctx):
        if ctx not in cond:
            logits = rng.standard_normal(A) / max(temp, 1e-6)
            p = np.exp(logits - logits.max())
            p = p / p.sum()
            cond[ctx] = p
        return cond[ctx]

    for idx in product(range(A), repeat=N):
        prob = 1.0
        for i in range(N):
            ctx = tuple(idx[max(0, i - order):i])
            prob *= get_cond(ctx)[idx[i]]
        P[idx] = prob
    P /= P.sum()
    return P


def joint_entropy_subset(P, T):
    """H(X_T) in bits from the full joint table P by marginalising out
    the complement (the EXACT, brute-force oracle)."""
    N = P.ndim
    T = tuple(sorted(T))
    if len(T) == 0:
        return 0.0
    axes_to_sum = tuple(i for i in range(N) if i not in T)
    M = P.sum(axis=axes_to_sum) if axes_to_sum else P
    m = M.ravel()
    m = m[m > 0]
    return float(-(m * np.log2(m)).sum())


def total_entropy(P):
    p = P.ravel()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


# ======================================================================
# V1: submodularity + monotonicity + non-negativity, NO calibration.
# ======================================================================

def check_V1():
    print("=" * 70)
    print("V1: M(T)=H(X_T) is monotone, non-negative, SUBMODULAR for the")
    print("    discrete autoregressive field (Fujishige 1978) -- NO floor.")
    print("=" * 70)

    ok = True
    for trial, (N, A, order, temp) in enumerate(
        [(5, 2, 1, 1.0), (5, 3, 1, 0.7), (6, 2, 2, 1.3), (5, 2, 2, 0.5)]
    ):
        P = random_markov_field(N, A, order=order, temp=temp)
        V = list(range(N))
        H = {frozenset(): 0.0}
        for r in range(1, N + 1):
            for T in combinations(V, r):
                H[frozenset(T)] = joint_entropy_subset(P, T)

        # non-negativity
        nonneg = all(v >= -TOL for v in H.values())
        # monotonicity: H(T + i) >= H(T)
        mono = True
        for T in list(H):
            for i in V:
                if i in T:
                    continue
                if H[T | {i}] < H[T] - TOL:
                    mono = False
        # submodularity (diminishing returns):
        #   H(A+i) - H(A) >= H(B+i) - H(B)  for A subset B, i not in B
        submod = True
        worst = math.inf
        for A_ in list(H):
            for B_ in list(H):
                if not A_.issubset(B_):
                    continue
                for i in V:
                    if i in B_:
                        continue
                    gA = H[A_ | {i}] - H[A_]
                    gB = H[B_ | {i}] - H[B_]
                    worst = min(worst, gA - gB)
                    if gA - gB < -1e-7:
                        submod = False
        print(
            f"  trial {trial}: N={N} A={A} order={order} temp={temp}:  "
            f"nonneg={nonneg} monotone={mono} submodular={submod} "
            f"(min DR slack={worst:+.2e})"
        )
        ok = ok and nonneg and mono and submod

    print("  --> discrete entropy needs NO sigma^2>=1/(2 pi e) calibration")
    print("      (Remark 6.9d's necessary Gaussian floor simply vanishes).")
    print(f"  V1 {'PASS' if ok else 'FAIL'}")
    return ok


# ======================================================================
# V2: greedy (1-1/e) + additive-regret residual identity vs brute force.
# ======================================================================

def greedy_max_entropy(P, t):
    """Greedy: build reconstruction set T (|T|=t) maximising M(T)=H(X_T)."""
    N = P.ndim
    T = set()
    for _ in range(t):
        best_i, best_gain, best_val = None, -math.inf, None
        baseH = joint_entropy_subset(P, T)
        for i in range(N):
            if i in T:
                continue
            val = joint_entropy_subset(P, T | {i})
            gain = val - baseH
            if gain > best_gain:
                best_i, best_gain, best_val = i, gain, val
        T.add(best_i)
    return T


def brute_max_entropy(P, t):
    N = P.ndim
    best_T, best_val = None, -math.inf
    for T in combinations(range(N), t):
        val = joint_entropy_subset(P, T)
        if val > best_val:
            best_T, best_val = set(T), val
    return best_T, best_val


def check_V2():
    print("=" * 70)
    print("V2: greedy M(Tgreedy) >= (1-1/e) M(Topt); residual transfers")
    print("    ADDITIVELY  C_sel(greedy)-C_sel(opt) = M(opt)-M(greedy)")
    print("    (discrete mirror of Cor 6.9b eq 6.9b.1).  vs BRUTE FORCE.")
    print("=" * 70)

    ok = True
    factor = 1.0 - 1.0 / math.e
    for trial, (N, A, order, temp) in enumerate(
        [(6, 2, 1, 1.0), (7, 2, 1, 0.6), (6, 3, 1, 1.4), (7, 2, 2, 0.8)]
    ):
        P = random_markov_field(N, A, order=order, temp=temp)
        Hfull = total_entropy(P)
        t = N // 2  # reconstruct half, code half
        Tg = greedy_max_entropy(P, t)
        Mg = joint_entropy_subset(P, Tg)
        Topt, Mopt = brute_max_entropy(P, t)

        # C_sel(S) = H(X) - M(Sbar);  here Sbar = T (the reconstruction set)
        Cg = Hfull - Mg
        Copt = Hfull - Mopt
        add_identity = abs((Cg - Copt) - (Mopt - Mg)) < 1e-7
        ratio_ok = Mg >= factor * Mopt - 1e-7
        add_bound = (Cg - Copt) <= (1.0 / math.e) * Mopt + 1e-7

        print(
            f"  trial {trial}: N={N} A={A} order={order} t={t}:  "
            f"M_greedy={Mg:.4f}  M_opt={Mopt:.4f}  "
            f"ratio={Mg/max(Mopt,1e-12):.3f} (>= {factor:.3f}? {ratio_ok})"
        )
        print(
            f"            additive identity C_g-C_opt == M_opt-M_g: "
            f"{add_identity};  additive bound <= M_opt/e: {add_bound}"
        )
        ok = ok and add_identity and ratio_ok and add_bound

    # Honest caveat (mirrors Remark 6.9c): the (1-1/e) is a WORST-CASE
    # guarantee; on small instances greedy is near-optimal. Scan many random
    # fields and report the worst achieved ratio -- it stays well above
    # 1-1/e, confirming the bound holds without being tight at small N.
    worst_ratio = 1.0
    for seed in range(2000):
        rng = np.random.default_rng(seed)
        P = random_markov_field(7, 2, order=2, rng=rng, temp=rng.uniform(0.3, 2.0))
        Tg = greedy_max_entropy(P, 3)
        Mg = joint_entropy_subset(P, Tg)
        _, Mopt = brute_max_entropy(P, 3)
        if Mopt > 1e-6:
            worst_ratio = min(worst_ratio, Mg / Mopt)
    print(
        f"  worst greedy/opt over 2000 random fields = {worst_ratio:.4f} "
        f"(>= 1-1/e = {factor:.4f}? {worst_ratio >= factor - 1e-9}); "
        f"bound holds, not tight at small N (asymptotic worst case)."
    )
    ok = ok and (worst_ratio >= factor - 1e-9)

    print(f"  V2 {'PASS' if ok else 'FAIL'}")
    return ok


# ======================================================================
# V3: the oracle dichotomy.  Order-1 arbitrary subset via transfer matrix
#     == brute marginal entropy (exact).  Brute blows up; transfer-matrix
#     stays poly.  Contiguous bounded-order via forward algo == brute.
# ======================================================================

def markov1_params(A, rng=RNG, temp=1.0):
    """Return (pi0, Tmat) for an order-1 Markov chain: pi0 initial dist,
    Tmat[a,b] = P(next=b | cur=a)."""
    l0 = rng.standard_normal(A) / temp
    pi0 = np.exp(l0 - l0.max())
    pi0 /= pi0.sum()
    Tmat = np.zeros((A, A))
    for a in range(A):
        l = rng.standard_normal(A) / temp
        p = np.exp(l - l.max())
        Tmat[a] = p / p.sum()
    return pi0, Tmat


def markov1_joint_table(pi0, Tmat, N):
    A = len(pi0)
    P = np.zeros((A,) * N)
    for idx in product(range(A), repeat=N):
        prob = pi0[idx[0]]
        for i in range(1, N):
            prob *= Tmat[idx[i - 1], idx[i]]
        P[idx] = prob
    return P


def markov1_subset_marginal(pi0, Tmat, N, T):
    """EXACT marginal P(X_T) of an order-1 chain over an ARBITRARY subset T,
    computed by TRANSFER-MATRIX gap summation (poly-time, no exponential
    sum over the complement). Returns a table over the kept positions."""
    T = sorted(T)
    A = len(pi0)
    # distribution over the first kept position, after marginalising
    # positions 0..T[0]-1 by advancing pi0 through T[0] steps.
    cur = pi0.copy()
    for _ in range(T[0]):
        cur = cur @ Tmat                      # advance one step (sum out)
    # state distribution at position T[0] (marginal there): `cur`
    # Build the joint over kept positions incrementally.
    # joint[state_at_last_kept, kept_indices...] -- we keep a tensor whose
    # last axis is the current chain state at the last kept position.
    marg = cur.reshape(A)                      # P(X_{T[0]} = .)
    table = marg.copy()                        # shape (A,)
    last = T[0]
    for k in range(1, len(T)):
        gap = T[k] - last
        step = np.linalg.matrix_power(Tmat, gap)   # P(X_{T[k]} | X_{last})
        # table currently: shape (..., A) with last axis = state at `last`.
        # extend: new_table[..., a, b] = table[..., a] * step[a, b]
        table = table[..., :, None] * step[None, ...] if table.ndim == 1 \
            else np.einsum('...a,ab->...ab', table, step)
        last = T[k]
    return table   # shape A^{|T|}, indexed by kept positions in order


def entropy_of_table(tab):
    t = tab.ravel()
    t = t[t > 0]
    return float(-(t * np.log2(t)).sum())


def _entropy_vec(p):
    p = np.asarray(p, dtype=float).ravel()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def markov1_subset_entropy_polytime(pi0, Tmat, N, T):
    """GENUINELY poly-time H(X_T) for an order-1 chain over an ARBITRARY
    subset T, WITHOUT ever materialising the A^{|T|} joint table.

    The induced process on T = {t_1<...<t_r} is itself an order-1 Markov
    chain: the OBSERVED previous kept symbol X_{t_{k-1}} *is* the full state
    (order 1), so no hidden state / belief branching arises. By the chain
    rule,
        H(X_T) = H(X_{t_1}) + sum_{k>=2} H(X_{t_k} | X_{t_{k-1}})
               = H(X_{t_1}) + sum_{k>=2} [ H(X_{t_{k-1}},X_{t_k}) - H(X_{t_{k-1}}) ],
    and each pairwise marginal is an A x A table obtained from a gap-power of
    Tmat.  Total work O(r*A^2) + O(sum_k log(gap_k) * A^3) for the matrix
    powers -- polynomial in N and A, with NO exponential dependence on |T|.
    (Contrast markov1_subset_marginal, which builds the A^{|T|} table and is
    used only as an independent correctness oracle on small N.)"""
    T = sorted(T)
    A = len(pi0)
    # marginal state distribution at the first kept position
    cur = pi0.copy()
    for _ in range(T[0]):
        cur = cur @ Tmat
    H = _entropy_vec(cur)                       # H(X_{t_1})
    last = T[0]
    for k in range(1, len(T)):
        gap = T[k] - last
        step = np.linalg.matrix_power(Tmat, gap)   # P(X_{t_k} | X_{last})
        joint = cur[:, None] * step                # A x A pairwise marginal
        nxt = cur @ step                           # marginal at t_k
        # H(X_{t_k} | X_{last}) = H(X_{last},X_{t_k}) - H(X_{last})
        H += _entropy_vec(joint) - _entropy_vec(cur)
        cur = nxt
        last = T[k]
    return H


def order2_params(A, rng=RNG, temp=1.0):
    """Explicit order-2 chain: pair_init[a,b]=P(X_0=a,X_1=b) and
    T2[a,b,c]=P(X_i=c | X_{i-2}=a, X_{i-1}=b)."""
    l = rng.standard_normal((A, A)) / temp
    pe = np.exp(l - l.max()); pair_init = pe / pe.sum()
    T2 = np.zeros((A, A, A))
    for a in range(A):
        for b in range(A):
            ll = rng.standard_normal(A) / temp
            p = np.exp(ll - ll.max()); T2[a, b] = p / p.sum()
    return pair_init, T2


def order2_joint_table(pair_init, T2, N):
    A = pair_init.shape[0]
    P = np.zeros((A,) * N)
    for idx in product(range(A), repeat=N):
        prob = pair_init[idx[0], idx[1]]
        for i in range(2, N):
            prob *= T2[idx[i - 2], idx[i - 1], idx[i]]
        P[idx] = prob
    return P


def order2_contiguous_block_entropy_forward(pair_init, T2, N, a, b):
    """GENUINE forward-algorithm H(X_{a..b}) for a CONTIGUOUS block of an
    order-2 chain, O((b-a) A^3), never building the A^{b-a+1} table.

    Marginalising the prefix X_{<a} leaves the block an order-2 process whose
    pair-marginal pi_i(x,y)=P(X_{i-1},X_i) evolves by pi_{i+1}(y,z)=sum_x
    pi_i(x,y) T2[x,y,z].  By the chain rule, with all conditionals order-2,
        H(X_{a..b}) = H(X_a,X_{a+1}) + sum_{j=a+2}^b H(X_j | X_{j-2},X_{j-1}),
        H(X_j|X_{j-2},X_{j-1}) = sum_{x,y} pi_{j-1}(x,y) * H( T2[x,y,:] ).
    Each pi is a single A x A table advanced by one matrix contraction --
    poly in block length, NO exponential dependence.  (Contiguity is what
    keeps the conditioning set inside the block; an ARBITRARY order-2 subset
    with gaps would condition X_{t_k} on a HIDDEN pair, branching the belief.)"""
    A = pair_init.shape[0]
    # advance the pair-marginal from positions (0,1) to (a, a+1):
    pi = pair_init.copy()                       # P(X_0, X_1)
    for _ in range(a):                          # shift the window forward by a
        pi = np.einsum('xy,xyz->yz', pi, T2)    # P(X_{k},X_{k+1}) from P(X_{k-1},X_k)
    # now pi = P(X_a, X_{a+1}); H(X_a,X_{a+1}) starts the chain-rule sum
    H = _entropy_vec(pi)
    rowH = np.array([[_entropy_vec(T2[x, y]) for y in range(A)] for x in range(A)])
    for _j in range(a + 2, b + 1):
        # H(X_j | X_{j-2}, X_{j-1}) with pi = P(X_{j-2}, X_{j-1})
        H += float((pi * rowH).sum())
        pi = np.einsum('xy,xyz->yz', pi, T2)    # advance to P(X_{j-1}, X_j)
    return H


def check_V3():
    print("=" * 70)
    print("V3: the ORACLE dichotomy (no Gaussian analogue).")
    print("    (a) order-1 arbitrary subset: transfer-matrix marginal entropy")
    print("        == brute marginal entropy (EXACT); poly-time.")
    print("    (b) contiguous bounded-order block: forward algo == brute.")
    print("    (c) brute marginalisation is EXPONENTIAL (the bottleneck).")
    print("=" * 70)

    ok = True

    # (a) order-1, arbitrary subsets: transfer-matrix == brute == poly-time
    #     chain-rule, all exact.  The poly-time path NEVER builds the A^{|T|}
    #     table -- it is the operative O(r A^2) algorithm; the table and brute
    #     methods are independent correctness oracles on small N.
    for trial, (N, A, temp) in enumerate([(7, 2, 1.0), (6, 3, 0.7), (8, 2, 1.5)]):
        pi0, Tmat = markov1_params(A, temp=temp)
        P = markov1_joint_table(pi0, Tmat, N)
        max_err = 0.0
        max_err_poly = 0.0
        tested = 0
        for r in range(1, N + 1):
            for T in combinations(range(N), r):
                H_brute = joint_entropy_subset(P, T)
                tab = markov1_subset_marginal(pi0, Tmat, N, T)
                H_tm = entropy_of_table(tab)
                H_poly = markov1_subset_entropy_polytime(pi0, Tmat, N, T)
                max_err = max(max_err, abs(H_brute - H_tm))
                max_err_poly = max(max_err_poly, abs(H_brute - H_poly))
                tested += 1
        passed = max_err < 1e-7 and max_err_poly < 1e-7
        print(
            f"  (a) trial {trial}: N={N} A={A}:  tested all {tested} subsets, "
            f"max |H_brute - H_transfermatrix| = {max_err:.2e}, "
            f"max |H_brute - H_polytime(O(rA^2))| = {max_err_poly:.2e}  "
            f"{'OK' if passed else 'FAIL'}"
        )
        ok = ok and passed

    # (b) contiguous bounded-order block: GENUINE forward algorithm
    #     O((b-a)A^3) == brute, exact; never builds the A^{b-a+1} table.
    ctg_ok = True
    max_err_fwd = 0.0
    for trial, (N, A, temp) in enumerate([(7, 2, 0.9), (6, 3, 1.2)]):
        pair_init, T2 = order2_params(A, temp=temp)
        P = order2_joint_table(pair_init, T2, N)
        for a in range(N - 1):
            for b in range(a + 1, N):
                h_fwd = order2_contiguous_block_entropy_forward(pair_init, T2, N, a, b)
                h_brute = joint_entropy_subset(P, list(range(a, b + 1)))
                max_err_fwd = max(max_err_fwd, abs(h_fwd - h_brute))
                if abs(h_fwd - h_brute) > 1e-7:
                    ctg_ok = False
    print(f"  (b) order-2 CONTIGUOUS forward algo vs brute: "
          f"max err = {max_err_fwd:.2e}  {'OK' if ctg_ok else 'FAIL'}")
    ok = ok and ctg_ok

    # (b') the order-1 METHOD does not extend to noncontiguous subsets: an
    #     ARBITRARY order-2 subset with a gap induces dependence beyond order 1
    #     -- X_{t3} NOT indep of X_{t1} given X_{t2} -- so the order-1 chain-rule
    #     trick of (a) cannot apply.  HONEST SCOPE: this shows the (a) method is
    #     order-1-SPECIFIC and that (b)'s contiguity is a SUFFICIENT condition;
    #     it does NOT prove hardness.  Whether noncontiguous bounded-order subset
    #     entropy is poly-time or #P-hard is left OPEN (it sits between the poly
    #     endpoints (a)/(b) and the #P-hard general case (c)).
    pair_init, T2 = order2_params(2, temp=0.8)
    Pfull = order2_joint_table(pair_init, T2, 6)
    # subset T = {0,2,4} (gaps): test conditional independence X_4 ⟂ X_0 | X_2
    def cond_mi(P, i, j, k):
        # I(X_i ; X_j | X_k) from the joint table P
        Pijk = P.sum(axis=tuple(ax for ax in range(P.ndim) if ax not in (i, j, k)))
        # reorder axes to (i,j,k)
        order = np.argsort([i, j, k])
        Pijk = np.transpose(Pijk, axes=order)  # now axes sorted asc; relabel
        # compute I(A;B|C) over a 3-index table indexed by sorted(i,j,k)
        idx = {v: r for r, v in enumerate(sorted([i, j, k]))}
        A_, B_, C_ = idx[i], idx[j], idx[k]
        Pc = Pijk.sum(axis=tuple(x for x in range(3) if x != C_))
        mi = 0.0
        for ci in range(Pijk.shape[C_]):
            sl = [slice(None)] * 3; sl[C_] = ci
            Pab_c = Pijk[tuple(sl)]
            pc = Pab_c.sum()
            if pc <= 0:
                continue
            Pab = Pab_c / pc
            Pa = Pab.sum(axis=1 if A_ < B_ else 0)
            Pb = Pab.sum(axis=0 if A_ < B_ else 1)
            for ai in range(Pab.shape[0]):
                for bi in range(Pab.shape[1]):
                    p = Pab[ai, bi]
                    if p > 0:
                        marg = (Pa[ai] * Pb[bi])
                        mi += pc * p * math.log2(p / marg) if marg > 0 else 0.0
        return mi
    cmi = cond_mi(Pfull, 0, 4, 2)   # I(X_0 ; X_4 | X_2)
    order1_method_fails = cmi > 1e-6
    print(f"  (b') arbitrary subset {{0,2,4}}: I(X_0;X_4|X_2) = {cmi:.4f} bits "
          f"> 0 ({order1_method_fails}) => induced process is NOT order-1; the")
    print(f"       order-1 trick of (a) is order-1-SPECIFIC (fails here).")
    print(f"       Contiguity in (b) is SUFFICIENT, not shown necessary --")
    print(f"       noncontiguous bounded-order poly-vs-#P-hard is OPEN.")
    ok = ok and order1_method_fails

    # (c) exhibit the exponential cost of brute marginalisation: the number
    #     of complement-configurations summed grows as A^{N-|T|}.
    print("  (c) brute marginal cost A^{N-|T|} (exponential in N-|T|):")
    A = 4
    for N in [5, 10, 15, 20]:
        Tsz = max(1, N // 2)
        brute_terms = A ** (N - Tsz)
        # transfer-matrix order-1 cost: O(N * A^3) regardless of subset.
        tm_cost = N * A ** 3
        print(
            f"      N={N:2d}, |T|={Tsz:2d}: brute sums {brute_terms:>12,} "
            f"configs vs transfer-matrix ~{tm_cost:,} ops "
            f"(ratio {brute_terms / tm_cost:.1e}x)"
        )
    print("      => for general/high-order AR the ORACLE is a SEPARATE #P-hard")
    print("         obstacle (V5) ON TOP of NP-hard selection (V4); the")
    print("         Gaussian case has neither (Sigma_T read in O(N^2)).")

    print(f"  V3 {'PASS' if ok else 'FAIL'}")
    return ok


# ======================================================================
# V4: hardness TRANSFERS -- discrete max-entropy subset generalises
#     MAX-COVERAGE (NP-hard, Feige 1998). Reduction witness on a small
#     instance: the max-entropy subset == the max-coverage subset.
# ======================================================================

def build_coverage_field(sets, U):
    """Construct the EXPLICIT discrete field whose FREE per-position max-entropy
    selection is max-coverage -- faithful to Lemma 6.9 (selection ranges over
    POSITIONS, free choice of any size-s subset, exactly as MESP does), not a
    candidate-restricted variant.

    Latent: U independent fair coins E_1..E_U ~ Bern(1/2)  (NOT selectable).
    Positions (the ground set V the encoder selects from): one per candidate
    set C_i, with value X_i := (E_e)_{e in C_i}  -- the tuple of coins it
    covers (a deterministic function of the latent coins).

    Then for ANY subset S of positions,
        H(X_S) = H( (E_e)_{e in union_{i in S} C_i} ) = | union_{i in S} C_i |
    because the coins are independent fair bits and X_S reveals exactly the
    coins in the covered union (each contributing one bit, none twice).  So
        max_{|S|=s} H(X_S)  ==  max-coverage with s sets,
    a FREE selection over positions -- NP-hard (Feige 1998).  The value oracle
    H(X_S)=|union| is POLY-time here, so the SELECTION is hard even with an
    easy oracle: V4 (selection NP-hard) and V5 (oracle #P-hard) are genuinely
    two SEPARATE obstacles, not one relocated.

    Returns the entropy set-function H_field(S) over position-subsets S."""
    sets = [set(c) for c in sets]

    def H_field(S):
        cov = set()
        for i in S:
            cov |= sets[i]
        return float(len(cov))      # H(X_S) in bits (distinct covered coins)

    return H_field


def check_V4():
    print("=" * 70)
    print("V4: hardness TRANSFERS -- FREE per-position max-entropy subset")
    print("    selection generalises MAX-COVERAGE (Feige 1998 NP-hard).")
    print("    Faithful to Lemma 6.9: selection ranges over POSITIONS X_i,")
    print("    each X_i = tuple of independent coins it covers; H(X_S)=|union|.")
    print("=" * 70)

    # small set-cover instance.
    U = 6                      # universe of latent coins {0..5}
    sets = [
        {0, 1, 2},
        {2, 3},
        {3, 4, 5},
        {0, 5},
        {1, 4},
    ]
    s = 2                      # FREE choice of s positions (candidates-as-X_i)

    H_field = build_coverage_field(sets, U)

    # EXPLICIT cross-check that H_field is a genuine joint-entropy of an
    # independent-coin field: build the real uniform joint over {0,1}^U,
    # define X_i = tuple of covered coins, and confirm H(X_S) computed from
    # the joint table equals |union| for every subset.
    coins_joint = np.full((2,) * U, 1.0 / (2 ** U))   # uniform on {0,1}^U
    def H_from_joint(S):
        cov = sorted(set().union(*[set(sets[i]) for i in S])) if S else []
        if not cov:
            return 0.0
        axes_out = tuple(e for e in range(U) if e not in cov)
        M = coins_joint.sum(axis=axes_out) if axes_out else coins_joint
        m = M.ravel(); m = m[m > 0]
        return float(-(m * np.log2(m)).sum())
    field_consistent = all(
        abs(H_field(S) - H_from_joint(S)) < 1e-9
        for r in range(0, len(sets) + 1)
        for S in combinations(range(len(sets)), r)
    )

    # FREE max-entropy over s positions == max-coverage.
    best_cov, best_choice = -1.0, None
    for S in combinations(range(len(sets)), s):
        c = H_field(S)
        if c > best_cov:
            best_cov, best_choice = c, S

    # greedy on the entropy set-function reproduces the (1-1/e) guarantee.
    def greedy_field(s):
        S, cov = [], set()
        for _ in range(s):
            best_i, best_gain = None, -1.0
            for i in range(len(sets)):
                if i in S:
                    continue
                gain = len(cov | set(sets[i])) - len(cov)
                if gain > best_gain:
                    best_i, best_gain = i, gain
            S.append(best_i); cov |= set(sets[best_i])
        return S, float(len(cov))

    g_choice, g_cov = greedy_field(s)
    factor = 1.0 - 1.0 / math.e
    ratio_ok = g_cov >= factor * best_cov - 1e-9

    print(f"  set system: U={U}, sets={[sorted(c) for c in sets]}, s={s}")
    print(f"  field-vs-joint entropy consistent on all subsets: {field_consistent}")
    print(f"  optimal FREE-position max H(X_S): {best_cov:.0f} bits, "
          f"positions {best_choice}")
    print(f"  greedy max H(X_S): {g_cov:.0f} bits, positions {tuple(g_choice)}")
    print(f"  greedy/opt = {g_cov/best_cov:.3f} >= 1-1/e = {factor:.3f}? {ratio_ok}")
    print("  => FREE per-position max-entropy selection IS max-coverage,")
    print("     NP-hard (Feige 1998); discreteness does NOT soften §6.9")
    print("     hardness. Oracle H(X_S)=|union| is poly here => selection")
    print("     hard even with an easy oracle (separate from V5's #P-oracle).")

    # sanity: {0,1,2} U {3,4,5} covers all 6 -> best should be 6 at s=2.
    full = H_field((0, 2))
    print(f"  (sanity: positions 0,2 cover {full:.0f} of {U} coins)")
    ok = ratio_ok and field_consistent and (abs(best_cov - full) < 1e-9)
    print(f"  V4 {'PASS' if ok else 'FAIL'}")
    return ok


# ======================================================================
# V5: the entropy ORACLE is #P-hard -- explicit reduction marginal
#     PROBABILITY -> entropy oracle (patches inference != oracle gap).
# ======================================================================

def _hbin(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def _hbin_inv_lower(y):
    """Inverse of binary entropy on [0, 1/2] by bisection."""
    lo, hi = 0.0, 0.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if _hbin(mid) < y:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def check_V5():
    print("=" * 70)
    print("V5: entropy ORACLE is #P-hard (reduction marginal-prob -> oracle).")
    print("    Y = E AND fair-bit B  =>  P(Y=1)=P(E)/2 in [0,1/2],")
    print("    H(Y)=h(P(E)/2) injective => oracle recovers #P-hard P(E).")
    print("=" * 70)

    # h strictly increasing on [0,1/2] (injectivity).
    ps = np.linspace(0.0, 0.5, 2000)
    hs = np.array([_hbin(p) for p in ps])
    strictly_inc = bool(np.all(np.diff(hs) > -1e-15))

    # exact recovery of an arbitrary #P-hard marginal P(E) from H(Y) alone.
    max_err = 0.0
    for pE in np.linspace(0.0, 1.0, 201):
        pY = pE / 2.0          # P(Y=1)
        HY = _hbin(pY)         # the only thing the entropy oracle returns
        pE_rec = 2.0 * _hbin_inv_lower(HY)
        max_err = max(max_err, abs(pE - pE_rec))

    print(f"  h strictly increasing on [0,1/2] (injective): {strictly_inc}")
    print(f"  max |P(E) - recovered-from-oracle| = {max_err:.2e}")
    print("  => computing H(Y) solves the #P-hard marginal P(E); the entropy")
    print("     oracle for general AR fields is #P-hard. The poly-time")
    print("     transfer-matrix oracle (V3) is special to low-order/contiguous.")
    ok = strictly_inc and max_err < 1e-6
    print(f"  V5 {'PASS' if ok else 'FAIL'}")
    return ok


# ======================================================================

def main():
    print()
    print("#" * 70)
    print("# Remark 7.15a verification: transformer-predictor support")
    print("# selection -- §6.9 Gaussian-MESP extends to the discrete AR field")
    print("#" * 70)
    print()

    r1 = check_V1()
    print()
    r2 = check_V2()
    print()
    r3 = check_V3()
    print()
    r4 = check_V4()
    print()
    r5 = check_V5()
    print()

    print("=" * 70)
    allok = r1 and r2 and r3 and r4 and r5
    print(f"V1 submodular(no floor)      : {'PASS' if r1 else 'FAIL'}")
    print(f"V2 greedy(1-1/e)+additive    : {'PASS' if r2 else 'FAIL'}")
    print(f"V3 oracle dichotomy          : {'PASS' if r3 else 'FAIL'}")
    print(f"V4 hardness transfers        : {'PASS' if r4 else 'FAIL'}")
    print(f"V5 oracle #P-hard reduction  : {'PASS' if r5 else 'FAIL'}")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
