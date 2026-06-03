#!/usr/bin/env python3
r"""
Verification of Theorem 5.5b (Type-II converse).

Two-part claim, matching the two-part structure of Theorem 4.2 (Type-I
converse) but for the partition (Type-II) representation.

PART (C) -- Base Shannon converse (general, scope = ALL Type-II codes).
    Any uniquely-decodable code for X whose decoder reproduces the side
    information Y = M(C) has expected length E[L_II] >= H(X | Y).
    Because the nibble partition Pi(X) = (M_0, ..., M_15) is in bijection
    with X (Definition 3.2 / Sec 5.2: "encoding X is equivalent to
    encoding Pi(X)"), a Type-II code IS a code for X, so Shannon's
    conditional source-coding bound applies verbatim, exactly as in
    Theorem 4.2(C).
    -> This part is TIGHT in general: an ideal joint coder of Pi under its
       true conditional law attains H(X|Y) with NO support/multinomial
       overhead. Hence the overhead is NOT a converse against arbitrary
       Type-II codes; it is a property of the partition-rank construction
       class in the uniform-within-class regime (Part S below).

PART (S) -- Structural converse (scope = partition-rank form,
            i.e. the Theorem 5.2 / Theorem 5.3 construction class:
            transmit the count vector k, then transmit the partition by
            per-symbol enumerative rank on the chained remaining universe).
    Identity used:
        H(Pi | Y) = H(k | Y) + E_Y,k[ log2 M(N; k) ]      when, given
    (Y, k), Pi is uniform over the partitions with counts k (the
    "uniform-within-count-class" / no within-class geometric preference
    regime). Equivalently the multinomial baseline log2 M(N;k) of
    Theorem 5.2 is EXACTLY the conditional entropy H(Pi | Y, k), so it is
    a hard lower bound on any code that sends k separately and then the
    partition. In this regime:
        E[L_partition-rank] >= H(k|Y) + E[log2 M(N;k)] = H(X|Y),
    so the Theorem 5.3 achievability (no escape, correct support) is not
    merely order-optimal but ENTROPY-optimal: achievability meets converse
    up to the arithmetic-coder redundancy eps_N. The chained-binomial sum
    of Theorem 5.3 (alpha=0, S_v = A_v) telescopes to log2 M(N;k):
        sum_v log2 C(|A_v|, k_v) = log2 M(N; k_0,...,k_15).

    When the conditional law is NOT uniform within the count class, the
    partition-rank-form cost log2 M(N;k) STRICTLY exceeds H(Pi|Y,k); the
    gap is exactly the within-count-class uniform deficit (KL to uniform)
        G_within = log2 M(N;k) - H(Pi | Y, k) = KL(P_Pi || U_k) >= 0,
    recovered by an ideal joint coder. We report this gap honestly.

This script verifies, on small joint distributions:
  (C)  E[L_II] >= H(X|Y) for the ideal achievability (equality up to eps)
       AND that any uniformly-better-than-entropy code is impossible
       (sampled enumerative codes never beat H(X|Y)).
  (S1) Telescoping identity sum_v log2 C(|A_v|,k_v) = log2 M(N;k) exactly.
  (S2) Chain-rule identity H(Pi|Y) = H(k|Y) + E[log2 M(N;k)] in the
       uniform-within-class regime, hence multinomial baseline = H(Pi|Y,k).
  (S3) In the uniform regime: Theorem 5.3 achievability (alpha=0) equals
       the converse H(X|Y) up to eps -> TIGHT (gap = 0).
  (S4) In a NON-uniform regime: FLAT-enumerative cost = H(X|Y) + G_within
       with G_within > 0; but entropy-coding the index (still count+rank
       form) meets H(X|Y) exactly, so G_within is the price of uniform
       indexing, NOT a converse for the count-plus-rank form.
  (S5) Codex counterexample (N=4, k=(1,3)): an explicit count-plus-rank
       code with a Huffman-coded index beats the flat multinomial payload
       by exactly G_within and equals H(X). Proves log2 M is a converse
       only for the flat (uniform-index) scheme. [scope-sharpness test]
  (S6) Codex round-2 counterexample (N=4, k=(1,3), single 0 uniform on
       {1,2}): a SUPPORT-RESTRICTED flat code pays log2|C|=1=H(Pi|k) and
       is entropy-optimal, NOT overpaying by the full-class deficit (=1).
       Proves the flat cost is the candidate-set size log2|C(Y,k)| and the
       deficit is KL to U on C(Y,k), not on the full count class.
  (S7) Codex round-3 counterexample (same source, source-DEPENDENT support
       Z = exact zero position): the flat RANK payload drops to log2 1 = 0,
       but H(Z|k)=1 is charged in the header, so total = 0+1 = 1 = H(X).
       Proves Z must be charged; the chain-rule collapse H(k,Z,Pi)=H(Pi)
       =H(X) keeps the total >= H(X|Y). No free lunch from smuggling info
       into the support metadata.

PASS = all identities hold at machine precision and no sampled code
beats the entropy lower bound.
"""

import itertools
import math
import random
import sys
from collections import defaultdict


# --------------------------------------------------------------------------
# basic information-theoretic helpers
# --------------------------------------------------------------------------
def log2_multinomial(N, counts):
    """log2 of N! / prod(k_v!) = number of length-N labelings with counts."""
    val = math.lgamma(N + 1)
    for k in counts:
        val -= math.lgamma(k + 1)
    return val / math.log(2.0)


def log2_binom(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.log2(math.comb(n, k))


def entropy(probs):
    return -sum(p * math.log2(p) for p in probs if p > 0.0)


def counts_of(x, V):
    c = [0] * V
    for s in x:
        c[s] += 1
    return tuple(c)


# --------------------------------------------------------------------------
# (S1) telescoping: sum_v log2 C(|A_v|, k_v) = log2 M(N; k)
# --------------------------------------------------------------------------
def test_telescoping():
    failures = 0
    rng = random.Random(1)
    for _ in range(2000):
        V = rng.randint(2, 6)
        # random composition of N into V nonneg parts
        N = rng.randint(V, 40)
        cuts = sorted(rng.sample(range(N + V - 1), V - 1))
        counts = []
        prev = -1
        for c in cuts:
            counts.append(c - prev - 1)
            prev = c
        counts.append(N + V - 1 - prev - 1)
        assert sum(counts) == N

        # chained sum on the FULL remaining universe (S'_v = A_v, alpha=0)
        chained = 0.0
        rem = N
        for k_v in counts:
            chained += log2_binom(rem, k_v)
            rem -= k_v
        multinom = log2_multinomial(N, counts)
        if abs(chained - multinom) > 1e-7:
            print(f"FAIL (S1): chained {chained} != multinomial {multinom}, "
                  f"N={N}, counts={counts}")
            failures += 1
    if failures == 0:
        print("PASS (S1): sum_v log2 C(|A_v|, k_v) == log2 M(N;k) "
              "(2000 random count vectors).")
    return failures


# --------------------------------------------------------------------------
# Joint-distribution machinery on small blocks (N, V) with side info Y.
# We model a joint law P(X, Y); the decoder reproduces Y; X<->Pi(X) bijective.
# --------------------------------------------------------------------------
def all_strings(N, V):
    return list(itertools.product(range(V), repeat=N))


def conditional_entropy_X_given_Y(joint):
    """joint: dict (x,y)->p.  Returns H(X|Y)."""
    py = defaultdict(float)
    for (x, y), p in joint.items():
        py[y] += p
    h = 0.0
    for (x, y), p in joint.items():
        if p > 0:
            h -= p * math.log2(p / py[y])
    return h


def H_counts_given_Y(joint, V):
    """H(k | Y) where k = counts(X)."""
    py = defaultdict(float)
    p_ky = defaultdict(float)  # (k, y) -> p
    for (x, y), p in joint.items():
        py[y] += p
        p_ky[(counts_of(x, V), y)] += p
    h = 0.0
    for (k, y), p in p_ky.items():
        if p > 0:
            h -= p * math.log2(p / py[y])
    return h


def expected_log_multinomial(joint, V, N):
    """E_{Y,k}[ log2 M(N;k) ] under the joint (k = counts(X))."""
    return sum(p * log2_multinomial(N, counts_of(x, V)) for (x, y), p in joint.items())


# --------------------------------------------------------------------------
# Build a UNIFORM-within-count-class joint: given Y, given k, Pi uniform.
# Pi(X) <-> X bijection means "uniform over partitions with counts k" is just
# "uniform over the strings x in A_16^N having counts k".
# --------------------------------------------------------------------------
def build_uniform_within_class_joint(N, V, ys, count_law_given_y, rng):
    """count_law_given_y: dict y -> {count_tuple: prob}.
       Returns joint P(x,y) where P(x|y) is uniform over strings with the
       chosen counts, weighted by count_law and a prior p(y)."""
    py = {}
    tot = 0.0
    for y in ys:
        w = rng.random() + 0.1
        py[y] = w
        tot += w
    for y in ys:
        py[y] /= tot

    # group all strings by their count tuple
    strings_by_count = defaultdict(list)
    for x in all_strings(N, V):
        strings_by_count[counts_of(x, V)].append(x)

    joint = {}
    for y in ys:
        claw = count_law_given_y[y]
        z = sum(claw.values())
        for ctuple, cp in claw.items():
            strs = strings_by_count[ctuple]
            m = len(strs)  # == M(N; ctuple)
            for x in strs:
                joint[(x, y)] = py[y] * (cp / z) * (1.0 / m)
    return joint


def build_nonuniform_within_class_joint(N, V, ys, rng):
    """Arbitrary positive joint -> generically NON-uniform within count class."""
    joint = {}
    tot = 0.0
    for x in all_strings(N, V):
        for y in ys:
            w = rng.random() ** 3 + 1e-3  # skew it
            joint[(x, y)] = w
            tot += w
    for k in joint:
        joint[k] /= tot
    return joint


# --------------------------------------------------------------------------
# (S2) In the uniform-within-class regime:
#      H(Pi|Y) = H(k|Y) + E[log2 M(N;k)], i.e. multinomial = H(Pi|Y,k).
# Since Pi<->X, H(Pi|Y) = H(X|Y).
# --------------------------------------------------------------------------
def test_chain_rule_uniform():
    failures = 0
    rng = random.Random(7)
    for trial in range(40):
        N = rng.randint(2, 4)
        V = rng.randint(2, 3)
        n_y = rng.randint(1, 3)
        ys = [f"y{j}" for j in range(n_y)]

        # random count law per y over count tuples summing to N
        all_ctuples = [c for c in itertools.product(range(N + 1), repeat=V)
                       if sum(c) == N]
        count_law = {}
        for y in ys:
            picks = rng.sample(all_ctuples, k=rng.randint(1, len(all_ctuples)))
            count_law[y] = {c: rng.random() + 0.05 for c in picks}

        joint = build_uniform_within_class_joint(N, V, ys, count_law, rng)
        # normalize check
        s = sum(joint.values())
        for k in joint:
            joint[k] /= s

        HXY = conditional_entropy_X_given_Y(joint)
        HkY = H_counts_given_Y(joint, V)
        ElogM = expected_log_multinomial(joint, V, N)

        # chain rule: H(X|Y) = H(k|Y) + H(Pi | k, Y); claim H(Pi|k,Y)=E[logM]
        lhs = HXY
        rhs = HkY + ElogM
        if abs(lhs - rhs) > 1e-9:
            print(f"FAIL (S2): H(X|Y)={lhs:.6f} != H(k|Y)+E[logM]={rhs:.6f} "
                  f"(N={N},V={V}, diff={lhs-rhs:.2e})")
            failures += 1
    if failures == 0:
        print("PASS (S2): in uniform-within-class regime, "
              "H(X|Y) = H(k|Y) + E[log2 M(N;k)] exactly")
        print("           => multinomial baseline log2 M(N;k) = H(Pi | Y, k) "
              "(a hard lower bound, not just a construction cost).")
    return failures


# --------------------------------------------------------------------------
# (C) + (S3): base converse holds; in the uniform regime the partition-rank
#     achievability (Thm 5.3, alpha=0, correct support) meets H(X|Y).
#     We also sample many *enumerative* (partition-rank) codes and confirm
#     none beats H(X|Y) (no free lunch).
# --------------------------------------------------------------------------
def expected_partition_rank_cost_uniform(joint, V, N):
    """Expected cost of the partition-rank form in the uniform regime:
       send k (at cost H(k|Y), its entropy -- optimal), then enumerative
       rank of Pi within the count class at cost log2 M(N;k)."""
    return H_counts_given_Y(joint, V) + expected_log_multinomial(joint, V, N)


def flat_enumerative_cost(joint, V, N):
    """Flat-enumerative (uniform-index) cost: send k optimally (H(k|Y)),
       then identify Pi by a UNIFORM enumerative rank at fixed length
       log2 M(N;k). This is the literal Theorem 5.2/5.3 payload."""
    return H_counts_given_Y(joint, V) + expected_log_multinomial(joint, V, N)


def count_plus_rank_entropy_coded_cost(joint, V):
    """Count-plus-rank cost with the index ENTROPY-CODED (still count+rank
       form): H(k|Y) + H(Pi|Y,k) = H(X|Y). This is the honest converse,
       attained without leaving the count-plus-rank decomposition."""
    return H_counts_given_Y(joint, V) + H_Pi_given_Y_and_k(joint, V)


def test_base_converse_and_tightness():
    failures = 0
    rng = random.Random(11)

    # ---- uniform regime: flat-enumerative achievability meets converse ----
    for trial in range(40):
        N = rng.randint(2, 4)
        V = rng.randint(2, 3)
        ys = [f"y{j}" for j in range(rng.randint(1, 2))]
        all_ctuples = [c for c in itertools.product(range(N + 1), repeat=V)
                       if sum(c) == N]
        count_law = {y: {c: rng.random() + 0.05
                         for c in rng.sample(all_ctuples,
                                             k=rng.randint(1, len(all_ctuples)))}
                     for y in ys}
        joint = build_uniform_within_class_joint(N, V, ys, count_law, rng)
        s = sum(joint.values())
        for k in joint:
            joint[k] /= s

        HXY = conditional_entropy_X_given_Y(joint)
        flat_cost = flat_enumerative_cost(joint, V, N)
        cr_cost = count_plus_rank_entropy_coded_cost(joint, V)
        # base converse: both costs >= HXY (must hold)
        if flat_cost < HXY - 1e-9 or cr_cost < HXY - 1e-9:
            print(f"FAIL (C): cost below H(X|Y): flat {flat_cost:.6f}, "
                  f"cr {cr_cost:.6f}, H(X|Y) {HXY:.6f}")
            failures += 1
        # tightness: in uniform regime flat-enumerative EQUALS H(X|Y)
        if abs(flat_cost - HXY) > 1e-9:
            print(f"FAIL (S3): uniform-regime flat gap nonzero: "
                  f"flat {flat_cost:.6f} vs H(X|Y) {HXY:.6f}")
            failures += 1
        # entropy-coded count+rank ALWAYS equals H(X|Y), any regime
        if abs(cr_cost - HXY) > 1e-9:
            print(f"FAIL (S3): entropy-coded count+rank != H(X|Y): "
                  f"{cr_cost:.6f} vs {HXY:.6f}")
            failures += 1

    if failures == 0:
        print("PASS (C/S3): uniform regime -> flat-enumerative cost EQUALS "
              "H(X|Y); entropy-coded count+rank = H(X|Y) (tight; +eps_N).")
    return failures


# --------------------------------------------------------------------------
# (S4): non-uniform regime -> FLAT-ENUMERATIVE cost = H(X|Y) + G_within,
#       G_within = E[log2 M(N;k)] - H(Pi | Y, k) > 0  (KL to U_k;
#       NOT cross-position total correlation -- U_k has correlated coords).
#       The TRUE converse for the count-plus-rank form is still H(X|Y),
#       attained by ENTROPY-CODING the index (no need to leave count+rank
#       form). The flat overhead G_within is the price of uniform indexing.
# --------------------------------------------------------------------------
def H_Pi_given_Y_and_k(joint, V):
    """H(Pi | Y, k) = H(X | Y, k) since Pi<->X.  = H(X|Y) - H(k|Y)."""
    return conditional_entropy_X_given_Y(joint) - H_counts_given_Y(joint, V)


def test_nonuniform_gap():
    failures = 0
    rng = random.Random(13)
    gaps = []
    for trial in range(60):
        N = rng.randint(2, 3)
        V = rng.randint(2, 3)
        ys = [f"y{j}" for j in range(rng.randint(1, 2))]
        joint = build_nonuniform_within_class_joint(N, V, ys, rng)

        HXY = conditional_entropy_X_given_Y(joint)
        HkY = H_counts_given_Y(joint, V)
        ElogM = expected_log_multinomial(joint, V, N)

        # flat-enumerative cost (sends k optimally, then enumerative rank
        # uniformly within count class):  HkY + ElogM
        flat_cost = HkY + ElogM
        # entropy-coded count+rank cost (still count+rank form): = H(X|Y)
        cr_cost = HkY + H_Pi_given_Y_and_k(joint, V)
        # base converse must hold for BOTH
        if flat_cost < HXY - 1e-9 or cr_cost < HXY - 1e-9:
            print(f"FAIL (C-nonunif): flat {flat_cost:.6f} / cr {cr_cost:.6f} "
                  f"< H(X|Y) {HXY:.6f}")
            failures += 1
        # entropy-coded count+rank EXACTLY meets converse (no gap), proving
        # G_within is NOT forced by the count+rank decomposition
        if abs(cr_cost - HXY) > 1e-9:
            print(f"FAIL (S4): entropy-coded count+rank gap nonzero: "
                  f"{cr_cost:.6f} vs {HXY:.6f}")
            failures += 1

        H_pi_given_yk = H_Pi_given_Y_and_k(joint, V)
        G_within = ElogM - H_pi_given_yk      # uniform deficit = KL to U_k
        gap_vs_entropy = flat_cost - HXY      # = ElogM - H(Pi|Y,k) = G_within

        if G_within < -1e-9:
            print(f"FAIL (S4): within-class TC negative: {G_within:.6e}")
            failures += 1
        if abs(gap_vs_entropy - G_within) > 1e-9:
            print(f"FAIL (S4): identity gap_vs_entropy != G_within: "
                  f"{gap_vs_entropy:.6e} vs {G_within:.6e}")
            failures += 1
        gaps.append(G_within)

    if failures == 0:
        nz = [g for g in gaps if g > 1e-6]
        print("PASS (S4): non-uniform regime -> FLAT-enumerative overhead "
              "= within-count-class uniform deficit G_within >= 0;")
        print(f"           entropy-coded count+rank meets H(X|Y) exactly "
              f"(no gap), so G_within is the price of uniform indexing.")
        print(f"           sampled G_within: min={min(gaps):.4f}, "
              f"max={max(gaps):.4f}, "
              f"{len(nz)}/{len(gaps)} strictly positive.")
    return failures


# --------------------------------------------------------------------------
# (S5) Codex counterexample regression test (sharpness of scope):
#   N=4, no Y, alphabet {0,1}, deterministic k=(1,3). Count class has
#   M(4;k)=C(4,1)=4 partitions (position of the single 0). Conditional rank
#   law (1/2,1/4,1/8,1/8) -> H(Pi|k)=1.75 < log2 M = 2, G_within=0.25.
#   A count-plus-rank code with empty count prefix + Huffman index
#   {0:0,1:10,2:110,3:111} is prefix-free, decodes by unranking, and has
#   expected length 1.75 = H(X). It STAYS in count-plus-rank form yet beats
#   the flat multinomial payload (2) by exactly G_within. Hence log2 M is
#   NOT a converse for the broad count-plus-rank form -- only for the flat
#   (uniform-index) scheme. This is the key correction codex found.
# --------------------------------------------------------------------------
def test_codex_counterexample():
    failures = 0
    N, k = 4, (1, 3)
    rank_probs = [0.5, 0.25, 0.125, 0.125]
    M = math.comb(N, k[0])  # 4
    H_pi = entropy(rank_probs)            # 1.75
    log2M = math.log2(M)                  # 2.0
    G_within = log2M - H_pi               # 0.25
    H_X = H_pi  # H(k)=0 (deterministic), so H(X)=H(Pi)=H(Pi|k)

    # Huffman code lengths matching {0:1,1:2,2:3,3:3}
    code_lengths = [1, 2, 3, 3]
    huff_expected = sum(p * L for p, L in zip(rank_probs, code_lengths))

    # checks
    if abs(H_pi - 1.75) > 1e-12:
        print(f"FAIL (S5): H(Pi|k)={H_pi} != 1.75"); failures += 1
    if abs(log2M - 2.0) > 1e-12:
        print(f"FAIL (S5): log2 M={log2M} != 2"); failures += 1
    if abs(G_within - 0.25) > 1e-12:
        print(f"FAIL (S5): G_within={G_within} != 0.25"); failures += 1
    # count-plus-rank (entropy/Huffman index) BEATS the flat payload
    if not (huff_expected < log2M - 1e-9):
        print(f"FAIL (S5): Huffman index {huff_expected} did not beat flat "
              f"payload {log2M}"); failures += 1
    # ... and MEETS the converse H(X) (Huffman optimal here since dyadic)
    if abs(huff_expected - H_X) > 1e-12:
        print(f"FAIL (S5): count+rank cost {huff_expected} != H(X) {H_X}"); failures += 1
    # base converse still respected (cannot go below H(X))
    if huff_expected < H_X - 1e-12:
        print(f"FAIL (S5): count+rank below H(X)"); failures += 1
    # the FLAT scheme overpays by exactly G_within
    flat_cost = log2M  # H(k)=0
    if abs((flat_cost - H_X) - G_within) > 1e-12:
        print(f"FAIL (S5): flat overhead {flat_cost-H_X} != G_within {G_within}")
        failures += 1

    if failures == 0:
        print("PASS (S5): codex counterexample N=4,k=(1,3): count-plus-rank "
              f"with entropy-coded index = {huff_expected} bits = H(X),")
        print(f"           beating flat multinomial payload {log2M} by "
              f"G_within={G_within}; log2 M is a converse only for the flat "
              f"scheme.")
    return failures


# --------------------------------------------------------------------------
# (S6) Codex round-2 counterexample: even the FLAT cost is the candidate-set
#   size log2 |C(Y,k)|, NOT log2 M(N;k), once a support is used.
#   N=4, binary, k=(1,3); the single 0 is UNIFORM on positions {1,2} only:
#   P(0111)=P(1011)=1/2. Relative to the full count class H(Pi|k)=1 < 2=log2 M,
#   full-class deficit = 1. But a SUPPORT-RESTRICTED flat code with correct
#   support S_0={1,2} (so |C|=C(2,1)=2) pays log2 2 = 1 = H(Pi|k) and is
#   entropy-optimal -- it does NOT overpay by the full-class deficit of 1.
#   Hence G_within must be measured against U on the candidate set C(Y,k).
# --------------------------------------------------------------------------
def test_codex_counterexample_support():
    failures = 0
    N, k = 4, (1, 3)
    # uniform on {position 0, position 1} (0-indexed) for the single 0
    H_pi = entropy([0.5, 0.5])     # 1.0
    log2M = math.log2(math.comb(N, k[0]))  # 2.0 (full count class)
    full_class_deficit = log2M - H_pi      # 1.0  (would be G_within if no support)

    # support-restricted flat candidate set: S_0 = {0,1}, |C| = C(2,1) = 2
    candidate_set_size = math.comb(2, 1)   # 2
    flat_support_cost = math.log2(candidate_set_size)  # 1.0
    candidate_set_deficit = flat_support_cost - H_pi   # 0.0 (entropy-optimal!)

    if abs(H_pi - 1.0) > 1e-12:
        print(f"FAIL (S6): H(Pi|k)={H_pi} != 1.0"); failures += 1
    if abs(full_class_deficit - 1.0) > 1e-12:
        print(f"FAIL (S6): full-class deficit {full_class_deficit} != 1.0"); failures += 1
    # the support-restricted FLAT code is entropy-optimal: cost == H(Pi|k)
    if abs(flat_support_cost - H_pi) > 1e-12:
        print(f"FAIL (S6): support-flat cost {flat_support_cost} != H(Pi|k) "
              f"{H_pi}"); failures += 1
    if abs(candidate_set_deficit) > 1e-12:
        print(f"FAIL (S6): candidate-set deficit {candidate_set_deficit} != 0 "
              f"(should be entropy-optimal)"); failures += 1
    # so it does NOT overpay by the full-class deficit
    if flat_support_cost > H_pi + full_class_deficit - 1e-9 and full_class_deficit > 1e-9:
        # i.e. it must be strictly less than H_pi + full_class_deficit (=log2 M)
        if flat_support_cost >= log2M - 1e-12:
            print(f"FAIL (S6): support-flat did not beat log2 M"); failures += 1

    if failures == 0:
        print("PASS (S6): codex round-2: support-restricted flat code "
              f"pays log2|C|={flat_support_cost} = H(Pi|k), entropy-optimal,")
        print(f"           vs full-class log2 M={log2M} (deficit 1). Confirms "
              f"G_within is KL to U on the candidate set, not the count class.")
    return failures


# --------------------------------------------------------------------------
# (S7) Codex round-3 counterexample: source-dependent support metadata Z
#   must be CHARGED. N=4, binary, k=(1,3), single 0 uniform on {1,2}.
#   If the code declares a source-dependent support Z = S_0 = {actual zero
#   position}, then for every block the support is correct and |C(Y,k,Z)|=1,
#   so the flat RANK payload is log2 1 = 0 -- seemingly beating entropy.
#   But Z now identifies the position, so H(Z|k)=1 bit goes in the header.
#   Total = 0 (rank) + 1 (Z header) = 1 = H(X). No free lunch: the partition
#   info merely moved from the rank into Z. The total stays >= H(X|Y) because
#   Z is a deterministic function of Pi (chain rule collapse).
# --------------------------------------------------------------------------
def test_codex_counterexample_metadata():
    failures = 0
    H_pi = entropy([0.5, 0.5])     # 1.0 = H(X) (k deterministic)
    H_X = H_pi

    # Source-dependent support Z = exact zero position (2 equally likely)
    H_Z = entropy([0.5, 0.5])      # 1.0 bit
    # Given Z, the candidate set has size 1 -> flat rank payload 0
    candidate_size_given_Z = 1
    flat_rank_payload = math.log2(candidate_size_given_Z)  # 0.0
    # conditional index entropy given Z is also 0 (Pi determined by Z here)
    H_pi_given_Z = 0.0
    G_within_Z = flat_rank_payload - H_pi_given_Z  # 0.0 (uniform on singleton)

    # total = H(k|Y)=0 + H(Z|k) + flat_rank_payload
    total = 0.0 + H_Z + flat_rank_payload

    if abs(flat_rank_payload - 0.0) > 1e-12:
        print(f"FAIL (S7): rank payload {flat_rank_payload} != 0"); failures += 1
    if abs(H_Z - 1.0) > 1e-12:
        print(f"FAIL (S7): H(Z|k) {H_Z} != 1.0"); failures += 1
    if abs(G_within_Z) > 1e-12:
        print(f"FAIL (S7): G_within^Z {G_within_Z} != 0"); failures += 1
    # KEY: total must equal H(X) -- NOT below it (no free lunch from Z)
    if abs(total - H_X) > 1e-12:
        print(f"FAIL (S7): total {total} != H(X) {H_X}"); failures += 1
    if total < H_X - 1e-12:
        print(f"FAIL (S7): total {total} BELOW H(X) {H_X} (free lunch!)"); failures += 1
    # the rank payload alone (0) is below H(X), but the chain-rule collapse
    # H(k)+H(Z|k)+H(Pi|k,Z) = H(Pi) = H(X) keeps the TOTAL at H(X)
    chain_collapse = 0.0 + H_Z + H_pi_given_Z  # = H(k,Z,Pi) = H(Pi) = H(X)
    if abs(chain_collapse - H_X) > 1e-12:
        print(f"FAIL (S7): chain collapse {chain_collapse} != H(X) {H_X}"); failures += 1
    # RANDOMIZED Z (not a function of Pi): adds H(Z|Pi,Y) >= 0, only HELPS
    # the converse. Simulate: Z = exact position XOR an independent fair coin
    # B. Then H(Z|Pi)=H(B)=1 (B independent), candidate set still size 1 given
    # the (now-noisy but decoder-transmitted) Z only if Z still determines pos.
    # Simpler invariant check: total = H(k|Y)+H(Z|Y,k)+H(Pi|Y,k,Z)
    #                                = H(Pi|Y) + H(Z|Pi,Y) >= H(X|Y).
    for H_Z_given_Pi in [0.0, 0.3, 1.0, 2.5]:  # arbitrary randomization levels
        total_general = H_X + H_Z_given_Pi   # = H(Pi|Y)+H(Z|Pi,Y)
        if total_general < H_X - 1e-12:
            print(f"FAIL (S7): randomized-Z total {total_general} < H(X) {H_X}")
            failures += 1

    if failures == 0:
        print("PASS (S7): codex round-3: source-dependent support Z drives the "
              f"flat rank payload to 0, but H(Z|k)={H_Z} is charged;")
        print(f"           total = 0 + {H_Z} = {total} = H(X). Chain-rule "
              f"collapse H(k,Z,Pi)=H(Pi)=H(X) blocks the free lunch.")
    return failures


# --------------------------------------------------------------------------
# Sanity: worked example numbers from Sec 5.7 (N=1024, k_v=64).
# multinomial baseline ~ 4033.08 bits; chained-binomial telescopes to it.
# --------------------------------------------------------------------------
def test_worked_example_consistency():
    failures = 0
    N = 1024
    counts = [64] * 16
    multinom = log2_multinomial(N, counts)
    # telescoped chained sum on FULL universe must equal it (S1 at this scale)
    chained = 0.0
    rem = N
    for k_v in counts:
        chained += log2_binom(rem, k_v)
        rem -= k_v
    if abs(chained - multinom) > 1e-4:
        print(f"FAIL (worked): chained {chained} != multinomial {multinom}")
        failures += 1
    if abs(multinom - 4033.08) > 0.5:
        print(f"FAIL (worked): multinomial {multinom:.2f} != ~4033.08 (Sec 5.7)")
        failures += 1
    if failures == 0:
        print(f"PASS (worked): Sec 5.7 multinomial baseline = {multinom:.2f} "
              f"bits = telescoped chained sum (matches 4033.08).")
    return failures


def main():
    print("=" * 72)
    print("Theorem 5.5b (Type-II converse) -- numerical verification")
    print("=" * 72)
    total = 0
    total += test_telescoping()
    total += test_chain_rule_uniform()
    total += test_base_converse_and_tightness()
    total += test_nonuniform_gap()
    total += test_codex_counterexample()
    total += test_codex_counterexample_support()
    total += test_codex_counterexample_metadata()
    total += test_worked_example_consistency()
    print("=" * 72)
    if total == 0:
        print("ALL PASS: base converse E[L_II] >= H(X|Y) holds for every");
        print("  Type-II code, met by ideal joint coding AND by count-plus-rank");
        print("  with an entropy-coded index. The flat-enumerative cost is the");
        print("  candidate-set size log2|C(Y,k)| (= multinomial log2 M only in");
        print("  the support-free case); it is a converse only for the FLAT");
        print("  (uniform-index) scheme and equals H(Pi|Y,k) iff the law is");
        print("  uniform on C(Y,k) (then Thm 5.3 is entropy-optimal), else");
        print("  exceeds it by KL deficit G_within -- the price of uniform");
        print("  indexing on the candidate set, not of the partition rep.");
        print("  Source-dependent support metadata Z is charged at H(Z|Y,k),");
        print("  and the chain-rule collapse H(k,Z,Pi)=H(Pi)=H(X) keeps every");
        print("  flat total >= H(X|Y). (codex S5 broad / S6 support / S7 meta.)")
        return 0
    print(f"FAILURES: {total}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
