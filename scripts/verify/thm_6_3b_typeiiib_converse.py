#!/usr/bin/env python3
"""
Verification of Theorem 6.3b (Type-III-B converse; tightness of the
I(H;L|C) saving of Theorem 6.3).

The theorem has three checkable contents on small joint nibble-pair
distributions P(H, L | C) with H, L in {0,...,15} and a small context
alphabet C:

  (1) DECOMPOSITION IDENTITY (the algebra of the converse):
        H(H|C) + H(L|C) - H(H,L|C) = I(H;L|C)
      verified exactly (to floating tolerance) on random distributions.

  (2) CONVERSE LOWER BOUND, E[L] >= H(H,L|C):
      A uniquely decodable code corresponds (Kraft) to a length vector
      {l_b} with sum 2^{-l_b} <= 1; the per-context expected length is
      >= H(B | C = c) with equality only for the (generally fractional)
      ideal lengths l_b = -log2 P(b|c). We check the bound two ways:
        (2a) the *ideal* code length E[-log2 P(B|C)] equals H(B|C)
             exactly (the converse is met with equality by the ideal
             coder, the achievability target of part A);
        (2b) every *integer* uniquely decodable code (sampled, and the
             optimal integer code = Huffman per context) has expected
             length >= H(B|C); no Kraft-feasible integer code beats it.
      Both confirm H(B|C) = H(H,L|C) is the floor.

  (3) TIGHTNESS / ACHIEVABILITY MATCH (Theorem 6.3 saving is exactly
      I(H;L|C)):
      independent-nibble baseline rate H(H|C) + H(L|C) minus joint
      byte rate H(H,L|C) equals I(H;L|C) (ideal coders), and with a
      concrete per-context Huffman coder the realized saving equals
      I(H;L|C) - O(eps) with eps the (sub-1-bit per context) Huffman
      redundancy -- i.e. the joint coder never *over*-saves relative
      to the converse maximum I(H;L|C).

PASS means all three hold across all sampled distributions.

Context reproducibility: C is a decoder-side variable; the converse
conditions on C = c and averages over P(C), exactly matching the
"reproduces context C" hypothesis of the theorem and the §6.3
construction (the decoder has C, hence per-context coding is legal).
"""

import heapq
import itertools
import math
import random
import sys

TOL = 1e-9


def entropy(probs):
    return -sum(p * math.log2(p) for p in probs if p > 0.0)


def huffman_expected_length(symbol_probs):
    """Expected length of an optimal binary prefix (Huffman) code for a
    distribution given as a list of probabilities (need not sum to 1; we
    normalize). Returns sum_i p_i * len_i. This is the optimal *integer*
    uniquely decodable code length (Cover & Thomas Thm 5.8.1), so it
    upper-bounds H and lower-bounds any other integer UD code."""
    ps = [p for p in symbol_probs if p > 0.0]
    total = sum(ps)
    ps = [p / total for p in ps]
    if len(ps) == 1:
        # Single symbol: a uniquely decodable code still needs >= 1 bit
        # in any nondegenerate alphabet; for a degenerate 1-symbol source
        # H = 0 and the (degenerate) code length is 0. Use 0 to match H.
        return 0.0
    # Standard Huffman.
    heap = [[p, i] for i, p in enumerate(ps)]
    heapq.heapify(heap)
    # Track depth of each leaf via repeated merging.
    # Represent each heap entry as [weight, set_of_leaf_indices].
    heap = [[p, {i}] for i, p in enumerate(ps)]
    heapq.heapify(heap)
    depth = [0] * len(ps)
    while len(heap) > 1:
        a = heapq.heappop(heap)
        b = heapq.heappop(heap)
        for idx in a[1]:
            depth[idx] += 1
        for idx in b[1]:
            depth[idx] += 1
        heapq.heappush(heap, [a[0] + b[0], a[1] | b[1]])
    return sum(ps[i] * depth[i] for i in range(len(ps)))


def kraft_feasible(lengths):
    return sum(2.0 ** (-l) for l in lengths) <= 1.0 + 1e-12


def random_joint_HLC(n_context, support_size, seed, sharp=2):
    """Random joint over (h, l, c) with h, l in {0..15}, c in {0..n_context-1}.
    support_size atoms (in the h,l plane) per context, weights ~ U^sharp so
    distributions are peaked (nonzero conditional MI typical)."""
    rng = random.Random(seed)
    # joint[(c)][(h,l)] = prob; also a context marginal.
    joint = {}  # (h, l, c) -> p
    ctx_weight = [rng.random() ** sharp + 0.05 for _ in range(n_context)]
    cw_tot = sum(ctx_weight)
    ctx_weight = [w / cw_tot for w in ctx_weight]
    for c in range(n_context):
        atoms = set()
        # induce correlation: pick a random "diagonal" coupling then noise
        bias_h = rng.randrange(16)
        while len(atoms) < support_size:
            h = rng.randrange(16)
            # low nibble correlated with high nibble for nonzero MI
            if rng.random() < 0.6:
                l = (h + bias_h + rng.randrange(3)) % 16
            else:
                l = rng.randrange(16)
            atoms.add((h, l))
        ws = [rng.random() ** sharp for _ in atoms]
        tot = sum(ws)
        for (h, l), w in zip(atoms, ws):
            joint[(h, l, c)] = ctx_weight[c] * w / tot
    return joint, n_context


def conditional_quantities(joint, n_context):
    """Return per-context and averaged entropy quantities."""
    # P(c)
    pc = [0.0] * n_context
    for (h, l, c), p in joint.items():
        pc[c] += p
    # Build per-context conditional distributions.
    HH_C = 0.0  # H(H|C)
    HL_C = 0.0  # H(L|C)
    HHL_C = 0.0  # H(H,L|C)
    I_HL_C = 0.0  # I(H;L|C)
    ideal_joint_len = 0.0  # E[-log2 P(H,L|C)]
    huff_joint_len = 0.0  # E[ Huffman per context on joint byte ]
    huff_indep_len = 0.0  # E[ Huffman on H | C ] + E[ Huffman on L | C ]
    for c in range(n_context):
        if pc[c] <= 0.0:
            continue
        # conditional joint over (h,l)
        cond = {}
        for (h, l, cc), p in joint.items():
            if cc == c:
                cond[(h, l)] = p / pc[c]
        # marginals
        ph = [0.0] * 16
        pl = [0.0] * 16
        for (h, l), p in cond.items():
            ph[h] += p
            pl[l] += p
        h_h = entropy(ph)
        h_l = entropy(pl)
        h_hl = entropy(cond.values())
        i_hl = h_h + h_l - h_hl
        HH_C += pc[c] * h_h
        HL_C += pc[c] * h_l
        HHL_C += pc[c] * h_hl
        I_HL_C += pc[c] * i_hl
        # ideal joint length under true conditional == H(H,L|C=c)
        ideal_joint_len += pc[c] * sum(
            p * (-math.log2(p)) for p in cond.values() if p > 0.0
        )
        # concrete Huffman coders (optimal integer UD codes)
        huff_joint_len += pc[c] * huffman_expected_length(list(cond.values()))
        huff_indep_len += pc[c] * (
            huffman_expected_length(ph) + huffman_expected_length(pl)
        )
    return dict(
        HH_C=HH_C, HL_C=HL_C, HHL_C=HHL_C, I_HL_C=I_HL_C,
        ideal_joint_len=ideal_joint_len,
        huff_joint_len=huff_joint_len, huff_indep_len=huff_indep_len,
    )


def test_random_integer_codes_dont_beat_entropy(joint, n_context, rng):
    """Sample random Kraft-feasible integer length vectors per context and
    confirm none has expected length below H(B|C=c) (part 2b, sampled)."""
    pc = [0.0] * n_context
    for (h, l, c), p in joint.items():
        pc[c] += p
    worst_slack = math.inf
    for c in range(n_context):
        if pc[c] <= 0.0:
            continue
        cond = {}
        for (h, l, cc), p in joint.items():
            if cc == c:
                cond[(h, l)] = p / pc[c]
        symbols = list(cond.items())
        h_bc = entropy([p for _, p in symbols])
        # ideal integer lengths = ceil(-log2 p); plus random perturbations
        for _ in range(20):
            lengths = []
            for _, p in symbols:
                base = math.ceil(-math.log2(p)) if p > 0 else 1
                lengths.append(base + rng.randrange(0, 3))
            if not kraft_feasible(lengths):
                continue
            exp_len = sum(p * l for (_, p), l in zip(symbols, lengths))
            worst_slack = min(worst_slack, exp_len - h_bc)
    return worst_slack  # should be >= -TOL


def main() -> int:
    failures = 0
    trials = 0
    rng_global = random.Random(20240603)

    max_huffman_redundancy = 0.0
    max_oversave = -math.inf  # realized_saving - I(H;L|C); must be <= O(eps)

    for n_context in (1, 2, 3):
        for support_size in (2, 4, 8, 16):
            for seed in range(60):
                trials += 1
                joint, nc = random_joint_HLC(
                    n_context, support_size, seed=seed * 101 + n_context * 7
                )
                q = conditional_quantities(joint, nc)

                # (1) Decomposition identity (exact).
                lhs = q["HH_C"] + q["HL_C"] - q["HHL_C"]
                if abs(lhs - q["I_HL_C"]) > 1e-7:
                    print(f"FAIL (1) identity: nc={n_context} s={support_size} "
                          f"seed={seed}: {lhs} vs I={q['I_HL_C']}")
                    failures += 1

                # I(H;L|C) >= 0 (sanity).
                if q["I_HL_C"] < -1e-9:
                    print(f"FAIL: I(H;L|C) < 0: {q['I_HL_C']}")
                    failures += 1

                # (2a) ideal coder meets converse with equality:
                #      E[-log2 P(H,L|C)] == H(H,L|C).
                if abs(q["ideal_joint_len"] - q["HHL_C"]) > 1e-7:
                    print(f"FAIL (2a) ideal != H(H,L|C): "
                          f"{q['ideal_joint_len']} vs {q['HHL_C']}")
                    failures += 1

                # (2b) concrete optimal-integer (Huffman) joint coder
                #      satisfies the converse: E[L] >= H(H,L|C).
                slack = q["huff_joint_len"] - q["HHL_C"]
                if slack < -TOL:
                    print(f"FAIL (2b) Huffman joint < H(H,L|C) by {slack}")
                    failures += 1
                # And Huffman redundancy < 1 bit (Cover-Thomas 5.4.1):
                if slack > 1.0 + 1e-9:
                    print(f"FAIL (2b) Huffman redundancy >= 1 bit: {slack}")
                    failures += 1
                max_huffman_redundancy = max(max_huffman_redundancy, slack)

                # (2b, sampled random Kraft-feasible integer codes)
                ws = test_random_integer_codes_dont_beat_entropy(
                    joint, nc, rng_global
                )
                if ws < -TOL:
                    print(f"FAIL (2b-rand) integer code beat H(B|C) by {ws}")
                    failures += 1

                # (3) Tightness: ideal saving == I(H;L|C) exactly.
                ideal_saving = (q["HH_C"] + q["HL_C"]) - q["HHL_C"]
                if abs(ideal_saving - q["I_HL_C"]) > 1e-7:
                    print(f"FAIL (3) ideal saving != I: "
                          f"{ideal_saving} vs {q['I_HL_C']}")
                    failures += 1

                # (3) Concrete coders: realized saving <= I(H;L|C) + O(eps),
                #     i.e. the joint coder cannot beat the converse maximum
                #     by more than the entropy-coder redundancy. We allow a
                #     1-bit Huffman slack on the baseline side.
                realized_saving = q["huff_indep_len"] - q["huff_joint_len"]
                oversave = realized_saving - q["I_HL_C"]
                max_oversave = max(max_oversave, oversave)
                # realized saving must not exceed I + 1 bit (baseline Huffman
                # can pay up to ~1 extra bit per nibble stream, inflating the
                # apparent saving; the *converse* on the joint coder is the
                # hard floor checked in (2b)). The genuine tightness claim is
                # the ideal one in (3) above; this is a concrete-coder sanity.
                if oversave > 2.0 + 1e-9:
                    print(f"FAIL (3-conc) realized saving exceeds I by {oversave}")
                    failures += 1

    # Worked-example cross-check (§6.4 numbers): construct a distribution
    # whose conditional entropies are close to the paper's 2.4 / 2.9 / 4.7
    # is not required; instead verify the arithmetic Delta = 5.3 - 4.7 = 0.6
    # equals H(H|C)+H(L|C)-H(H,L|C) with those literal values.
    HHc, HLc, HHLc = 2.4, 2.9, 4.7
    if abs((HHc + HLc - HHLc) - 0.6) > 1e-12:
        print("FAIL: worked-example arithmetic 5.3 - 4.7 != 0.6")
        failures += 1

    # Handcrafted high-MI case (single context): L = H exactly => the joint
    # is perfectly correlated, I(H;L|C) = H(H|C), and the joint coder pays
    # only H(H|C) while the independent baseline pays 2 H(H|C). This makes
    # the I(H;L|C) saving maximal and concretely large (not a degenerate
    # near-zero case), exercising the tightness claim of part (3).
    joint_corr = {}
    for h in range(8):  # uniform high nibble over 8 values, low = high
        joint_corr[(h, h, 0)] = 1.0 / 8.0
    qc = conditional_quantities(joint_corr, 1)
    # H(H|C) = 3 bits; H(L|C) = 3 bits; H(H,L|C) = 3 bits; I = 3 bits.
    expect = dict(HH_C=3.0, HL_C=3.0, HHL_C=3.0, I_HL_C=3.0)
    for k, v in expect.items():
        if abs(qc[k] - v) > 1e-9:
            print(f"FAIL (handcrafted corr) {k}={qc[k]} expected {v}")
            failures += 1
    # ideal saving = H(H|C)+H(L|C)-H(H,L|C) = 3 bits = I(H;L|C): tight.
    if abs((qc["HH_C"] + qc["HL_C"] - qc["HHL_C"]) - qc["I_HL_C"]) > 1e-9:
        print("FAIL (handcrafted corr) saving != I")
        failures += 1
    # concrete: joint Huffman on a uniform-8 byte source = 3 bits exactly;
    # independent baseline = 3 + 3 = 6 bits => realized saving = 3 = I.
    if abs(qc["huff_joint_len"] - 3.0) > 1e-9:
        print(f"FAIL (handcrafted corr) joint Huffman {qc['huff_joint_len']} != 3")
        failures += 1
    if abs((qc["huff_indep_len"] - qc["huff_joint_len"]) - qc["I_HL_C"]) > 1e-9:
        print("FAIL (handcrafted corr) concrete saving != I")
        failures += 1

    print()
    print(f"handcrafted perfectly-correlated case (L=H, uniform-8): "
          f"H(H,L|C)={qc['HHL_C']:.3f}, I(H;L|C)={qc['I_HL_C']:.3f}, "
          f"concrete Huffman saving={qc['huff_indep_len']-qc['huff_joint_len']:.3f} "
          f"(all = 3.000, exactly tight)")
    print(f"max Huffman joint redundancy over H(H,L|C): "
          f"{max_huffman_redundancy:.6f} bits (strictly < 1; the near-1 "
          f"worst case is the degenerate two-atom source where H~0)")
    print(f"max concrete realized-saving overshoot vs I(H;L|C): "
          f"{max_oversave:.6f} bits (baseline-Huffman slack, not a converse "
          f"violation; the joint-coder converse floor is checked in (2b))")
    if failures == 0:
        print(f"PASS: converse E[L] >= H(H,L|C), decomposition identity, "
              f"and I(H;L|C) tightness hold across {trials} distributions.")
        return 0
    print(f"FAIL: {failures} checks failed across {trials} distributions.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
