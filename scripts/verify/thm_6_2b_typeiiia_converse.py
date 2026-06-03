#!/usr/bin/env python3
"""
Verify Theorem 6.2b: Type-III-A absolute-rate converse and optimality conditions.

Claims checked numerically on small explicit joint distributions:

  (6.2b.1) E[L_IIIA] >= H(X|Y)                          -- Shannon floor
  (6.2b.2) R_tile := sum_t E_{C_t}[ min_m E[L_m(t)|C_t] ]    (CAUSAL-ADAPTIVE)
                  >= sum_t H(X_t|C_t)
                  >= sum_t H(X_t|X_<t,Y) = H(X|Y)
            (a) first ineq  equality iff per-context cost-min mode expressive
            (b) second ineq equality iff mode contexts causally sufficient
            (c) last step is the chain rule, EXACT
  (6.2b.3) context deficit = sum_t [H(X_t|C_t) - H(X_t|X_<t,Y)] >= 0;
           for BLOCK-INDEPENDENT modes (C_t = Y) this equals the inter-tile
           total correlation TC(X_1..X_T | Y).

KEY HONEST FINDING (Remark 6.2b): the gap is the *context deficit*, NOT an
inherent per-tile penalty. On a dependent-tile source:
  * block-independent modes (C_t = Y) overpay by exactly TC(.|Y) > 0;
  * causal-context modes (C_t = (Y, X_<t)) -- legitimate because the §6.2
    mode set includes the Type-III-C dictionary and Type-I streaming, both
    cross-tile -- close the gap EXACTLY to H(X|Y) (chain rule).
So TC(.|Y) is the worst case, removable by autoregressive conditioning.

REGRESSION TEST (Source C, Remark 6.2b'): the OPERATIVE benchmark is the
causal-adaptive floor (ii), NOT the global-commit cost (iii)
sum_t min_m E[L_m(t)]. An earlier draft used (iii) and was refuted: with
context-dependent mode selection the operative scheme legally beats (iii)
flag-free. Source C is the refuting example -- global-commit OVER-estimates
(3.1875) while the causal-adaptive selector attains the floor H(X)=2.75.
The three benchmarks obey (i) hindsight <= (ii) causal-adaptive <= (iii)
global-commit; only (ii) is the operative rate and only (ii) is a valid floor.

PASS iff all relations hold to tolerance, the block-independent floor on the
dependent source exhibits the correctly-sized TC gap, the causal-context floor
attains H(X|Y), and Source C exhibits (i)<=(ii)<=(iii) with (ii)==H(X)<(iii).
"""
import math

TOL = 1e-9


def H(p):
    """Shannon entropy (bits) of a flat iterable of probabilities."""
    return -sum(x * math.log2(x) for x in p if x > 0.0)


def marginal(joint, axes, shape):
    """Marginalize a flat joint over the kept `axes` (tuple of dim indices)."""
    out = {}
    for idx, pr in joint.items():
        key = tuple(idx[a] for a in axes)
        out[key] = out.get(key, 0.0) + pr
    return out


def cond_entropy(joint, X_axes, cond_axes, shape):
    """H(X_axes | cond_axes) = H(X,cond) - H(cond) for a flat-dict joint."""
    XC = marginal(joint, tuple(X_axes) + tuple(cond_axes), shape)
    C = marginal(joint, tuple(cond_axes), shape)
    return H(XC.values()) - H(C.values())


def total_correlation_given_Y(joint, tile_axes, Y_axis, shape):
    """sum_t H(X_t|Y) - H(X_1..X_T | Y)."""
    sum_marg = sum(cond_entropy(joint, [t], [Y_axis], shape) for t in tile_axes)
    joint_cond = cond_entropy(joint, list(tile_axes), [Y_axis], shape)
    return sum_marg - joint_cond, sum_marg, joint_cond


def oracle_blockindep(joint, tile_axes, Y_axis, shape):
    """
    BLOCK-INDEPENDENT modes (context C_t = Y only): best mode for tile t
    achieves E[L_m(t)] = H(X_t | Y). Oracle = sum_t H(X_t|Y).
    """
    return sum(cond_entropy(joint, [t], [Y_axis], shape) for t in tile_axes)


def oracle_causal(joint, tile_axes, Y_axis, shape):
    """
    CAUSAL-CONTEXT modes (context C_t = (Y, X_<t)): best mode for tile t
    achieves E[L_m(t)] = H(X_t | X_<t, Y). Oracle = sum_t H(X_t | X_<t, Y),
    which by the chain rule equals H(X_1..X_T | Y) exactly.
    """
    s = 0.0
    for i, t in enumerate(tile_axes):
        cond_axes = list(tile_axes[:i]) + [Y_axis]   # X_<t and Y
        s += cond_entropy(joint, [t], cond_axes, shape)
    return s


def make_joint(probs_by_tuple):
    s = sum(probs_by_tuple.values())
    return {k: v / s for k, v in probs_by_tuple.items()}


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    return cond


def main():
    ok = True
    shape = (2, 2, 2)
    th = {0: 0.2, 1: 0.7}
    pY = {0: 0.5, 1: 0.5}
    print("Theorem 6.2b: Type-III-A absolute-rate converse")
    print("=" * 64)

    # ----------------------------------------------------------------
    # Source A: tiles conditionally INDEPENDENT given Y.
    #   given Y=y, X_t ~ Bern(theta_y) i.i.d.  => TC(.|Y)=0, both oracles == floor.
    # ----------------------------------------------------------------
    print("\nSource A (conditionally independent tiles given Y):")
    jointA = {}
    for y in (0, 1):
        for x1 in (0, 1):
            for x2 in (0, 1):
                p1 = th[y] if x1 == 1 else 1 - th[y]
                p2 = th[y] if x2 == 1 else 1 - th[y]
                jointA[(x1, x2, y)] = pY[y] * p1 * p2
    jointA = make_joint(jointA)
    tcA, sum_margA, jcA = total_correlation_given_Y(jointA, (0, 1), 2, shape)
    HxyA = cond_entropy(jointA, [0, 1], [2], shape)
    bi_A = oracle_blockindep(jointA, (0, 1), 2, shape)
    ca_A = oracle_causal(jointA, (0, 1), 2, shape)
    print(f"    H(X|Y)                  = {HxyA:.6f}")
    print(f"    TC(X1,X2|Y)             = {tcA:.6e}")
    print(f"    oracle block-independent= {bi_A:.6f}")
    print(f"    oracle causal-context   = {ca_A:.6f}")
    ok &= check("A: TC == 0 (cond. independent)", abs(tcA) < TOL)
    ok &= check("A: block-indep oracle == H(X|Y)", abs(bi_A - HxyA) < TOL)
    ok &= check("A: causal oracle == H(X|Y)", abs(ca_A - HxyA) < TOL)

    # ----------------------------------------------------------------
    # Source B: tiles DEPENDENT given Y (X1 == X2 given Y).
    #   given Y=y, Z ~ Bern(theta_y), X1 = X2 = Z.
    #   block-independent: oracle = 2 H(Z|Y), gap = TC = H(Z|Y) > 0.
    #   causal-context:    oracle = H(X1|Y) + H(X2|X1,Y) = H(Z|Y) + 0 = H(X|Y).
    # ----------------------------------------------------------------
    print("\nSource B (perfectly dependent tiles given Y, X1=X2):")
    jointB = {}
    for y in (0, 1):
        for z in (0, 1):
            p = th[y] if z == 1 else 1 - th[y]
            jointB[(z, z, y)] = pY[y] * p
    jointB = make_joint(jointB)
    tcB, sum_margB, jcB = total_correlation_given_Y(jointB, (0, 1), 2, shape)
    HxyB = cond_entropy(jointB, [0, 1], [2], shape)
    bi_B = oracle_blockindep(jointB, (0, 1), 2, shape)
    ca_B = oracle_causal(jointB, (0, 1), 2, shape)
    HzgivenY = sum(pY[y] * H([th[y], 1 - th[y]]) for y in (0, 1))
    print(f"    H(X|Y)                  = {HxyB:.6f}   (= H(Z|Y) = {HzgivenY:.6f})")
    print(f"    TC(X1,X2|Y)             = {tcB:.6f}")
    print(f"    oracle block-independent= {bi_B:.6f}   gap = {bi_B - HxyB:.6f}")
    print(f"    oracle causal-context   = {ca_B:.6f}   gap = {ca_B - HxyB:.6f}")
    # block-independent: pays the full inter-tile total correlation
    ok &= check("B: TC > 0 (dependent tiles)", tcB > 1e-3)
    ok &= check("B: TC == H(Z|Y)", abs(tcB - HzgivenY) < TOL)
    ok &= check("B: block-indep oracle == 2*H(Z|Y) (pays twice)", abs(bi_B - 2 * HzgivenY) < TOL)
    ok &= check("B: block-indep gap == TC (= context deficit at C_t=Y)",
                abs((bi_B - HxyB) - tcB) < TOL)
    ok &= check("B: block-indep strictly above floor", bi_B - HxyB > 1e-3)
    # causal-context: closes the gap EXACTLY (the corrected finding)
    ok &= check("B: causal oracle == H(X|Y) (gap closed by causal conditioning)",
                abs(ca_B - HxyB) < TOL)
    ok &= check("B: causal context deficit == 0", abs(ca_B - HxyB) < TOL)

    # ----------------------------------------------------------------
    # (6.2b.2) chain + (c) chain-rule identity sum_t H(X_t|X_<t,Y) == H(X|Y).
    # ----------------------------------------------------------------
    print("\nChain (6.2b.2) and chain-rule identity (6.2b.2c):")
    for nm, j, bi, ca, hxy, smarg in (
        ("A", jointA, bi_A, ca_A, HxyA, sum_margA),
        ("B", jointB, bi_B, ca_B, HxyB, sum_margB),
    ):
        ok &= check(f"{nm}: block-indep oracle >= sum_t H(X_t|Y)", bi >= smarg - TOL)
        ok &= check(f"{nm}: sum_t H(X_t|Y) >= H(X|Y)", smarg >= hxy - TOL)
        ok &= check(f"{nm}: causal oracle >= H(X|Y) (>= floor)", ca >= hxy - TOL)
        ok &= check(f"{nm}: chain rule sum_t H(X_t|X_<t,Y) == H(X|Y)", abs(ca - hxy) < TOL)

    # ----------------------------------------------------------------
    # Mode-inexpressiveness gap (first inequality strict): a mode coding X_t
    # under the WRONG law q pays cross-entropy = H + KL > H(X_t|C_t).
    # ----------------------------------------------------------------
    print("\nMode-inexpressiveness gap (first inequality of 6.2b.2 strict):")
    q = 0.5
    ce = 0.0
    h_true = 0.0
    for y in (0, 1):
        for x1 in (0, 1):
            pxy = sum(jointA.get((x1, x2, y), 0.0) for x2 in (0, 1))
            if pxy <= 0:
                continue
            pcond = pxy / pY[y]
            ce += pY[y] * pcond * (-math.log2(q if x1 == 1 else 1 - q))
            h_true += pY[y] * pcond * (-math.log2(pcond))
    print(f"    H(X_1|Y)              = {h_true:.6f}")
    print(f"    E[L | wrong mode q]   = {ce:.6f}")
    ok &= check("wrong-law mode pays strictly above H(X_1|Y) (gap = avg KL > 0)",
                ce > h_true + 1e-3)

    # ----------------------------------------------------------------
    # Source C (Remark 6.2b' regression test): the OPERATIVE benchmark is the
    # causal-adaptive floor (ii), NOT global-commit (iii). This is the example
    # that refuted the global-commit draft.
    #   T=2, Y empty. X1 ~ Bern(1/2), 1-bit code.
    #   X2 in {a,b,c,d}: P(.|X1=0)=(1/2,1/4,1/8,1/8), P(.|X1=1)=(1/8,1/8,1/2,1/4).
    #   Two prefix modes for tile 2: A lengths (1,2,3,3), B lengths (3,3,1,2).
    # ----------------------------------------------------------------
    print("\nSource C (Remark 6.2b' -- global-commit is NOT a floor):")
    p0 = [1 / 2, 1 / 4, 1 / 8, 1 / 8]
    p1 = [1 / 8, 1 / 8, 1 / 2, 1 / 4]
    A = [1, 2, 3, 3]
    B = [3, 3, 1, 2]
    half = 0.5
    HX = 1.0 + half * H(p0) + half * H(p1)          # H(X1) + H(X2|X1)
    EA0 = sum(l * p for l, p in zip(A, p0)); EA1 = sum(l * p for l, p in zip(A, p1))
    EB0 = sum(l * p for l, p in zip(B, p0)); EB1 = sum(l * p for l, p in zip(B, p1))
    # (i) hindsight: per-realization-of-X2 best mode (chosen after seeing X2)
    hindsight = 1.0 + half * sum(min(a, b) * p for a, b, p in zip(A, B, p0)) \
                    + half * sum(min(a, b) * p for a, b, p in zip(A, B, p1))
    # (ii) causal-adaptive: best mode per realization of C_2 = X1 (the operative selector)
    causal = 1.0 + half * min(EA0, EB0) + half * min(EA1, EB1)
    # (iii) global-commit: one mode fixed for the whole source
    glob = 1.0 + min(half * EA0 + half * EA1, half * EB0 + half * EB1)
    print(f"    H(X)                     = {HX:.4f}")
    print(f"    (i)   hindsight          = {hindsight:.4f}")
    print(f"    (ii)  causal-adaptive    = {causal:.4f}   (operative)")
    print(f"    (iii) global-commit      = {glob:.4f}   (NOT a floor)")
    ok &= check("C: ordering (i) <= (ii) <= (iii)", hindsight <= causal + TOL <= glob + TOL)
    ok &= check("C: causal-adaptive (ii) == H(X) (operative attains floor)", abs(causal - HX) < TOL)
    ok &= check("C: global-commit (iii) STRICTLY above H(X) (over-estimate, not a floor)",
                glob > HX + 1e-3)
    ok &= check("C: hindsight (i) <= H(X) (sub-floor, non-UD without flags)", hindsight <= HX + TOL)

    print("\n" + "=" * 64)
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
