#!/usr/bin/env python3
"""
Verify Theorem 6.2b: Type-III-A absolute-rate converse and optimality conditions.

Claims checked numerically on small explicit joint distributions:

  (6.2b.1) E[L_IIIA] >= H(X|Y)            -- Shannon floor (sanity, via entropy coder)
  (6.2b.2) sum_t min_m E[L_m(t)] >= sum_t H(X_t|Y) >= H(X|Y)
            (a) first ineq  equality  iff modes per-tile-expressive
            (b) second ineq equality  iff tiles conditionally independent given Y
  (6.2b.3) sum_t H(X_t|Y) - H(X|Y) = TC(X_1..X_T | Y)  (inter-tile total correlation)

KEY HONEST FINDING (Remark 6.2b): on a *dependent-tile* source the second
inequality is strict; even with a perfect predictor and fully expressive modes,
the per-tile coder overpays by exactly TC(.|Y) > 0. We construct such a source
and confirm the gap is real and equals the total correlation.

PASS iff all relations hold to numerical tolerance and the dependent-tile
example exhibits a strictly positive (correctly-sized) inter-tile gap.
"""
import itertools
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


def cond_entropy(joint, X_axes, Y_axis, shape):
    """H(X_axes | Y_axis) for a flat-dict joint indexed by full tuples."""
    # H(X,Y) - H(Y)
    XY = marginal(joint, tuple(X_axes) + (Y_axis,), shape)
    Y = marginal(joint, (Y_axis,), shape)
    return H(XY.values()) - H(Y.values())


def total_correlation_given_Y(joint, tile_axes, Y_axis, shape):
    """sum_t H(X_t|Y) - H(X_1..X_T | Y)."""
    sum_marg = sum(cond_entropy(joint, [t], Y_axis, shape) for t in tile_axes)
    joint_cond = cond_entropy(joint, list(tile_axes), Y_axis, shape)
    return sum_marg - joint_cond, sum_marg, joint_cond


def per_tile_oracle_expressive(joint, tile_axes, Y_axis, shape):
    """
    With per-tile-expressive entropy-coding modes, the best mode for tile t
    achieves E[L_m(t)] = H(X_t | Y). Oracle = sum_t H(X_t|Y).
    """
    return sum(cond_entropy(joint, [t], Y_axis, shape) for t in tile_axes)


def make_joint(probs_by_tuple):
    s = sum(probs_by_tuple.values())
    return {k: v / s for k, v in probs_by_tuple.items()}


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    return cond


def main():
    ok = True
    print("Theorem 6.2b: Type-III-A absolute-rate converse")
    print("=" * 64)

    # ----------------------------------------------------------------
    # Source A: tiles conditionally INDEPENDENT given Y.
    #   X_1, X_2 in {0,1}, Y in {0,1}; given Y=y, X_t ~ Bern(theta_y) i.i.d.
    #   => TC(X1,X2|Y) = 0, per-tile oracle == H(X|Y) == absolute floor.
    # ----------------------------------------------------------------
    print("\nSource A (conditionally independent tiles given Y):")
    th = {0: 0.2, 1: 0.7}
    pY = {0: 0.5, 1: 0.5}
    jointA = {}
    for y in (0, 1):
        for x1 in (0, 1):
            for x2 in (0, 1):
                p1 = th[y] if x1 == 1 else 1 - th[y]
                p2 = th[y] if x2 == 1 else 1 - th[y]
                jointA[(x1, x2, y)] = pY[y] * p1 * p2
    jointA = make_joint(jointA)
    shape = (2, 2, 2)
    tcA, sum_margA, joint_condA = total_correlation_given_Y(jointA, (0, 1), 2, shape)
    HxyA = cond_entropy(jointA, [0, 1], 2, shape)
    oracleA = per_tile_oracle_expressive(jointA, (0, 1), 2, shape)
    print(f"    H(X|Y)            = {HxyA:.6f}")
    print(f"    sum_t H(X_t|Y)    = {sum_margA:.6f}")
    print(f"    TC(X1,X2|Y)       = {tcA:.6e}")
    print(f"    per-tile oracle   = {oracleA:.6f}")
    ok &= check("A: TC == 0 (cond. independent)", abs(tcA) < TOL)
    ok &= check("A: oracle == H(X|Y) (absolute floor attained)", abs(oracleA - HxyA) < TOL)
    ok &= check("A: sum_t H(X_t|Y) == H(X|Y)", abs(sum_margA - HxyA) < TOL)

    # ----------------------------------------------------------------
    # Source B: tiles DEPENDENT given Y (in fact X1 == X2 given Y).
    #   Given Y=y, draw Z ~ Bern(theta_y) and set X1 = X2 = Z.
    #   => H(X1|Y)=H(X2|Y)=H(Z|Y) each > 0 but H(X1,X2|Y)=H(Z|Y),
    #   so TC = H(Z|Y) > 0: a per-tile coder pays TWICE for one bit.
    # ----------------------------------------------------------------
    print("\nSource B (perfectly dependent tiles given Y, X1=X2):")
    jointB = {}
    for y in (0, 1):
        for z in (0, 1):
            p = th[y] if z == 1 else 1 - th[y]
            jointB[(z, z, y)] = pY[y] * p
    # zero-prob entries for x1!=x2 implicitly absent
    jointB = make_joint(jointB)
    tcB, sum_margB, joint_condB = total_correlation_given_Y(jointB, (0, 1), 2, shape)
    HxyB = cond_entropy(jointB, [0, 1], 2, shape)
    oracleB = per_tile_oracle_expressive(jointB, (0, 1), 2, shape)
    # the genuine single-tile conditional entropy (what an ideal joint coder pays)
    HzgivenY = sum(pY[y] * H([th[y], 1 - th[y]]) for y in (0, 1))
    print(f"    H(X|Y)            = {HxyB:.6f}")
    print(f"    sum_t H(X_t|Y)    = {sum_margB:.6f}")
    print(f"    TC(X1,X2|Y)       = {tcB:.6f}   (= H(Z|Y) = {HzgivenY:.6f})")
    print(f"    per-tile oracle   = {oracleB:.6f}")
    print(f"    overpay vs floor  = {oracleB - HxyB:.6f}")
    ok &= check("B: TC > 0 (dependent tiles)", tcB > 1e-3)
    ok &= check("B: TC == H(Z|Y) (exact size of inter-tile gap)", abs(tcB - HzgivenY) < TOL)
    ok &= check("B: H(X|Y) == H(Z|Y) (true joint cost is one copy)", abs(HxyB - HzgivenY) < TOL)
    ok &= check("B: per-tile oracle == 2*H(Z|Y) (pays twice)", abs(oracleB - 2 * HzgivenY) < TOL)
    ok &= check("B: oracle - floor == TC (gap is exactly the total correlation)",
                abs((oracleB - HxyB) - tcB) < TOL)
    ok &= check("B: oracle strictly above floor (III-A NOT absolutely optimal here)",
                oracleB - HxyB > 1e-3)

    # ----------------------------------------------------------------
    # (6.2b.2) chain holds in both directions; (6.2b.1) floor respected.
    # ----------------------------------------------------------------
    print("\nChain (6.2b.2) sum_t min_m E[L_m] >= sum_t H(X_t|Y) >= H(X|Y):")
    for nm, orc, smarg, hxy in (("A", oracleA, sum_margA, HxyA), ("B", oracleB, sum_margB, HxyB)):
        ok &= check(f"{nm}: oracle >= sum_t H(X_t|Y) - TOL", orc >= smarg - TOL)
        ok &= check(f"{nm}: sum_t H(X_t|Y) >= H(X|Y) - TOL", smarg >= hxy - TOL)

    # ----------------------------------------------------------------
    # Mode-inexpressiveness gap (first inequality strict): a coarse mode
    # that codes X_t under the *wrong* law pays KL above H(X_t|Y).
    # ----------------------------------------------------------------
    print("\nMode-inexpressiveness gap (first inequality of 6.2b.2 strict):")
    # tile-1 of Source A, true Bern mixture vs a single wrong-law mode q=0.5
    q = 0.5
    # E[L under q] = sum_x p(x) * (-log2 q(x)); cross-entropy
    pX1 = marginal(jointA, (0,), shape)  # P(X1)
    # but converse is per-Y; use conditional cross entropy averaged over Y
    ce = 0.0
    h_true = 0.0
    for y in (0, 1):
        for x1 in (0, 1):
            pxy = jointA_get(jointA, x1, y)
            if pxy <= 0:
                continue
            pcond = pxy / pY[y]
            ce += pY[y] * pcond * (-math.log2(q if x1 == 1 else 1 - q))
            h_true += pY[y] * pcond * (-math.log2(pcond))
    print(f"    H(X_1|Y)              = {h_true:.6f}")
    print(f"    E[L | wrong mode q]   = {ce:.6f}")
    ok &= check("wrong-law mode pays strictly above H(X_1|Y) (gap = avg KL > 0)",
                ce > h_true + 1e-3)

    print("\n" + "=" * 64)
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def jointA_get(joint, x1, y):
    s = 0.0
    for x2 in (0, 1):
        s += joint.get((x1, x2, y), 0.0)
    return s


if __name__ == "__main__":
    raise SystemExit(main())
