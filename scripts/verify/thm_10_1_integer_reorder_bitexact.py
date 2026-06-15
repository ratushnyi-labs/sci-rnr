#!/usr/bin/env python3
r"""
Verification of the MATH CORE of Theorem 10.1 (BUG-007-B):
"Bit-exact reproduction under deterministic integer inference".

WHAT THEOREM 10.1 ACTUALLY CLAIMS (section 10.5 of the main paper).
  Let M be an inference graph with L layers, layer l an integer
  matrix-vector product of fan-in d_l followed by a non-linearity. If
    (i)   all weights and inputs are integers;
    (ii)  every layer-l accumulator has signed bit-width b_l satisfying
              b_l >= ceil(log2(d_l * W_l * I_{l-1} + 1)) + 1,
          where W_l bounds |weight entry| and I_{l-1} bounds |input|
          (the +1 inside the log accounts for the asymmetric signed
          two's-complement range [-2^{b-1}, 2^{b-1}-1], so the positive
          extremum +d_l*W_l*I_{l-1} is representable);
    (iii) all non-linearities are deterministic integer->integer maps
          (e.g. lookup tables);
    (iv)  all reduction sums follow a fixed deterministic order pinned
          by the model description;
  THEN no accumulator overflows and M(C) is BIT-IDENTICAL on any two
  conforming implementations, for every input C.

It is a SUFFICIENT-CONDITION / specification theorem. It does NOT claim
floating point is reproducible -- on the contrary, float non-associativity
is the thing the integer model avoids. The end-to-end cross-platform
demonstration on real hardware (a reference codec + the engineering-spec
section 12 conformance suite across targets) is deferred to section 13.3 of
the main paper and is NOT what this probe establishes.

WHAT THIS PROBE ESTABLISHES (the math core only; no real LM, no real HW):
  (1) PER-LAYER NO-OVERFLOW BIT-WIDTH BOUND. For integer inputs/weights in
      the stated ranges and fan-in K, the accumulator width
          b = ceil(log2(K * W * I + 1)) + 1
      keeps |sum| < 2^{b-1} ALWAYS -- verified against the saturating
      worst case (every addend at the extreme, same sign) and against
      random inputs. We confirm the worst-case sum is representable and
      that b-1 (one bit narrower) is NOT enough at the worst case.
  (2) REDUCTION-ORDER INVARIANCE. Under that width, the integer sum over
      MANY random permutations of the addends, and tree- vs sequential-
      reduction, is bit-identical -- because integer addition (within
      range) is exactly associative and commutative. We then chain a few
      integer "layers" (integer matmul + a deterministic integer
      non-linearity / requantization) and show the full forward pass is
      bit-identical across independently reordered reductions.
  (3) NEGATIVE CONTROL. An UNDERSIZED accumulator (width < the bound)
      breaks the guarantee, in the two overflow conventions that real
      integer ALUs use (Theorem 10.6 condition (c) pins exactly this
      "saturation-vs-wrap behaviour"):
        - WRAP (mod 2^b): addition mod 2^b is an abelian group, so the
          result stays reorder-INVARIANT but is WRONG (!= the exact value
          a conforming-width accumulator would give) -> the bound is
          necessary for CORRECTNESS.
        - SATURATE (clamp to range): NON-associative, so the result becomes
          REORDER-DEPENDENT -> the integer analogue of IEEE-754 float
          non-associativity, the failure hypothesis (ii) precludes by
          keeping every partial sum in range.
      We exhibit both, including a concrete hand-checkable saturating
      instance and a full chained forward pass that diverges across
      reorderings under an undersized saturating accumulator.

Pure Python integers are unbounded, so "the accumulator" under hypothesis
(ii) is modeled as Python int (exact). The negative control models a
*hardware* accumulator of insufficient width by reducing into width b after
every partial add, under each of the two conventions above -- exactly what
an overflowing fixed-width integer ALU does.

Emits PASS/FAIL per check and "OVERALL -> PASS" on success (CI contract).
"""

import math
import random
import sys

SEED = 20260615
RNG = random.Random(SEED)


def line(ok, label, detail):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {label}: {detail}")
    return bool(ok)


# --------------------------------------------------------------------------
# Theorem 10.1 hypothesis (ii): the per-layer no-overflow accumulator width.
# --------------------------------------------------------------------------
def accumulator_bits(K, W, I):
    """Theorem 10.1 (ii): b >= ceil(log2(K*W*I + 1)) + 1.

    K = fan-in (d_l), W = bound on |weight|, I = bound on |input|.
    The worst-case |partial sum| is K*W*I; +1 makes the positive extremum
    representable in the asymmetric signed range [-2^{b-1}, 2^{b-1}-1]."""
    worst = K * W * I
    return math.ceil(math.log2(worst + 1)) + 1


def representable_signed(value, b):
    """True iff value fits in signed two's-complement of width b:
    -2^{b-1} <= value <= 2^{b-1} - 1."""
    return -(1 << (b - 1)) <= value <= (1 << (b - 1)) - 1


def wrap_signed(value, b):
    """Two's-complement wraparound into width b (models a fixed-width ALU)."""
    m = 1 << b
    v = value % m
    if v >= (1 << (b - 1)):
        v -= m
    return v


# --------------------------------------------------------------------------
# Integer reductions: exact (unbounded) vs. fixed-width-wrapping.
# --------------------------------------------------------------------------
def _shuffle(seq):
    """Return a randomly reordered copy (does not mutate the argument)."""
    out = list(seq)
    RNG.shuffle(out)
    return out


def reduce_sequential(addends):
    """Exact left-to-right integer sum (Python int = unbounded accumulator)."""
    acc = 0
    for a in addends:
        acc += a
    return acc


def reduce_tree(addends):
    """Exact pairwise/tree integer reduction."""
    vals = list(addends)
    if not vals:
        return 0
    while len(vals) > 1:
        nxt = []
        for i in range(0, len(vals) - 1, 2):
            nxt.append(vals[i] + vals[i + 1])
        if len(vals) % 2 == 1:
            nxt.append(vals[-1])
        vals = nxt
    return vals[0]


def reduce_sequential_wrapping(addends, b):
    """Left-to-right sum in a fixed-width-b accumulator (wraps each add).

    NOTE: addition modulo 2^b is itself an abelian-group operation, so a
    pure-wrap accumulator stays reorder-INVARIANT -- it just returns the
    WRONG (wrapped) value once the bound is violated. Wrap therefore
    demonstrates the bound is necessary for CORRECTNESS, not for order-
    independence. The order-DEPENDENT failure mode is saturation, below."""
    acc = 0
    for a in addends:
        acc = wrap_signed(acc + a, b)
    return acc


def sat_signed(value, b):
    """Saturating clamp into signed width b: clamp to [-2^{b-1}, 2^{b-1}-1]."""
    lo, hi = -(1 << (b - 1)), (1 << (b - 1)) - 1
    if value < lo:
        return lo
    if value > hi:
        return hi
    return value


def reduce_sequential_saturating(addends, b):
    """Left-to-right sum in a fixed-width-b SATURATING accumulator.

    Saturation is NON-associative: clamp(clamp(a+b)+c) != clamp(clamp(a+c)+b)
    in general once the running partial sum exceeds the range. This is the
    integer analogue of IEEE-754 float non-associativity, and it is exactly
    what Theorem 10.1's hypothesis (ii) precludes by keeping every partial
    sum within range. Theorem 10.6 condition (c) pins 'saturation-vs-wrap
    behaviour at each arithmetic operation' for the same reason."""
    acc = 0
    for a in addends:
        acc = sat_signed(acc + a, b)
    return acc


# --------------------------------------------------------------------------
# (1) Per-layer no-overflow bit-width bound.
# --------------------------------------------------------------------------
def check_bound():
    print("\n(1) Per-layer no-overflow bit-width bound "
          "b = ceil(log2(K*W*I + 1)) + 1:")
    ok = True

    # Stated numerical sanity checks from section 10.5.
    #   int8:  I=W=2^7,  d=2^12  -> worst sum 2^26, b=32 suffices comfortably.
    #   int16: I=W=2^15, d=2^12  -> worst sum 2^42, need b>=44 (b=32 fails).
    b_int8 = accumulator_bits(2**12, 2**7, 2**7)
    ok &= line(
        2**12 * 2**7 * 2**7 == 2**26 and b_int8 <= 32,
        "int8 example (d=4096, W=I=128): worst sum 2^26",
        f"worst=2^26, bound b={b_int8} <= 32 (paper: b=32 suffices)",
    )
    b_int16 = accumulator_bits(2**12, 2**15, 2**15)
    ok &= line(
        2**12 * 2**15 * 2**15 == 2**42 and b_int16 >= 44 and b_int16 > 32,
        "int16 example (d=4096, W=I=2^15): worst sum 2^42",
        f"worst=2^42, bound b={b_int16} >= 44 > 32 (paper: b=32 overflows)",
    )
    # Degenerate scalar int16 case d=1 fits b=32 (paper's stated edge case);
    # any fan-in d>=2 under int16 bounds exceeds 2^31 and needs b>=33.
    b_scalar16 = accumulator_bits(1, 2**15, 2**15)
    b_fanin2_16 = accumulator_bits(2, 2**15, 2**15)
    ok &= line(
        b_scalar16 <= 32 and b_fanin2_16 >= 33,
        "int16 edge: scalar d=1 fits b=32; d=2 needs b>=33",
        f"b(d=1)={b_scalar16} <= 32; b(d=2)={b_fanin2_16} >= 33",
    )

    # Worst-case saturating inputs: all addends at the extreme, same sign.
    # The summed magnitude must be representable at b and NOT at b-1.
    print("  worst-case saturating + random inputs across (K,W,I) cases:")
    cases = [
        (1, 1, 1), (2, 1, 1), (4096, 127, 127), (4096, 2**15, 2**15),
        (16, 7, 31), (1000, 255, 255), (3, 2, 2), (777, 1023, 4095),
    ]
    for (K, W, I) in cases:
        b = accumulator_bits(K, W, I)
        worst = K * W * I  # all addends +W*I (or all -W*I)

        # bound holds: +worst and -worst both representable at b
        rep_pos = representable_signed(worst, b)
        rep_neg = representable_signed(-worst, b)
        # bound is tight-ish: b-1 cannot represent the positive worst case
        # (the asymmetric range means -worst may still fit at b-1, but the
        #  positive extremum is what the +1-inside-log guards, so we test +worst)
        narrow_fails = not representable_signed(worst, b - 1)

        ok &= line(
            rep_pos and rep_neg and narrow_fails,
            f"K={K},W={W},I={I}: b={b}",
            f"+/-worst={worst} representable at b={b}, "
            f"NOT at b-1={b - 1} => bound necessary & sufficient",
        )

        # random inputs at this (K,W,I): exact sum never overflows width b
        bad = 0
        for _ in range(200):
            addends = [RNG.randint(-W, W) * RNG.randint(-I, I) for _ in range(K)]
            # clamp to the per-term magnitude bound W*I (product of bounds)
            s = reduce_sequential(addends)
            if not representable_signed(s, b):
                bad += 1
        ok &= line(
            bad == 0,
            f"K={K},W={W},I={I}: 200 random sums fit width b={b}",
            f"{200 - bad}/200 within signed range [-2^{b-1}, 2^{b-1}-1]",
        )
    return ok


# --------------------------------------------------------------------------
# (2) Reduction-order invariance under the no-overflow width.
# --------------------------------------------------------------------------
def check_reorder_invariance():
    print("\n(2) Reduction-order invariance under the no-overflow width:")
    ok = True

    print("  single reduction: sequential == tree == many permutations:")
    for trial in range(40):
        K = RNG.randint(2, 4096)
        W = RNG.choice([1, 7, 127, 255, 2**15])
        I = RNG.choice([1, 7, 127, 255, 2**15])
        b = accumulator_bits(K, W, I)
        addends = [RNG.randint(-W, W) * RNG.randint(-I, I) for _ in range(K)]

        ref = reduce_sequential(addends)
        # tree reduction
        tree = reduce_tree(addends)
        # many random permutations
        perm_ok = True
        for _ in range(25):
            shuffled = addends[:]
            RNG.shuffle(shuffled)
            if reduce_sequential(shuffled) != ref:
                perm_ok = False
                break
        # all within width b (sanity: we are in the no-overflow regime)
        within = representable_signed(ref, b)
        good = (tree == ref) and perm_ok and within
        if not good or trial < 3:
            ok &= line(
                good,
                f"trial {trial}: K={K},W={W},I={I},b={b}",
                f"seq=tree={tree == ref}, 25 perms identical={perm_ok}, "
                f"sum={ref} fits b={b}={within}",
            )
        else:
            ok &= good
    print(f"    ... 40 single-reduction trials checked "
          f"(seq==tree==perms, all within width)")

    # Chain a few integer layers: matmul + deterministic integer non-linearity.
    print("  chained integer forward pass: bit-identical across reorderings:")

    def int_relu_requantize(vec, shift):
        """Deterministic integer->integer non-linearity (hypothesis iii):
        clamp at 0 (ReLU) then arithmetic-right-shift requantize. Pure
        integer, no rounding ambiguity, identical on any implementation."""
        return [max(0, v) >> shift for v in vec]

    def forward_pass(x, weights, shift, reducer):
        """weights: list of layers; layer = matrix (rows = out, cols = in).
        reducer: how the per-output-neuron dot product is summed."""
        cur = x
        for W_mat in weights:
            out = []
            for row in W_mat:
                products = [w * a for w, a in zip(row, cur)]
                out.append(reducer(products))
            cur = int_relu_requantize(out, shift)
        return cur

    pass_ok = True
    for trial in range(30):
        d_in = RNG.randint(2, 64)
        depth = RNG.randint(2, 4)
        Wb, Ib = 127, 127  # int8-style bounds
        shift = 7          # requantize so magnitudes stay bounded
        # build layers; keep widths modest so the exact-int forward pass is cheap
        weights = []
        d = d_in
        for _ in range(depth):
            d_out = RNG.randint(2, 64)
            mat = [[RNG.randint(-Wb, Wb) for _ in range(d)] for _ in range(d_out)]
            weights.append(mat)
            d = d_out
        x = [RNG.randint(0, Ib) for _ in range(d_in)]

        # reference: strict sequential reduction (the pinned order, hyp. iv)
        ref = forward_pass(x, weights, shift, reduce_sequential)

        # independently reordered reductions must agree bit-for-bit
        def permuted_reducer(products):
            p = products[:]
            RNG.shuffle(p)
            return reduce_sequential(p)

        out_tree = forward_pass(x, weights, shift, reduce_tree)
        out_perm = forward_pass(x, weights, shift, permuted_reducer)
        good = (out_tree == ref) and (out_perm == ref)
        pass_ok &= good
        if not good or trial < 2:
            ok &= line(
                good,
                f"forward-pass trial {trial}: d_in={d_in}, depth={depth}",
                f"tree==seq:{out_tree == ref}, perm==seq:{out_perm == ref}, "
                f"out_len={len(ref)}",
            )
    ok &= pass_ok
    print(f"    ... 30 chained-forward-pass trials checked "
          f"(tree & permuted reductions bit-identical to pinned order)")
    return ok


# --------------------------------------------------------------------------
# (3) NEGATIVE CONTROL: an undersized accumulator breaks the guarantee.
#
# Two distinct hardware overflow conventions (Theorem 10.6 condition (c)
# pins exactly this 'saturation-vs-wrap behaviour'):
#   - WRAP (mod 2^b): addition mod 2^b is an abelian group, so the result
#     stays reorder-INVARIANT but is WRONG (!= exact) once the bound is
#     violated -> the bound is necessary for CORRECTNESS.
#   - SATURATE (clamp to range): NON-associative, so the result becomes
#     reorder-DEPENDENT -> the integer analogue of float non-associativity,
#     the failure Theorem 10.1 hypothesis (ii) precludes.
# A conforming engine must therefore EITHER size the accumulator per the
# bound (no overflow ever triggers, both conventions agree with exact and
# with each other) OR pin the convention; the bound makes the convention
# irrelevant, which is the cleaner guarantee.
# --------------------------------------------------------------------------
def check_negative_control():
    print("\n(3) Negative control: UNDERSIZED accumulator breaks the "
          "guarantee (bound is necessary):")
    ok = True

    # (3a) WRAP: reorder-invariant but WRONG value (bound necessary for
    #      correctness; wrap alone does NOT exhibit reorder divergence).
    print("  (3a) wrap mode: reorder-INVARIANT but WRONG once bound violated:")
    wrong_count = 0
    invariant_count = 0
    trials = 200
    used = 0
    for _ in range(trials):
        K = RNG.randint(8, 512)
        W = RNG.choice([127, 255, 2**15])
        I = RNG.choice([127, 255, 2**15])
        b = accumulator_bits(K, W, I)
        b_bad = max(4, b - RNG.randint(2, 8))
        addends = [RNG.randint(-W, W) * RNG.randint(-I, I) for _ in range(K)]
        exact = reduce_sequential(addends)
        if representable_signed(exact, b_bad):
            continue  # no overflow in this trial; not adversarial
        used += 1
        # wrap result is the same for every ordering (mod-2^b is a group) ...
        results = {reduce_sequential_wrapping(_shuffle(addends), b_bad)
                   for _ in range(30)}
        if len(results) == 1:
            invariant_count += 1
        # ... but it is WRONG: differs from the exact (conforming-width) value
        if next(iter(results)) != exact:
            wrong_count += 1
    ok &= line(
        used > 0 and invariant_count == used,
        "wrap mode stays reorder-invariant (mod 2^b is an abelian group)",
        f"{invariant_count}/{used} adversarial trials: one wrapped value "
        f"across 30 orderings",
    )
    ok &= line(
        used > 0 and wrong_count == used,
        "wrap mode is WRONG vs exact => bound necessary for correctness",
        f"{wrong_count}/{used} adversarial trials: wrapped != exact sum",
    )

    # (3b) SATURATE: reorder-DEPENDENT (non-associative) once bound violated.
    print("  (3b) saturate mode: reorder-DEPENDENT once bound violated:")
    found_divergence = 0
    used_sat = 0
    for _ in range(trials):
        K = RNG.randint(8, 512)
        W = RNG.choice([127, 255, 2**15])
        I = RNG.choice([127, 255, 2**15])
        b = accumulator_bits(K, W, I)
        b_bad = max(4, b - RNG.randint(2, 8))
        addends = [RNG.randint(-W, W) * RNG.randint(-I, I) for _ in range(K)]
        exact = reduce_sequential(addends)
        if representable_signed(exact, b_bad):
            continue
        used_sat += 1
        results = {reduce_sequential_saturating(_shuffle(addends), b_bad)
                   for _ in range(30)}
        if len(results) > 1:
            found_divergence += 1
    ok &= line(
        used_sat > 0 and found_divergence > 0,
        "saturating undersized reduction is reorder-DEPENDENT",
        f"{found_divergence} of {used_sat} adversarial trials produced "
        f">1 distinct saturated value across orderings",
    )

    # Concrete, hand-checkable saturating instance (signed width 5: [-16, 15]).
    print("  concrete saturating instance (b_bad = 5, range [-16, 15]):")
    b_bad = 5
    # exact sum = 13, but partial sums overshoot +15 and clamp differently by
    # order: [15, 14, -16] vs [15, -16, 14] vs [-16, 15, 14].
    base = [15, 14, -16]
    exact = sum(base)  # 13, fits the range -- so the *correct* answer is 13
    s1 = reduce_sequential_saturating([15, 14, -16], b_bad)  # 15 ->sat(29)=15 ->sat(-1)=-1
    s2 = reduce_sequential_saturating([15, -16, 14], b_bad)  # 15 ->sat(-1)=-1 ->sat(13)=13
    s3 = reduce_sequential_saturating([-16, 15, 14], b_bad)  # -16 ->sat(-1)=-1 ->sat(13)=13
    distinct = len({s1, s2, s3})
    ok &= line(
        exact == 13 and distinct >= 2 and 13 in {s1, s2, s3} and -1 in {s1, s2, s3},
        "concrete saturating reorder divergence",
        f"exact={exact} (in range); saturated paths "
        f"{{[15,14,-16]->{s1}, [15,-16,14]->{s2}, [-16,15,14]->{s3}}} "
        f"=> {distinct} distinct values (only some equal the exact 13)",
    )

    # Full forward pass with an undersized SATURATING accumulator: diverges.
    print("  chained forward pass, undersized saturating accumulator: diverges:")

    def int_relu_requantize(vec, shift):
        return [max(0, v) >> shift for v in vec]

    def forward_pass_saturating(x, weights, shift, b_acc, reducer_order):
        cur = x
        for W_mat in weights:
            out = []
            for row in W_mat:
                products = [w * a for w, a in zip(row, cur)]
                out.append(reduce_sequential_saturating(reducer_order(products), b_acc))
            cur = int_relu_requantize(out, shift)
        return cur

    diverged = 0
    fp_trials = 60
    for _ in range(fp_trials):
        d_in = RNG.randint(8, 48)
        depth = RNG.randint(2, 3)
        Wb = 127
        shift = 3  # under-shift so magnitudes stay large and saturation triggers
        weights = []
        d = d_in
        for _ in range(depth):
            d_out = RNG.randint(8, 48)
            mat = [[RNG.randint(-Wb, Wb) for _ in range(d)] for _ in range(d_out)]
            weights.append(mat)
            d = d_out
        x = [RNG.randint(0, 127) for _ in range(d_in)]
        # deliberately undersized for these int8 magnitudes & fan-in (~22-24 bits)
        b_bad = 12
        out_seq = forward_pass_saturating(x, weights, shift, b_bad, lambda p: p[:])
        out_perm = forward_pass_saturating(x, weights, shift, b_bad, _shuffle)
        if out_seq != out_perm:
            diverged += 1
    ok &= line(
        diverged > 0,
        "undersized saturating forward pass is reorder-dependent",
        f"{diverged} of {fp_trials} trials: permuted reduction != pinned order "
        f"(b_bad=12 << required width)",
    )
    return ok


# --------------------------------------------------------------------------
# Anti-strawman guard: confirm the POSITIVE result is real, not vacuous.
# --------------------------------------------------------------------------
def check_nonvacuous():
    print("\n(4) Anti-vacuity: the same data that diverges when undersized "
          "is bit-identical at the correct width:")
    ok = True
    matched = 0
    trials = 50
    for _ in range(trials):
        K = RNG.randint(8, 512)
        W = RNG.choice([127, 255])
        I = RNG.choice([127, 255])
        b = accumulator_bits(K, W, I)
        addends = [RNG.randint(-W, W) * RNG.randint(-I, I) for _ in range(K)]
        exact = reduce_sequential(addends)
        # correct width: exact (unbounded) == wrapping-at-b (no wrap occurs)
        wide = reduce_sequential_wrapping(addends, b)
        # and any permutation agrees
        shuffled = addends[:]
        RNG.shuffle(shuffled)
        wide_perm = reduce_sequential_wrapping(shuffled, b)
        if wide == exact == wide_perm:
            matched += 1
    ok &= line(
        matched == trials,
        "correct width: wrapping == exact == permuted (no divergence)",
        f"{matched}/{trials} trials bit-identical at the Theorem 10.1 bound",
    )
    return ok


def main() -> int:
    print("Verification: Theorem 10.1 math core -- integer reduction-order "
          "determinism (BUG-007-B)")
    print("=" * 74)
    print("MATH CORE ONLY. Demonstrates integer reorder-invariance under the")
    print("per-layer no-overflow bit-width bound; NOT a real LM, NOT real")
    print("hardware, NOT a claim that float is reproducible. The end-to-end")
    print("cross-platform demonstration is deferred to section 13.3 (external).")

    results = {
        "(1) per-layer no-overflow bit-width bound": check_bound(),
        "(2) reduction-order invariance (sum + forward pass)": check_reorder_invariance(),
        "(3) negative control: undersized accumulator (wrap WRONG / saturate DIVERGES)": check_negative_control(),
        "(4) anti-vacuity: correct width is bit-identical": check_nonvacuous(),
    }

    print("\n" + "=" * 74)
    print("Summary:")
    all_ok = True
    for name, ok in results.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        all_ok &= ok

    print()
    if all_ok:
        print("OVERALL -> PASS")
        return 0
    print("OVERALL -> FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
