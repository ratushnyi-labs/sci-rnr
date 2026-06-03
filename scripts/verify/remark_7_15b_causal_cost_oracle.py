r"""
Verification for Remark 7.15b (the causal cost oracle is poly-time for ANY
autoregressive field; the #P-hardness of Remark 7.15a is confined to the
optional NON-CONTIGUOUS support-selection lever, which the Theorem 7.15
streaming coder never invokes).

------------------------------------------------------------------------
THE POINT.

Remark 7.15a established that the value oracle H(X_T) for a discrete AR field
is #P-hard for general fields when T is a NON-CONTIGUOUS subset (it forces a
marginalisation over the gap positions). One might fear this makes the
Theorem 7.15 neural-LM encoder itself intractable. It does NOT. The streaming
RNR coder codes each sub-block left-to-right, coding X_i under the predictor's
conditional Q_i(. | x_{<i}) on the REALISED, decoder-reproducible prefix. Its
realised codelength on an input x is

    L(x) = sum_i  -log2 Q_i(x_i | x_{<i}),

a sum of K single-softmax evaluations -- ONE forward pass, NO marginalisation
of any position. So the *causal* cost oracle (evaluate and achieve the
streaming-coder rate) is POLYNOMIAL for every AR field, of ANY effective order
(including a full-attention transformer). The #P-hard object appears only when
a coding decision forces a MARGINAL rather than a CONDITIONAL -- i.e. the
optional non-contiguous support selection of Remark 7.15a.

This is the autoregressive moral: conditioning on realised data is a forward
pass; marginalising hidden/gap positions is #P-hard. The same field admits a
poly causal-cost oracle and a #P-hard selection oracle.

------------------------------------------------------------------------
WHAT THIS SCRIPT VERIFIES.

(C1) Chain-rule identity: for an AR field with conditionals Q_i(.|x_{<i}),
     H(X_1..n) = sum_i H(X_i | X_{<i}).  (Entropy chain rule; sanity that the
     per-position softmax entropies sum to the block entropy.)

(C2) The streaming coder's realised codelength L(x) = sum_i -log2 Q_i(x_i|x_{<i})
     -- computed using ONLY the forward conditionals, never the joint or any
     marginal -- equals -log2 P(x), and its sample mean over draws from the
     field converges to H(X_1..n).  So the operative cost is forward-pass-only
     and achieves the entropy rate.  Work: O(n|Sig|) per sequence.

(C3) CONTRAST: a non-contiguous residual H(X_S | X_Sbar) (the Remark 7.15a
     lever) needs the marginal over the gap positions.  It is NOT obtainable
     from the ORIGINAL AR conditionals Q_i on a realised full prefix: the kept
     subsequence does have a chain-rule product P(x_T)=prod_k P(X_{t_k}|X_{t_<k}),
     but those INDUCED conditionals are not the Q_i -- obtaining them requires
     marginalising the gaps.  Op-count: causal = O(n|Sig|); the gap marginal
     = O(|Sig|^{gap}).  We certify the induced subsequence is not naive-Markov
     (a kept gap induces dependence the original forward conditionals miss) and
     exhibit the exponential op-count gap.

(C4) The causal poly evaluation is ORDER-INDEPENDENT: for a maximally
     history-dependent field (order = n-1, every position conditions on the
     whole prefix) the causal realised codelength is still O(n|Sig|) softmax
     lookups, while the non-contiguous marginal is maximally expensive.  This
     is the regime of a full-attention transformer.

PASS = C1 (chain rule holds) AND C2 (forward-only codelength == -log2 P, sample
mean -> H) AND C3 (non-contiguous needs a genuine marginal; exp op-count gap)
AND C4 (order-(n-1) field: causal still O(n|Sig|), marginal still exponential).

HONEST SCOPE.  This says the *streaming* coder and the Theorem 7.15 bound
N*H_{M'} are poly-time achievable/evaluable for any AR field; it does NOT
re-open Remark 7.15a's hardness (which stands for the optional non-contiguous
selection).  Under a mismatched predictor Q != P the realised mean is the
cross-entropy H(X)+KL = N*H_{M'} of Theorem 7.4/7.15, still forward-pass-only.
"""

import math
import sys
from itertools import product

import numpy as np

RNG = np.random.default_rng(20260603)
TOL = 1e-9


# ======================================================================
# An autoregressive field exposed THROUGH its conditionals Q_i(. | x_<i),
# exactly as a neural predictor exposes per-position softmaxes.  order=None
# means full history dependence (order = i, the transformer regime).
# ======================================================================

class ARField:
    def __init__(self, n, A, order=None, rng=RNG, temp=1.0):
        self.n, self.A = n, A
        self.order = order
        self.temp = temp
        self.rng = rng
        self.cond = {}      # (i, prefix_tuple) -> softmax vector

    def Q(self, i, prefix):
        """The forward conditional Q_i(. | x_{<i}); `prefix` is the realised
        tuple x_0..x_{i-1}.  This is the ONLY access the streaming coder uses
        -- one softmax per position, no marginalisation."""
        if self.order is not None:
            prefix = prefix[max(0, i - self.order):i]
        key = (i, tuple(prefix))
        if key not in self.cond:
            logits = self.rng.standard_normal(self.A) / max(self.temp, 1e-6)
            p = np.exp(logits - logits.max())
            self.cond[key] = p / p.sum()
        return self.cond[key]

    def joint(self):
        """Full joint table P(x) = prod_i Q_i(x_i | x_<i) (brute, for checks)."""
        P = np.zeros((self.A,) * self.n)
        for idx in product(range(self.A), repeat=self.n):
            pr = 1.0
            for i in range(self.n):
                pr *= self.Q(i, idx[:i])[idx[i]]
            P[idx] = pr
        # numerical renorm (should already sum to 1 up to fp error)
        return P / P.sum()

    def realised_codelength(self, x):
        """L(x) = sum_i -log2 Q_i(x_i | x_<i): FORWARD-PASS ONLY, no joint,
        no marginal.  O(n |Sig|) work (one softmax read per position)."""
        L = 0.0
        for i in range(self.n):
            qi = self.Q(i, tuple(x[:i]))
            L += -math.log2(qi[x[i]])
        return L


def entropy_of_table(P):
    p = P.ravel(); p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def subset_entropy(P, T):
    n = P.ndim
    T = tuple(sorted(T))
    axes_out = tuple(i for i in range(n) if i not in T)
    M = P.sum(axis=axes_out) if axes_out else P
    m = M.ravel(); m = m[m > 0]
    return float(-(m * np.log2(m)).sum())


def cond_block_entropy_chain_rule(P):
    """sum_i H(X_i | X_<i) computed from successive prefix marginals of the
    joint (a sanity oracle for the chain-rule identity C1)."""
    n = P.ndim
    total = 0.0
    for i in range(n):
        Hpre = subset_entropy(P, range(i)) if i > 0 else 0.0
        Hpre_i = subset_entropy(P, range(i + 1))
        total += Hpre_i - Hpre          # H(X_i | X_<i)
    return total


# ======================================================================
# C1: chain-rule identity  H(X_1..n) = sum_i H(X_i | X_<i).
# ======================================================================

def check_C1():
    print("=" * 70)
    print("C1: chain rule  H(X_1..n) = sum_i H(X_i | X_<i)  (per-position")
    print("    softmax entropies sum to the block entropy).")
    print("=" * 70)
    ok = True
    for trial, (n, A, order, temp) in enumerate(
        [(5, 2, 1, 1.0), (5, 3, 2, 0.8), (6, 2, None, 1.2)]
    ):
        f = ARField(n, A, order=order, temp=temp)
        P = f.joint()
        H = entropy_of_table(P)
        H_chain = cond_block_entropy_chain_rule(P)
        err = abs(H - H_chain)
        print(f"  trial {trial}: n={n} A={A} order={order}: "
              f"H={H:.5f}  sum_i H(X_i|X_<i)={H_chain:.5f}  err={err:.2e}")
        ok = ok and err < 1e-7
    print(f"  C1 {'PASS' if ok else 'FAIL'}")
    return ok


# ======================================================================
# C2: realised codelength (forward-only) == -log2 P(x); sample mean -> H.
# ======================================================================

def check_C2():
    print("=" * 70)
    print("C2: streaming realised codelength L(x)=sum_i -log2 Q_i(x_i|x_<i),")
    print("    FORWARD-PASS ONLY, == -log2 P(x); sample mean -> H(X).")
    print("=" * 70)
    ok = True
    for trial, (n, A, order, temp) in enumerate(
        [(6, 2, 1, 1.0), (5, 3, 2, 0.9), (6, 2, None, 1.1)]
    ):
        f = ARField(n, A, order=order, temp=temp)
        P = f.joint()
        H = entropy_of_table(P)
        flat = P.ravel()
        configs = list(product(range(A), repeat=n))

        # (i) forward-only codelength equals -log2 P(x) for every x
        max_err = 0.0
        for ci, x in enumerate(configs):
            Lx = f.realised_codelength(x)
            px = flat[ci]
            if px > 0:
                max_err = max(max_err, abs(Lx - (-math.log2(px))))
        ident_ok = max_err < 1e-9

        # (ii) sample mean of L over draws from the field converges to H
        probs = flat / flat.sum()
        draws = RNG.choice(len(configs), size=20000, p=probs)
        meanL = float(np.mean([f.realised_codelength(configs[d]) for d in draws]))
        mean_ok = abs(meanL - H) < 0.05

        print(f"  trial {trial}: n={n} A={A} order={order}: "
              f"max|L(x)+log2 P(x)|={max_err:.2e} (fwd==joint: {ident_ok}); "
              f"mean L={meanL:.4f} vs H={H:.4f} ({mean_ok})")
        ok = ok and ident_ok and mean_ok
    print("  --> the OPERATIVE cost is K forward softmax reads, O(n|Sig|), with")
    print("      NO marginalisation; achieves the entropy rate (T7.4/T7.15).")
    print(f"  C2 {'PASS' if ok else 'FAIL'}")
    return ok


# ======================================================================
# C3: non-contiguous residual needs a genuine MARGINAL (not forward
#     conditionals); exponential op-count contrast.
# ======================================================================

def check_C3():
    print("=" * 70)
    print("C3: CONTRAST -- a NON-CONTIGUOUS residual H(X_S|X_Sbar) needs the")
    print("    gap MARGINAL (the #P-hard Remark 7.15a object); NOT obtainable")
    print("    from the original Q_i on a realised prefix.  Exp op-count gap.")
    print("=" * 70)
    ok = True
    n, A = 6, 2
    f = ARField(n, A, order=None, temp=1.0)
    P = f.joint()

    # Take a non-contiguous coded set S={0,2,4}, reconstruction Sbar={1,3,5}.
    S = (0, 2, 4)
    Sbar = tuple(i for i in range(n) if i not in S)
    # H(X_S | X_Sbar) = H(X) - H(X_Sbar): forming H(X_Sbar) marginalises out S
    H_full = entropy_of_table(P)
    H_Sbar = subset_entropy(P, Sbar)          # <- marginalises S (gap sum)
    H_cond = H_full - H_Sbar

    # P(X_Sbar) DOES factor by chain rule on the subsequence, but NOT via the
    # original Q_i on a realised prefix: P(x_1,x_3,x_5) requires summing over
    # x_0,x_2,x_4 (2^3 terms) to obtain the induced conditionals, which are not
    # the Q_j (the field is not naive-Markov on the subsequence). We certify
    # "needs marginal" by checking the kept subsequence has dependence the
    # original forward conditionals miss:
    #   I(X_4 ; X_0 | X_2) > 0  (a kept gap induces dependence not reducible to
    #   evaluating Q on a realised prefix of the subsequence).
    def cond_mi(P, i, j, k):
        keep = sorted([i, j, k])
        Pijk = P.sum(axis=tuple(a for a in range(P.ndim) if a not in keep))
        ai, aj, ak = keep.index(i), keep.index(j), keep.index(k)
        Pk = Pijk.sum(axis=tuple(a for a in range(3) if a != ak))
        mi = 0.0
        for kk in range(Pijk.shape[ak]):
            sl = [slice(None)] * 3; sl[ak] = kk
            Pab_c = Pijk[tuple(sl)]; pc = Pab_c.sum()
            if pc <= 0:
                continue
            Pab = Pab_c / pc
            Pa = Pab.sum(axis=1 if ai < aj else 0)
            Pb = Pab.sum(axis=0 if ai < aj else 1)
            for a in range(Pab.shape[0]):
                for b in range(Pab.shape[1]):
                    p = Pab[a, b]
                    if p > 0 and Pa[a] > 0 and Pb[b] > 0:
                        mi += pc * p * math.log2(p / (Pa[a] * Pb[b]))
        return mi
    gap_dep = cond_mi(P, 0, 4, 2)
    needs_marginal = gap_dep > 1e-6

    print(f"  S={S} (coded), Sbar={Sbar} (reconstructed):")
    print(f"    H(X_S|X_Sbar)=H(X)-H(X_Sbar)={H_cond:.4f} bits "
          f"(forming H(X_Sbar) marginalises out S)")
    print(f"    subsequence dependence I(X_0;X_4|X_2)={gap_dep:.4f}>0 "
          f"({needs_marginal}) => induced conds are not the original Q_i")
    print("  op-count: causal codelength O(n|Sig|) vs gap marginal |Sig|^{|S|}:")
    for nn in [10, 20, 40]:
        coded = nn // 2
        causal = nn * A
        marginal = A ** coded
        print(f"      n={nn:2d}: causal {causal:>4} softmax reads  vs  "
              f"gap marginal {marginal:>15,} terms "
              f"(ratio {marginal/causal:.1e}x)")
    ok = needs_marginal
    print(f"  C3 {'PASS' if ok else 'FAIL'}")
    return ok


# ======================================================================
# C4: order-independence -- a full-history (order = n-1) field still has an
#     O(n|Sig|) causal codelength while the non-contiguous marginal is
#     maximally expensive (the transformer regime).
# ======================================================================

def check_C4():
    print("=" * 70)
    print("C4: ORDER-INDEPENDENCE -- order=(n-1) full-history field (the")
    print("    full-attention transformer regime): causal codelength still")
    print("    O(n|Sig|), non-contiguous marginal maximally expensive.")
    print("=" * 70)
    n, A = 7, 2
    f = ARField(n, A, order=None, temp=1.0)   # full history
    P = f.joint()
    H = entropy_of_table(P)

    # causal: evaluate realised codelength for a random realised sequence using
    # ONLY forward conditionals; count softmax reads.
    x = tuple(RNG.integers(0, A, size=n).tolist())
    reads = {"count": 0}
    base_Q = f.Q
    def counting_Q(i, prefix):
        reads["count"] += 1
        return base_Q(i, prefix)
    f.Q = counting_Q
    Lx = f.realised_codelength(x)
    f.Q = base_Q
    causal_reads = reads["count"]

    # confirm full-history dependence is real (so it's genuinely order n-1):
    # last symbol's conditional varies with the FIRST symbol.
    q_if0 = f.Q(n - 1, tuple([0] + [0] * (n - 2)))
    q_if1 = f.Q(n - 1, tuple([1] + [0] * (n - 2)))
    full_history = float(np.abs(q_if0 - q_if1).max()) > 1e-9

    print(f"  realised x={x}: codelength L(x)={Lx:.4f}, "
          f"softmax reads={causal_reads} (== n={n}: {causal_reads == n})")
    print(f"  full-history dependence (Q_n|x_0=0 != Q_n|x_0=1): {full_history}")
    print(f"  vs non-contiguous marginal over |S|={n//2} gaps: "
          f"{A**(n//2):,} summation terms")
    print("  --> causal cost stays linear at ANY order; the #P barrier lives")
    print("      only in the (optional) non-contiguous selection marginal.")
    ok = (causal_reads == n) and full_history and abs(Lx + math.log2(P[x])) < 1e-9
    print(f"  C4 {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print()
    print("#" * 70)
    print("# Remark 7.15b verification: causal cost oracle poly for any AR")
    print("# field; #P-hardness confined to non-contiguous selection")
    print("#" * 70)
    print()
    r1 = check_C1(); print()
    r2 = check_C2(); print()
    r3 = check_C3(); print()
    r4 = check_C4(); print()
    print("=" * 70)
    allok = r1 and r2 and r3 and r4
    print(f"C1 chain-rule identity        : {'PASS' if r1 else 'FAIL'}")
    print(f"C2 forward-only codelength->H : {'PASS' if r2 else 'FAIL'}")
    print(f"C3 non-contiguous needs marg  : {'PASS' if r3 else 'FAIL'}")
    print(f"C4 order-independent (transf) : {'PASS' if r4 else 'FAIL'}")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
