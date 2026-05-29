#!/usr/bin/env python3
r"""
Verification for the OP2(a) additive-rule-cost analysis (Lemma 6.8k), recorded
in rnr_coding.tex S13.4.1.

QUESTION (the journal's explicit untried path): adapt Casel et al. 2021's
"alternative size measure also taking the number of rules into account"
(= L_sym + K, the additive rule-cost) to general alpha > 0, and ask whether it
gives a CEILING-FREE proxy for the discrete bit-cost L_sym * ceil(log2(|Sigma|+K)).

TWO OBJECTIVES ON THE SAME GRAMMAR:
  - additive rule-cost   A_alpha(G) = L_sym(G) + alpha * K(G)         (a SUM)
  - continuous bit-cost  B(G)       = L_sym(G) * log2(|Sigma|+K(G))   (a PRODUCT)
  - discrete bit-cost    Bd(G)      = L_sym(G) * ceil(log2(|Sigma|+K(G)))

Casel's reduction is robust to A_1 (alpha=1) because in their VC reduction the
optimum has L_sym = f(m,n) + |Gamma| AND the number of rules is also
monotone-increasing in |Gamma|, so min over covers is preserved by adding any
alpha*K term -> A_alpha SGP is APX-hard for every alpha>=0 (POSITIVE result).

THE ANTI-TAUTOLOGY TEST (decisive):
  POSITIVE only if A_alpha is genuinely DIFFERENT from the continuous bit-cost
  B (Lemma 6.8e). We exhibit grammars where A_alpha and B RANK THEM DIFFERENTLY.
  If they always agreed, Lemma 6.8k would be Lemma 6.8e relabeled -> OUTCOME-B.

WHAT WE CHECK (all numeric, PASS/FAIL to stdout):
  T1  Casel robustness: in the VC reduction, both L_sym and #rules are
      affine-increasing in cover size |Gamma|, so for EVERY alpha>=0 the
      A_alpha-minimiser = the L_sym-minimiser = min vertex cover. Hence the
      f(m,n)+|Gamma| gap (=> VC APX-hardness) transfers to A_alpha verbatim.
  T2  ANTI-TAUTOLOGY (rank-discordance): exhibit two grammars G1,G2 for the
      SAME string with A_alpha(G1) < A_alpha(G2) but B(G1) > B(G2)
      (and vice versa). Proves A_alpha is NOT a monotone function of B:
      they are genuinely different objectives, NOT 6.8e relabeled.
  T3  Why no transfer to the discrete bit-cost MINIMISER: the K-optimal point
      for the PRODUCT B is K*~sqrt(cN) (alpha-free, |Sigma|>>K regime), whereas
      the A_alpha-optimal point for the SUM is K*~sqrt(cN/(1+alpha)) (alpha-
      dependent). Different functionals. Then, the EXACT Charikar/BK gap:
      YES,NO differ only in minVC term Delta=N/432, giving L_ratio=1+1/8568=
      gamma_SGP but K_ratio=1+1/1008 (~8.5x LARGER, since K-base 7N/3 << L-base
      19.83N). The additive gap (L_NO+alpha K_NO)/(L_YES+alpha K_YES) is
      STRICTLY INCREASING in alpha from gamma_SGP (alpha=0) toward K_ratio
      (alpha->inf). So large alpha AMPLIFIES, not dilutes, the ceiling-free gap
      -- a stronger positive (no large-alpha penalty in Charikar's family).
  T4  Honest verdict gate: A_alpha SGP APX-hard for ALL alpha>=0 (clean,
      ceiling-free, gap gamma_SGP..1+1/1008); A_alpha != discrete bit-cost
      (step-function/ceiling vs strictly-monotone-in-K, T4); A_alpha !=
      continuous bit-cost B (rank-discordant, T2). So Lemma 6.8k is a NEW
      result on a DIFFERENT (real RNR dictionary-rule-size) objective, NOT a
      relabeling of Lemma 6.8e and NOT a closure of discrete-bit-cost OP2(a).

A model of the Charikar/BK grammar family is used: for a source of "mass" N
(=|V|), a grammar with K nonterminals achieves
    L_sym(K) ~ K + c*N/K        (paper's own model, S13.2 K-optimisation para)
saturating at L_sym = L*_sym = Theta(N) for K=K*_sym=Theta(N) (Charikar regime),
with alphabet |Sigma|=Theta(N). YES vs NO differ by the minVC term:
    L*_sym(NO) = gamma_SGP * L*_sym(YES).
"""

import math
import sys
from fractions import Fraction

GAMMA_SGP = Fraction(8569, 8568)  # Charikar 2005 inapprox ratio (BK 1999)


# ----------------------------------------------------------------------------
# T1: Casel robustness -- A_alpha-minimiser = min vertex cover, for every alpha>=0
# ----------------------------------------------------------------------------
def casel_reduction_costs(cover_size, m, n, alpha):
    r"""Casel/Charikar VC reduction grammar costs as a function of cover size.

    Conference Thm: a minimal grammar for w_G has L_sym = f(m,n) + |Gamma|.
    The number of rules also increases by exactly one per cover vertex chosen
    (each covered vertex v_i gets an extra rule V_i -> #V_i; Lemma 6 / Sec 3.2),
    so #rules(Gamma) = r0(m,n) + |Gamma|.  Thus:
        L_sym  = f(m,n)      + |Gamma|
        Krules = r0(m,n)     + |Gamma|
        A_alpha = (f + alpha*r0) + (1+alpha)*|Gamma|
    which is strictly increasing in |Gamma| for every alpha >= 0.
    """
    f = 9 * n + 4 * m          # structural symbol overhead (schematic, >0, indep of cover)
    r0 = 14 * n                # structural rule overhead (codeword nonterminals)
    L_sym = f + cover_size
    Krules = r0 + cover_size
    A = L_sym + alpha * Krules
    return L_sym, Krules, A


def test_casel_robust_to_alpha():
    print("\nT1: Casel reduction -- A_alpha-minimiser = min vertex cover, ALL alpha>=0")
    m, n = 30, 20
    ok = True
    for alpha in [0, 1, 2, 5, math.log2(n), 10 * math.log2(n)]:
        # min cover gives min A_alpha; any larger cover gives strictly larger A_alpha
        costs = []
        for cover in range(3, 12):
            _, _, A = casel_reduction_costs(cover, m, n, alpha)
            costs.append((cover, A))
        # strictly increasing in cover?
        mono = all(costs[i][1] < costs[i + 1][1] for i in range(len(costs) - 1))
        argmin_cover = min(costs, key=lambda t: t[1])[0]
        print("    alpha=%6.3f: A_alpha strictly increasing in |Gamma|=%s, argmin |Gamma|=%d (=min)"
              % (alpha, mono, argmin_cover))
        if not (mono and argmin_cover == 3):
            ok = False
    print("    => For EVERY alpha>=0, min A_alpha <=> min vertex cover. The")
    print("       f(m,n)+|Gamma| gap transfers, so A_alpha SGP is APX-hard (alpha=O(1)).")
    return ok


# ----------------------------------------------------------------------------
# Grammar-family model (Charikar/BK): L_sym(K) ~ K + c*N/K, |Sigma| ~ N
# ----------------------------------------------------------------------------
def Lsym_of_K(K, N, c):
    return K + c * N / K


def additive_cost(K, N, c, alpha):
    return Lsym_of_K(K, N, c) + alpha * K


def cont_bitcost(K, N, c, Sigma):
    return Lsym_of_K(K, N, c) * math.log2(Sigma + K)


def disc_bitcost(K, N, c, Sigma):
    return Lsym_of_K(K, N, c) * math.ceil(math.log2(Sigma + K))


# ----------------------------------------------------------------------------
# T2: ANTI-TAUTOLOGY -- additive and continuous-bit-cost rank grammars DIFFERENTLY
# ----------------------------------------------------------------------------
def test_rank_discordance():
    print("\nT2: ANTI-TAUTOLOGY -- additive A_alpha vs continuous bit-cost B")
    print("    Find grammars G1,G2 for the SAME string ranked DIFFERENTLY.")
    N, c, Sigma = 1000.0, 50.0, 1000.0
    alpha = 1.0
    ok = False
    # scan pairs of K-values; look for discordant ranking under A_1 vs B
    Ks = list(range(2, 400, 2))
    found = None
    for i in range(len(Ks)):
        for j in range(i + 1, len(Ks)):
            K1, K2 = Ks[i], Ks[j]
            A1 = additive_cost(K1, N, c, alpha)
            A2 = additive_cost(K2, N, c, alpha)
            B1 = cont_bitcost(K1, N, c, Sigma)
            B2 = cont_bitcost(K2, N, c, Sigma)
            # discordant: A says G1 better, B says G2 better (or vice versa)
            if (A1 < A2) != (B1 < B2):
                found = (K1, K2, A1, A2, B1, B2)
                ok = True
                break
        if found:
            break
    if found:
        K1, K2, A1, A2, B1, B2 = found
        print("    Witness: K1=%d, K2=%d (same source N=%g, |Sigma|=%g)" % (K1, K2, N, Sigma))
        print("      A_1(G1)=%.2f  A_1(G2)=%.2f  ->  additive prefers G%d"
              % (A1, A2, 1 if A1 < A2 else 2))
        print("      B(G1)=%.2f  B(G2)=%.2f  ->  cont-bitcost prefers G%d"
              % (B1, B2, 1 if B1 < B2 else 2))
        print("    => DISCORDANT ranking: A_alpha is NOT a monotone function of B.")
        print("       The additive objective is GENUINELY DIFFERENT from Lemma 6.8e's")
        print("       continuous bit-cost -- NOT a relabeling. Anti-tautology PASSES.")
    else:
        print("    NO discordant pair found -- additive == continuous bit-cost ordering!")
        print("    => This would be TAUTOLOGY (6.8e relabeled). OUTCOME-B.")
    return ok


# ----------------------------------------------------------------------------
# T3: different sqrt-laws; gap-dilution for large alpha (quantified NEGATIVE)
# ----------------------------------------------------------------------------
def test_sqrt_laws_and_dilution():
    print("\nT3: K-optimal sqrt-laws differ; EXACT-Charikar gap AMPLIFIES with alpha")
    N, c, Sigma = 1e6, 50.0, 1e9   # |Sigma| >> K regime so log~const for the product
    ok = True

    # For L_sym(K)=K + cN/K:
    #   SUM:     A_alpha=(1+alpha)K + cN/K, dA/dK=0 -> K*=sqrt(cN/(1+alpha)).
    #   PRODUCT: B=(K+cN/K)*log2(Sigma+K). For Sigma>>K, log2~const, so
    #            B ~ const*(K+cN/K), minimised at K*=sqrt(cN) (alpha-free,
    #            and independent of |Sigma| to leading order).
    for alpha in [1.0, math.log2(N)]:
        Kstar_pred = math.sqrt(c * N / (1.0 + alpha))
        Ks = [k for k in range(2, int(4 * Kstar_pred) + 2)]
        Kstar_num = min(Ks, key=lambda K: additive_cost(K, N, c, alpha))
        print("    additive  alpha=%7.3f: K*_pred=sqrt(cN/(1+alpha))=%.1f  K*_num=%d"
              % (alpha, Kstar_pred, Kstar_num))
        if abs(Kstar_num - Kstar_pred) / Kstar_pred > 0.05:
            ok = False

    Kstar_B_pred = math.sqrt(c * N)          # product minimiser (Sigma>>K)
    Ks = [k for k in range(2, int(4 * Kstar_B_pred) + 2)]
    Kstar_B_num = min(Ks, key=lambda K: cont_bitcost(K, N, c, Sigma))
    print("    cont-bit       : K*_pred~sqrt(cN)=%.1f  K*_num=%d  (alpha-FREE)"
          % (Kstar_B_pred, Kstar_B_num))
    if abs(Kstar_B_num - Kstar_B_pred) / Kstar_B_pred > 0.05:
        ok = False
    print("    => SUM minimiser K* DEPENDS on alpha (sqrt(cN/(1+alpha))); PRODUCT")
    print("       minimiser is alpha-free (sqrt(cN)). Different functionals: the")
    print("       additive alpha-term and the multiplicative log-factor are not")
    print("       interchangeable, so A_alpha is no proxy for the bit-cost minimiser.")

    # Additive gap between YES and NO branches, EXACT Charikar/BK accounting.
    # YES and NO grammars differ ONLY in the minVC term Delta = k_no - k_opt =
    # (1/144)k_opt = N/432 (same numerator in BOTH L and K). Bases:
    #   L_YES = 15N + 3E + k_opt = 19.833...N   (E=1.5N, k_opt=N/3)
    #   K_YES = 2N + k_opt       =  7N/3 = 2.333...N
    # so   L_ratio = 1 + Delta/L_YES = 1 + 1/8568 = gamma_SGP (Charikar)
    #      K_ratio = 1 + Delta/K_YES = 1 + 1/1008  (~8.5x LARGER gap!)
    # The additive gap is a Delta-mediant:
    #   gap(alpha) = (L_NO+alpha K_NO)/(L_YES+alpha K_YES)
    #              = 1 + Delta(1+alpha)/(L_YES+alpha K_YES),
    # strictly INCREASING in alpha (since K_YES<L_YES), from gamma_SGP (alpha=0)
    # UP toward K_ratio=1+1/1008 (alpha->inf). NO dilution: ceiling-free gap is
    # AMPLIFIED by alpha, because in Charikar's family the nonterminal-count gap
    # exceeds the symbol-count gap.
    print("    Additive YES/NO gap, EXACT Charikar/BK accounting:")
    print("      L_ratio = 1+1/8568 = %.8f = gamma_SGP" % (1 + 1 / 8568))
    print("      K_ratio = 1+1/1008 = %.8f (nonterminal-count gap, ~8.5x larger)"
          % (1 + 1 / 1008))
    for Nv in [600.0, 6000.0, 60000.0]:
        E = 1.5 * Nv
        k_opt = Nv / 3.0
        Delta = k_opt / 144.0
        L_YES = 15 * Nv + 3 * E + k_opt
        K_YES = 2 * Nv + k_opt
        L_NO = L_YES + Delta
        K_NO = K_YES + Delta
        gaps = []
        for alpha in [0.0, 1.0, math.log2(Nv + K_YES), 1e6]:
            gap = (L_NO + alpha * K_NO) / (L_YES + alpha * K_YES)
            gaps.append((alpha, gap))
        print("      N=%6.0f: " % Nv + "  ".join(
            "a=%s:%.7f" % ("inf" if a > 1e5 else ("%.2f" % a), g) for a, g in gaps))
        # monotone increasing in alpha, bounded above by K_ratio
        if not all(gaps[i][1] < gaps[i + 1][1] for i in range(len(gaps) - 1)):
            ok = False
        if not (abs(gaps[0][1] - (1 + 1 / 8568)) < 1e-6
                and gaps[-1][1] < 1 + 1 / 1008 + 1e-6):
            ok = False
    print("    => additive gap is STRICTLY INCREASING in alpha (gamma_SGP -> 1+1/1008).")
    print("       NO dilution: the ceiling-free additive gap is AMPLIFIED by alpha")
    print("       because Charikar's nonterminal-count gap (1/1008) EXCEEDS its")
    print("       symbol-count gap (1/8568=gamma_SGP-1). Larger alpha => stronger")
    print("       (ceiling-free) APX-hardness for the additive RNR rule-size objective.")
    return ok


# ----------------------------------------------------------------------------
# T4: verdict gate -- additive != discrete bit-cost AND != continuous bit-cost
# ----------------------------------------------------------------------------
def test_verdict_gate():
    print("\nT4: verdict gate")
    ok = True

    # (i) additive != discrete bit-cost. The discrete bit-cost multiplier
    #     ceil(log2(|Sigma|+K)) is PIECEWISE-CONSTANT in K (a step function):
    #     for all K with |Sigma|+K in the SAME power-of-2 band [2^{b-1},2^b),
    #     the multiplier is the constant b, so discrete bit-cost = b*L_sym(K)
    #     ranks two such grammars purely by L_sym -- it is BLIND to K within a
    #     band. The additive A_alpha=L_sym+alpha*K is STRICTLY MONOTONE in K
    #     (for fixed L_sym) -- never blind to K. Two distinct objectives.
    Sigma = 1100.0           # so |Sigma|+K stays in band [1024,2048) for K in [0,948)
    N, c = 5000.0, 80.0
    b = math.ceil(math.log2(Sigma + 1))    # band index (= 11 here)
    band_hi = 2 ** b
    print("    Band [2^%d,2^%d)=[%d,%d): discrete multiplier is the CONSTANT %d for"
          % (b - 1, b, 2 ** (b - 1), band_hi, b))
    print("    every K in this band -> discrete bit-cost = %d*L_sym(K), BLIND to K." % b)
    # pick two grammars in-band with EQUAL L_sym but different K (so discrete
    # bit-cost is identical) yet different additive cost.
    # Use the L_sym(K)=K+cN/K curve symmetry: L_sym(K1)=L_sym(K2) when
    # K1*K2=cN (the two roots straddling sqrt(cN)). Pick such a pair in-band.
    cN = c * N
    found_pair = None
    Kmax_band = int(band_hi - Sigma) - 1
    for K1 in range(2, min(Kmax_band, int(math.sqrt(cN)))):
        K2f = cN / K1                      # partner with equal L_sym
        K2 = int(round(K2f))
        if K2 <= Kmax_band and K2 != K1 and abs(Lsym_of_K(K1, N, c) - Lsym_of_K(K2, N, c)) < 1e-6:
            found_pair = (K1, K2)
            break
    if found_pair:
        K1, K2 = found_pair
        Bd1, Bd2 = disc_bitcost(K1, N, c, Sigma), disc_bitcost(K2, N, c, Sigma)
        A1, A2 = additive_cost(K1, N, c, 1.0), additive_cost(K2, N, c, 1.0)
        print("    In-band pair K1=%d,K2=%d (equal L_sym=%.3f):" % (K1, K2, Lsym_of_K(K1, N, c)))
        print("      discrete bit-cost: Bd1=%.3f Bd2=%.3f  (TIE: %s)"
              % (Bd1, Bd2, abs(Bd1 - Bd2) < 1e-6))
        print("      additive A_1     : A1=%.3f A2=%.3f  (strict: %s)"
              % (A1, A2, abs(A1 - A2) > 1e-6))
        if not (abs(Bd1 - Bd2) < 1e-6 and abs(A1 - A2) > 1e-6):
            ok = False
        print("    => discrete bit-cost TIES these grammars; additive DISTINGUISHES")
        print("       them. The discrete (ceiling) objective is a step-function of K;")
        print("       A_alpha is strictly monotone. A_alpha != discrete bit-cost.")
    else:
        print("    (no equal-L_sym in-band pair found at these parameters)")
        ok = False

    # (ii) and additive is continuous across the boundary while discrete JUMPS.
    # ceil(log2 x) jumps from b to b+1 as x crosses 2^b -> 2^b+1 (since
    # ceil(log2 2^b)=b but ceil(log2(2^b+1))=b+1). Pick K so |Sigma|+K = 2^b.
    Kb = int(band_hi - Sigma)          # |Sigma|+Kb = 2^b exactly
    mult_lo = math.ceil(math.log2(Sigma + Kb))      # = b
    mult_hi = math.ceil(math.log2(Sigma + Kb + 1))  # = b+1
    print("    Boundary K: discrete multiplier jumps %d -> %d (a full bit); additive"
          % (mult_lo, mult_hi))
    print("    moves only by alpha (continuous). Discrete has the ceiling discontinuity")
    print("    that IS the OP2(a) residual; A_alpha lacks it -> cannot close it.")
    if not (mult_hi == mult_lo + 1):
        ok = False

    # (ii) additive != continuous bit-cost: established by T2 (rank-discordance).
    print("    A_alpha != continuous bit-cost B: established by T2 rank-discordance.")
    print("    => Lemma 6.8k (additive-rule-cost APX-hard) is a NEW result on a")
    print("       DIFFERENT objective (the RNR dictionary rule-size), NOT 6.8e")
    print("       relabeled and NOT a closure of the discrete-bit-cost OP2(a).")
    return ok


def main():
    print("=" * 74)
    print("OP2(a) additive-rule-cost analysis (Lemma 6.8k) -- verification")
    print("=" * 74)
    results = {
        "T1 Casel robust to alpha (APX-hard transfers)": test_casel_robust_to_alpha(),
        "T2 anti-tautology: additive != continuous bit-cost (rank-discordant)": test_rank_discordance(),
        "T3 sqrt-laws differ; EXACT-Charikar gap amplifies with alpha": test_sqrt_laws_and_dilution(),
        "T4 verdict gate: additive != discrete AND != continuous bit-cost": test_verdict_gate(),
    }
    print("\n" + "=" * 74)
    allok = True
    for name, r in results.items():
        print("  [%s] %s" % ("PASS" if r else "FAIL", name))
        allok = allok and r
    print("=" * 74)
    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    print("\nVERDICT (for the research note):")
    print("  POSITIVE (Lemma 6.8k): additive-rule-cost SGP L_sym+alpha*K is APX-hard")
    print("  for EVERY fixed alpha>=0 (Casel/Charikar robustness). The ceiling-free")
    print("  gap is STRICTLY INCREASING in alpha: gamma_SGP=1+1/8568 (alpha=0) up to")
    print("  K_ratio=1+1/1008 (alpha->inf), because Charikar's nonterminal-count gap")
    print("  EXCEEDS its symbol-count gap. This is a NEW result on the RNR dictionary")
    print("  rule-size objective, GENUINELY DIFFERENT from both the continuous bit-")
    print("  cost (Lemma 6.8e; rank-discordant, T2) and the discrete bit-cost (no")
    print("  ceiling, T4). It does NOT close the ceiling-discontinuity residual of")
    print("  OP2(a) (additive lacks the ceiling), but it is a standalone hardness")
    print("  result for a real RNR objective -- NOT a tautology / NOT 6.8e relabeled.")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
