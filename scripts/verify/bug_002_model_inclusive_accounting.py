#!/usr/bin/env python3
r"""
Verification of the model-inclusive vs. amortized accounting (BUG-002-B/C):
the consolidated table in section 10.2 and the break-even caveat of
Theorem 7.15, together with the Theorem 7.21 total-cost optimum.

INTERNAL-CONSISTENCY ONLY. This probe recomputes the paper's *own*
arithmetic and asserts it matches the stated figures. It does NOT measure
the empirical H_M' = 0.664 bpb datum (that is the Deletang et al. 2024
input constant, verified externally -- BUG-002-D). It treats the paper's
stated constants as given and checks the derived quantities:

  (1) self-contained single-archive size vs. amortized bitstream size, and
      the ~140x single-archive expansion (negative savings);
  (2) the break-even volume V* = 140*8/(8 - 0.664) ~ 152 GB from the
      identity V* * (log2|Sigma| - H_M') = L(M) * log2|Sigma|, plus the
      general identity for arbitrary (L(M), |Sigma|, H_M');
  (3) the headline bitstream figures N*H_M' ~ 79 MB and 12.0x = 8/0.664;
  (4) Theorem 7.21's TotalCost optimum L*(N) = (N C alpha)^{1/(1+alpha)},
      the sub-linear N^{1/(1+alpha)} overhead scaling, and the
      monotone "net savings positive once amortized V > V*" check.

Stated constants (source: Theorem 7.15 caveat / worked example, Theorem 7.21):
  H_M'        = 0.664 bpb          (Deletang et al. 2024 ICLR Table 1)
  |Sigma|     = 256  => log2|Sigma| = 8 bits/byte
  L(M)        = 140 GB             (Chinchilla-70B, fp16)
  N (enwik9)  = 1 GB
  alpha       in {0.34, 0.5, 1.0}  (Kaplan 2020 / Hoffmann 2022)

Emits PASS/FAIL per check and "OVERALL -> PASS" on success (CI contract).
"""

import math
import sys

# ---- stated constants (the paper's source-of-truth inputs) --------------
H_MP = 0.664            # bpb, Deletang et al. 2024 Table 1 (T7.15)
SIGMA_BITS = 8.0        # log2|Sigma|, |Sigma| = 256 (T7.15)
MODEL_GB = 140.0        # L(M), Chinchilla-70B fp16 (T7.15 caveat)
N_GB = 1.0              # enwik9 source size (T7.15 worked example)

GB_BITS = 8.0 * 1e9     # 1 GB (decimal) in bits

# tolerances
ATOL_RATIO = 0.05       # "~12.0x", "~140x" tolerances
ATOL_VSTAR_GB = 1.0     # "~152 GB"
RTOL_OPT = 1e-9         # calculus optimum (closed form)


def line(ok, label, detail):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {label}: {detail}")
    return ok


# ---- (1) self-contained vs amortized; ~140x expansion -------------------
def check_self_contained_vs_amortized():
    print("\n(1) Self-contained (model-inclusive) vs. amortized (bitstream-only):")
    ok = True

    # amortized bitstream-only archive at N = 1 GB
    bitstream_bits = N_GB * GB_BITS * (H_MP / SIGMA_BITS)   # N * H_M'
    bitstream_MB = bitstream_bits / 8 / 1e6                 # decimal MB
    bitstream_MiB = bitstream_bits / 8 / 1024**2
    # paper writes "~79 MB"; N*H_M' = 0.664 Gbit / 8 = 83 MB decimal = 79.15 MiB
    ok &= line(
        abs(bitstream_MiB - 79.0) < 1.0,
        "amortized bitstream N*H_M'",
        f"{bitstream_MiB:.2f} MiB ({bitstream_MB:.1f} MB decimal); paper ~79 MB",
    )

    # self-contained single archive: model + bitstream charged to one 1 GB source
    self_contained_GB = MODEL_GB + N_GB * (H_MP / SIGMA_BITS)  # 140 + 0.083
    expansion = self_contained_GB / N_GB
    # paper: "expands the input by a factor of roughly 140x"
    ok &= line(
        abs(expansion - 140.0) < MODEL_GB * ATOL_RATIO,  # within ~7 of 140
        "single-archive expansion factor",
        f"{expansion:.1f}x (self-contained {self_contained_GB:.2f} GB / {N_GB:.0f} GB); paper ~140x",
    )

    # negative savings: self-contained archive strictly LARGER than raw source
    ok &= line(
        self_contained_GB > N_GB,
        "single-archive net savings sign",
        f"self-contained {self_contained_GB:.2f} GB > raw {N_GB:.0f} GB => negative savings (expansion)",
    )
    return ok


# ---- (2) break-even V* identity -----------------------------------------
def breakeven_volume(model_size, sigma_bits, h_mp):
    """V* * (sigma_bits - H_M') = model_size * sigma_bits."""
    return model_size * sigma_bits / (sigma_bits - h_mp)


def check_breakeven():
    print("\n(2) Break-even volume V* from V*(log2|Sigma| - H_M') = L(M) log2|Sigma|:")
    ok = True

    vstar = breakeven_volume(MODEL_GB, SIGMA_BITS, H_MP)
    # paper: V* ~ 140*8/(8 - 0.664) ~ 152 GB
    explicit = 140.0 * 8 / (8 - 0.664)
    ok &= line(
        abs(vstar - 152.0) < ATOL_VSTAR_GB,
        "break-even V*",
        f"{vstar:.2f} GB; explicit 140*8/(8-0.664)={explicit:.2f}; paper ~152 GB",
    )
    ok &= line(
        abs(vstar - explicit) < 1e-9,
        "V* matches the explicit arithmetic in the caveat",
        f"{vstar:.6f} == {explicit:.6f}",
    )

    # general identity holds for arbitrary inputs (no new claim; structural check)
    print("  general identity V*(log2|Sigma| - H_M') = L(M) log2|Sigma| (arbitrary inputs):")
    cases = [(140.0, 8.0, 0.664), (2.0, 8.0, 1.0), (70.0, 8.0, 0.5), (10.0, 4.0, 1.2)]
    for (Lm, sb, h) in cases:
        v = breakeven_volume(Lm, sb, h)
        lhs = v * (sb - h)
        rhs = Lm * sb
        ok &= line(
            abs(lhs - rhs) < 1e-9 * max(1.0, abs(rhs)),
            f"identity L(M)={Lm}GB, log2|S|={sb}, H={h}",
            f"V*={v:.3f}GB; LHS={lhs:.4f} == RHS={rhs:.4f}",
        )

    # break-even must be sensible: above the raw source size (else no expansion)
    ok &= line(
        vstar > N_GB,
        "V* > single-archive N (single docs do NOT amortize)",
        f"V*={vstar:.1f} GB >> N={N_GB:.0f} GB",
    )
    return ok


# ---- (3) headline bitstream ratio ---------------------------------------
def check_headline_ratio():
    print("\n(3) Headline bitstream compression ratio:")
    ok = True
    ratio = SIGMA_BITS / H_MP    # 8 / 0.664
    ok &= line(
        abs(ratio - 12.0) < ATOL_RATIO,
        "bitstream ratio log2|Sigma| / H_M'",
        f"{ratio:.3f}x; paper ~12.0x",
    )
    return ok


# ---- (4) Theorem 7.21 TotalCost optimum + monotone net-savings ----------
def optimal_L(N, C, alpha):
    """L*(N) = (N C alpha)^{1/(1+alpha)} (T7.21)."""
    return (N * C * alpha) ** (1.0 / (1.0 + alpha))


def total_cost(L, N, h_X, C, alpha):
    """TotalCost = L + N*h(X) + N*C*L^{-alpha} (+ O(sqrt N), dropped here)."""
    return L + N * h_X + N * C * (L ** (-alpha))


def check_t721_optimum():
    print("\n(4) Theorem 7.21 TotalCost optimum L*(N) = (N C alpha)^{1/(1+alpha)}:")
    ok = True
    h_X = 0.8
    C = 10.0
    for alpha in (0.34, 0.5, 1.0):
        exponent = 1.0 / (1.0 + alpha)
        for N in (1e6, 1e9, 1e12):
            Lopt = optimal_L(N, C, alpha)

            # (a) closed-form vs. numerical minimum agreement
            TC_opt = total_cost(Lopt, N, h_X, C, alpha)
            # numerical first-derivative root: d/dL = 1 - N C alpha L^{-(1+alpha)} = 0
            deriv = 1.0 - N * C * alpha * Lopt ** (-(alpha + 1.0))
            ok &= line(
                abs(deriv) < 1e-6 * max(1.0, Lopt),
                f"alpha={alpha}, N={N:.0e}: dTotalCost/dL = 0 at L*",
                f"L*={Lopt:.3e}, deriv={deriv:.2e}",
            )

            # (b) sub-linear: L*/N -> 0, exponent < 1
            ok &= line(
                Lopt < N and exponent < 1.0,
                f"alpha={alpha}, N={N:.0e}: sub-linear L* (exp={exponent:.4f})",
                f"L*/N = {Lopt / N:.2e}",
            )

            # (c) overhead = (1 + 1/alpha) L* = Theta(N^{1/(1+alpha)})
            overhead = TC_opt - N * h_X
            predicted = (1.0 + 1.0 / alpha) * Lopt
            ok &= line(
                abs(overhead - predicted) < 1e-6 * predicted,
                f"alpha={alpha}, N={N:.0e}: overhead = (1+1/alpha)L*",
                f"overhead={overhead:.3e} == {predicted:.3e}",
            )

        # (d) scaling exponent of the overhead matches N^{1/(1+alpha)}:
        #     overhead(N2)/overhead(N1) == (N2/N1)^{1/(1+alpha)}
        N1, N2 = 1e9, 1e12
        o1 = (1.0 + 1.0 / alpha) * optimal_L(N1, C, alpha)
        o2 = (1.0 + 1.0 / alpha) * optimal_L(N2, C, alpha)
        observed = math.log(o2 / o1) / math.log(N2 / N1)
        ok &= line(
            abs(observed - exponent) < 1e-9,
            f"alpha={alpha}: overhead scaling exponent",
            f"measured {observed:.6f} == 1/(1+alpha) {exponent:.6f}",
        )
    return ok


# ---- (5) monotone net-savings crossing at V* ----------------------------
def check_monotone_crossing():
    print("\n(5) Net savings becomes positive once amortized volume V > V*:")
    ok = True
    vstar = breakeven_volume(MODEL_GB, SIGMA_BITS, H_MP)
    # self-contained archive over V GB of text shared by one fixed L(M):
    #   archive(V) = L(M) + V * (H_M' / log2|Sigma|)   [GB]
    #   raw(V)     = V                                  [GB]
    #   net saving = raw - archive ; positive iff V > V*
    def net_saving(V):
        archive = MODEL_GB + V * (H_MP / SIGMA_BITS)
        return V - archive
    for V in (10.0, 100.0, vstar - 5.0, vstar + 5.0, 500.0, 5000.0):
        ns = net_saving(V)
        expect_pos = V > vstar
        ok &= line(
            (ns > 0) == expect_pos,
            f"V={V:8.1f} GB",
            f"net saving {ns:+.2f} GB ({'positive' if ns > 0 else 'negative/zero'}); "
            f"V {'>' if V > vstar else '<='} V*={vstar:.1f}",
        )
    # exact crossing at V*: net saving ~ 0
    ok &= line(
        abs(net_saving(vstar)) < 1e-6,
        "net saving = 0 exactly at V = V*",
        f"net saving({vstar:.3f}) = {net_saving(vstar):.2e} GB",
    )
    return ok


def main() -> int:
    print("Verification: model-inclusive vs. amortized accounting (BUG-002)")
    print("=" * 70)
    print("INTERNAL-CONSISTENCY ONLY -- recomputes the paper's stated arithmetic;")
    print("does NOT measure the empirical 0.664 bpb datum (Deletang 2024 input).")

    results = {
        "self-contained vs amortized / ~140x expansion": check_self_contained_vs_amortized(),
        "break-even V* identity (~152 GB)": check_breakeven(),
        "headline bitstream ratio (~12.0x)": check_headline_ratio(),
        "T7.21 TotalCost optimum L*(N) & scaling": check_t721_optimum(),
        "monotone net-savings crossing at V*": check_monotone_crossing(),
    }

    print("\n" + "=" * 70)
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
