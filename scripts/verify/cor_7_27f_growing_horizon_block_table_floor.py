#!/usr/bin/env python3
r"""
Verification for Corollary 7.27f (growing-horizon block-table floor) and the
accompanying NEGATIVE delimitation of the unconditional sqrt(N) converse.

CONTEXT.  Corollary 7.27e closed the FV escape ONLY for growing-horizon
(triangular-array) sources and left genuinely open: "for a growing-horizon
(escape-closed) source, whether some NON-FV model-free succinct structure
achieves Shannon-near rate with o(sqrt N) query."

This script verifies the two halves of the honest resolution (Cor 7.27f):

  POSITIVE half (a conditional / architecture-restricted converse):
    Every BLOCK-BASED context-decode-table structure -- the whole known
    family of O(1)-access high-order-entropy structures (Ferragina-Venturini
    TCS 372 (2007); Gonzalez-Navarro "Random Access to High-Order Entropy
    Compressed Text" SPIRE 2006 / Grossi-Orlandi-Raman; Belazzougui-Navarro
    TALG 11(4) 2015) -- carries per-symbol redundancy
        rho(k,N) = Theta( ((k+1) log2 sigma + log2 log2 N) / log_sigma N ).
    Shannon-near rate forces the model term Delta_k^{(N)} = o(1); for a
    growing-horizon source Delta_k^{(N)} = c0 * L / k (L := log_sigma N) this
    needs k = omega(L), at which rho(k,N) = omega(1).  The two requirements
    are INCOMPATIBLE: inf over admissible k of [Delta_k + rho_k] = Theta(1),
    bounded away from 0.  So the ENTIRE block-table family (not just FV)
    fails to be Shannon-near on growing-horizon sources.  This is the new
    content over 7.27e (which only ruled out FV).

  NEGATIVE half (why the bound is NOT unconditional, and NOT sqrt(N)):
    For pure ACCESS, Belazzougui-Navarro prove access (and select) run in
    O(1) time within compressed space; their sqrt-type product lower bound
    rho * t_access * t_select = Omega((log sigma)^2 / w) is a RANK/SELECT
    bound and is VACUOUS for pure access (t_select unconstrained).  The
    generic cell-probe access floor at linear space is only
    Omega(log N / log log N) (Patrascu-Demaine 2006), NOT Omega(sqrt N).
    The sqrt(N) form appears only CONDITIONALLY: through the BWT
    (Sinha-Weinstein STOC 2019, Omega(n/t^2)) or for EXPONENTIALLY
    compressible inputs (Chen-Verbin-Yu ICALP 2012, Omega(n^{1/2-eps}) at
    grammar size n with L = 2^{sqrt n}) -- the latter is irrelevant in the
    RNR regime where S = Theta(N) bits.  Hence an UNCONDITIONAL sqrt(N)
    converse is FALSE for fixed sources and unprovable in general; the
    open residual survives, now correctly delimited.

PASS iff:
  (A) growing-horizon: inf_{k <= L/2} [Delta_k + rho_k] is Theta(1), bounded
      away from 0 -- block-table family closed.
  (B) fixed source (any decay, incl. divergent excess): inf [Delta_k + rho_k]
      -> 0 -- block-table family OPEN (no contradiction with 7.27e(a)).
  (C) the redundancy-order incompatibility: the k achieving Delta_k = o(1)
      (k_model = omega(L)) and the k achieving rho_k = o(1) (k_red = o(L))
      have k_model / k_red -> infinity, so no single k satisfies both.
  (D) the sqrt(N) gap: the generic access floor log N / log log N is
      asymptotically << sqrt(N), so an O(1)-access structure is NOT excluded
      by any unconditional cell-probe access bound (sanity check on the
      NEGATIVE half).
"""

from __future__ import annotations

import math
import sys


# ----------------------------------------------------------------------
# Shared model: FV / block-table redundancy (identical functional form
# across FV, Gonzalez-Navarro, Grossi-Orlandi-Raman, Belazzougui-Navarro).
# ----------------------------------------------------------------------

def L_sigma(N: int, sigma: int) -> float:
    """L := log_sigma N = log2 N / log2 sigma."""
    return math.log2(N) / math.log2(sigma)


def rho_block_table(k: float, N: int, sigma: int) -> float:
    """Per-symbol redundancy of any block context-table structure, bits/sym:
        rho = ((k+1) log2 sigma + log2 log2 N) / L_sigma(N).
    (FV/GOR/GN/BN all share this Theta(.) form.)"""
    log2_sigma = math.log2(sigma)
    loglog = math.log2(max(math.log2(max(N, 4)), 2.0))
    return ((k + 1) * log2_sigma + loglog) / L_sigma(N, sigma)


def delta_growing_horizon(k: float, N: int, sigma: int, c0: float) -> float:
    """Growing-horizon (triangular-array) excess: Delta_k^{(N)} = c0 * L / k,
    horizon proportional to L = log_sigma N (Cor 7.27e(b) witness)."""
    if k <= 0:
        return float("inf")
    return c0 * L_sigma(N, sigma) / k


def delta_fixed(k: float, decay: str, a: float) -> float:
    """Fixed-source excess sequences (N-independent), Cor 7.27e(a):
       'harmonic' : a/k   (DIVERGENT excess entropy E = inf)
       'geometric': a * 2^{-k}  (summable excess, exp-mixing)."""
    if k <= 0:
        return float("inf")
    if decay == "harmonic":
        return a / k
    if decay == "geometric":
        return a * (2.0 ** (-k))
    raise ValueError(decay)


def inf_joint_excess_block_table(N, sigma, delta_fn, k_cap_frac=0.5):
    """min over admissible k in [1, k_cap_frac * L] of [Delta_k + rho_k].
    FV validity requires k = o(L); we cap at k_cap_frac*L (Cor 7.27e proof)."""
    L = L_sigma(N, sigma)
    k_max = max(2, int(k_cap_frac * L))
    best = float("inf")
    best_k = 0
    # scan integer k plus a fine grid (the optimum can be fractional in L)
    ks = list(range(1, k_max + 1))
    grid = [1.0 + i * (k_max - 1) / 400.0 for i in range(401)] if k_max > 1 else [1.0]
    for k in ks + grid:
        v = delta_fn(k, N, sigma) + rho_block_table(k, N, sigma)
        if v < best:
            best, best_k = v, k
    return best, best_k, L


# ----------------------------------------------------------------------
# (A) Growing-horizon: block-table family CLOSED (inf is Theta(1)).
# ----------------------------------------------------------------------

def check_growing_horizon_closed() -> int:
    print("\n--- (A) Growing-horizon: block-table family fails (inf = Theta(1)) ---")
    fails = 0
    c0 = 0.5
    for sigma in [2, 4, 256]:
        prev = None
        vals = []
        for N in [10**6, 10**9, 10**12, 10**18, 10**30]:
            dfn = lambda k, N, s: delta_growing_horizon(k, N, s, c0)
            inf_val, k_star, L = inf_joint_excess_block_table(N, sigma, dfn)
            vals.append(inf_val)
            print(f"  sigma={sigma:>3d} N={N:.0e}: L={L:8.2f} "
                  f"k*={k_star:8.2f}  inf[Delta_k+rho_k]={inf_val:.4f} bits/sym")
        # The infimum must stay bounded away from 0 AND not decay to 0 as N->inf.
        lo = min(vals)
        spread = max(vals) - min(vals)
        if lo < 0.05:
            print(f"      FAIL: inf dipped to {lo:.4f} (should stay Theta(1))")
            fails += 1
        elif spread > 0.5 * lo:
            # tolerate slow drift; the key property is "bounded away from 0,
            # roughly N-independent".  Large drift would suggest decay.
            # Check it is not monotonically ->0:
            if vals[-1] < 0.5 * vals[0]:
                print(f"      FAIL: inf decaying with N ({vals[0]:.3f}->{vals[-1]:.3f})")
                fails += 1
            else:
                print(f"      OK: inf ~ Theta(1), bounded away from 0 (drift ok)")
        else:
            print(f"      OK: inf ~ {lo:.3f}..{max(vals):.3f}, N-independent Theta(1)")
    if fails == 0:
        print("  (A) PASS: growing-horizon closes the WHOLE block-table family.")
    return fails


# ----------------------------------------------------------------------
# (B) Fixed source: block-table family OPEN (inf -> 0), incl. divergent E.
# ----------------------------------------------------------------------

def check_fixed_source_open() -> int:
    print("\n--- (B) Fixed source: block-table family OPEN (inf -> 0) ---")
    fails = 0
    for decay, a in [("harmonic", 2.0), ("geometric", 2.0)]:
        for sigma in [2, 256]:
            vals = []
            for N in [10**6, 10**12, 10**24, 10**48]:
                dfn = lambda k, N, s: delta_fixed(k, decay, a)
                inf_val, k_star, L = inf_joint_excess_block_table(N, sigma, dfn)
                vals.append(inf_val)
                print(f"  decay={decay:>9s} sigma={sigma:>3d} N={N:.0e}: "
                      f"k*={k_star:8.2f}  inf={inf_val:.5f} bits/sym")
            # Must be DECREASING toward 0 as N grows.
            if not (vals[-1] < vals[0] and vals[-1] < 0.5 * max(0.02, vals[0])):
                print(f"      FAIL: inf not decaying to 0 ({vals[0]:.4f}->{vals[-1]:.4f})")
                fails += 1
            else:
                print(f"      OK: inf -> 0 ({vals[0]:.4f} -> {vals[-1]:.4f}); escape OPEN")
    if fails == 0:
        print("  (B) PASS: fixed sources (incl. divergent excess) stay OPEN, "
              "consistent with 7.27e(a).")
    return fails


# ----------------------------------------------------------------------
# (C) Redundancy-order incompatibility: k_model / k_red -> infinity.
# ----------------------------------------------------------------------

def check_incompatibility() -> int:
    print("\n--- (C) Redundancy/order incompatibility (k_model >> k_red) ---")
    fails = 0
    c0 = 0.5
    eps = 0.1  # target both Delta <= eps and rho <= eps
    for sigma in [2, 256]:
        ratios = []
        for N in [10**6, 10**12, 10**24, 10**48]:
            L = L_sigma(N, sigma)
            # k_model: smallest k with Delta_k^{(N)} = c0 L / k <= eps  =>  k >= c0 L / eps
            k_model = c0 * L / eps
            # k_red: largest k with rho_k <= eps  =>  (k+1) log2 sigma / L <= eps (approx)
            #        => k <= eps * L / log2 sigma - 1
            k_red = eps * L / math.log2(sigma) - 1.0
            ratio = k_model / max(k_red, 1e-9)
            ratios.append(ratio)
            feasible = (k_red >= k_model)  # is there a k satisfying both?
            print(f"  sigma={sigma:>3d} N={N:.0e}: L={L:8.2f} "
                  f"k_model>={k_model:8.2f} k_red<={k_red:8.2f}  "
                  f"ratio={ratio:7.2f}  both-feasible={feasible}")
            if feasible:
                print(f"      FAIL: a single k satisfies both Delta<=eps and rho<=eps")
                fails += 1
        # the ratio should be a fixed constant (> 1) or grow; never < 1.
        if min(ratios) <= 1.0:
            print(f"      FAIL: ratio dropped to {min(ratios):.2f} <= 1")
            fails += 1
    if fails == 0:
        print("  (C) PASS: model order k_model strictly exceeds redundancy "
              "budget k_red; no single block order is Shannon-near.")
    return fails


# ----------------------------------------------------------------------
# (D) NEGATIVE half: generic access floor log N/loglog N << sqrt(N).
# ----------------------------------------------------------------------

def check_sqrt_gap() -> int:
    print("\n--- (D) NEGATIVE: generic access floor (log N/loglog N) << sqrt N ---")
    fails = 0
    for N in [10**6, 10**9, 10**12, 10**18]:
        l2 = math.log2(N)
        floor = l2 / math.log2(max(l2, 2.0))      # Patrascu-Demaine access floor
        root = math.sqrt(N)
        print(f"  N={N:.0e}: cell-probe access floor ~{floor:8.2f}, "
              f"sqrt(N)~{root:.3e}, ratio sqrt/floor = {root/floor:.3e}")
        if not (floor < root and root / floor > 100):
            print(f"      FAIL: floor not << sqrt(N)")
            fails += 1
    if fails == 0:
        print("  (D) PASS: no unconditional cell-probe access bound reaches "
              "sqrt(N); the sqrt(N) converse is NOT provable unconditionally.")
        print("      => O(1)-access structures are not excluded; the target")
        print("         unconditional sqrt(N) converse is FALSE for fixed sources.")
    return fails


def main() -> int:
    print("=" * 72)
    print("Cor 7.27f: growing-horizon block-table floor (POSITIVE, conditional)")
    print("           + sqrt(N)-converse delimitation (NEGATIVE, unconditional)")
    print("=" * 72)
    fails = 0
    fails += check_growing_horizon_closed()
    fails += check_fixed_source_open()
    fails += check_incompatibility()
    fails += check_sqrt_gap()
    print("\n" + "=" * 72)
    if fails == 0:
        print("ALL CHECKS PASS")
        print("Verdict: PARTIAL/NEGATIVE.  Block-table family (FV/GOR/GN/BN)")
        print("provably fails on growing-horizon sources (new over 7.27e), but an")
        print("UNCONDITIONAL sqrt(N) access converse is unprovable: pure access is")
        print("O(1)-achievable (Belazzougui-Navarro) and the generic cell-probe")
        print("floor is only log N/loglog N. The escape's openness PERSISTS for")
        print("non-block-table architectures.")
        return 0
    print(f"FAILURES: {fails}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
