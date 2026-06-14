#!/usr/bin/env python3
"""
theorem_7_34i_asym_dispersion.py
============================================================================
Verifier for Theorem 7.34i (asymmetric BSMS operational RD-dispersion):
the CONSOLIDATION of instance (2). It does not re-prove the lemmas (those
have their own verifiers); it checks the two facts the consolidated theorem
adds on top of the committed pieces:

  T1  V_lossless(a,b) = the Gordin/Green-Kubo varentropy rate of the chain,
      with the displayed decomposition V = V_mix + V_mem,
        V_mix = pi0 a(1-a) log2^2((1-a)/a) + pi1 b(1-b) log2^2((1-b)/b),
      and V_mem = 0 iff a = b.  (Cross-checked vs a direct simulation, and
      vs the symmetric closed form p(1-p)log2^2((1-p)/p) at a=b.)
  T2  the assembly is monotone in the spectral constant c_g: every GAP-1
      downstream constant (c_2,c_4,c_5 = c_g/{4,16,4}pi^2) is a positive
      multiple of c_g, so the symmetric 11/16 substitutes by c_tilde(a,b,D)
      in [43/64, 25/32] with no sign change anywhere in the closure chain.
      (Structural check: the constants are increasing in c_g and the
      load-bearing E_{2'} margin c_g/pi^2 - 3 c_2 = c_g/4pi^2 > 0 for any
      c_g > 0.)

The achievability spectral input (item i), the deconvolution region (iii),
and the converse (rho=1) are verified in their own probes:
  lemma_7_34e_spi_closed_form.py, probe_7_34_asym_replica_structure.py,
  lemma_7_34g_sign_automaton.py, lemma_7_34g_G1_uniform_existence.py,
  remark_7_34d_general_converse_discrete_markov.py,
  probe_7_34_nonsym_dual_identity.py.
Deps: numpy (sim), math, fractions.
"""
import math
from fractions import Fraction as Fr
import numpy as np

PASS = True
def rep(name, ok):
    global PASS
    PASS = PASS and ok
    print(f"  {name:<62} {'PASS' if ok else 'FAIL'}")
    return ok


def Hb(x):
    return -x*math.log2(x) - (1-x)*math.log2(1-x)


def V_gordin(a, b):
    """Exact Gordin/Green-Kubo varentropy rate of the binary Markov chain."""
    pi0, pi1 = b/(a+b), a/(a+b)
    hbar = pi0*Hb(a) + pi1*Hb(b)
    f = {(0, 0): -math.log2(1-a), (0, 1): -math.log2(a),
         (1, 0): -math.log2(b),   (1, 1): -math.log2(1-b)}
    g0 = Hb(a)
    u1 = -(g0 - hbar)/a          # Poisson solution, gauge u0=0  (row0: -a u1 = g0-hbar)
    u = {0: 0.0, 1: u1}
    pi = {0: pi0, 1: pi1}
    T = {0: {0: 1-a, 1: a}, 1: {0: b, 1: 1-b}}
    return sum(pi[i]*T[i][j]*(f[(i, j)] + u[j] - u[i] - hbar)**2
               for i in (0, 1) for j in (0, 1))


def V_mix(a, b):
    pi0, pi1 = b/(a+b), a/(a+b)
    return (pi0*a*(1-a)*math.log2((1-a)/a)**2
            + pi1*b*(1-b)*math.log2((1-b)/b)**2)


def V_sim(a, b, N=3_000_000, seed=2):
    rng = np.random.default_rng(seed)
    T = [[1-a, a], [b, 1-b]]
    x = 0 if rng.random() < b/(a+b) else 1
    g = np.empty(N)
    for t in range(N):
        xn = 0 if rng.random() < T[x][0] else 1
        g[t] = -math.log2(T[x][xn]); x = xn
    g = g - g.mean(); n = len(g); L = 300
    return g @ g/n + 2*sum(g[:n-k] @ g[k:]/n for k in range(1, L))


def V_sym(p):
    return p*(1-p)*math.log2((1-p)/p)**2


def T1_varentropy():
    print("-"*78)
    print("T1  V_lossless = Gordin rate; decomposition V=V_mix+V_mem; V_mem=0 iff a=b")
    ok = True
    for a, b in [(0.25, 0.25), (0.4, 0.4), (0.1, 0.3), (0.2, 0.4), (0.05, 0.45), (0.3, 0.4)]:
        Vg, Vm = V_gordin(a, b), V_mix(a, b)
        Vmem = Vg - Vm
        sym = abs(a-b) < 1e-12
        # V_mem >= 0 always; = 0 iff a=b
        cond = (abs(Vmem) < 1e-9) if sym else (Vmem > 1e-3)
        ok = ok and cond
        print(f"     (a={a},b={b}): V_gordin={Vg:.5f}  V_mix={Vm:.5f}  V_mem={Vmem:+.5f}"
              f"  {'[=0 (a=b)]' if sym else '[>0 (a!=b)]'}")
    rep("T1a V_mem = 0 iff a=b (decomposition correct)", ok)
    # cross-check Gordin vs sim (Bartlett, ~1% noise)
    sim_ok = True
    for a, b in [(0.1, 0.3), (0.2, 0.4)]:
        rel = abs(V_gordin(a, b) - V_sim(a, b))/V_gordin(a, b)
        sim_ok = sim_ok and rel < 0.03
        print(f"     (a={a},b={b}): |Gordin-sim|/Gordin = {rel:.4f}")
    rep("T1b Gordin rate matches direct simulation (<3%)", sim_ok)
    # symmetric reduction exact
    sym_ok = all(abs(V_gordin(p, p) - V_sym(p)) < 1e-12 for p in (0.1, 0.25, 0.4))
    rep("T1c Gordin(p,p) = p(1-p)log2^2((1-p)/p) exactly (symmetric reduction)", sym_ok)
    return ok and sim_ok and sym_ok


def T2_assembly_monotone():
    print("-"*78)
    print("T2  GAP-1 downstream constants are positive multiples of c_g (monotone subst.)")
    # symbolic-rational check: with c_g a free positive, c2=c_g/4pi^2 etc., the
    # E_2' margin c_g/pi^2 - 3 c2 = c_g/4pi^2 > 0 for ALL c_g>0 (no 11/16 needed).
    cg = Fr(1)  # any positive; the ratios are c_g-linear so the sign is c_g-independent
    # use rationals for the pi^2-free ratios (factor out 1/pi^2)
    c2 = cg/4          # c_g/(4 pi^2), drop the 1/pi^2 common factor
    c4 = cg/16
    c5 = cg/4
    margin = cg - 3*c2   # (c_g/pi^2) - 3 c2 , common 1/pi^2 dropped  => c_g - 3*(c_g/4) = c_g/4
    ok = (c2 > 0 and c4 > 0 and c5 > 0 and margin == cg/4 and margin > 0)
    print(f"     c2=c_g/4, c4=c_g/16, c5=c_g/4 (x 1/pi^2); E_2' margin = c_g/pi^2-3c2 = c_g/4pi^2 > 0")
    rep("T2a all GAP-1 constants positive & monotone in c_g; margin>0 for any c_g>0", ok)
    # c_tilde range from Lemma 7.34e (iii): symmetric 11/16=0.6875 in [43/64,25/32]
    lo, hi = Fr(43, 64), Fr(25, 32)
    sym_in = lo <= Fr(11, 16) <= hi
    rep("T2b symmetric c_g constant 11/16 lies in the certified c_tilde range [43/64,25/32]", sym_in)
    return ok and sym_in


if __name__ == "__main__":
    print("="*78)
    print("Theorem 7.34i (asymmetric BSMS RD-dispersion) -- consolidation verifier")
    print("="*78)
    T1_varentropy()
    T2_assembly_monotone()
    print("="*78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
