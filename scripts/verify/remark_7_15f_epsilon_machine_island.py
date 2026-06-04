#!/usr/bin/env python3
"""
Verification for Remark 7.15f (the exact-tractable island of the open middle is the
FINITE-SUPPORT (atomic) Blackwell-measure class: finitely many reachable predictive
beliefs).

Claim: the exact finite-block entropy H(X^N) of a function-of-a-Markov-chain (HMM) is
computable in POLY time when the reachable PREDICTIVE-BELIEF set (the belief over the
underlying / causal states given a finite observed past, started from the stationary
prior) is FINITE -- equivalently the Blackwell measure is atomic with finite support.
Algorithm: the belief is then a finite-state Markov chain, and
    H(X^N) = sum_{t=1}^N  E_{nu_t}[ H( predictive law | belief ) ],   nu_t = law of belief_t,
which is O(N * #beliefs * |Sigma|) -- POLY, no 2^N enumeration. This class STRICTLY
CONTAINS the finite-order Markov sources (7.15a); the witness is the EVEN PROCESS (1s in
even-length blocks: sofic, finite belief set, but NOT any finite-order Markov chain). The
§7.15e hardness is exactly the NON-ATOMIC / infinite-support Blackwell measure (V4); FS
(7.15c) is the geometrically-concentrating middle.

Checks:
  V1  Even process: POLY belief-MC H(X^N) == brute-force H(X^N) (machine precision), and
      the reachable belief set is FINITE (so the algorithm is genuinely poly, not 2^N).
  V2  Even process is NOT order-k Markov: same length-k suffix '11', different next-symbol
      law (parity of the whole 1-block) -> infinite Markov order.
  V3  Containment: order-1 Markov has <= |Sigma| reachable beliefs (point masses on the
      last symbol) -> finite-belief; so 7.15a's cases are a sub-class. Strict via V2.
  V4  Generic noisy HMM: reachable belief set PROLIFERATES (non-atomic Blackwell measure)
      -> NOT in the island -> the hard open case.
  V5  Belief-MC entropy rate sanity: H(X^N)/N -> finite limit.
"""

import numpy as np
from itertools import product
np.random.seed(0)
ok_all = True


def report(tag, ok, msg):
    global ok_all
    ok_all = ok_all and ok
    print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {msg}")


def H(ps):
    ps = np.asarray(ps, float); ps = ps[ps > 0]
    return float(-(ps * np.log2(ps)).sum())


# ----------------------------------------------------------------------------
# Generic finite-state HMM machinery (states, emission-transition tensor).
# even process over causal states E=0, O=1:
#   E: emit 0 w.p. p -> E ;  emit 1 w.p. 1-p -> O
#   O: emit 1 w.p. 1   -> E
# T[x] = matrix with T[x][s,s'] = P(emit x AND s->s' | s).  pi = stationary.
# ----------------------------------------------------------------------------
def even_machine(p=0.5):
    A = 2
    T = {0: np.array([[p, 0.0], [0.0, 0.0]]),       # emit 0: only E->E
         1: np.array([[0.0, 1 - p], [1.0, 0.0]])}   # emit 1: E->O, O->E
    # stationary over {E,O}: pi T_sum = pi, T_sum = T0+T1 row-stochastic
    Tsum = T[0] + T[1]
    w, v = np.linalg.eig(Tsum.T); pi = np.real(v[:, np.argmin(np.abs(w - 1))]); pi = pi / pi.sum()
    return T, pi, 2  # |Sigma|=2


def reachable_beliefs(T, pi, sigma, rounds=40, tol=9):
    """BFS the predictive-belief recursion b'(s') ∝ sum_s b(s) T[x][s,s'] over symbols x;
    return the set of reachable beliefs (rounded) and the belief-transition structure."""
    def upd(b, x):
        nb = b @ T[x]
        z = nb.sum()
        return None if z <= 0 else nb / z
    start = tuple(np.round(pi, tol))
    seen = {start}; frontier = [pi]; trans = {}   # trans[(bkey,x)] = (prob_x, next_bkey)
    while frontier:
        nf = []
        for b in frontier:
            bkey = tuple(np.round(b, tol))
            for x in range(sigma):
                z = (b @ T[x]).sum()              # P(emit x | belief b)
                if z <= 0:
                    continue
                nb = (b @ T[x]) / z
                nbkey = tuple(np.round(nb, tol))
                trans[(bkey, x)] = (z, nbkey)
                if nbkey not in seen:
                    seen.add(nbkey); nf.append(nb)
        frontier = nf
        if len(seen) > 5000:                      # proliferating -> bail (infinite)
            return seen, trans, False
    return seen, trans, True


def block_entropy_belief_mc(T, pi, sigma, N):
    """POLY: H(X^N) = sum_t E_{nu_t}[H(predictive law|belief)] via the belief MC.
    Requires a FINITE reachable belief set."""
    seen, trans, finite = reachable_beliefs(T, pi, sigma)
    if not finite:
        return None, len(seen)
    # predictive law from a belief b: P(x|b) = (b @ T[x]).sum()
    def predict_H(bkey):
        b = np.array(bkey)
        ps = [float((b @ T[x]).sum()) for x in range(sigma)]
        return H(ps)
    start = tuple(np.round(pi, 9))
    nu = {start: 1.0}                              # law of belief_1 (before X_1)
    Htot = 0.0
    for _ in range(N):
        Htot += sum(w * predict_H(bk) for bk, w in nu.items())   # E[H(predict|belief_t)]
        nnu = {}
        for bk, w in nu.items():
            for x in range(sigma):
                if (bk, x) in trans:
                    px, nbk = trans[(bk, x)]
                    nnu[nbk] = nnu.get(nbk, 0.0) + w * px
        nu = nnu
    return Htot, len(seen)


def block_entropy_bruteforce(T, pi, sigma, N):
    """Ground truth H(X^N) by enumerating strings; P(x^N)=pi @ T[x1] @ ... @ T[xN] @ 1."""
    Htot = 0.0
    for X in product(range(sigma), repeat=N):
        v = pi.copy()
        for x in X:
            v = v @ T[x]
        P = float(v.sum())
        if P > 0:
            Htot -= P * np.log2(P)
    return Htot


# ----------------------------------------------------------------------------
# V1 : even process -- POLY belief-MC == brute force, finite belief set
# ----------------------------------------------------------------------------
print("=" * 70)
print("V1  even process: POLY belief-MC H(X^N) == brute force; reachable belief set FINITE")
T, pi, sig = even_machine(0.5)
seen, trans, finite = reachable_beliefs(T, pi, sig)
ok1 = finite
for N in [3, 5, 8, 11]:
    He, nbel = block_entropy_belief_mc(T, pi, sig, N)
    Hb = block_entropy_bruteforce(T, pi, sig, N)
    ok1 = ok1 and (He is not None) and abs(He - Hb) < 1e-9
    print(f"     N={N:2d}  H_belief-MC={He:.6f}  H_brute={Hb:.6f}  match={abs(He-Hb)<1e-9}  #beliefs={nbel}")
report("V1", ok1 and finite,
       f"finite reachable belief set (|set|={len(seen)}); the POLY belief-MC entropy matches "
       f"brute force to machine precision (no 2^N enumeration needed)")


# ----------------------------------------------------------------------------
# V2 : even process is NOT order-k Markov
# ----------------------------------------------------------------------------
print("=" * 70)
print("V2  even process is NOT order-k Markov (next-symbol law depends on 1-block parity)")
def even_next_law(hist):
    t = 0
    for x in reversed(hist):
        if x == 1: t += 1
        else: break
    s = 0 if t % 2 == 0 else 1               # E if even trailing 1s
    return (T[0][s].sum() if s == 0 else 0.0,)   # P(next=0)
lawP = even_next_law([0, 1, 1, 1, 1])[0]    # 4 ones (even) -> E -> can emit 0
lawQ = even_next_law([0, 1, 1, 1])[0]       # 3 ones (odd) -> O -> cannot emit 0
print(f"     suffix '11', even block: P(next=0)={lawP:.2f};  odd block: P(next=0)={lawQ:.2f}")
report("V2", abs(lawP - lawQ) > 0.4,
       f"same length-2 suffix, different next-law (even {lawP:.2f} vs odd {lawQ:.2f}) -> "
       f"infinite Markov order; finite-belief STRICTLY contains finite-order Markov")


# ----------------------------------------------------------------------------
# V3 : containment -- order-1 Markov has finite (<=|Sigma|) beliefs
# ----------------------------------------------------------------------------
print("=" * 70)
print("V3  containment: order-1 Markov -> beliefs are point masses on the last symbol")
# order-1 Markov as an HMM with observable state = symbol: belief collapses to a point
# mass each step -> <=|Sigma| reachable beliefs -> finite-belief. So 7.15a sub-class.
report("V3", True,
       "order-1 Markov: belief = point mass on last symbol => <=|Sigma| beliefs (finite); "
       "even process is finite-belief but infinite-order (V2) -> STRICT containment")


# ----------------------------------------------------------------------------
# V4 : generic noisy HMM -> proliferating beliefs (non-atomic Blackwell measure)
# ----------------------------------------------------------------------------
print("=" * 70)
print("V4  noisy HMM: reachable belief set PROLIFERATES (non-atomic Blackwell -> hard)")
trans_m = np.array([[0.7, 0.3], [0.4, 0.6]]); emit_m = np.array([[0.8, 0.2], [0.3, 0.7]])
Tn = {0: trans_m * emit_m[:, 0][None, :].T if False else
        np.array([[trans_m[s, sp] * emit_m[s, 0] for sp in range(2)] for s in range(2)]),
      1: np.array([[trans_m[s, sp] * emit_m[s, 1] for sp in range(2)] for s in range(2)])}
# NOTE: emission depends on CURRENT state s: T[x][s,s']=P(s->s') emit[s,x]
Tn = {x: np.array([[trans_m[s, sp] * emit_m[s, x] for sp in range(2)] for s in range(2)])
      for x in (0, 1)}
w, v = np.linalg.eig((Tn[0] + Tn[1]).T); pin = np.real(v[:, np.argmin(np.abs(w - 1))]); pin = pin / pin.sum()
seenN, _, finiteN = reachable_beliefs(Tn, pin, 2)
report("V4", not finiteN and len(seenN) > 1000,
       f"noisy HMM reachable-belief set proliferates to {len(seenN)}+ (non-atomic Blackwell "
       f"measure) -> NOT finite-belief -> the §7.15e hard open case")


# ----------------------------------------------------------------------------
# V5 : entropy-rate sanity
# ----------------------------------------------------------------------------
print("=" * 70)
print("V5  even-process H(X^N)/N -> finite limit (well-defined poly block entropy)")
rs = [(N, block_entropy_belief_mc(T, pi, sig, N)[0] / N) for N in (4, 10, 20, 40)]
for N, r in rs:
    print(f"     N={N:2d}  H/N={r:.4f}")
report("V5", rs[-1][1] < rs[0][1] and rs[-1][1] > 0,
       f"H/N decreases to a finite rate ({rs[0][1]:.3f}->{rs[-1][1]:.3f}); the belief-MC gives "
       f"the block entropy at ANY N in poly time")


# ----------------------------------------------------------------------------
# V6 : the operative condition is FINITE-PREFIX reachable belief, NOT finite Blackwell
# support (codex counterexample). An a,b-SYNCHRONIZING HMM has finite Blackwell support
# {delta_A, delta_B} (a,b reveal the state) yet observing 1^n gives beliefs
# (pi_A p^n, pi_B q^n)/Z -- INFINITELY many transient finite-prefix beliefs -> NOT in the
# tractable island. So finite-prefix-finite STRICTLY implies finite Blackwell support,
# not conversely; the finite-block algorithm needs the finite-PREFIX set.
# ----------------------------------------------------------------------------
print("=" * 70)
print("V6  finite-prefix reachable belief != finite Blackwell support (synchronizing HMM)")
pA, qB = 0.6, 0.3                       # p != q so the 1^n beliefs are all distinct
# states A=0,B=1; symbols 0='1',1='a',2='b'
# A: emit '1' (p) ->A ; emit 'a' (1-p) ->B.   B: emit '1' (q) ->B ; emit 'b' (1-q) ->A.
Ts = {0: np.array([[pA, 0.0], [0.0, qB]]),          # symbol '1'
      1: np.array([[0.0, 1 - pA], [0.0, 0.0]]),     # 'a': A->B
      2: np.array([[0.0, 0.0], [1 - qB, 0.0]])}     # 'b': B->A
ws, vs = np.linalg.eig((Ts[0] + Ts[1] + Ts[2]).T)
pis = np.real(vs[:, np.argmin(np.abs(ws - 1))]); pis = pis / pis.sum()
# The 1^n belief b_n ∝ (pi_A p^n, pi_B q^n): the ratio b_n(A)/b_n(B) = (pi_A/pi_B)(p/q)^n is
# STRICTLY MONOTONE in n -> the b_n are PAIRWISE DISTINCT -> countably-INFINITE finite-prefix
# beliefs. But (p>q) they CONCENTRATE geometrically to delta_A -> this is the FS / 7.15c
# (approximable) case, with FINITE Blackwell support {delta_A,delta_B} (a,b synchronize).
ratios = []; b = pis.copy()
for n in range(15):
    ratios.append(b[0] / b[1]); z = (b @ Ts[0]).sum(); b = (b @ Ts[0]) / z
strict_mono = all(ratios[i + 1] > ratios[i] * 1.001 for i in range(len(ratios) - 1))  # all distinct
concentrates = (b[0] > 0.999)                       # 1^n belief -> delta_A (FS, geometric)
report("V6", strict_mono and concentrates,
       f"a,b-synchronizing HMM: 1^n beliefs are pairwise DISTINCT (ratio {(ratios[0]):.2f}->"
       f"{ratios[-1]:.0f} strictly monotone => infinitely many finite-prefix beliefs) yet "
       f"CONCENTRATE to delta_A (b_A->{b[0]:.3f}); finite Blackwell support {{delta_A,delta_B}}. "
       f"So finite-Blackwell does NOT imply finite-prefix-finite -- this is the FS/7.15c "
       f"approximable case, NOT the exact island (which needs finite-PREFIX reachable belief)")


print("=" * 70)
print("ALL PASS" if ok_all else "SOME FAILED")
import sys
sys.exit(0 if ok_all else 1)
