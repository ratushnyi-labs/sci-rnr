#!/usr/bin/env python3
"""
lemma_7_34b_polymer_gas.py
============================================================================
Lemma 7.34b'' (the ELLIPTIC DOMAIN-WALL POLYMER GAS).  The quenched shell
polynomial of the source around a word x,
    S_x(eta) := E_{X'}[eta^{d_H(x,X')}] / P_X(x) = sum_e eta^{|e|} P(x XOR e)/P(x),
is the grand partition function of a 1-D HARD-CORE polymer gas on x:
  - polymers = maximal runs ("flip-blocks") I=[a,b] of flipped sites,
  - activity  z_I(eta) = eta^{|I|} * prod_{edge bonds j of I} rho^{+-1},
       rho=p/(1-p); exponent +1 if the bond of x at the block edge is a STAY,
       -1 if a SWITCH; word-boundary edges contribute nothing,
  - hard core: distinct blocks separated by >= 1 unflipped site.
This holds because -log P_X is affine in the switch count and flipping a block
toggles exactly its two edge bonds.  (Via Lemma 7.34b', the posterior char.
function of the radius K=d_H(x,Y*) is Phi(s;x)=(C_s/C_0)^n S_x(eta_s).)

  V1  block transfer-matrix DP == brute-force 2^n sum, multiple words
      (random source-typical + pathological all-zeros / alternating), to 1e-9.
"""
import math
import numpy as np
from itertools import product

def S_brute(x, n, p, eta):
    def sw(z): return int(np.sum(z[1:] != z[:-1]))
    s0 = sw(x); r = p/(1-p); tot = 0.0
    for e in product([0,1], repeat=n):
        e = np.array(e)
        tot += (eta**int(e.sum())) * r**(sw(x ^ e) - s0)
    return tot

def S_blocks(x, n, p, eta):
    ratio = p/(1-p)
    def act(a, b):
        z = eta**(b-a+1)
        if a >= 1:      z *= ratio if (x[a]==x[a-1]) else 1.0/ratio
        if b+1 <= n-1:  z *= ratio if (x[b+1]==x[b]) else 1.0/ratio
        return z
    memo = [None]*(n+2); memo[n] = 1.0
    if n+1 <= n+1: memo[n+1] = 1.0
    for i in range(n-1, -1, -1):
        tot = memo[i+1]
        for b in range(i, n):
            nxt = memo[b+2] if b+2 <= n else 1.0
            tot += act(i, b)*nxt
        memo[i] = tot
    return memo[0]

if __name__ == "__main__":
    print("="*78)
    print("Lemma 7.34b'': S_x(eta) = hard-core domain-wall polymer gas (block DP == brute)")
    print("="*78)
    rng = np.random.default_rng(1)
    ok = True
    for p in (0.25, 0.4):
        n = 10
        words = []
        for _ in range(3):
            words.append((np.cumsum(rng.random(n) < p) % 2).astype(int))
        words.append(np.zeros(n, dtype=int))                       # pathological
        words.append(np.array([i % 2 for i in range(n)]))          # alternating
        worst = 0.0
        for x in words:
            for eta in (0.05, 0.1, -0.1, 0.2, -0.2, 0.3, 0.0):
                a = S_brute(x, n, p, eta); b = S_blocks(x, n, p, eta)
                d = abs(a-b)/max(abs(a), 1e-300); worst = max(worst, d)
                ok &= d < 1e-9
        print(f"  p={p}, n={n}: 5 words (incl. all-zeros, alternating) x 7 eta values; "
              f"worst rel.diff={worst:.1e}")
    print("="*78)
    print(f"RESULT: {'ALL PASS' if ok else 'FAIL'} -- the polymer (flip-block) representation is exact:")
    print("  activities eta^|I| * rho^{+-1 per edge bond}, hard-core separation, on the")
    print("  ELLIPTIC source word. The quenched lattice object lives on this gas.")
