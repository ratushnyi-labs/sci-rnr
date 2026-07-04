#!/usr/bin/env python3
r"""
lemma_7_34_a5_negside_certificates.py
============================================================================
BUG-009-D / Route B mid-band program, A=5 REPLICATION: the FIVE
NEGATIVE-SIDE certificates CERTIFIED for A=5 on the entire rationalized
Gray superdomain -- all at the root box, exact rational arithmetic,
closed-face towers included.

TARGETS: with psi(x) := -chi(-x) (monic; roots = -lambda_i(Q)), certify
    psi^(k)(Lambda) > 0   for k = 0,1,2,3,4,
equivalently (-1)^{k+1} chi^(k)(-Lambda) > 0, at Lambda = L/Cr,
L = 1 - (3/2) D(1-D) t.  Since psi's coefficients are x^5 + e1 x^4 +
e2 x^3 + e3 x^2 + e4 x + e5 (the e_k of Q with ALL PLUS signs), the
t-domain small-piece assembly of lemma_7_34_a5_chi_certificates.py applies
verbatim with the (-1)^k factors dropped.

WHY THESE ARE NEEDED: Budan-Fourier at +Lambda (certs C1-C4 of the chi
script) bounds only real eigenvalues ABOVE +Lambda.  The spectral radius
also needs no real eigenvalue BELOW -Lambda.  These five certificates
close that side: zero sign variations of (psi(Lambda), ..., psi''''(Lambda),
1)  =>  psi has no real root in [Lambda, inf)  =>  Q has no real
eigenvalue in (-inf, -Lambda].

EXACT SIGN BOOKKEEPING: the shared P3 gates (chat = p(p-1)B, CrN = -B,
CrD = -(5D-4)^2, chat*CrN = p(1-p)B^2, B > 0 on the claimed domain, with
B = (5D-4)^2 + 2tD^2(1-D)(4-D)) make sign(target) == sign(psi^(k)(Lambda))
EXACT; the dps-80 numeric chains are cross-checks and no sign flip is
tolerated (the A=4 upgrade over the committed A=3 negside script, kept
verbatim).

RESULT: all five certify AT THE ROOT BOX on sigma,theta in [0,1], t in
[0,2] (the full domain: all p in (0,4/5), all D in (0,Dbar], all angles),
and all face towers close (numbers from the recorded run of this script):
  N4 (k=4): target 367 terms,   strip (1-sg)^1, core 360,   min ~ 5.68
  N3 (k=3): target 2268 terms,  strip (1-sg)^2, core 2146,  min ~ 3.52e-2
  N2 (k=2): target 6784 terms,  strip (1-sg)^3, core 6496,  min ~ 1.36e-4
  N1 (k=1): target 15075 terms, strip (1-sg)^4, core 14283, min ~ 3.51e-7
  N0 (k=0): target 28281 terms, strip (1-sg)^5, core 26976, min ~ 6.77e-10
(zero-face grading (1-sg)^k on the negative side, exactly as at A=3/A=4;
term counts for k <= 3 are IDENTICAL to the A=3/A=4 negside counts -- the
assembly's monomial support is A-independent there -- and the k=4 count
367 matches A=4.)  The minima are ~15-30x HEALTHIER than the A=4 run
(N0: 6.77e-10 vs 2.32e-11; N1: 3.51e-7 vs 2.36e-8).  NOTE: as at A=3/A=4,
k=0 has NO pocket -- the negative-side bound holds on the FULL superdomain
including beyond-Gray, in contrast to the positive-side cert5.

COMBINED SPECTRAL STATE for A=5 after this lemma (see the chi script's
header for the positive side):
  * real eigenvalues > +Lambda: chi', ..., chi'''' > 0 certified
    full-domain (chi script); chi(Lambda) > 0 full-domain is FALSE (A=5
    real-Perron pocket, scout witness) -- cert5-on-Gray via the TWO-PIECE
    restricted-curve splice (curveM = 1 - (3/4)sg^2 on (0,3/4] +
    curveP = (9-5sg^2)/10 on [7/10,1), companion script
    lemma_7_34_a5_cert5_on_gray.py; a single quadratic curve is DEAD at
    A=5, scout G4 corridor);
  * real eigenvalues < -Lambda: EXCLUDED FULL-DOMAIN (this lemma);
  * realness on the superdomain: FALSE at A=5 (new vs A=3/A=4) -- the
    s = pi complexification island lies INSIDE the strip (W < 0 for
    sg ~ [0.58,0.77], bottom theta ~ 0.7275 at sg = 0.64, t = 2, scouts
    2b/4b), so the eventual realness certificate must be CURVE-RESTRICTED
    (the island clears curveM by ~+0.024, curveP by ~+0.028).
  => once cert5-on-Gray + curve-restricted interior-angle realness are
    certified, rho(Q) < Lambda on Gray follows, i.e. R_5(s) >= 3/2 on
    all of Gray.

CHECKS (default run, ~6 min; shared P* gates rerun here, then one block
per k): dps-80 numeric sign chain (90 pts, no flip allowed) + dps-80 value
chain against eigenvalue-built psi^(k)(Lambda); zero-face strip; EXACT
root-box Bernstein (Fraction tensor engine); face towers theta=1 / t=2 +
edge Sturm, closing {sigma in (0,1)} x {theta in (0,1]} x {t in (0,2]}.

Deps: sympy, mpmath; imports the pipeline from
lemma_7_34_a5_chi_certificates.py (same directory).
Python: /Users/para/.venvs/rnr/bin/python.
"""
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lemma_7_34_a5_chi_certificates import (   # noqa: E402
    build_pieces, run_target)

mp.mp.dps = 40
PASS = True


def rep(name, ok):
    global PASS
    PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    t00 = time.time()
    print("=" * 78)
    print("A=5 negative-side certificates psi^(k)(Lambda) > 0, k=0..4, full")
    print("rationalized Gray superdomain -- exact root-box Bernstein positivity")
    print("+ closed-face towers")
    print("=" * 78)
    PC = build_pieces(rep, "N*")
    if PC is None:
        print("=" * 78)
        print("OVERALL -> FAIL")
        return
    for j in range(4, -1, -1):
        print("-" * 78)
        print(f"N{j}  psi^({j})(Lambda) > 0")
        run_target(j, PC, neg=True, name=f"N{j}", repfn=rep, seed=29)
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}   "
          f"({time.time() - t00:.0f}s)")
    print("psi^(k)(Lambda) > 0, k=0..4, on the whole A=5 superdomain => Q has")
    print("NO real eigenvalue <= -Lambda anywhere (negative side closed, full")
    print("domain, no pocket).  With cert5-on-Gray (two-piece splice) +")
    print("curve-restricted interior-angle realness this yields rho(Q) <")
    print("Lambda on Gray, i.e. R_5(s) >= 3/2.")


if __name__ == "__main__":
    main()
