#!/usr/bin/env python3
r"""
lemma_7_34_a4_negside_certificates.py
============================================================================
BUG-009-D / Route B mid-band program, A=4 REPLICATION: the FIVE
NEGATIVE-SIDE certificates CERTIFIED for A=4 on the entire rationalized
Gray superdomain -- all at the root box, exact rational arithmetic,
closed-face towers included.

TARGETS: with psi(x) := -chi(-x) (monic; roots = -lambda_i(Q)), certify
    psi^(k)(Lambda) > 0   for k = 0,1,2,3,4,
equivalently (-1)^{k+1} chi^(k)(-Lambda) > 0, at Lambda = L/Cr,
L = 1 - (3/2) D(1-D) t.  Since psi's coefficients are x^5 + e1 x^4 +
e2 x^3 + e3 x^2 + e4 x + e5 (the e_k of Q with ALL PLUS signs), the
t-domain small-piece assembly of lemma_7_34_a4_chi_certificates.py applies
verbatim with the (-1)^k factors dropped.

WHY THESE ARE NEEDED: Budan-Fourier at +Lambda (certs C1-C4 of the chi
script) bounds only real eigenvalues ABOVE +Lambda.  The spectral radius
also needs no real eigenvalue BELOW -Lambda.  These five certificates
close that side: zero sign variations of (psi(Lambda), ..., psi''''(Lambda),
1)  =>  psi has no real root in [Lambda, inf)  =>  Q has no real
eigenvalue in (-inf, -Lambda].

EXACT SIGN BOOKKEEPING: the shared P3 gates (chat = p(p-1)B, CrN = -B,
CrD = -(4D-3)^2, chat*CrN = p(1-p)B^2, B > 0 on the claimed domain) make
sign(target) == sign(psi^(k)(Lambda)) EXACT; the dps-80 numeric chains are
cross-checks and no sign flip is tolerated (upgrade over the committed A=3
negside script, which fixed the global sign numerically with an auto-flip).

RESULT: all five certify AT THE ROOT BOX on sigma,theta in [0,1], t in
[0,2] (the full domain: all p in (0,3/4), all D in (0,Dbar], all angles),
and all face towers close (numbers from the recorded run of this script):
  N4 (k=4): target 367 terms,   strip (1-sg)^1, core 360,   min ~ 2.38
  N3 (k=3): target 2268 terms,  strip (1-sg)^2, core 2146,  min ~ 9.13e-3
  N2 (k=2): target 6784 terms,  strip (1-sg)^3, core 6496,  min ~ 1.80e-5
  N1 (k=1): target 15075 terms, strip (1-sg)^4, core 14283, min ~ 2.36e-8
  N0 (k=0): target 28281 terms, strip (1-sg)^5, core 26976, min ~ 2.32e-11
(zero-face grading (1-sg)^k on the negative side, exactly as at A=3; term
counts for k <= 3 are IDENTICAL to the A=3 negside counts -- the assembly's
monomial support is A-independent there -- and k=4 grows 347 -> 367 with
E1's support.)  NOTE: as at A=3, k=0 has NO pocket -- the negative-side
bound holds on the FULL superdomain including beyond-Gray, in contrast to
the positive-side cert5.

COMBINED SPECTRAL STATE for A=4 after this lemma (see the chi script's
header for the positive side):
  * real eigenvalues > +Lambda: chi', ..., chi'''' > 0 certified
    full-domain (chi script); chi(Lambda) > 0 full-domain is FALSE (A=4
    real-Perron pocket, scout witness) -- cert5-on-Gray via a restricted
    curve (candidate thetatilde4 = 1 - sigma^2/2) is the open task;
  * real eigenvalues < -Lambda: EXCLUDED FULL-DOMAIN (this lemma);
  * realness on the superdomain: numeric-grade only (scout 4: 500-point
    scan, max |Im|/rho ~ 1.6e-28), and the A=4 collision margin is much
    tighter than A=3's (min collision-surface/Dbar ~ 1.005 at sigma ~ 0.7,
    t = 2, vs 1.35 at A=3) -- an interior-angle realness certificate at
    A=4 must confront this thin margin.
  => once cert5-on-Gray + interior-angle realness are certified,
    rho(Q) < Lambda on Gray follows, i.e. R_4(s) >= 3/2 on all of Gray.

CHECKS (default run, ~3 min; shared P* gates rerun here, then one block
per k): dps-80 numeric sign chain (90 pts, no flip allowed) + dps-80 value
chain against eigenvalue-built psi^(k)(Lambda); zero-face strip; EXACT
root-box Bernstein (Fraction tensor engine); face towers theta=1 / t=2 +
edge Sturm, closing {sigma in (0,1)} x {theta in (0,1]} x {t in (0,2]}.

Deps: sympy, mpmath; imports the pipeline from
lemma_7_34_a4_chi_certificates.py (same directory).
Python: /Users/para/.venvs/rnr/bin/python.
"""
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lemma_7_34_a4_chi_certificates import (   # noqa: E402
    build_pieces, run_target)

mp.mp.dps = 30
PASS = True


def rep(name, ok):
    global PASS
    PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    t00 = time.time()
    print("=" * 78)
    print("A=4 negative-side certificates psi^(k)(Lambda) > 0, k=0..4, full")
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
    print("psi^(k)(Lambda) > 0, k=0..4, on the whole A=4 superdomain => Q has")
    print("NO real eigenvalue <= -Lambda anywhere (negative side closed, full")
    print("domain, no pocket).  With cert5-on-Gray + interior-angle realness")
    print("this yields rho(Q) < Lambda on Gray, i.e. R_4(s) >= 3/2.")


if __name__ == "__main__":
    main()
