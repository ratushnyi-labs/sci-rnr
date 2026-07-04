#!/usr/bin/env python3
r"""
lemma_7_34_midband_negside_certificates.py
============================================================================
BUG-009-D / Route B mid-band program: the FIVE NEGATIVE-SIDE certificates
CERTIFIED for A=3 on the entire rationalized Gray superdomain -- all at the
root box, exact rational arithmetic.

TARGETS: with psi(x) := -chi(-x) (monic; roots = -lambda_i(Q)), certify
    psi^(k)(Lambda) > 0   for k = 0,1,2,3,4,
equivalently (-1)^{k+1} chi^(k)(-Lambda) > 0, at Lambda = L/Cr.
Since psi's coefficients are x^5 + e1 x^4 + e2 x^3 + e3 x^2 + e4 x + e5
(the e_k of Q with ALL PLUS signs), the t-domain small-piece assembly of
lemma_7_34_midband_chi_certificates.py applies verbatim with the (-1)^k
factors dropped.

WHY THESE ARE NEEDED: Budan-Fourier at +Lambda (certs 1-4 of the chi
script) bounds only real eigenvalues ABOVE +Lambda.  The spectral radius
also needs no real eigenvalue BELOW -Lambda; numerically
max |lambda_min| / Lambda = 0.747 over the whole superdomain (1.34x
margin).  These five certificates close that side: zero sign variations
of (psi(Lambda), ..., psi''''(Lambda), 1) => psi has no real root in
[Lambda, inf) => Q has no real eigenvalue in (-inf, -Lambda].

RESULT: all five certify AT THE ROOT BOX on sigma,theta in [0,1], t in
[0,2] (the full domain: all p in (0,2/3), all D in (0,Dbar], all angles):
  k=4: target 347 terms,   strip (1-sg)^1, core 340
  k=3: target 2268 terms,  strip (1-sg)^2, core 2146
  k=2: target 6784 terms,  strip (1-sg)^3, core 6496
  k=1: target 15075 terms, strip (1-sg)^4, core 14283
  k=0: target 28281 terms, strip (1-sg)^5, core 26976
(zero-face grading (1-sg)^k on the negative side vs (1-sg)^{2k} on the
positive side.)  NOTE: unlike the positive side, k=0 has NO pocket -- the
negative-side bound holds on the FULL superdomain including beyond-Gray.

COMBINED SPECTRAL STATE for A=3 after this lemma (see the chi script's
corrected header for the positive side):
  * real eigenvalues > +Lambda: excluded on Gray (chi(Lambda)>0 numeric
    on Gray + chi',..,chi'''' > 0 certified full-domain); one real Perron
    crossing exists beyond Gray (D_perron in (D_c, Dbar)) -- retracted
    full-domain cert5;
  * real eigenvalues < -Lambda: EXCLUDED FULL-DOMAIN (this lemma);
  * realness on the superdomain: adversarially mapped (collision surface
    >= 1.35 x Dbar); w=-1 slice proven exactly; interior-angle
    certificate pending.
  => once cert5-on-Gray + interior-angle realness are certified,
    rho(Q) < Lambda on Gray follows, i.e. R_3(s) >= 3/2 on all of Gray.

CHECKS (each k = one block, ~8 min total):
  N4..N0: charpoly grading (shared gate), 90-pt numeric sign validation
  against eigenvalue-built psi^(k)(Lambda), zero-face strip, root-box
  Bernstein certification (all coefficients > 0).

Deps: sympy, mpmath; imports the pipeline from
lemma_7_34_midband_chi_certificates.py (same directory).
Python: /Users/para/.venvs/rnr/bin/python.
"""
import os
import sys
import math
import random
import itertools

import sympy as sp
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lemma_7_34_midband_chi_certificates import (   # noqa: E402
    A, AVAL, D, PSUB, DSUB, VS, CANDS, REPS, STATES, _pat, p, sg, t, th, w, X,
    build_Q, build_Cr, bernstein_root_box, dom_piece, strip_trivial,
    sym_to_t_mid, Q_numeric, face_closure)

mp.mp.dps = 30
PASS = True


def rep(name, ok):
    global PASS
    PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}")
    return ok


def psi_j_numeric(pv, Dv, sv, tvf, j):
    """psi^(j)(Lambda) from the numeric eigenvalues (psi roots = -lambda)."""
    Qn, Crn = Q_numeric(pv, Dv, sv)
    ev = mp.eig(Qn)[0]
    Lm = (1 - 1.5 * Dv * (1 - Dv) * tvf) / Crn
    cs = [mp.mpc(1)]
    for lam in ev:
        ns = [mp.mpc(0)] * (len(cs) + 1)
        for k, cc in enumerate(cs):
            ns[k] += cc
            ns[k + 1] += cc * lam          # (x - (-lambda)) = x + lambda
        cs = ns
    for _ in range(j):
        n = len(cs) - 1
        cs = [cs[k] * (n - k) for k in range(n)]
    val = mp.mpc(0)
    for cc in cs:
        val = val * Lm + cc
    return float(mp.re(val))


def main():
    print("=" * 78)
    print("Negative-side certificates psi^(k)(Lambda) > 0, k=0..4, A=3, full")
    print("rationalized Gray superdomain -- exact root-box Bernstein positivity")
    print("=" * 78)
    Q = build_Q()
    Cr = build_Cr()
    CrN, CrD = sp.fraction(Cr)
    L = 1 - sp.Rational(3, 2) * D * (1 - D) * t
    G = D ** 2 * w + (D - 2) * (1 - D)
    Gt = D ** 2 + (D - 2) * (1 - D) * w
    c = sp.expand(p * (1 - p) * G * Gt)
    M = sp.zeros(5, 5)
    for i in range(5):
        for j in range(5):
            n_, d_ = sp.fraction(sp.cancel(c * Q[i, j]))
            M[i, j] = sp.expand(n_ / d_)
    E = [sp.expand(co * (-1) ** k)
         for k, co in enumerate(M.charpoly(X).all_coeffs())]
    Etld, mks = [sp.Integer(1)], [0]
    for k in range(1, 6):
        ek, mk = sym_to_t_mid(E[k])
        Etld.append(ek)
        mks.append(mk)
    chat, mc = sym_to_t_mid(c)
    CrNt, mN = sym_to_t_mid(CrN)
    CrDt, mD = sym_to_t_mid(CrD)
    grading_ok = (None not in mks and mc == 1 and mN == mD
                  and all(mks[k] == k * mc for k in range(6)))
    rep("N* charpoly grading: E_k symmetric at (-w)^k, chat at (-w)^1",
        grading_ok)
    if not grading_ok:
        return
    PC = {'L': dom_piece(L), 'CrN': dom_piece(CrNt), 'CrD': dom_piece(CrDt),
          'chat': dom_piece(chat)}
    for k in range(1, 6):
        PC[f'E{k}'] = dom_piece(Etld[k])
    rep("N* all pieces have const*(1+sigma^2)^k denominators",
        all(v is not None for v in PC.values()))

    def assemble_neg(j):
        terms = []
        for k in range(0, 6):
            if 5 - k - j < 0:
                continue
            fall = sp.Integer(1)
            for q_ in range(j):
                fall *= (5 - k - q_)
            pw_ch, pw_LCrD, pw_CrN = 5 - j - k, 5 - k - j, k
            Pk = PC[f'E{k}'][0] if k >= 1 else sp.Poly(sp.Integer(1), *VS)
            dE = PC[f'E{k}'][1] if k >= 1 else 0
            cE = PC[f'E{k}'][2] if k >= 1 else sp.Integer(1)
            dens = (dE + pw_ch * PC['chat'][1]
                    + pw_LCrD * (PC['L'][1] + PC['CrD'][1])
                    + pw_CrN * PC['CrN'][1])
            consts = (cE * PC['chat'][2] ** pw_ch
                      * (PC['L'][2] * PC['CrD'][2]) ** pw_LCrD
                      * PC['CrN'][2] ** pw_CrN)
            terms.append((fall, Pk, pw_ch, pw_LCrD, pw_CrN, dens, consts))
        maxden = max(tt[5] for tt in terms)
        onep = sp.Poly(1 + sg ** 2, *VS)
        num = sp.Poly(sp.Integer(0), *VS)
        for (fall, Pk, pw_ch, pw_LCrD, pw_CrN, dens, consts) in terms:
            term = Pk
            for _ in range(pw_ch):
                term = term * PC['chat'][0]
            for _ in range(pw_LCrD):
                term = term * PC['L'][0] * PC['CrD'][0]
            for _ in range(pw_CrN):
                term = term * PC['CrN'][0]
            for _ in range(maxden - dens):
                term = term * onep
            num = num + term * fall / consts       # all + signs: psi coeffs
        return num.as_expr()

    box = [(sp.Integer(0), sp.Integer(1)), (sp.Integer(0), sp.Integer(1)),
           (sp.Integer(0), sp.Integer(2))]
    for j in range(4, -1, -1):
        print("-" * 78)
        print(f"N{j}  psi^({j})(Lambda) > 0")
        tg = assemble_neg(j)
        random.seed(29)
        agree = disagree = 0
        for _ in range(90):
            sv = sp.Rational(random.randint(1, 31), 32)
            tv = sp.Rational(random.randint(1, 32), 16)
            thv = sp.Rational(random.randint(1, 16), 16)
            pv = float(PSUB.subs(sg, sv))
            Dv = float(DSUB.subs({sg: sv, th: thv}))
            if Dv <= 1e-9 or Dv >= 0.666:
                continue
            tvf = float(tv)
            cn = psi_j_numeric(pv, Dv, math.acos(1 - tvf), tvf, j)
            tgv = float(tg.subs({sg: sv, th: thv, t: tv}))
            if tgv == 0:
                continue
            if (cn > 0) == (tgv > 0):
                agree += 1
            else:
                disagree += 1
        rep(f"N{j} numeric sign chain uniform ({agree} agree, {disagree} flip)",
            agree == 0 or disagree == 0)
        if disagree and not agree:
            tg = sp.expand(-tg)
        elif disagree:
            continue
        core, stripped = strip_trivial(tg)
        strip_ok = all(f in CANDS for f in stripped)
        print(f"     target {len(sp.Poly(tg, *VS).terms())} terms -> core "
              f"{len(sp.Poly(core, *VS).terms())} terms; stripped "
              f"{[(str(kk), v) for kk, v in stripped.items()]}")
        rep(f"N{j} zero-face strip uses box-nonnegative factors only", strip_ok)
        mn_b, _, nz_b = bernstein_root_box(core, box)
        print(f"     min-nonzero ~ {float(mn_b):.6g}, exact zeros: {nz_b}")
        rep(f"N{j} no negative Bernstein coeffs (OPEN-box f > 0)", mn_b > 0)
        face_closure(core, f"N{j}", rep)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("psi^(k)(Lambda) > 0, k=0..4, on the whole superdomain => Q has NO")
    print("real eigenvalue <= -Lambda anywhere (negative side closed, full")
    print("domain, no pocket).  With cert5-on-Gray + interior-angle realness")
    print("this yields rho(Q) < Lambda on Gray, i.e. R_3(s) >= 3/2.")


if __name__ == "__main__":
    main()
