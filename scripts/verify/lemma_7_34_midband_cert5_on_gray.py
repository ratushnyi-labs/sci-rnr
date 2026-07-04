#!/usr/bin/env python3
r"""
lemma_7_34_midband_cert5_on_gray.py
============================================================================
CERT5-ON-GRAY (BUG-009-D / Route B, A=3): exact certificate that the k=0
charpoly sign  chi(Lambda) > 0  holds on a GRAY-CONTAINING subdomain of the
rationalized superdomain -- the certificate that is FALSE on the full
superdomain (real Perron pocket D_perron(p) in (D_c, Dbar), exact witness
chi(Lambda) = -0.00073 at p = 8/15, D = Dbar, w = -1; reproduced as check P5).

STATEMENT (all exact, A = 3).  Rationalized domain
    p = (2/3)(1 - sg^2),   Dbar = (2/3)(1-sg)^2/(1+sg^2),   D = th*Dbar,
    t = 1 - cos s in (0, 2],  sg in (0, 1),  th in (0, 1].
Separating curve (rational, degree 2):
    thetatilde(sg) = 1 - (5/16) sg^2.
 (A)  chi(Lambda) > 0  for all  sg in (0, 25/32],  th in (0, thetatilde(sg)],
      t in (0, 2].
 (B)  D_c(p) < thetatilde(sg) * Dbar(p)  for all  sg in (0, 25/32]
      (equivalently p in [133/512, 2/3) ~ [0.2598, 2/3), since
      p(25/32) = (2/3)(1 - 625/1024) = 133/512), where D_c is the ternary
      Gray threshold = smallest positive
      root of the committed quartic Q(D,p) of Lemma 7.34j
      (scripts/verify/probe_7_34j_ternary_gray.py).
 =>   the Gray region {0 < D <= D_c(p)} x {angles s in (0, pi]} for
      p in [0.2598, 2/3) lies strictly inside the region where (A) certifies
      chi(Lambda) > 0.

SCOPE SPLIT (explicit, honest grades).  (A)+(B) cover sg in (0, 25/32],
i.e. p in [133/512, 2/3) (p/pmax >= 399/1024 ~ 0.3896), at FULL CONTINUUM
grade -- this is the certified region of this lemma.  Below that the
committed Schur small/moderate-p certificate reaches p/pmax <= 0.40, i.e.
sg >= sqrt(3/5) ~ 0.7746, and the interval arithmetic overlaps (exactly:
3/5 < 625/1024, check P14) -- BUT the Schur lemma is by its own header
DISCRETE in p (grid p/pmax in {0.05,...,0.40}) and grid-grade in (D,s),
so coverage of p < 133/512 is discrete-p/grid-grade evidence, NOT a
continuum certificate.  Whole-Gray continuum coverage below p = 133/512
awaits the Schur grid-to-interval upgrade (known residual).

SIGN LINKAGE (exact, no numerics needed).  With the committed t-domain
assembly (lemma_7_34_midband_chi_certificates.py, run_certs_45):
    target := assemble(0) = (1+sg^2)^K * chat_t^5 * CrN_t^5 * chi(Lambda)
on the rationalized domain, where (check P3, exact expansions)
    chat_t = p(p-1) * B(D,t),      CrN_t = -B(D,t),
    B(D,t) = 2 D^2 t (1-D)(2-D) + (3D-2)^2  >  0   for D < 2/3
    (sum of nonnegative terms, second strictly positive),
hence  chat_t^5 CrN_t^5 = [p(1-p)]^5 B^10 > 0  on the claim region and
    sign(target) = sign(chi(Lambda))            (P4 = 100-pt numeric guard).
The identity itself is the charpoly expansion chi = sum (-1)^k e_k Lam^{5-k}
with e_k = Etld_k / chat^k, Lam = L*CrD_t/CrN_t, valid once the grading
gates P1 hold (E_k symmetric at (-w)^k, chat at (-w)^1, mN = mD).

CERTIFICATE (A) STRUCTURE.  target = (1-sg)^7 * th * t * core / lcm  (exact
division, P6; all three stripped factors > 0 on the claim region).  Substitute
    sg = (25/32) x,   th = u * (16384 - 3125 x^2)/16384,   t = 2 tau
and clear denominators (integer tensor, exact cross-check P7).  Then on
[0,1]^3 in (x, u, tau):
    P8   ALL scaled Bernstein coefficients >= 0, not all 0
         =>  core_u > 0 on the OPEN box (0,1)^3;
    P9   faces u=1, tau=1, x=1: 2-var Bernstein no-negatives + nonzero
         =>  core_u > 0 on the three open faces;
    P10  edges (u=1,tau=1), (x=1,tau=1), (x=1,u=1): 1-var no-negatives +
         nonzero => positive on the open edges; corner (1,1,1): exact
         positive value.
    Union = (0,1]^3 exactly = the claim region of (A).  (x=0 <=> p = 2/3,
    u=0 <=> D=0, tau=0 <=> s=0 are excluded limits of the Gray family.)
    A zero-coefficient count is reported: zeros are allowed by the open-box
    +closure logic; the boundary zero at (x=0, u=1) (the degenerate corner
    p -> 2/3, D -> 2/3, where chi -> 0 like sg^10) is exactly why NO closed-
    box certificate can exist and the faces/edges tower is used instead.

CERTIFICATE (B) (exact Sturm).  Qq(0,p) = 16 p^2 (p^2+1) > 0 (identity
P11a) and Qq(thetatilde*Dbar, p) < 0 for sg in (0, 25/32]:  after clearing
the positive denominator 20736 (1+sg^2)^2 and stripping sg^3, the degree-17
integer polynomial has NO roots in [0, 25/32] (Sturm) and is negative there
(P11).  A continuous function positive at D=0 and negative at
D = thetatilde*Dbar has a root below the curve; the smallest positive root
of Qq is D_c (committed Lemma 7.34j; only the SAFE direction D_c <= that
root is used -- 7.34j V3/V4 identify the root as the alternating-word
threshold, an upper bound for the all-word threshold D_c).

HONEST NOTES.
 *  The gap theta_perron - theta_c does NOT survive at sg -> 0 in theta
    units (both -> 1); the pocket itself is absent for p/pmax > 0.92
    (sg < ~0.28), and thetatilde(0) = 1 rides through the pinch; the
    binding constraint near sg -> 0 is theta_c < thetatilde, which holds
    with margin (theta_c = 1 - sg^2(1+o(1)), thetatilde = 1 - (5/16)sg^2).
 *  At sg -> 1 (p -> 0) the clearance theta_x - theta_c -> 0; that end is
    NOT covered here -- it is covered by the committed Schur certificate
    (overlap at [sqrt(3/5), 25/32], check P14).
 *  chi(Lambda) > 0 is the k=0 member of the Budan-Fourier tower; the
    negative-side certificates chi^(k)(-Lambda) and interior-angle realness
    remain separate open items (unchanged by this script).

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.
Runtime: ~10-15 min.  Machinery adapted from the committed
scripts/verify/lemma_7_34_midband_chi_certificates.py (build_Q, build_Cr,
Chebyshev t-domain assembly) and probe_7_34j_ternary_gray.py (quartic).
"""
import itertools
import math
import random
import time
from fractions import Fraction as Fr
from math import comb

import sympy as sp
import mpmath as mp

mp.mp.dps = 30
T0 = time.time()
PASS = True


def rep(name, ok):
    global PASS
    PASS = PASS and bool(ok)
    print(f"  {name:<68} {'PASS' if ok else 'FAIL'}  [{time.time()-T0:.0f}s]",
          flush=True)
    return ok


# ---------------------------------------------------------------------------
# symbols and the rationalized domain (A = 3)
# ---------------------------------------------------------------------------
p, D, w, t = sp.symbols('p D w t')
sg, th, U, X = sp.symbols('sigma theta u x')
AVAL = 3
A = sp.Integer(AVAL)
PSUB = sp.Rational(2, 3) * (1 - sg ** 2)
DBAR = sp.Rational(2, 3) * (1 - sg) ** 2 / (1 + sg ** 2)
DSUB = th * DBAR
VS = (sg, th, t)
THT = 1 - sp.Rational(5, 16) * sg ** 2          # the separating curve
SGMAX = sp.Rational(25, 32)

# ---------------------------------------------------------------------------
# replica pattern-quotient machinery (committed pipeline, A=3)
# ---------------------------------------------------------------------------
STATES = list(itertools.product(range(AVAL), repeat=3))


def _pat(tr):
    x_, u_, up = tr
    if x_ == u_ == up: return 0
    if x_ == u_ and u_ != up: return 1
    if x_ == up and u_ != up: return 2
    if u_ == up and x_ != u_: return 3
    return 4


REPS = [None] * 5
for _k, _tr in enumerate(STATES):
    if REPS[_pat(_tr)] is None:
        REPS[_pat(_tr)] = _k

eth = D / ((A - 1) * (1 - D))


def _etaf(W):
    E = eth * W
    return ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))


def build_Q():
    eta, etb = _etaf(w), _etaf(1 / w)
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(AVAL)]
         for i in range(AVAL)]
    Q = sp.zeros(5, 5)
    for a_ in range(5):
        x_, xp, xq = STATES[REPS[a_]]
        for (y, yp, yq) in STATES:
            term = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x_][y])
            if yp != y: term *= eta
            if yq != y: term *= etb
            Q[a_, _pat((y, yp, yq))] += term
    return Q


def build_Cr():
    def Cf(W):
        E = eth * W
        return ((A - 1) * D * E + D - (A - 1)) / (A * D - (A - 1))
    return sp.cancel(Cf(w) * Cf(1 / w) / (Cf(1) ** 2))


UU = sp.symbols('uu')


def laurent_sym_to_t(expr, K=60):
    """w-symmetric Laurent polynomial -> t-polynomial via Chebyshev in
    uu = w + 1/w = 2 - 2t on |w| = 1; None if asymmetric."""
    ex = sp.expand(expr)
    npoly = sp.Poly(sp.expand(ex * w ** K), w)
    co = {}
    for (md,), c in npoly.terms():
        co[md - K] = co.get(md - K, 0) + c
    ks = sorted(co.keys())
    maxk = max(abs(k) for k in ks) if ks else 0
    C = {0: sp.Integer(2), 1: UU}
    for k2 in range(2, maxk + 1):
        C[k2] = sp.expand(UU * C[k2 - 1] - C[k2 - 2])
    out = sp.Integer(0)
    done = set()
    for k2 in ks:
        if k2 in done: continue
        if k2 == 0:
            out += co[0]; done.add(0)
        else:
            kk = abs(k2)
            if sp.simplify(co.get(kk, 0) - co.get(-kk, 0)) != 0:
                return None
            out += co.get(kk, 0) * C[kk]
            done.add(kk); done.add(-kk)
    return sp.expand(out.subs(UU, 2 - 2 * t))


def sym_to_t_mid(expr, K=60):
    ex = sp.expand(expr)
    poly = sp.Poly(sp.expand(ex * w ** K), w)
    degs = [md - K for (md,), cc in poly.terms()]
    mid = sp.Rational(min(degs) + max(degs), 2)
    if int(mid) != mid:
        return None, None
    r = laurent_sym_to_t(sp.cancel(ex / (-w) ** int(mid)), K)
    return r, (int(mid) if r is not None else None)


def dom_piece(ex):
    """(Poly on (sg,th,t), (1+sg^2)-exponent, numeric const) with denominator
    exactly const*(1+sg^2)^k after the domain substitution, else None."""
    n_, d_ = sp.fraction(sp.cancel(sp.together(ex.subs({p: PSUB, D: DSUB}))))
    f = sp.factor_list(d_)
    const, expo = f[0], 0
    for b, e in f[1]:
        if sp.simplify(b - (1 + sg ** 2)) == 0:
            expo += e
        elif b.is_number:
            const *= b ** e
        else:
            return None
    return sp.Poly(sp.expand(n_), *VS), expo, const


def Q_numeric(pv, Dv, sv):
    thv = mp.log(Dv / ((AVAL - 1) * (1 - Dv)))
    E = mp.e ** (thv + 1j * sv)
    et = ((AVAL - 1) * Dv * E + Dv - (AVAL - 1) * E) / \
        ((AVAL - 1) * Dv * E + Dv - (AVAL - 1))
    Cn = ((AVAL - 1) * Dv * E + Dv - (AVAL - 1)) / (AVAL * Dv - (AVAL - 1))
    C0 = ((AVAL - 1) * Dv * mp.e ** thv + Dv - (AVAL - 1)) / \
        (AVAL * Dv - (AVAL - 1))
    Crn = abs(Cn / C0) ** 2
    Tn = [[(1 - pv) if i == j else pv / (AVAL - 1) for j in range(AVAL)]
          for i in range(AVAL)]
    etb = mp.conj(et)
    Qn = mp.zeros(5, 5)
    for a2 in range(5):
        x_, xp, xq = STATES[REPS[a2]]
        for (y, yp, yq) in STATES:
            Qn[a2, _pat((y, yp, yq))] += (Tn[xp][yp] * Tn[xq][yq] / Tn[x_][y]
                                          * (et if yp != y else 1)
                                          * (etb if yq != y else 1))
    return Qn, Crn


def chi_numeric(pv, Dv, sv, tvf):
    """chi(Lambda) from the numeric eigenvalues of Q."""
    Qn, Crn = Q_numeric(pv, Dv, sv)
    ev = mp.eig(Qn)[0]
    Lm = (1 - 1.5 * Dv * (1 - Dv) * tvf) / Crn
    val = mp.mpf(1)
    out = mp.mpc(1)
    for lam in ev:
        out *= (Lm - lam)
    return float(mp.re(out))


# ---------------------------------------------------------------------------
# the committed ternary Gray quartic (Lemma 7.34j)
# ---------------------------------------------------------------------------
Qquartic = (9 * p ** 2 * (3 * p - 4) ** 2 * D ** 4
            - 4 * p * (3 * p - 4) * (27 * p ** 2 - 30 * p - 4) * D ** 3
            + (441 * p ** 4 - 864 * p ** 3 + 432 * p ** 2 + 64) * D ** 2
            - 8 * (27 * p ** 4 - 33 * p ** 3 + 36 * p ** 2 - 32 * p + 16) * D
            + 16 * p ** 2 * (p ** 2 + 1))


def Dc3_numeric(pv):
    rts = sp.Poly(Qquartic.subs({p: sp.Float(pv, 25)}), D).nroots(n=25)
    cand = sorted(float(sp.re(r)) for r in rts
                  if abs(float(sp.im(r))) < 1e-15
                  and 0 < float(sp.re(r)) < 2.0 / 3.0)
    return cand[0] if cand else None


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("CERT5-ON-GRAY: chi(Lambda) > 0 on th <= thetatilde(sg) = 1-(5/16)sg^2,")
    print("sg in (0, 25/32], all angles; + D_c < thetatilde*Dbar (exact Sturm).")
    print("A = 3.  Exact rational/integer arithmetic end to end.")
    print("=" * 78)

    # ---- build the t-domain pieces (committed assembly) -------------------
    print("building Q, Cr, charpoly pieces ...", flush=True)
    Q = build_Q()
    Cr = build_Cr()
    CrN, CrD = sp.fraction(Cr)
    L = 1 - sp.Rational(3, 2) * D * (1 - D) * t
    G = D ** 2 * w + (D - 2) * (1 - D)
    Gt = D ** 2 + (D - 2) * (1 - D) * w
    c = sp.expand(p * (1 - p) * G * Gt)
    XX = sp.symbols('XX')
    M = sp.zeros(5, 5)
    for i in range(5):
        for j in range(5):
            n_, d_ = sp.fraction(sp.cancel(c * Q[i, j]))
            M[i, j] = sp.expand(n_ / d_)
    E = [sp.expand(co * (-1) ** k)
         for k, co in enumerate(M.charpoly(XX).all_coeffs())]
    Etld, mks = [sp.Integer(1)], [0]
    for k in range(1, 6):
        ek, mk = sym_to_t_mid(E[k])
        Etld.append(ek); mks.append(mk)
    chat, mc = sym_to_t_mid(c)
    CrNt, mN = sym_to_t_mid(CrN)
    CrDt, mD = sym_to_t_mid(CrD)

    print("-" * 78)
    print("P1  charpoly grading gates (exact)")
    rep("P1 E_k symmetric at (-w)^k, chat at (-w)^1, mN == mD",
        (None not in mks and mc == 1 and mN == mD
         and all(mks[k] == k * mc for k in range(6))))

    PC = {'L': dom_piece(L), 'CrN': dom_piece(CrNt), 'CrD': dom_piece(CrDt),
          'chat': dom_piece(chat)}
    for k in range(1, 6):
        PC[f'E{k}'] = dom_piece(Etld[k])
    print("-" * 78)
    print("P2  domain pieces")
    rep("P2 all pieces have const*(1+sigma^2)^k denominators",
        all(v is not None for v in PC.values()))

    # ---- P3: exact multiplier sign structure ------------------------------
    print("-" * 78)
    print("P3  multiplier identities (exact expansions)")
    Bpoly = 2 * D ** 2 * t * (1 - D) * (2 - D) + (3 * D - 2) ** 2
    ok3 = (sp.expand(chat - p * (p - 1) * Bpoly) == 0)
    ok3 &= (sp.expand(CrNt + Bpoly) == 0)
    ok3 &= (sp.expand(CrDt + (3 * D - 2) ** 2) == 0)
    rep("P3 chat_t = p(p-1)B, CrN_t = -B, B = 2D^2t(1-D)(2-D)+(3D-2)^2", ok3)
    rep("P3b => chat_t^5 CrN_t^5 = [p(1-p)]^5 B^10 > 0 for 0<p<1, D<2/3",
        ok3)   # algebraic consequence; recorded as a gate

    # ---- assemble target = (1+sg^2)^K * chat^5 CrN^5 * chi ---------------
    print("-" * 78)
    print("P4  target assembly + numeric sign guard")

    def assemble(j):
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
            terms.append((k, fall, Pk, pw_ch, pw_LCrD, pw_CrN, dens, consts))
        maxden = max(tt[6] for tt in terms)
        onep = sp.Poly(1 + sg ** 2, *VS)
        num = sp.Poly(sp.Integer(0), *VS)
        for (k, fall, Pk, pw_ch, pw_LCrD, pw_CrN, dens, consts) in terms:
            term = Pk
            for _ in range(pw_ch):
                term = term * PC['chat'][0]
            for _ in range(pw_LCrD):
                term = term * PC['L'][0] * PC['CrD'][0]
            for _ in range(pw_CrN):
                term = term * PC['CrN'][0]
            for _ in range(maxden - dens):
                term = term * onep
            num = num + term * sp.Rational((-1) ** k, 1) * fall / consts
        return num

    Ptg = assemble(0)
    tg = Ptg.as_expr()
    print(f"     target: {len(Ptg.terms())} terms, degs "
          f"{[Ptg.degree(v) for v in VS]}", flush=True)

    # NOTE: the target has huge inter-term cancellation near the degenerate
    # corner (true values ~1e-17 of term size); double precision FLIPS signs
    # there (caught in development: 12/140 spurious flips with a float
    # lambdify).  Evaluate the target with mpmath at 60 digits.
    fgt = sp.lambdify(VS, tg, modules="mpmath")
    random.seed(23)
    agree = disagree = npts = 0
    with mp.workdps(60):
        for _ in range(140):
            sv = sp.Rational(random.randint(1, 31), 32)
            tv = sp.Rational(random.randint(1, 32), 16)
            thv = sp.Rational(random.randint(1, 16), 16)
            pv = float(PSUB.subs(sg, sv))
            Dv = float(DSUB.subs({sg: sv, th: thv}))
            if Dv <= 1e-9 or Dv >= 0.666:
                continue
            tvf = float(tv)
            cn = chi_numeric(pv, Dv, math.acos(1 - tvf), tvf)
            tgv = fgt(mp.mpf(sv.p) / sv.q, mp.mpf(thv.p) / thv.q,
                      mp.mpf(tv.p) / tv.q)
            if tgv == 0 or abs(cn) < 1e-24:
                continue
            npts += 1
            if (cn > 0) == (tgv > 0):
                agree += 1
            else:
                disagree += 1
    rep(f"P4 numeric sign chain sign(target)=sign(chi), {npts} pts, "
        f"{disagree} flips", npts >= 80 and disagree == 0)

    # ---- P5: pocket witness (why the restriction is necessary) ------------
    print("-" * 78)
    print("P5  pocket witness: target < 0 at (sg=559/1250, th=1, t=2) [exact]")
    wit = tg.subs({sg: sp.Rational(559, 1250), th: 1, t: 2})
    print(f"     target = {float(wit):.6f}")
    rep("P5 exact target < 0 at the pocket point (p ~ 8/15, D = Dbar, s=pi)",
        wit < 0)

    # ---- P6: strip trivial factors, integerize ----------------------------
    print("-" * 78)
    print("P6  strip (1-sg)^7 * th * t and integerize (exact divisions)")
    corePoly = Ptg
    mults = {}
    for f, expect in (((1 - sg), 7), (th, 1), (t, 1), (sg, 0), ((1 + sg), 0),
                      ((1 + sg ** 2), 0)):
        fP = sp.Poly(f, *VS)
        while True:
            q_, r_ = sp.div(corePoly, fP)
            if r_.is_zero and not q_.is_zero:
                corePoly = q_
                mults[f] = mults.get(f, 0) + 1
            else:
                break
    ok6 = (mults.get(1 - sg, 0) == 7 and mults.get(th, 0) == 1
           and mults.get(t, 0) == 1 and len(mults) == 3)
    print(f"     multiplicities: {[(str(k), v) for k, v in mults.items()]}; "
          f"core {len(corePoly.terms())} terms, degs "
          f"{[corePoly.degree(v) for v in VS]}")
    rep("P6 target = (1-sg)^7 * th * t * core exactly (all factors > 0 on "
        "claim region)", ok6)
    dens = [sp.fraction(cc)[1] for cc in corePoly.coeffs()]
    lc = sp.ilcm(*[int(d_) for d_ in dens])
    Tcore = {}
    intOK = True
    for mon, cc in zip(corePoly.monoms(), corePoly.coeffs()):
        val = cc * lc
        if not (val.is_rational and sp.Integer(val) == val):
            intOK = False
            break
        Tcore[mon] = int(val)
    rep(f"P6b core * lcm({lc}) has integer coefficients", intOK)

    # ---- P7: u-substitution to (x,u,tau) on [0,1]^3, exact cross-check ----
    print("-" * 78)
    print("P7  substitution sg=(25/32)x, th=u(16384-3125x^2)/16384, t=2tau")
    DS = corePoly.degree(sg)
    DTH = corePoly.degree(th)
    DT = corePoly.degree(t)
    DX = DS + 2 * DTH
    POW = {0: {0: 1}}
    for b in range(1, DTH + 1):
        prev = POW[b - 1]
        cur = {}
        for jj, cf in prev.items():
            cur[jj] = cur.get(jj, 0) + 16384 * cf
            cur[jj + 1] = cur.get(jj + 1, 0) - 3125 * cf
        POW[b] = cur
    P25 = [25 ** a for a in range(DS + 1)]
    P32 = [32 ** a for a in range(DS + 1)]
    P2l = [2 ** cc for cc in range(DT + 1)]
    P16384 = [16384 ** b for b in range(DTH + 1)]
    # core(sg,th,t)*32^DS*16384^DTH with the substitution; the 32-power
    # compensates ONLY the (25/32)^a sigma-scaling (x^{2j} terms from the
    # thetatilde expansion carry no 32-denominator -- a wrong 32^(DS-a-2j)
    # here would move the curve into the pocket; caught in development by
    # exactly the P7 cross-check below).
    TU = {}
    for (a, b, cc), cf in Tcore.items():
        base = cf * P25[a] * P32[DS - a] * P2l[cc] * P16384[DTH - b]
        for jj, pc in POW[b].items():
            key = (a + 2 * jj, b, cc)
            TU[key] = TU.get(key, 0) + base * pc
    TU = {k: v for k, v in TU.items() if v != 0}
    print(f"     integer tensor: {len(TU)} nonzero, degs ({DX},{DTH},{DT})")
    SCALE = Fr(32) ** DS * Fr(16384) ** DTH
    random.seed(11)
    ok7 = True
    for _ in range(12):
        xv = Fr(random.randint(0, 16), 16)
        uv = Fr(random.randint(0, 16), 16)
        tv = Fr(random.randint(0, 16), 16)
        sgv = Fr(25, 32) * xv
        thv = uv * (1 - Fr(5, 16) * sgv ** 2)
        ttv = 2 * tv
        lhs = sum(Fr(cf) * xv ** a2 * uv ** b2 * tv ** c2
                  for (a2, b2, c2), cf in TU.items())
        rhs = SCALE * sum(Fr(cf) * sgv ** a2 * thv ** b2 * ttv ** c2
                          for (a2, b2, c2), cf in Tcore.items())
        ok7 &= (lhs == rhs)
    rep("P7 exact Fraction cross-check core_u == scaled core at 12 points", ok7)
    # P7b: conclusive polynomial-identity check by Schwartz-Zippel over F_q,
    # q = 2^61 - 1 (prime).  The difference polynomial has total degree
    # <= 223, so a nonzero difference survives 20 uniform F_q points with
    # probability <= (223/q)^20 ~ 1e-370.
    qP = (1 << 61) - 1
    inv32 = pow(32, qP - 2, qP)
    inv16384 = pow(16384, qP - 2, qP)
    scale_q = pow(32, DS, qP) * pow(16384, DTH, qP) % qP
    random.seed(101)
    ok7b = True
    for _ in range(20):
        xq = random.randrange(qP)
        uq = random.randrange(qP)
        tq = random.randrange(qP)
        lhs = 0
        pxa = {a2: pow(xq, a2, qP) for a2 in set(k[0] for k in TU)}
        pub = {b2: pow(uq, b2, qP) for b2 in set(k[1] for k in TU)}
        ptc = {c2: pow(tq, c2, qP) for c2 in set(k[2] for k in TU)}
        for (a2, b2, c2), cf in TU.items():
            lhs = (lhs + cf * pxa[a2] % qP * pub[b2] % qP * ptc[c2]) % qP
        sgq = 25 * xq % qP * inv32 % qP
        thq = uq * ((16384 - 3125 * xq * xq) % qP) % qP * inv16384 % qP
        ttq = 2 * tq % qP
        psa = {a2: pow(sgq, a2, qP) for a2 in set(k[0] for k in Tcore)}
        psb = {b2: pow(thq, b2, qP) for b2 in set(k[1] for k in Tcore)}
        psc = {c2: pow(ttq, c2, qP) for c2 in set(k[2] for k in Tcore)}
        rhs = 0
        for (a2, b2, c2), cf in Tcore.items():
            rhs = (rhs + cf * psa[a2] % qP * psb[b2] % qP * psc[c2]) % qP
        ok7b &= (lhs % qP == rhs * scale_q % qP)
    rep("P7b Schwartz-Zippel identity check mod 2^61-1, 20 uniform points",
        ok7b)

    # ---- P8: root-box scaled-Bernstein, no negatives ----------------------
    print("-" * 78)
    print("P8  root box [0,1]^3: scaled Bernstein coefficients")

    def bern_axis(T, degs, ax):
        d_ = degs[ax]
        out = {}
        for mon, cval in T.items():
            a2 = mon[ax]
            lm = list(mon)
            for b2 in range(a2, d_ + 1):
                f = comb(d_ - a2, b2 - a2)
                lm[ax] = b2
                key = tuple(lm)
                out[key] = out.get(key, 0) + f * cval
        return {k: v for k, v in out.items() if v != 0}

    degs3 = [DX, DTH, DT]
    BB = dict(TU)
    for ax in (2, 1, 0):
        BB = bern_axis(BB, degs3, ax)
    total = (DX + 1) * (DTH + 1) * (DT + 1)
    negs = sum(1 for v in BB.values() if v < 0)
    nz = total - len(BB)
    print(f"     nonzero {len(BB)}, exact zeros {nz}, negatives {negs}")
    rep("P8 NO negative Bernstein coefficients (=> core_u > 0 on OPEN box)",
        negs == 0 and len(BB) > 0)

    # ---- P9: face certificates (direct, 2-var) ----------------------------
    print("-" * 78)
    print("P9  faces u=1, tau=1, x=1 (2-var Bernstein no-negatives + nonzero)")

    def collapse(T, ax, degs):
        out = {}
        for mon, cf in T.items():
            key = tuple(m for i, m in enumerate(mon) if i != ax)
            out[key] = out.get(key, 0) + cf     # value at var=1: sum coeffs
        return {k: v for k, v in out.items() if v != 0}, \
            [d_ for i, d_ in enumerate(degs) if i != ax]

    ok9 = True
    for ax, nm in ((1, "u=1"), (2, "tau=1"), (0, "x=1")):
        Tf, degf = collapse(TU, ax, degs3)
        Bf = dict(Tf)
        for ax2 in range(len(degf) - 1, -1, -1):
            Bf = bern_axis(Bf, degf, ax2)
        negf = sum(1 for v in Bf.values() if v < 0)
        okf = (negf == 0 and len(Bf) > 0)
        ok9 &= okf
        rep(f"P9 face {nm}: no negatives ({negf}), nonzero ({len(Bf)}) "
            f"=> > 0 on open face", okf)

    # ---- P10: edges and corner --------------------------------------------
    print("-" * 78)
    print("P10 edges (1-var no-negatives + nonzero) and the corner value")
    ok10 = True
    for keep, nm in ((0, "(u=1,tau=1), poly in x"),
                     (1, "(x=1,tau=1), poly in u"),
                     (2, "(x=1,u=1), poly in tau")):
        co = {}
        for mon, cf in TU.items():
            co[mon[keep]] = co.get(mon[keep], 0) + cf
        d_ = degs3[keep]
        # 1-var scaled Bernstein
        bb = [0] * (d_ + 1)
        for a2, cf in co.items():
            if cf == 0:
                continue
            for b2 in range(a2, d_ + 1):
                bb[b2] += comb(d_ - a2, b2 - a2) * cf
        negb = sum(1 for v in bb if v < 0)
        oke = (negb == 0 and any(v > 0 for v in bb))
        ok10 &= oke
        rep(f"P10 edge {nm}: no negatives => > 0 on open edge", oke)
    corner = sum(TU.values())
    rep("P10d corner (x,u,tau)=(1,1,1): exact value > 0", corner > 0)
    print("     [union: open box + 3 open faces + 3 open edges + corner")
    print("      = (0,1]^3 exactly = claim region of (A)]")

    # ---- P11: certificate (B), exact Sturm --------------------------------
    print("-" * 78)
    print("P11 certificate (B): D_c < thetatilde * Dbar on sg in (0, 25/32]")
    id0 = sp.expand(Qquartic.subs(D, 0) - 16 * p ** 2 * (p ** 2 + 1))
    rep("P11a Qq(0,p) == 16 p^2 (p^2+1)  (> 0 for p > 0)", id0 == 0)
    Xcurve = THT * DBAR
    expr = Qquartic.subs({D: Xcurve, p: PSUB})
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    denf = sp.factor_list(den)
    okden = all(b.is_number or sp.simplify(b - (1 + sg ** 2)) == 0
                for b, e in denf[1]) and denf[0] > 0
    rep("P11b cleared denominator = 20736 (1+sg^2)^2  (positive)",
        okden and sp.simplify(den - 20736 * (1 + sg ** 2) ** 2) == 0)
    Rp = sp.Poly(sp.expand(num), sg)
    mindeg = min(m[0] for m in Rp.monoms())
    Rcore = sp.Poly(sp.expand(num / sg ** mindeg), sg)
    nroots = Rcore.count_roots(0, SGMAX)
    v0 = Rcore.eval(0)
    vh = Rcore.eval(sp.Rational(1, 2))
    ve = Rcore.eval(SGMAX)
    print(f"     stripped sg^{mindeg}; deg {Rcore.degree()}; Sturm roots in "
          f"[0,25/32]: {nroots}; values (0,1/2,25/32) signs: "
          f"{sp.sign(v0)},{sp.sign(vh)},{sp.sign(ve)}")
    rep("P11c numerator core has NO roots in [0,25/32] and is negative",
        nroots == 0 and v0 < 0 and vh < 0 and ve < 0 and mindeg == 3)
    rep("P11d => Qq(thetatilde*Dbar) < 0 on (0,25/32]; with P11a and the "
        "first-root logic, D_c < thetatilde*Dbar",
        nroots == 0 and v0 < 0)

    # ---- P12: exact anchor ------------------------------------------------
    print("-" * 78)
    print("P12 exact anchor: chi(Lambda) at sg=1/2, th=thetatilde(1/2)=59/64,")
    print("    w=-1 (rational point on the curve) vs core_u sign")
    pv = sp.Rational(1, 2)
    Dv = sp.Rational(59, 64) * sp.Rational(2, 15)
    ethv = Dv / ((A - 1) * (1 - Dv))

    def etaf_w(Wv):
        E_ = ethv * Wv
        return ((A - 1) * Dv * E_ + Dv - (A - 1) * E_) / \
            ((A - 1) * Dv * E_ + Dv - (A - 1))

    ew = sp.nsimplify(etaf_w(-1))
    Tw = [[(1 - pv) if i == j else pv / (A - 1) for j in range(AVAL)]
          for i in range(AVAL)]
    Qw = sp.zeros(5, 5)
    for a_ in range(5):
        x_, xp, xq = STATES[REPS[a_]]
        for (y, yp, yq) in STATES:
            term = Tw[xp][yp] * Tw[xq][yq] / Tw[x_][y]
            if yp != y: term *= ew
            if yq != y: term *= ew
            Qw[a_, _pat((y, yp, yq))] += term

    def Cfw(Wv):
        E_ = ethv * Wv
        return ((A - 1) * Dv * E_ + Dv - (A - 1)) / (A * Dv - (A - 1))

    Crw = sp.nsimplify(Cfw(-1) * Cfw(-1) / (Cfw(1) ** 2))
    Lw = 1 - sp.Rational(3, 2) * Dv * (1 - Dv) * 2
    chiw = sp.nsimplify(((Lw / Crw) * sp.eye(5) - Qw).det())
    vc = sum(Fr(cf) * Fr(16, 25) ** a2 for (a2, b2, c2), cf in TU.items())
    print(f"     chi(Lambda) = {float(chiw):.6e} (exact rational), "
          f"core_u(16/25,1,1) sign {'+' if vc > 0 else '-'}")
    rep("P12 exact chi > 0 at the anchor and sign matches core_u",
        chiw > 0 and vc > 0)

    # ---- P13: numeric separation scan (consistency) ------------------------
    print("-" * 78)
    print("P13 numeric separation: theta_c < thetatilde < theta_perron (or no")
    print("    pocket) at 25 grid points (mpmath consistency scan)")
    ok13 = True
    for k in range(1, 26):
        sgv = float(SGMAX) * k / 25.0
        pv2 = (2.0 / 3.0) * (1.0 - sgv * sgv)
        Db = (2.0 / 3.0) * (1.0 - sgv) ** 2 / (1.0 + sgv * sgv)
        dth = 1.0 - (5.0 / 16.0) * sgv * sgv
        dc = Dc3_numeric(pv2)
        ok13 &= (dc is not None and dc / Db < dth)
        chc = chi_numeric(pv2, dth * Db, math.pi, 2.0)
        ok13 &= (chc > 0)
    rep("P13 curve strictly above theta_c and chi > 0 on the curve at t=2",
        ok13)

    # ---- P14: scope split --------------------------------------------------
    print("-" * 78)
    print("P14 scope split vs the committed Schur certificate")
    rep("P14 sqrt(3/5) < 25/32 exactly (3/5 < 625/1024): overlap with the "
        "Schur p <= 0.40 pmax region", sp.Rational(3, 5) < SGMAX ** 2)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("CLAIM CERTIFIED (A): chi(Lambda) > 0 for sg in (0,25/32],")
    print("th in (0, 1-(5/16)sg^2], t in (0,2]   [exact Bernstein tower]")
    print("CLAIM CERTIFIED (B): D_c(p) < (1-(5/16)sg^2) * Dbar(p) there")
    print("[exact Sturm; D_c = smallest positive quartic root, Lemma 7.34j]")
    print("=> Gray region for p in [p(25/32), 2/3) = [0.2598.., 2/3) is")
    print("   strictly inside the certified chi > 0 region; p <= 0.40 pmax")
    print("   is covered by the committed Schur certificate (overlap")
    print("   sg in [sqrt(3/5), 25/32]).")
    return PASS


if __name__ == "__main__":
    main()
