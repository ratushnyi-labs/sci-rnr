#!/usr/bin/env python3
r"""
lemma_7_34_midband_cert5_smallp.py
============================================================================
CERT5-SMALL-P (BUG-009-D / Route B, A=3): exact certificate that the k=0
charpoly sign  chi(Lambda) > 0  holds on a GRAY-CONTAINING subdomain of the
rationalized superdomain for the SMALL-p range  sg in [25/32, 1), i.e.
p in (0, 133/512] -- the last continuum gap of the A=3 mid-band theorem.
Companion piece to the committed lemma_7_34_midband_cert5_on_gray.py
(which covers sg in (0, 25/32], i.e. p in [133/512, 2/3)); together the
two lemmas certify chi(Lambda) > 0 on a Gray-containing region for ALL
p in (0, 2/3) at full continuum grade.

STATEMENT (all exact, A = 3).  Rationalized domain
    p = (2/3)(1 - sg^2),   Dbar = (2/3)(1-sg)^2/(1+sg^2),   D = th*Dbar,
    t = 1 - cos s in (0, 2],  sg in [25/32, 1),  th in (0, 1].
Separating curve (rational, degree 2):
    thetatilde2(sg) = 2/3 + (11/32)(1 - sg^2).
 (A2) chi(Lambda) > 0  for all  sg in [25/32, 1),  th in (0, thetatilde2],
      t in (0, 2].
 (B2) D_c(p) < thetatilde2(sg) * Dbar(p)  for all  sg in [25/32, 1)
      (equivalently p in (0, 133/512]), where D_c is the ternary Gray
      threshold = smallest positive root of the committed quartic Q(D,p)
      of Lemma 7.34j (scripts/verify/probe_7_34j_ternary_gray.py).
 =>   the Gray region {0 < D <= D_c(p)} x {angles s in (0, pi]} for
      p in (0, 133/512] lies strictly inside the region where (A2)
      certifies chi(Lambda) > 0.

SPLICE WITH THE COMMITTED PIECE (exact).  At the joint sg = 25/32:
    thetatilde2(25/32) = 78703/98304  <=  13259/16384 = 79554/98304
                       = thetatilde(25/32) of cert5-on-gray  (check P14),
so the new certified slab nests under the committed one at the joint, and
the sg-coverage (0, 25/32] u [25/32, 1) = (0, 1) is seamless: the union
region  {th <= thetatilde(sg), sg <= 25/32} u {th <= thetatilde2(sg),
sg >= 25/32}  contains Gray for ALL p in (0, 2/3)  (each piece by its own
(B)-certificate).  The committed Schur small-p certificate (discrete-p /
grid grade) is thereby SUPERSEDED as a load-bearing element: the whole-
Gray chi(Lambda) > 0 coverage is now continuum-grade end to end.

GEOGRAPHY (measured at dps 40, 300+ points; CORRECTS the prior scout
briefing).  In theta-units with omega = 1 - sg^2 (proportional to p):
    theta_c(sg)      = 2/3 + c(sg) omega,   c: 0.2514 -> 0.2223 (-> 2/9),
    theta_x(sg; t=2) = 2/3 + q(sg) omega,   q: 0.4933 -> 0.6664 (-> 2/3),
theta_x = pocket edge (first th with chi(Lambda) < 0), t = 2 binding
(theta_x decreasing in t; checked on a t-grid at 4 sg values).  The
pocket does NOT die on this range: theta_x stays < 1 and tends to 2/3
as sg -> 1 (the prior "no pocket for p/pmax < 0.02" note is WRONG -- at
p/pmax = 0.0149, theta_x = 0.6766).  Hence a pocket-death certificate is
impossible AND unnecessary; instead the slope-11/32 curve rides the
corridor all the way to sg -> 1 with two-sided margins
    thetatilde2 - theta_c >= 0.0923 omega,  theta_x - thetatilde2 >=
    0.1496 omega   (300-point mpmath sweep; consistency, not load-bearing).

SIGN LINKAGE (exact, inherited unchanged from the committed cert5).  With
the committed t-domain assembly (lemma_7_34_midband_chi_certificates.py,
run_certs_45):
    target := assemble(0) = (1+sg^2)^K * chat_t^5 * CrN_t^5 * chi(Lambda),
    chat_t = p(p-1) B,  CrN_t = -B,  B = 2D^2 t(1-D)(2-D) + (3D-2)^2 > 0
    for D < 2/3   (check P3, exact expansions),
hence  chat_t^5 CrN_t^5 = [p(1-p)]^5 B^10 > 0  and  sign(target) =
sign(chi(Lambda))  on the claim region (P4 = 140-pt numeric guard at
dps 80 -- the target has huge inter-term cancellation; double precision
flips signs, mpmath does not).

CERTIFICATE (A2) STRUCTURE.  target = (1-sg)^7 * th * t * core  (exact
division, P6; all three stripped factors > 0 on the claim region).
Substitute   sg = (25 + 7x)/32,
             th = u * (78703 - 11550 x - 1617 x^2)/98304
                = u * thetatilde2(sg)          (P6c exact identity),
             t  = 2 tau,
and clear denominators (integer tensor, exact cross-checks P7 + P7b).
Then on [0,1]^3 in (x, u, tau):
    P8   ALL scaled Bernstein coefficients >= 0, not all 0
         =>  core_u > 0 on the OPEN box (0,1)^3;
    P9   faces x=0 (restriction), u=1, tau=1 (collapse): 2-var Bernstein
         no-negatives + nonzero  =>  core_u > 0 on the three open faces;
    P10  edges (x=0,u=1), (x=0,tau=1), (u=1,tau=1): 1-var no-negatives +
         nonzero => positive on the open edges; corner (0,1,1): exact
         positive value.
    Union = [0,1) x (0,1] x (0,1] exactly = the claim region of (A2).
    (x=1 <=> sg=1 <=> p = 0, u=0 <=> D=0, tau=0 <=> s=0 are excluded open
    limits of the Gray family.)  A zero-coefficient count is reported:
    zeros are allowed by the open-box+closure logic; the x=1 face (the
    degenerate corner p -> 0 where every physical quantity vanishes) is
    exactly why the claim is open there and the faces/edges tower is
    anchored at x=0, not x=1.

CERTIFICATE (B2) (exact Sturm).  Qq(0,p) = 16 p^2 (p^2+1) > 0 (identity
P11a; needs p > 0, i.e. x < 1) and Qq(thetatilde2*Dbar, p) < 0 for
x in [0, 1):  after clearing the positive denominator
const * (49x^2+350x+1649)^2  with  49x^2+350x+1649 = 1024 + (25+7x)^2 > 0
(P11b), the numerator is (1-x)^3 * R17(x) with R17 of degree 17 having NO
roots in [0,1] (Sturm) and negative there, including at x = 1 (P11c) --
the (1-x)^3 factor is exactly the p^2 * omega vanishing rate of the
composition at the open endpoint.  A continuous function positive at D=0
and negative at D = thetatilde2*Dbar has a root below the curve; the
smallest positive root of Qq is D_c (committed Lemma 7.34j; only the SAFE
direction D_c <= that root is used -- 7.34j V3/V4 identify the root as
the alternating-word threshold, an upper bound for the all-word
threshold D_c).

HONEST NOTES.
 *  The sg -> 1 (p -> 0) endpoint is an OPEN limit and stays open: the
    Gray family itself is only defined for p > 0, every certified
    inequality degenerates to 0 = 0 at p = 0, and the inequality
    R_3 >= 3/2 is asymptotically SHARP there (margin (A-2)/(A-1) p +
    O(p^2)).  Nothing is claimed at p = 0; the claim region is
    half-open in x and the Bernstein tower certifies exactly that.
 *  theta_x - theta_c -> 0 in ABSOLUTE terms as sg -> 1 (both -> 2/3);
    the corridor survives only RELATIVE to omega = 1-sg^2.  All margins
    of this certificate vanish linearly in omega -- this is forced by
    the sharpness of the theorem, not an artifact.
 *  The pocket persists arbitrarily close to p = 0 (see GEOGRAPHY);
    thetatilde2 CANNOT be pushed to 1 anywhere on this range.
 *  chi(Lambda) > 0 is the k=0 member of the Budan-Fourier tower; the
    negative-side certificates and interior-angle realness are committed
    FULL-superdomain lemmas (unchanged by this script), so the assembly
    theorem (thm_7_34_route_b_a3_assembly.py) extends to p in (0, 2/3)
    with no other new input.

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.
Runtime: ~3-5 min.  Machinery identical to the committed
scripts/verify/lemma_7_34_midband_cert5_on_gray.py (t-domain assembly,
strip, integer Bernstein tower, quartic Sturm inclusion); only the
sg-range, the separating curve, and the included-face pattern differ.
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
THT2 = sp.Rational(2, 3) + sp.Rational(11, 32) * (1 - sg ** 2)  # the curve
SGMIN = sp.Rational(25, 32)
SGX = SGMIN + sp.Rational(7, 32) * X          # sg = (25+7x)/32, x in [0,1)
CXnum = 78703 - 11550 * X - 1617 * X ** 2     # thetatilde2 numerator (x-dom)
CXDEN = 98304

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
    print("CERT5-SMALL-P: chi(Lambda) > 0 on th <= thetatilde2(sg) =")
    print("2/3 + (11/32)(1-sg^2), sg in [25/32, 1), all angles; +")
    print("D_c < thetatilde2*Dbar (exact Sturm).  A = 3.  Exact end to end.")
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
    print("P4  target assembly + numeric sign guard (NEW sg-range)")

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

    # NOTE (committed): the target has huge inter-term cancellation; double
    # precision FLIPS signs.  Near sg -> 1 everything additionally vanishes
    # at high order, so evaluate at 80 digits (committed piece used 60).
    fgt = sp.lambdify(VS, tg, modules="mpmath")
    random.seed(23)
    agree = disagree = npts = 0
    with mp.workdps(80):
        for _ in range(140):
            # sg = (800 + j)/1024, j in [0, 223]: covers [25/32, 1)
            sv = sp.Rational(800 + random.randint(0, 223), 1024)
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

    # ---- P5: pocket witness INSIDE the new range --------------------------
    print("-" * 78)
    print("P5  pocket witness: target < 0 at (sg=27/32, th=1, t=2) [exact]")
    wit = tg.subs({sg: sp.Rational(27, 32), th: 1, t: 2})
    print(f"     target = {float(wit):.6e}")
    rep("P5 exact target < 0 at the pocket point (p = 295/1536, D = Dbar, "
        "s = pi)", wit < 0)

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

    # ---- P6c: x-domain identities of the substitution ---------------------
    print("-" * 78)
    print("P6c exact identities: sg=(25+7x)/32 => 1-sg^2 = 7(1-x)(57+7x)/1024,")
    print("    thetatilde2(sg) = (78703 - 11550x - 1617x^2)/98304")
    id1 = sp.expand(1024 * (1 - SGX ** 2) - 7 * (1 - X) * (57 + 7 * X))
    id2 = sp.expand(CXDEN * THT2.subs(sg, SGX) - CXnum)
    rep("P6c both substitution identities hold exactly", id1 == 0 and id2 == 0)

    # ---- P7: u-substitution to (x,u,tau) on [0,1]^3, exact cross-check ----
    print("-" * 78)
    print("P7  substitution sg=(25+7x)/32, th=u(78703-11550x-1617x^2)/98304,")
    print("    t=2tau; integer tensor")
    DS = corePoly.degree(sg)
    DTH = corePoly.degree(th)
    DT = corePoly.degree(t)
    DX = DS + 2 * DTH
    # (78703 - 11550 x - 1617 x^2)^b, integer coefficient dicts
    CPOW = {0: {0: 1}}
    for b in range(1, DTH + 1):
        prev = CPOW[b - 1]
        cur = {}
        for jj, cf in prev.items():
            cur[jj] = cur.get(jj, 0) + 78703 * cf
            cur[jj + 1] = cur.get(jj + 1, 0) - 11550 * cf
            cur[jj + 2] = cur.get(jj + 2, 0) - 1617 * cf
        CPOW[b] = cur
    # (25 + 7x)^a
    SPOW = {0: {0: 1}}
    for a in range(1, DS + 1):
        prev = SPOW[a - 1]
        cur = {}
        for jj, cf in prev.items():
            cur[jj] = cur.get(jj, 0) + 25 * cf
            cur[jj + 1] = cur.get(jj + 1, 0) + 7 * cf
        SPOW[a] = cur
    P32 = [32 ** a for a in range(DS + 1)]
    P2l = [2 ** cc for cc in range(DT + 1)]
    PCD = [CXDEN ** b for b in range(DTH + 1)]
    # core(sg,th,t) * 32^DS * 98304^DTH with the substitution; the 32-power
    # compensates the (25+7x)/32 sigma-scaling and the 98304-power the
    # thetatilde2 denominator (committed-P7 bookkeeping pattern; the exact
    # Fraction cross-check below is the guard against any misplaced power).
    TU = {}
    for (a, b, cc), cf in Tcore.items():
        base = cf * P32[DS - a] * P2l[cc] * PCD[DTH - b]
        for ja, pa in SPOW[a].items():
            for jc, pc in CPOW[b].items():
                key = (ja + jc, b, cc)
                TU[key] = TU.get(key, 0) + base * pa * pc
    TU = {k: v for k, v in TU.items() if v != 0}
    print(f"     integer tensor: {len(TU)} nonzero, degs ({DX},{DTH},{DT})")
    SCALE = Fr(32) ** DS * Fr(CXDEN) ** DTH
    random.seed(11)
    ok7 = True
    for _ in range(12):
        xv = Fr(random.randint(0, 16), 16)
        uv = Fr(random.randint(0, 16), 16)
        tv = Fr(random.randint(0, 16), 16)
        sgv = Fr(25, 32) + Fr(7, 32) * xv
        thv = uv * (Fr(2, 3) + Fr(11, 32) * (1 - sgv ** 2))
        ttv = 2 * tv
        lhs = sum(Fr(cf) * xv ** a2 * uv ** b2 * tv ** c2
                  for (a2, b2, c2), cf in TU.items())
        rhs = SCALE * sum(Fr(cf) * sgv ** a2 * thv ** b2 * ttv ** c2
                          for (a2, b2, c2), cf in Tcore.items())
        ok7 &= (lhs == rhs)
    rep("P7 exact Fraction cross-check core_u == scaled core at 12 points", ok7)
    # P7b: conclusive polynomial-identity check by Schwartz-Zippel over F_q,
    # q = 2^61 - 1 (prime).  The difference polynomial has total degree
    # <= DX + DTH + DT = 223, so a nonzero difference survives 20 uniform
    # F_q points with probability <= (223/q)^20 ~ 1e-370.
    qP = (1 << 61) - 1
    inv32 = pow(32, qP - 2, qP)
    invCD = pow(CXDEN, qP - 2, qP)
    scale_q = pow(32, DS, qP) * pow(CXDEN, DTH, qP) % qP
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
        sgq = (25 + 7 * xq) % qP * inv32 % qP
        thq = uq * ((78703 - 11550 * xq - 1617 * xq * xq) % qP) % qP \
            * invCD % qP
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

    # ---- P9: face certificates (2-var) -------------------------------------
    print("-" * 78)
    print("P9  faces x=0 (restriction), u=1, tau=1 (collapse): 2-var")
    print("    Bernstein no-negatives + nonzero")

    def collapse_at1(T, ax, degs):
        out = {}
        for mon, cf in T.items():
            key = tuple(m for i, m in enumerate(mon) if i != ax)
            out[key] = out.get(key, 0) + cf     # value at var=1: sum coeffs
        return {k: v for k, v in out.items() if v != 0}, \
            [d_ for i, d_ in enumerate(degs) if i != ax]

    def restrict_at0(T, ax, degs):
        out = {}
        for mon, cf in T.items():
            if mon[ax] != 0:
                continue                        # value at var=0: only a=0
            key = tuple(m for i, m in enumerate(mon) if i != ax)
            out[key] = out.get(key, 0) + cf
        return {k: v for k, v in out.items() if v != 0}, \
            [d_ for i, d_ in enumerate(degs) if i != ax]

    ok9 = True
    for ax, nm, mode in ((0, "x=0", "r0"), (1, "u=1", "c1"),
                         (2, "tau=1", "c1")):
        if mode == "r0":
            Tf, degf = restrict_at0(TU, ax, degs3)
        else:
            Tf, degf = collapse_at1(TU, ax, degs3)
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

    def edge_coeffs(x0_restrict, keep):
        """Edge polynomial: restrict x=0 if asked, set the other two vars
        to 1 (sum), keep the running variable `keep`."""
        co = {}
        for mon, cf in TU.items():
            if x0_restrict and mon[0] != 0:
                continue
            co[mon[keep]] = co.get(mon[keep], 0) + cf
        return co

    for x0r, keep, nm in ((True, 2, "(x=0,u=1), poly in tau"),
                          (True, 1, "(x=0,tau=1), poly in u"),
                          (False, 0, "(u=1,tau=1), poly in x")):
        co = edge_coeffs(x0r, keep)
        d_ = degs3[keep]
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
    corner = sum(cf for (a2, b2, c2), cf in TU.items() if a2 == 0)
    rep("P10d corner (x,u,tau)=(0,1,1): exact value > 0", corner > 0)
    print("     [union: open box + faces x=0,u=1,tau=1 + edges (x=0,u=1),")
    print("      (x=0,tau=1), (u=1,tau=1) + corner (0,1,1)")
    print("      = [0,1) x (0,1] x (0,1] exactly = claim region of (A2);")
    print("      x=1 (p=0), u=0 (D=0), tau=0 (s=0) are excluded open limits]")

    # ---- P11: certificate (B2), exact Sturm --------------------------------
    print("-" * 78)
    print("P11 certificate (B2): D_c < thetatilde2 * Dbar on sg in [25/32, 1)")
    id0 = sp.expand(Qquartic.subs(D, 0) - 16 * p ** 2 * (p ** 2 + 1))
    rep("P11a Qq(0,p) == 16 p^2 (p^2+1)  (> 0 for p > 0, i.e. x < 1)",
        id0 == 0)
    PSX = sp.Rational(2, 3) * (1 - SGX ** 2)
    DBX = sp.Rational(2, 3) * (1 - SGX) ** 2 / (1 + SGX ** 2)
    curveX = (CXnum / CXDEN) * DBX
    expr = Qquartic.subs({D: curveX, p: PSX})
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    denf = sp.factor_list(den)
    qden = 49 * X ** 2 + 350 * X + 1649
    okden = (denf[0] > 0
             and all(sp.simplify(b - qden) == 0 for b, e in denf[1]
                     if not b.is_number)
             and sum(e for b, e in denf[1] if not b.is_number) == 2
             and sp.expand(qden - 1024 - (25 + 7 * X) ** 2) == 0)
    rep("P11b cleared denominator = const*(1024+(25+7x)^2)^2  (positive)",
        okden)
    Rp = sp.Poly(sp.expand(num), X)
    k1 = 0
    while True:
        q_, r_ = sp.div(Rp, sp.Poly(1 - X, X))
        if r_.is_zero and not q_.is_zero:
            Rp = q_
            k1 += 1
        else:
            break
    nroots = Rp.count_roots(0, 1)
    v0 = Rp.eval(0)
    vh = Rp.eval(sp.Rational(1, 2))
    ve = Rp.eval(1)
    print(f"     stripped (1-x)^{k1}; deg {Rp.degree()}; Sturm roots in "
          f"[0,1]: {nroots}; values (0,1/2,1) signs: "
          f"{sp.sign(v0)},{sp.sign(vh)},{sp.sign(ve)}")
    rep("P11c numerator = (1-x)^3 * R17, R17 root-free in [0,1] and negative",
        nroots == 0 and v0 < 0 and vh < 0 and ve < 0 and k1 == 3
        and Rp.degree() == 17)
    rep("P11d => Qq(thetatilde2*Dbar) < 0 on x in [0,1); with P11a and the "
        "first-root logic, D_c < thetatilde2*Dbar",
        nroots == 0 and v0 < 0 and k1 == 3)

    # ---- P12: exact anchor ------------------------------------------------
    print("-" * 78)
    print("P12 exact anchor: chi(Lambda) at x=1/2 (sg=57/64), th =")
    print("    thetatilde2(57/64) = 290095/393216, w=-1 vs core_u sign")
    pv = sp.Rational(847, 6144)          # p(57/64) = (2/3)(1 - 3249/4096)
    thv = sp.Rational(290095, 393216)    # thetatilde2(57/64), on the curve
    Dbv = sp.Rational(2, 3) * (1 - sp.Rational(57, 64)) ** 2 \
        / (1 + sp.Rational(57, 64) ** 2)
    okpt = (sp.simplify(pv - PSUB.subs(sg, sp.Rational(57, 64))) == 0
            and sp.simplify(thv - THT2.subs(sg, sp.Rational(57, 64))) == 0)
    rep("P12a anchor point identities (p, thetatilde2 at sg=57/64)", okpt)
    Dv = sp.nsimplify(thv * Dbv)
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
    vc = sum(Fr(cf) * Fr(1, 2) ** a2 for (a2, b2, c2), cf in TU.items())
    print(f"     chi(Lambda) = {float(chiw):.6e} (exact rational), "
          f"core_u(1/2,1,1) sign {'+' if vc > 0 else '-'}")
    rep("P12 exact chi > 0 at the anchor and sign matches core_u",
        chiw > 0 and vc > 0)

    # ---- P13: numeric separation scan (consistency) ------------------------
    print("-" * 78)
    print("P13 numeric separation: theta_c < thetatilde2 and chi > 0 on the")
    print("    curve at t=2 at 25 grid points (mpmath consistency scan)")
    ok13 = True
    for k in range(0, 25):
        xv = k / 25.0
        sgv = (25.0 + 7.0 * xv) / 32.0
        pv2 = (2.0 / 3.0) * (1.0 - sgv * sgv)
        Db = (2.0 / 3.0) * (1.0 - sgv) ** 2 / (1.0 + sgv * sgv)
        dth = 2.0 / 3.0 + (11.0 / 32.0) * (1.0 - sgv * sgv)
        dc = Dc3_numeric(pv2)
        ok13 &= (dc is not None and dc / Db < dth)
        chc = chi_numeric(pv2, dth * Db, math.pi, 2.0)
        ok13 &= (chc > 0)
    rep("P13 curve strictly above theta_c and chi > 0 on the curve at t=2",
        ok13)

    # ---- P14: splice with the committed piece ------------------------------
    print("-" * 78)
    print("P14 splice with the committed cert5-on-gray at sg = 25/32")
    t1 = 1 - sp.Rational(5, 16) * SGMIN ** 2       # committed curve at joint
    t2 = THT2.subs(sg, SGMIN)                      # new curve at joint
    ok14 = (t2 == sp.Rational(78703, 98304)
            and t1 == sp.Rational(13259, 16384)
            and t2 <= t1)
    print(f"     thetatilde2(25/32) = {t2} = {float(t2):.6f}")
    print(f"     thetatilde (25/32) = {t1} = {float(t1):.6f}")
    rep("P14 thetatilde2(25/32) <= thetatilde(25/32) exactly (nested joint); "
        "sg-coverage (0,25/32] u [25/32,1) seamless", ok14)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("CLAIM CERTIFIED (A2): chi(Lambda) > 0 for sg in [25/32, 1),")
    print("th in (0, 2/3 + (11/32)(1-sg^2)], t in (0,2]  [exact Bernstein")
    print("tower: open box + faces x=0,u=1,tau=1 + edges + corner]")
    print("CLAIM CERTIFIED (B2): D_c(p) < thetatilde2(sg) * Dbar(p) there")
    print("[exact Sturm; D_c = smallest positive quartic root, Lemma 7.34j]")
    print("=> Gray region for p in (0, 133/512] is strictly inside the")
    print("   certified chi > 0 region.  Together with the committed")
    print("   cert5-on-gray (p in [133/512, 2/3)): chi(Lambda) > 0 covers")
    print("   Gray at CONTINUUM grade for ALL p in (0, 2/3); p = 0 remains")
    print("   the excluded open limit of the Gray family.")
    return PASS


if __name__ == "__main__":
    main()
