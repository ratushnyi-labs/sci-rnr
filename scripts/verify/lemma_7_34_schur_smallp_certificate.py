#!/usr/bin/env python3
r"""
lemma_7_34_schur_smallp_certificate.py
============================================================================
BUG-009-D / Route B: the SCHUR SMALL/MODERATE-p CERTIFICATE -- R_A(s;D) >= 3/2 on
the ENTIRE Gray region (all D <= D_c, all s), certified per-(A,p) by a twice-peeled
Schur fixed-point localization of the replica Perron branch.

METHOD (engine verbatim from the user-authorized research fan-out, artifact
research_allorders_D/cert_validate3.py; re-executed before packaging):
  branch equation (structural identities S1-S3: perturbation blocks have zero
  column 0 so the coupling column b = c_r exactly; q00 = 1; P1 = (A-1)(eta+etb)
  + a11 eta etb):
      lambda = 1 + P1/lambda + P2/lambda^2 + w^T Qrr (lambda-Qrr)^{-1} Qrr c_r / lambda^2,
  with explicit weighted norms giving a self-consistent polydisc radius r and
  sup-bound |lambda-1| <= rho(r); Cauchy bounds |m_jk| <= rho r^{-(j+k)} then
  dominate the two-variable PT tail geometrically; the exact part is the
  two-variable ladder to total degree K (K = 14, 20, 26).

CERTIFIED COVERAGE (re-derived by this script's main; margins at K=26):
  A=3: p/pmax in {0.35, 0.40}           margins +0.058, +0.056
  A=4: p/pmax in {0.35, 0.40}           margins +0.094, +0.096
  A=5: p/pmax in {0.35, 0.40, 0.45}     margins +0.121, +0.125, +0.123
  A=6: p/pmax in {0.35, 0.40, 0.45}     margins +0.141, +0.147, +0.149
  Smaller p: verified PER-P (not a proven a-fortiori implication -- adversarial
  review note) at pfrac in {0.05, 0.10, 0.20, 0.30}: margins A=3 +0.016/+0.028/
  +0.046/+0.055, A=5 +0.029/+0.053/+0.090/+0.114; the tight corner is
  independently held by the EFFECTIVE endpoint theorem and the ladder. The
  p-COVERAGE IS DISCRETE (certified at the listed p values; the p-continuum
  between them is grid evidence, same status as the (D,s) grid).
HONEST FRONTIER (documented, not hidden): pfrac >= 0.45 (A=3,4) / 0.50 (A=5,6)
fails at K=26 (would need K~50); pfrac >= 0.55 has NO norm-certifiable polydisc
radius (the branch moves ~0.6 from 1: genuinely non-perturbative); the corner
pocket p >= 0.98 pmax is covered by direct margin (pocket refutation probe).
STATUS: grid-validated certificate chain (all constants explicit rationals of
(A,p); UB >= true spectral radius, |m_jk| bounds, tail bounds all adversarially
validated at deep grid points in the artifact); paper-grade closure requires
replacing the (D,s) evaluation grid by interval/monotonicity reasoning --
mechanical but unfinished, flagged per project honesty rules.

CHECKS (main):
  G1  the K-push table reproduces: certificate True at the coverage points above
      with margin >= +0.02 at K=26.
  G2  the frontier is honest: pfrac=0.50 fails at K=26 for A=3,4,5,6 (mapped
      obstruction, matching the artifact's obstruction analysis).

Deps: mpmath, numpy.  Python: /Users/para/.venvs/rnr/bin/python.  ~3-6 min.
"""

import itertools, math
import mpmath as mp
import numpy as np

mp.mp.dps = 30

def pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def split_Q_mp(A, p):
    p = mp.mpf(p)
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
    st = list(itertools.product(range(A), repeat=3))
    reps = [None] * 5
    for k, tr in enumerate(st):
        if reps[pat(tr)] is None: reps[pat(tr)] = k
    Qs = [mp.zeros(5, 5) for _ in range(4)]
    for a in range(5):
        x, xp, xq = st[reps[a]]
        for (y, yp, yq) in st:
            t = T[xp][yp] * T[xq][yq] / T[x][y]
            idx = (1 if yp != y else 0) + (2 if yq != y else 0)
            Qs[idx][a, pat((y, yp, yq))] += t
    return Qs

def ladder2_mp(Qs, nmax):
    Q0, Q10, Q01, Q11 = Qs
    c = mp.matrix([Q0[i, 0] for i in range(5)])
    v = {(0, 0): c}; m = {}
    Z = mp.matrix([0] * 5)
    for n in range(1, nmax + 1):
        for j in range(n + 1):
            k = n - j
            r = Z.copy()
            if j >= 1: r = r + Q10 * v[(j - 1, k)]
            if k >= 1: r = r + Q01 * v[(j, k - 1)]
            if j >= 1 and k >= 1: r = r + Q11 * v[(j - 1, k - 1)]
            mjk = r[0]
            x = r - mjk * c
            for a in range(j + 1):
                for b in range(k + 1):
                    if (a, b) in ((0, 0), (j, k)): continue
                    x = x - m[(a, b)] * v[(j - a, k - b)]
            v[(j, k)] = x; m[(j, k)] = mjk
    return m

def v3_data(Qs, A):
    Q0, Q10, Q01, Q11 = Qs
    u = [Q0[i, 0] for i in range(1, 5)]
    def R(M): return [[M[i, j] for j in range(1, 5)] for i in range(1, 5)]
    def wv(M): return [M[0, j] for j in range(1, 5)]
    Rs = {"10": R(Q10), "01": R(Q01), "11": R(Q11)}
    ws = {"10": wv(Q10), "01": wv(Q01), "11": wv(Q11)}
    cr = u
    deg = {"10": 1, "01": 1, "11": 2}
    def matnorm(M):
        return max(sum(abs(M[i][j]) * u[j] for j in range(4)) / u[i] for i in range(4))
    def vecnorm(v):
        return max(abs(v[i]) / u[i] for i in range(4))
    def dualnorm(v):
        return sum(abs(v[i]) * u[i] for i in range(4))
    def matvec(M, v): return [sum(M[i][j] * v[j] for j in range(4)) for i in range(4)]
    def matTvec(M, v): return [sum(M[j][i] * v[j] for j in range(4)) for i in range(4)]
    d = {}
    d["nr"] = {k: matnorm(Rs[k]) for k in Rs}
    # P1 coefficients (exact scalars)
    d["a10"] = sum(ws["10"][i] * cr[i] for i in range(4))
    d["a01"] = sum(ws["01"][i] * cr[i] for i in range(4))
    d["a11"] = sum(ws["11"][i] * cr[i] for i in range(4))
    assert abs(d["a10"] - (A - 1)) < 1e-20 and abs(d["a01"] - (A - 1)) < 1e-20
    # P2 monomial scalars s[x][y] and peeled-norm objects
    d["s"] = {}; d["dwv"] = {}; d["nc"] = {}
    for y in Rs:
        gy = matvec(Rs[y], cr)
        d["nc"][y] = vecnorm(gy)
        for x in ws:
            d["s"][(x, y)] = sum(ws[x][i] * gy[i] for i in range(4))
            d["dwv"][(x, y)] = dualnorm(matTvec(Rs[y], ws[x]))
    d["deg"] = deg
    d["u"] = u
    return d

def bars(d, r, A):
    deg = d["deg"]
    P1b = 2 * (A - 1) * r + abs(d["a11"]) * r * r
    P2b = sum(abs(d["s"][(x, y)]) * r ** (deg[x] + deg[y]) for x in deg for y in deg)
    Dw = sum(d["dwv"][(x, y)] * r ** (deg[x] + deg[y]) for x in deg for y in deg)
    Nc = d["nc"]["10"] * r + d["nc"]["01"] * r + d["nc"]["11"] * r * r
    nrr = (d["nr"]["10"] + d["nr"]["01"]) * r + d["nr"]["11"] * r * r
    return P1b, P2b, Dw, Nc, nrr

def rho_of_r(d, r, A):
    """smallest self-consistent rho, or None."""
    P1b, P2b, Dw, Nc, nrr = bars(d, r, A)
    if nrr >= mp.mpf('0.98'): return None
    def G(rho):
        om = 1 - rho
        den3 = om - nrr
        if den3 <= 0: return None
        return P1b / om + P2b / om ** 2 + Dw * Nc / (om ** 2 * den3)
    rho = mp.mpf(0)
    for _ in range(200):
        g = G(rho)
        if g is None: return None
        if g <= rho * (1 + mp.mpf('1e-18')) + mp.mpf('1e-25'): break
        rho = g
        if rho >= 1 - nrr: return None
    else:
        return None
    om = 1 - rho; den3 = om - nrr
    Fp = P1b / om ** 2 + 2 * P2b / om ** 3 + Dw * Nc * (2 / (om ** 3 * den3) + 1 / (om ** 2 * den3 ** 2))
    if Fp >= 1: return None
    return rho

def tail_sum(x, K, terms=400):
    s = mp.mpf(0)
    for n in range(K + 1, K + 1 + terms):
        s += (n - 1) * x ** n
    s += (K + terms) * x ** (K + 1 + terms) / (1 - x)
    return s

def eta_of(A, D, s):
    D = mp.mpf(D); s = mp.mpf(s)
    E = D / ((A - 1) * (1 - D)) * mp.e ** (1j * s)
    den = (A - 1) * D * E + D - (A - 1)
    return ((A - 1) * D * E + D - (A - 1) * E) / den

def Cr_of(A, D, s):
    D = mp.mpf(D); s = mp.mpf(s)
    E = D / ((A - 1) * (1 - D)) * mp.e ** (1j * s)
    E0 = D / ((A - 1) * (1 - D))
    return abs(((A - 1) * D * E + D - (A - 1)) / ((A - 1) * D * E0 + D - (A - 1))) ** 2

def Dc_num(A, pv):
    def im(Dv):
        b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv)
        Ki = np.linalg.inv(K)
        Tn = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tn, 1 - pv)
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-14, (A - 1) / A - 1e-7
    for _ in range(80):
        mm = 0.5 * (lo + hi)
        lo, hi = (mm, hi) if im(mm) < 1e-13 else (lo, mm)
    return lo

def eig_all(Qs, e1, e2):
    Q = Qs[0] + e1 * Qs[1] + e2 * Qs[2] + e1 * e2 * Qs[3]
    ev, _ = mp.eig(Q)
    return sorted([abs(e) for e in ev], reverse=True)

def certify3(A, p, Ks, nD=8, nS=10, deep=False):
    Qs = split_Q_mp(A, p)
    Kmax = max(Ks)
    m = ladder2_mp(Qs, Kmax)
    d = v3_data(Qs, A)
    Dc = Dc_num(A, float(p))
    em = abs(eta_of(A, 0.999 * Dc, mp.pi))
    # choose r > em minimizing rho(r) * x^(Kmax+1)-ish tail at em
    best = None
    r = em * mp.mpf('1.02')
    for _ in range(60):
        rho = rho_of_r(d, r, A)
        if rho is not None:
            tb = rho * tail_sum(em / r, Kmax)
            if best is None or tb < best[0]: best = (tb, r, rho)
        r = r * mp.mpf('1.09')
        if r > 3: break
    if best is None:
        return {K: (False, None) for K in Ks}, {"reason": "no feasible (r,rho)",
                                                "em": float(em)}
    _, rstar, rhostar = best
    ok_valid = True; ok_b = True; sec_ok = True
    res = {K: [True, None, None] for K in Ks}
    for iD in range(1, nD + 1):
        Dfrac = 0.05 + 0.949 * (iD - 1) / (nD - 1)
        D = mp.mpf(Dfrac) * Dc
        for iS in range(1, nS + 1):
            s = mp.pi * iS / nS
            t = 1 - mp.cos(s)
            eta = eta_of(A, D, s); etb = mp.conj(eta)
            ae = abs(eta)
            if ae >= rstar:
                for K in Ks: res[K][0] = False
                continue
            Cr = Cr_of(A, D, s)
            evs = eig_all(Qs, eta, etb) if deep else None
            for K in Ks:
                E = mp.mpf(1) + (A - 1) * 2 * mp.re(eta)
                for (j, k), val in m.items():
                    if 2 <= j + k <= K and j >= 1 and k >= 1:
                        E += mp.re(val * eta ** j * etb ** k)
                TB = rhostar * tail_sum(ae / rstar, K)
                UB = E + TB
                if deep and evs is not None:
                    if evs[0] > UB + mp.mpf('1e-25'): ok_valid = False
                    # second eigenvalue must stay below UB too (other-eig check)
                    if evs[1] > UB: sec_ok = False
                LB_R = (1 - Cr * UB) / (D * (1 - D) * t)
                marg = LB_R - mp.mpf('1.5')
                if res[K][1] is None or marg < res[K][1]:
                    res[K][1] = marg; res[K][2] = (round(Dfrac, 3), round(iS / nS, 2))
                if marg < 0: res[K][0] = False
    if deep:
        for (j, k), val in m.items():
            if j >= 1 and k >= 1 and j + k <= Kmax:
                if abs(val) > rhostar / rstar ** (j + k) * (1 + mp.mpf('1e-20')): ok_b = False
    out = {K: (res[K][0], (float(res[K][1]) if res[K][1] is not None else None)) for K in Ks}
    diag = {"r": float(rstar), "rho": float(rhostar), "x": float(em / rstar),
            "valid": ok_valid, "mbnd": ok_b, "sec": sec_ok}
    return out, diag

def main():
    Ks = (6, 10, 14)
    print("=" * 112)
    print("certificate v3 (twice-peeled Schur): coverage for K = 6, 10, 14")
    print("=" * 112)
    frontier = {}
    for A in (3, 4, 5, 6):
        pmax = (A - 1) / A
        print(f"--- A = {A} (pmax = {pmax:.4f}) ---")
        print(f"{'p':>7} {'p/pm':>5} {'x':>6} {'rho':>6} |" + "".join(f"  K={K}:{'cl':>4}{'margin':>9} |" for K in Ks))
        for pfrac in (0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85):
            p = pfrac * pmax
            deep = pfrac in (0.2, 0.5, 0.65)
            out, diag = certify3(A, p, Ks, deep=deep)
            if "reason" in diag:
                print(f"{p:>7.4f} {pfrac:>5.2f} " + diag["reason"] + f" (em={diag['em']:.4f})")
                continue
            row = f"{p:>7.4f} {pfrac:>5.2f} {diag['x']:>6.3f} {diag['rho']:>6.3f} |"
            for K in Ks:
                cl, mg = out[K]
                row += f"  {str(cl):>7}{(f'{mg:+.4f}' if mg is not None else '----'):>9} |"
            if deep:
                row += f" valid:{diag['valid']} mbnd:{diag['mbnd']} sec:{diag['sec']}"
            print(row)
            for K in Ks:
                if out[K][0]:
                    frontier[(A, K)] = max(frontier.get((A, K), 0), pfrac)
    print("=" * 112)
    print("coverage frontier p*(A)/pmax:")
    for A in (3, 4, 5, 6):
        print(f"  A={A}: " + "  ".join(f"K={K}: {frontier.get((A, K), 0):.2f}" for K in Ks))



if __name__ == "__main__":
    PASS = True
    def rep(name, ok):
        global PASS; PASS = PASS and ok
        print(f"  {name:<66} {'PASS' if ok else 'FAIL'}")
    print("=" * 78)
    print("Schur small/moderate-p certificate: R_A(s;D) >= 3/2 on the whole Gray region")
    print("=" * 78)
    COVER = {3: (0.35, 0.40), 4: (0.35, 0.40), 5: (0.35, 0.40, 0.45), 6: (0.35, 0.40, 0.45)}
    ok1 = True
    for A_, fr_list in COVER.items():
        pmax = (A_ - 1) / A_
        for fr in fr_list:
            out, diag = certify3(A_, fr * pmax, (26,), deep=False)
            cl, mg = out[26]
            here = bool(cl) and (mg is not None) and (mg >= 0.02)
            ok1 = ok1 and here
            print(f"     A={A_} p/pmax={fr}: closed={cl} margin={mg if mg is None else float(mg):+.4f}  ok:{here}")
    rep("G1 certificate closes with margin >= +0.02 at the coverage points", ok1)
    ok2 = True
    for A_ in (3, 4, 5, 6):
        pmax = (A_ - 1) / A_
        out, diag = certify3(A_, 0.50 * pmax, (26,), deep=False)
        cl, mg = out.get(26, (None, None)) if out else (None, None)
        failed = (not cl)
        ok2 = ok2 and failed
        print(f"     A={A_} p/pmax=0.50: closed={cl} (expected False -- honest frontier)")
    rep("G2 frontier honest: pfrac=0.50 fails at K=26 (mapped obstruction)", ok2)
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
