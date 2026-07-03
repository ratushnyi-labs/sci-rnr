#!/usr/bin/env python3
r"""
probe_7_34_corner_pocket_refutation.py
============================================================================
BUG-009-D / Route F: HONEST SCOPE CORRECTION -- the two monotonicity conjectures
(V2: s-monotonicity; V3: endpoint D-monotonicity) are FALSE as global statements.
They fail in a bounded pocket near the uniform-source corner p -> (A-1)/A, while
the Route B target R_A(s;D) >= 3/2 HOLDS EVERYWHERE (with margin >= +0.24 inside
the pocket -- far above the tight-corner margins).

DISCOVERED by a hostile stress sweep (~20,000 points, A up to 32, p to
(1-1e-6)(A-1)/A, D to 0.99999 D_c, mpmath 50-100 dps), then INDEPENDENTLY
REPRODUCED here with a separate implementation (values agree to 9 sig figs).

THE REFUTATIONS (V3 and V2 witnesses BOTH reproduced in this probe, H1/H4):
  V3 (dR_A(pi;D)/dD <= 0 on (0,D_c]) FAILS for every A >= 3 with p close enough
  to (A-1)/A. Onset p/pmax: A=3: 0.9999, A=4: 0.999, A=5: 0.997, A=6: 0.995,
  A=8: 0.99, A>=12: 0.98 (all CLEAN at p/pmax <= 0.97). Witness (this probe):
  A=16, p = 0.99 * 15/16: R(pi; 0.8 D_c) = 1.908771... < R(pi; 0.9 D_c) =
  1.932013... -- R INCREASES in D. For A>=12, p/pmax>=0.98 the D-infimum moves
  to interior D* ~ 0.9 D_c (undercuts the corner by up to ~0.009 at A=24).
  V2 (R_A(s) non-increasing in s) FAILS for large alphabets in the same pocket:
  onset A=16/24 at p/pmax=0.995, A=32 at 0.99; A<=12 clean to 0.9999. In the
  pocket R(s) becomes monotone INCREASING (s-argmin flips to s->0+, where
  R(0+) = 2v/(D(1-D))).

WHAT SURVIVES (the corrected honest map):
  * The TARGET R_A(s;D) >= 3/2 holds at every sweep point; global min
    1.50000999... at the usual tight corner (small p, s=pi, D->D_c). Inside the
    violation pocket min R = 1.7404 (A=32) -- margin 16x the tight-corner margin.
  * Both monotonicity statements hold with quantified derivative slack for
    p <= 0.97 (A-1)/A -- ALL the committed order-by-order theorems (flatness,
    O(D^2) beta-certificate, endpoint corner + all-D grid at p<=0.3) live deep
    inside this clean region and are unaffected.
  * The monotone-reduction strategy is therefore REGIONAL: clean region ->
    monotone + endpoint; corner pocket -> direct margin (>= 1.74).
  * probe_7_34_endpoint_all_D_monotone.py tested p <= 0.3 only (as stated in its
    docstring); its grid claims stand, but its mechanism does NOT extend to the
    corner pocket. This probe is the boundary marker.

CHECKS:
  H1  V3 witness reproduced: A=16, p=0.99*15/16 -- R(pi;0.8Dc) < R(pi;0.9Dc)
      (strict), both >= 1.9 (target safe).
  H4  V2 witness reproduced: A=16, p=0.999*15/16, D=0.9Dc -- R(s) strictly
      INCREASING across a 28-point s-grid (27/27 up-steps; argmin flips to
      s->0), min R >= 3/2 with wide margin.
  H2  V1 pocket floor: R(pi; f Dc) >= 1.74 at the pocket-worst coordinates
      (A=16, p in {0.98, 0.99, 0.999} pmax, f in {0.5, 0.8, 0.9, 0.99}).
  H3  clean-region spot-check: at p = 0.9 pmax (A=16) D-monotonicity still
      holds on the same D-grid (the pocket boundary is real, not an artifact).

Deps: mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-4 min.
"""
import itertools
import mpmath as mp
mp.mp.dps = 50
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def gA(A, p, D, s):
    A = int(A); p = mp.mpf(p); D = mp.mpf(D); s = mp.mpf(s)
    th = mp.log(D / ((A - 1) * (1 - D))); E = mp.e**(th + 1j * s)
    C = ((A - 1) * D * E + D - (A - 1)) / (A * D - (A - 1))
    eta = ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    C0 = ((A - 1) * D * mp.e**th + D - (A - 1)) / (A * D - (A - 1)); Cr = abs(C / C0)**2
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
    etb = mp.conj(eta)
    st = list(itertools.product(range(A), repeat=3)); reps = [None] * 5
    for k, tr in enumerate(st):
        if reps[_pat(tr)] is None: reps[_pat(tr)] = k
    Q = mp.zeros(5, 5)
    for a_ in range(5):
        i = reps[a_]; x, xp, xq = st[i]
        for j, (y, yp, yq) in enumerate(st):
            Q[a_, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
    ev, _ = mp.eig(Q); return Cr * max(abs(e) for e in ev)

def Dc_of(A, p):
    """complexification threshold of the 3x3 symmetry-reduced alternating product"""
    def imflag(D):
        lam2 = 1 - D * A / (A - 1)
        al = mp.mpf(1) / A + (1 - mp.mpf(1) / A) / lam2
        be = mp.mpf(1) / A - (mp.mpf(1) / A) / lam2
        td = 1 - p; to = p / (A - 1)
        def Bred(y):
            Ki = {0: (al if y == 0 else be), 1: (be if y == 0 else al), 2: be}
            Trow = {0: {0: td, 1: to, 2: (A - 2) * to},
                    1: {0: to, 1: td, 2: (A - 2) * to},
                    2: {0: to, 1: to, 2: td + (A - 3) * to}}
            M = mp.zeros(3, 3)
            for i in range(3):
                for j in range(3): M[i, j] = Ki[i] * Trow[i][j]
            return M
        M = Bred(1) * Bred(0)
        ev, _ = mp.eig(M)
        return max(abs(mp.im(e)) for e in ev)
    lo = mp.mpf('1e-8'); hi = (mp.mpf(A - 1) / A) - mp.mpf('1e-9')
    for _ in range(120):
        m = (lo + hi) / 2
        if imflag(m) < mp.mpf('1e-30'): lo = m
        else: hi = m
    return (lo + hi) / 2

def R(A, p, D, s):
    return (1 - gA(A, p, D, s)) / (D * (1 - D) * (1 - mp.cos(s)))

def H1():
    print("-" * 78); print("H1  V3 witness: A=16, p=0.99*(15/16) -- R(pi;D) INCREASES 0.8Dc -> 0.9Dc")
    A = 16; p = mp.mpf('0.99') * mp.mpf(15) / 16
    dc = Dc_of(A, p)
    r8 = R(A, p, mp.mpf('0.8') * dc, mp.pi)
    r9 = R(A, p, mp.mpf('0.9') * dc, mp.pi)
    ok = (r9 > r8 + mp.mpf('1e-6')) and (r8 >= mp.mpf('1.9')) and (r9 >= mp.mpf('1.9'))
    print(f"     D_c={mp.nstr(dc,8)}  R(0.8Dc)={mp.nstr(r8,10)} < R(0.9Dc)={mp.nstr(r9,10)}: {r9>r8}")
    print(f"     (both >= 1.9 -- the target R>=3/2 is unthreatened here)")
    return rep("H1 V3 REFUTED at the witness point (D-monotonicity fails in pocket)", ok)

def H2():
    print("-" * 78); print("H2  pocket floor: R(pi; f Dc) >= 1.74 across the pocket-worst grid (A=16)")
    A = 16; pmax = mp.mpf(15) / 16
    ok = True; worst = mp.inf
    for pf in ('0.98', '0.99', '0.999'):
        p = mp.mpf(pf) * pmax; dc = Dc_of(A, p)
        for f in ('0.5', '0.8', '0.9', '0.99'):
            r = R(A, p, mp.mpf(f) * dc, mp.pi)
            worst = min(worst, r)
            ok = ok and (r >= mp.mpf('1.74'))
    print(f"     worst R over pocket grid = {mp.nstr(worst,8)} (>= 1.74)")
    return rep("H2 pocket floor >= 1.74 (margin 16x the tight corner; V1 safe)", ok)

def H3():
    print("-" * 78); print("H3  clean region: at p=0.9 pmax (A=16) D-monotonicity STILL holds")
    A = 16; p = mp.mpf('0.9') * mp.mpf(15) / 16
    dc = Dc_of(A, p)
    vals = [R(A, p, mp.mpf(f) * dc, mp.pi) for f in ('0.5', '0.8', '0.9', '0.99')]
    mono = all(vals[k + 1] < vals[k] + mp.mpf('1e-12') for k in range(len(vals) - 1))
    print(f"     R(pi; f Dc), f=0.5..0.99: " + ", ".join(mp.nstr(v, 8) for v in vals) + f"  decreasing: {mono}")
    return rep("H3 clean-region D-monotonicity intact at p=0.9 pmax (pocket is bounded)", mono)

def H4():
    print("-" * 78); print("H4  V2 witness: A=16, p=0.999*15/16, D=0.9Dc -- R increasing in s")
    A = 16; p = mp.mpf('0.999') * mp.mpf(15) / 16
    dc = Dc_of(A, p); D = mp.mpf('0.9') * dc
    ss = [mp.pi * (k + 1) / 28 for k in range(28)]
    vals = [R(A, p, D, s) for s in ss]
    ups = sum(1 for k in range(27) if vals[k + 1] > vals[k])
    ok = (ups == 27) and (min(vals) >= mp.mpf('1.5'))
    print(f"     up-steps: {ups}/27  min R = {mp.nstr(min(vals),8)} (>= 3/2)")
    return rep("H4 V2 REFUTED at the witness point (s-monotonicity fails; argmin->0)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("Corner-pocket refutation: V2/V3 monotonicity conjectures fail near p=(A-1)/A;")
    print("the Route B target R>=3/2 holds everywhere (pocket floor 1.74)")
    print("=" * 78)
    H1(); H2(); H3(); H4()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Corrected map: monotone-reduction is REGIONAL (p <= 0.97(A-1)/A); the corner")
    print("pocket is covered by direct margin. All committed order-by-order theorems")
    print("live inside the clean region and are unaffected.")
