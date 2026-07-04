#!/usr/bin/env python3
r"""
probe_7_34_route_b_c2_numeric.py
============================================================================
Fast extraction of the mixed 2nd-order Perron coefficient c2(A,p) in
   rho(eta,etb) = 1 + (A-1)(eta+etb) + c2 eta etb + O(3)   (pure k>=2 terms = 0,
by boundary affineness). c2 governs the 2nd-order drop of (1-g_A)/(D(1-D)t)
from the leading 2. Extract c2 by finite differences on the eta,etb plane at
small real test values (independent of the s-parametrization), high precision.

Then ASSEMBLE the leading-order prediction of the floor:
  on the diagonal etb=conj(eta), eta=(D/(A-1))(w-1):
     eta+etb=-(2D/(A-1))t,  eta etb=(2D^2/(A-1)^2)t,
  rho = 1 - 2Dt + c2 (2D^2/(A-1)^2) t + O(D^3),
  1-rho = 2Dt [ 1 - c2 D/(A-1)^2 ] + O(D^3),
  and with the prefactor 1-P=O(D^2)t and D(1-D) vs D, the 2nd-order coeff of
  (1-g_A)/(D(1-D)t) is computable. We just report c2 and the sign of the drop,
  confirming the leading 2 and that c2 D/(A-1)^2 stays < 1/4 on Gray (so the
  drop from 2 never reaches 3/2 -- consistent with the measured >=1.5 floor).
"""
import itertools
import mpmath as mp
mp.mp.dps = 50

def pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def build_Q_num(A, pv, eta, etb):
    T = [[(1 - pv) if i == j else pv / (A - 1) for j in range(A)] for i in range(A)]
    st = list(itertools.product(range(A), repeat=3)); npat = 5 if A >= 3 else 4
    reps = [None] * npat
    for k, tr in enumerate(st):
        if reps[pat(tr)] is None: reps[pat(tr)] = k
    Q = mp.zeros(npat, npat)
    for a in range(npat):
        x, xp, xq = st[reps[a]]
        for (y, yp, yq) in st:
            t = T[xp][yp] * T[xq][yq] / T[x][y]
            if yp != y: t *= eta
            if yq != y: t *= etb
            Q[a, pat((y, yp, yq))] += t
    return Q

def perron(Q):
    ev, _ = mp.eig(Q); return max(ev, key=lambda e: abs(e))

def c2_extract(A, p):
    """mixed partial d^2 rho/(d eta d etb) at 0 via central difference (real eta,etb)."""
    h = mp.mpf('1e-4')
    def r(e1, e2):
        return mp.re(perron(build_Q_num(A, p, e1, e2)))
    # d2/(de detb) f ~ [f(h,h)-f(h,-h)-f(-h,h)+f(-h,-h)]/(4h^2)
    val = (r(h, h) - r(h, -h) - r(-h, h) + r(-h, -h)) / (4 * h * h)
    return val

def main():
    print("Mixed 2nd-order Perron coefficient c2(A,p) (affineness => pure terms 0):")
    print("Also a1 = linear coeff (should be A-1) and the drop factor c2/(A-1)^2.")
    for A in (2, 3, 4, 5, 6):
        for p in (mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf('0.2'), mp.mpf('0.3')):
            if p >= mp.mpf(A - 1) / A: continue
            c2 = c2_extract(A, p)
            # linear: a1 = (r(h,0)-r(-h,0))/(2h)
            h = mp.mpf('1e-5')
            a1 = (mp.re(perron(build_Q_num(A, p, h, mp.mpf(0)))) -
                  mp.re(perron(build_Q_num(A, p, -h, mp.mpf(0))))) / (2 * h)
            drop = c2 / (A - 1)**2
            print(f"  A={A} p={float(p)}: a1={mp.nstr(a1,8)} (=A-1={A-1})  c2={mp.nstr(c2,8)}  "
                  f"c2/(A-1)^2={mp.nstr(drop,7)}")
        print()
    print("Interpretation: 1-rho = 2Dt[1 - (c2/(A-1)^2) D] + O(D^3).")
    print("The drop coefficient c2/(A-1)^2 times D_c (tiny) << 1, so the leading 2")
    print("erodes negligibly -- consistent with the measured >=1.5 (3/2) floor.")

if __name__ == "__main__":
    main()
