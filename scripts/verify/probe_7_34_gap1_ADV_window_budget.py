#!/usr/bin/env python3
"""
ADVERSARIAL (1b)+(4): window-size self-consistency and vertical-end budget.

The draft imposes several constraints on the fixed window half-width eps:
  (W1) eps <= delta_KP,            delta_KP >= r_KP/A          (KP window, real axis)
  (W2) eps <= r_KP/(8 A)           (strip containment |eta_z|<=r_KP, the "4A" route)
  (W3) eps <= c' vbar/C3           (zone-2 quadratic domination -log|Phi|>=v s^2/4)
and y_max = (1/2) r_KP/(4A) = r_KP/(8A).
Claim in draft: eps = Theta(1) is "fixed".  With r_KP ~ r_lb = p^2/(12 e (1-p)^2),
A = D(1-D)/(1-2D), check the ACTUAL numeric eps and whether the saddle |s0| is
interior (|s0| <= y_max) and the vertical budget closes.

Also ATTACK (4): the vertical exponent
   -1/4 vbar n eps^2 + C3 n |s0| + t|s0|
must be negative.  C3 n |s0| = C3 sqrt(n) |beta|/sqrt(vbar) GROWS like sqrt(n);
zone-2 is Theta(n).  Check the crossover n and whether it holds for accessible n
AND in the asymptotic limit.
"""
import math
import numpy as np


def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))


def A_of_D(D):
    return D * (1 - D) / (1 - 2 * D)


def r_lb(p):
    return p ** 2 / (12 * math.e * (1 - p) ** 2)


def vbar_closed(p, D):
    kappa2 = (1 - 2 * p) ** 2 / (p ** 2 * (1 - p) ** 2)
    Psi0 = kappa2 * D * (1 - D) / (1 - 2 * D) ** 2
    return D * (1 - D) * (1 - Psi0)


# C3: third-Cauchy constant.  From probe_7_34_gap1_saddle_strip.py measurement
# |Phi'''(0)|/n ~ 0.068 at p=0.4.  But the Cauchy bound C3 is a SUP over the
# half-window, typically a few x the central value.  We TEST with measured/n and
# a conservative C3 multiplier.
def measure_C3_over_n(p, D, n=200):
    # quick central third-derivative magnitude per site (annealed proxy: use a
    # typical shell word).  We reuse the recursion from the contour probe inline.
    import mpmath as mp
    mp.mp.dps = 40

    def eta_closed(s):
        w = mp.e ** (1j * s)
        return D * (1 - D) * (w - 1) / ((1 - D) ** 2 - D ** 2 * w)

    def C_ratio(s):
        w = mp.e ** (1j * s)
        return ((1 - D) ** 2 - D ** 2 * w) / (1 - 2 * D)

    rng = np.random.default_rng(7)
    x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
    x = [int(b) for b in x]
    rho = mp.mpf(p) / (1 - p)
    l = [mp.mpf(1)] * n
    r = [mp.mpf(1)] * n
    for i in range(1, n):
        l[i] = rho if x[i] == x[i - 1] else 1 / rho
    for b in range(0, n - 1):
        r[b] = rho if x[b + 1] == x[b] else 1 / rho

    def logS(eta):
        Fip1 = mp.mpf(1); Fip2 = mp.mpf(1); Gip1 = mp.mpf(0); logsc = mp.mpf(0)
        for i in range(n - 1, -1, -1):
            Gi = eta * (r[i] * Fip2 + Gip1)
            Fi = Fip1 + l[i] * Gi
            Fip2 = Fip1; Fip1 = Fi; Gip1 = Gi
            a = abs(Fip1)
            if a > 0:
                Fip1 /= a; Fip2 /= a; Gip1 /= a; logsc += mp.log(a)
        return mp.log(Fip1) + logsc

    def lP(s):
        return n * mp.log(C_ratio(s)) + logS(eta_closed(s))

    h = mp.mpf('1e-5')
    d3 = (lP(2 * h) - 2 * lP(h) + 2 * lP(-h) - lP(-2 * h)) / (2 * h ** 3)
    # also sample the third derivative AWAY from 0 (toward window edge) to estimate sup
    edge = mp.mpf('0.03')
    d3e = (lP(edge + 2 * h) - 2 * lP(edge + h) + 2 * lP(edge - h) - lP(edge - 2 * h)) / (2 * h ** 3)
    return float(abs(d3)) / n, float(abs(d3e)) / n


if __name__ == "__main__":
    print("=" * 80)
    print("ADV (1b): window-size self-consistency  eps = min(delta_KP, r_KP/8A, c' vbar/C3)")
    print("=" * 80)
    for p in (0.4, 0.25, 0.1, 0.05):
        for frac in (0.9, 0.99):
            D = frac * Dc(p)
            A = A_of_D(D)
            rkp = r_lb(p)  # conservative: r_KP >= r_lb
            vbar = vbar_closed(p, D)
            c3_central, c3_edge = measure_C3_over_n(p, D, n=160)
            C3 = max(c3_central, c3_edge) * 1.5  # conservative sup multiplier
            cprime = 1.0
            delta_KP = rkp / A
            eps_W2 = rkp / (8 * A)
            eps_W3 = cprime * vbar / C3
            y_max = rkp / (8 * A)
            eps = min(delta_KP, eps_W2, eps_W3)
            print(f"p={p} D={frac:.2f}Dc={D:.5f} | A={A:.4f} r_KP>={rkp:.5f} vbar={vbar:.4f} C3~{C3:.4f}")
            print(f"    delta_KP={delta_KP:.5f}  eps_W2(strip)={eps_W2:.5f}  "
                  f"eps_W3(zone2)={eps_W3:.5f}  -> eps={eps:.6f}  y_max={y_max:.6f}")
            # is eps Theta(1)?  Report which constraint binds.
            binder = ("strip(W2)" if eps == eps_W2 else
                      "zone2(W3)" if eps == eps_W3 else "KPwindow(W1)")
            print(f"    binding constraint: {binder};  eps is {'O(1)-ish' if eps>0.01 else 'SMALL (<0.01)'}")
            print()

    print("=" * 80)
    print("ADV (4): vertical-end budget  -1/4 vbar n eps^2 + C3 n|s0| + t|s0| < 0 ?")
    print("   ( |s0| = |beta|/sqrt(vbar n);  beta on shell ~ O(sqrt(log n)) )")
    print("=" * 80)
    p, D = 0.4, 0.9 * Dc(0.4)
    A = A_of_D(D); rkp = r_lb(p); vbar = vbar_closed(p, D)
    c3c, c3e = measure_C3_over_n(p, D, 160); C3 = max(c3c, c3e) * 1.5
    eps = min(rkp / A, rkp / (8 * A), vbar / C3)
    print(f"p={p} D=0.9Dc={D:.5f}: A={A:.4f} vbar={vbar:.4f} C3={C3:.4f} eps={eps:.6f}")
    for n in (256, 1024, 10**4, 10**6, 10**9, 10**12):
        beta = 2 * math.sqrt(math.log(n))   # shell extreme
        s0 = beta / math.sqrt(vbar * n)
        zone2 = -0.25 * vbar * n * eps ** 2
        travelC3 = C3 * n * s0
        travelt = math.floor(n * D) * s0
        total = zone2 + travelC3 + travelt
        print(f"  n={n:>13}: beta={beta:.2f} |s0|={s0:.3e} y_max={rkp/(8*A):.5f} "
              f"interior?{s0 <= rkp/(8*A)} | zone2={zone2:.2e} C3n|s0|={travelC3:.2e} "
              f"t|s0|={travelt:.2e} TOTAL={total:.2e} {'CLOSES' if total<0 else '*** OPEN ***'}")
