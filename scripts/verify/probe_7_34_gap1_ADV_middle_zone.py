#!/usr/bin/env python3
"""
ADVERSARIAL: the intermediate zone  eps_strip <= |s| <= eps_P1.
The strip-shift captures the central Gaussian on |s|<=eps_strip~0.013.
P1's tail bound |Phi|<=2e^{-cn} only holds for |s|>=eps_P1~0.1-0.4.
QUESTION: is the band [eps_strip, eps_P1] covered?  The committed tex claims
"no middle zone: the KP window is fixed while g(s)~1-as^2 lets the replica bound
reach s~n^{-1/2+kappa}, so the handoff overlaps with Theta(1) slack."

Test: does -log|Phi(s)| >= (vbar/4) n s^2 (zone-2 quadratic domination) hold on
the WHOLE band [eps_strip, eps_P1]?  If so, the Gaussian-tail bound
int_{eps_strip}^{eps_P1} e^{-vbar n s^2/4} ds is exponentially small in n at the
LOWER edge eps_strip:  e^{-vbar n eps_strip^2/4}.  Is n eps_strip^2 -> infinity?
  eps_strip ~ 0.013 FIXED => vbar n eps_strip^2/4 = Theta(n) -> infinity.  GOOD if
quadratic domination holds down to eps_strip.  CHECK the constant C3 n s <= v/4
condition: zone-2 requires s <= vbar/(4 C3).  Is eps_strip <= vbar/(4 C3)?
"""
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 40


def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))


def eta_closed(p, D, s):
    w = mp.e ** (1j * s)
    return D * (1 - D) * (w - 1) / ((1 - D) ** 2 - D ** 2 * w)


def C_ratio(p, D, s):
    w = mp.e ** (1j * s)
    return ((1 - D) ** 2 - D ** 2 * w) / (1 - 2 * D)


def vbar_closed(p, D):
    kappa2 = (1 - 2 * p) ** 2 / (p ** 2 * (1 - p) ** 2)
    Psi0 = kappa2 * D * (1 - D) / (1 - 2 * D) ** 2
    return D * (1 - D) * (1 - Psi0)


def make_logPhi(x, p, D):
    n = len(x); rho = mp.mpf(p) / (1 - p)
    l = [mp.mpf(1)] * n; r = [mp.mpf(1)] * n
    for i in range(1, n):
        l[i] = rho if x[i] == x[i - 1] else 1 / rho
    for b in range(0, n - 1):
        r[b] = rho if x[b + 1] == x[b] else 1 / rho

    def lP(s):
        eta = eta_closed(p, D, s)
        Fip1 = mp.mpf(1); Fip2 = mp.mpf(1); Gip1 = mp.mpf(0); logsc = mp.mpf(0)
        for i in range(n - 1, -1, -1):
            Gi = eta * (r[i] * Fip2 + Gip1)
            Fi = Fip1 + l[i] * Gi
            Fip2 = Fip1; Fip1 = Fi; Gip1 = Gi
            a = abs(Fip1)
            if a > 0:
                Fip1 /= a; Fip2 /= a; Gip1 /= a; logsc += mp.log(a)
        return n * mp.log(C_ratio(p, D, s)) + mp.log(Fip1) + logsc
    return lP


if __name__ == "__main__":
    p = 0.4; D = 0.9 * Dc(p)
    A = D * (1 - D) / (1 - 2 * D)
    r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
    vbar = vbar_closed(p, D)
    eps_strip = r_lb / (8 * A)
    # zone-2 condition s <= vbar/(4 C3); estimate C3~0.1
    C3 = 0.10
    s_zone2 = vbar / (4 * C3)
    print("=" * 80)
    print(f"p={p} D=0.9Dc={D:.5f}  vbar={vbar:.4f}  eps_strip={eps_strip:.5f}")
    print(f"zone-2 valid for s <= vbar/(4 C3) = {s_zone2:.5f}  (C3~{C3})")
    print(f"so quadratic domination -log|Phi|>=vbar n s^2/4 should hold on "
          f"[{eps_strip:.4f}, {s_zone2:.4f}]")
    print("=" * 80)
    print("\nCheck -log|Phi(s)|/n  vs the floor vbar s^2/4  across the band:")
    rng = np.random.default_rng(9)
    n = 600
    nw = 8
    s_band = np.linspace(eps_strip, 0.5, 12)
    acc = np.zeros(len(s_band))
    cnt = 0
    while cnt < nw:
        x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
        sw = int(np.sum(x[1:] != x[:-1]))
        if abs(sw - (n - 1) * p) > 2 * math.sqrt(n * p * (1 - p)):
            continue
        x = [int(b) for b in x]
        lP = make_logPhi(x, p, D)
        for k, s in enumerate(s_band):
            acc[k] += -float((lP(s)).real) / n
        cnt += 1
    acc /= nw
    print(f"n={n}, avg over {nw} shell words:")
    print(f"{'s':>8} {'-log|Phi|/n':>14} {'floor vbar s^2/4':>18} {'ratio':>8} {'dominated?':>10}")
    all_ok = True
    for k, s in enumerate(s_band):
        floor = vbar * s ** 2 / 4
        ratio = acc[k] / floor if floor > 0 else float('inf')
        ok = acc[k] >= floor * 0.999
        all_ok = all_ok and ok
        print(f"{s:>8.4f} {acc[k]:>14.6f} {floor:>18.6f} {ratio:>8.2f} {str(ok):>10}")
    print(f"\nQuadratic domination on the band [eps_strip, 0.5]: {'HOLDS' if all_ok else 'FAILS somewhere'}")
    print(f"At the lower edge eps_strip={eps_strip:.4f}:  vbar n eps_strip^2/4 grows like "
          f"{vbar*eps_strip**2/4:.2e} * n  -> infinity  (so the Gaussian-tail beyond")
    print(f"eps_strip IS exponentially small in n, IF domination holds down to eps_strip).")
