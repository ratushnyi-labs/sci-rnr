#!/usr/bin/env python3
"""
ADVERSARIAL: the TWO eps's.  The draft uses a single symbol eps for:
  (a) the STRIP half-width: eps <= r_KP/(8A) (so |eta_z|<=r_KP on the shifted line)
      -- numerically ~0.013 at p=0.4.
  (b) the P1 TAIL threshold eps_{P1} = min(delta_KP, c' vbar/C3) where |Phi|<=2e^{-cn}
      -- this is Theta(1), the FIXED tail/window split point.
At the strip edge (a), is |Phi| actually exponentially small (as Step 4(a) needs)?
Measure |Phi(s)| vs s to locate where -log|Phi| becomes Theta(n).
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
    eps_strip = r_lb / (8 * A)
    print("=" * 80)
    print(f"p={p} D=0.9Dc={D:.5f}  A={A:.4f} r_lb={r_lb:.5f}")
    print(f"strip half-width eps_strip = r_lb/(8A) = {eps_strip:.5f}")
    print(f"(the P1 tail threshold delta_KP = r_lb/A = {r_lb/A:.5f}, and the *measured*")
    print(f" KP window from tex is ~0.35-1.01 -- ORDERS larger than eps_strip)")
    print("=" * 80)
    print("\n -log|Phi(s)|/n  (per-site rate) across s for shell words, several n:")
    rng = np.random.default_rng(2)
    for n in (200, 400, 800):
        # average rate over a few shell words
        srates = {}
        nw = 6
        cnt = 0
        s_probe = [eps_strip, 0.03, 0.05, 0.1, 0.2, 0.4, 0.8, math.pi]
        acc = {s: 0.0 for s in s_probe}
        while cnt < nw:
            x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
            sw = int(np.sum(x[1:] != x[:-1]))
            if abs(sw - (n - 1) * p) > 2 * math.sqrt(n * p * (1 - p)):
                continue
            x = [int(b) for b in x]
            lP = make_logPhi(x, p, D)
            for s in s_probe:
                acc[s] += -float((lP(s)).real) / n
            cnt += 1
        print(f"  n={n}:")
        for s in s_probe:
            rate = acc[s] / nw
            tag = "  <- eps_strip" if abs(s - eps_strip) < 1e-9 else ""
            print(f"     s={s:.4f}: -log|Phi|/n = {rate:.5f}   (|Phi|=e^{{-{rate*n:.2f}}}){tag}")
    print("\nDIAGNOSIS: at s=eps_strip the per-site rate is ~O(s^2)=tiny, so |Phi(eps_strip)|~1,")
    print("NOT 2e^{-cn}.  The draft's Step-4(a) substitution |Phi(eps)|<=2e^{-cn} uses the")
    print("WRONG eps (it needs the Theta(1) tail threshold, but the STRIP forces eps~0.013).")
