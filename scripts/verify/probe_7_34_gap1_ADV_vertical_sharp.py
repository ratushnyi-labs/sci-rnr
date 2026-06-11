#!/usr/bin/env python3
"""
ADVERSARIAL (4-SHARP): the CORRECT vertical-end bound.
The draft bounds the vertical travel by |psi'|<=C3 n  ==> e^{C3 n|s0|}=e^{Theta(sqrt n)},
then ADDS e^{t|s0|}=e^{Theta(sqrt n)}.  This double-counts and is far too lossy.
The TRUTH: on the vertical segment s=eps+iy, the integrand modulus is
   |Phi(eps+iy)| e^{ty} = e^{Re psi(eps+iy)} e^{ty}.
Re psi(eps+iy) = Re psi(eps) - y Im psi'(eps) + O(y^2 |psi''|).
Im psi'(eps) ~ mu_x (the cumulant!), so
   Re psi(eps+iy)+ty ~ Re psi(eps) + y(t - mu_x) + O(y^2 v_x)
                     = Re psi(eps) - y v_x s0 + O(y^2 v_x).
Over y in [0,s0]:  max of -y v_x s0 is at y=s0 (if s0<0) or 0; |y v_x s0|<=v_x s0^2=beta^2=O(1).
So the travel factor is e^{O(beta^2)} = O(1), and the vertical end is
   ~ |s0| * |Phi(eps)| * e^{O(beta^2)} = O(1/sqrt n) * e^{-Theta(n)} * O(1) = e^{-Theta(n)}.
The Theta(sqrt n) in the draft is a LOSSY-bound artifact, not real.

We VERIFY by measuring the EXACT vertical integrand and its max, and showing
log10(vert_sup) ~ -Theta(n)  (linear in n), with NO sqrt(n) inflation.
"""
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 50


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


def cumulants(lP, h=mp.mpf('1e-6')):
    f0 = lP(0); fp = lP(h); fm = lP(-h)
    mu = float(((fp - fm) / (2 * h) / 1j).real)
    v = float((-(fp - 2 * f0 + fm) / (h * h)).real)
    return mu, v


if __name__ == "__main__":
    p = 0.4; D = 0.9 * Dc(p)
    A = D * (1 - D) / (1 - 2 * D)
    r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
    eps_strip = r_lb / (8 * A)
    print("=" * 86)
    print(f"SHARP vertical end:  log10[ sup_y |Phi(eps+iy)| e^(ty) ]  vs n   (should be ~ -Theta(n))")
    print(f"p={p} D=0.9Dc={D:.5f}  eps_strip={eps_strip:.5f}")
    print("=" * 86)
    rng = np.random.default_rng(13)
    rows = []
    for n in (200, 400, 800, 1600):
        t = math.floor(n * D)
        # shell word
        while True:
            x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
            sw = int(np.sum(x[1:] != x[:-1]))
            if abs(sw - (n - 1) * p) <= 2 * math.sqrt(n * p * (1 - p)):
                break
        x = [int(b) for b in x]
        lP = make_logPhi(x, p, D)
        mu, v = cumulants(lP)
        beta = (t - mu) / math.sqrt(v)
        s0 = (mu - t) / v
        # EXACT vertical integrand max over y in [0,s0]
        ys = np.linspace(0, s0, 41)
        logvert = []
        for y in ys:
            val = float((lP(eps_strip + 1j * y)).real) + t * float(y)  # log of integrand modulus
            logvert.append(val)
        logsup = max(logvert)          # natural log of the sup
        log10sup = logsup / math.log(10)
        # central scale for comparison: log10(1/sqrt(2 pi v))
        log10central = math.log(1 / math.sqrt(2 * math.pi * v)) / math.log(10)
        rows.append((n, beta, s0, logsup, log10sup, log10central, v * s0 * s0))
        print(f" n={n:>5} beta={beta:+.3f} s0={s0:+.5f} | v s0^2(=beta^2)={v*s0*s0:.4f}  "
              f"log10[vert sup]={log10sup:+.3f}  log10[central~1/sqrt(2piv)]={log10central:+.3f}  "
              f"=> vert/central = 10^{log10sup-log10central:+.2f}")
    print("-" * 86)
    # is log(vert sup) linear in -n?  fit slope
    ns = np.array([r[0] for r in rows]); ls = np.array([r[3] for r in rows])
    slope = np.polyfit(ns, ls, 1)[0]
    print(f"slope d[log(vert sup)]/dn = {slope:.5f}  (NEGATIVE & ~const => -Theta(n), NO sqrt(n) inflation)")
    print(f"v s0^2 = beta^2 stays O(1) across n: {[f'{r[6]:.3f}' for r in rows]}  -- the travel factor is O(1), not e^(Theta(sqrt n)).")
    print("\nVERDICT: the SHARP vertical bound is e^{-Theta(n)} at ALL n (the mu_x in Im psi'(eps)")
    print("cancels t exactly into v_x s0^2 = beta^2 = O(1)).  The draft's e^{C3 n|s0|}e^{t|s0|}")
    print("= e^{Theta(sqrt n)} is a LOSSY artifact; the conclusion (negligible) holds, the")
    print("written bound is just far from sharp.")
