#!/usr/bin/env python3
"""
ADVERSARIAL (decisive, large-n): does the eps_strip-window split become VALID
as n grows past the window-vs-peak crossover (~6.7e4 at p=0.4)?
Fast numpy logS recursion (vectorized over s-grid).  Compute:
  central |s|<=eps_strip  vs  tail eps_strip<|s|<=pi  vs  q_t exact.
If central -> q_t and tail -> 0 (relative), the architecture is asymptotically
sound (large-n0).  We push n to 1e5, 3e5, 1e6 at p=0.4.
"""
import math
import numpy as np


def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))


def eta_C(p, D, s):
    th0 = math.log(D / (1 - D))
    eb = np.exp(th0 + 1j * s)
    ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
    eta = (1 - ze) / (1 + ze)
    C = (1 + eb) * (1 + ze) / 2.0
    return eta, C / (1.0 / (1.0 - D))


def logS_vec(x, p, eta):
    n = len(x); rho = p / (1.0 - p)
    # build l,r as arrays
    same_prev = np.empty(n, dtype=bool); same_prev[0] = True
    same_prev[1:] = (x[1:] == x[:-1])
    l = np.where(same_prev, rho, 1.0 / rho); l[0] = 1.0
    same_next = np.empty(n, dtype=bool); same_next[-1] = True
    same_next[:-1] = (x[1:] == x[:-1])
    r = np.where(same_next, rho, 1.0 / rho); r[-1] = 1.0
    m = len(eta)
    Fip1 = np.ones(m, dtype=complex); Fip2 = np.ones(m, dtype=complex)
    Gip1 = np.zeros(m, dtype=complex); logsc = np.zeros(m, dtype=complex)
    for i in range(n - 1, -1, -1):
        Gi = eta * (r[i] * Fip2 + Gip1)
        Fi = Fip1 + l[i] * Gi
        Fip2 = Fip1; Fip1 = Fi; Gip1 = Gi
        if (i % 8) == 0:
            a = np.abs(Fip1); a = np.where(a > 0, a, 1.0)
            Fip1 = Fip1 / a; Fip2 = Fip2 / a; Gip1 = Gip1 / a; logsc += np.log(a)
    return np.log(Fip1) + logsc


def trap_part(Phi, sg, t, mask):
    integ = Phi * np.exp(-1j * sg * t)
    return np.real(np.trapezoid(integ[mask], sg[mask])) / (2 * math.pi)


if __name__ == "__main__":
    p = 0.4; D = 0.9 * Dc(p)
    A = D * (1 - D) / (1 - 2 * D)
    r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
    eps_strip = r_lb / (8 * A)
    kappa2 = (1 - 2 * p) ** 2 / (p ** 2 * (1 - p) ** 2)
    vbar = D * (1 - D) * (1 - kappa2 * D * (1 - D) / (1 - 2 * D) ** 2)
    print("=" * 92)
    print(f"LARGE-N split validity.  p={p} D={D:.5f} eps_strip={eps_strip:.5f} vbar={vbar:.4f}")
    print(f"window-vs-peak crossover n ~ 1/(vbar eps_strip^2) = {1/(vbar*eps_strip**2):.2e}")
    print("=" * 92)
    rng = np.random.default_rng(77)
    for n in (10**4, 10**5, 3 * 10**5, 10**6):
        t = math.floor(n * D)
        # one shell word
        while True:
            x = (np.cumsum(rng.random(n) < p) % 2).astype(np.int8)
            sw = int(np.sum(x[1:] != x[:-1]))
            if abs(sw - (n - 1) * p) <= 2 * math.sqrt(n * p * (1 - p)):
                break
        # peak width ~ 1/sqrt(vbar n); grid: dense in [-eps_strip, eps_strip], then to pi
        pw = 1.0 / math.sqrt(vbar * n)
        # fine grid covering many peak widths AND the strip edge
        s_in = np.linspace(-eps_strip, eps_strip, 4001)
        s_out = np.concatenate([np.linspace(-math.pi, -eps_strip, 3000, endpoint=False),
                                np.linspace(eps_strip, math.pi, 3000)[1:]])
        sg = np.sort(np.concatenate([s_in, s_out]))
        eta_s, Cr = eta_C(p, D, sg)
        logPhi = n * np.log(Cr) + logS_vec(x, p, eta_s)
        Phi = np.exp(logPhi)
        qt = trap_part(Phi, sg, t, np.ones(len(sg), bool))
        mask_in = np.abs(sg) <= eps_strip + 1e-15
        central = trap_part(Phi, sg, t, mask_in)
        tail = trap_part(Phi, sg, t, ~mask_in)
        # cumulants by small-s finite diff
        sh = np.array([0.0008, 0.0016, -0.0008, -0.0016])
        eh, Ch = eta_C(p, D, sh)
        lp = n * np.log(Ch) + logS_vec(x, p, eh)
        mu = (8 * lp[0].imag - lp[1].imag) / (6 * 0.0008)
        vh = -2 * lp[0].real / 0.0008 ** 2; v2h = -2 * lp[1].real / 0.0016 ** 2
        v = (4 * vh - v2h) / 3
        beta = (t - mu) / math.sqrt(v)
        peakw = 1.0 / math.sqrt(v)
        print(f"\n n={n:>8} t={t} beta={beta:+.3f}  peak width 1/sqrt(v)={peakw:.5f}  "
              f"eps_strip/peakw = {eps_strip/peakw:.2f} (#widths in window)")
        print(f"   q_t exact      = {qt:.6e}")
        print(f"   central(<=eps) = {central:.6e}   (central/q_t = {central/qt:.4f})")
        print(f"   tail(>eps)     = {tail:.6e}   (tail/q_t   = {tail/qt:.4f})")
        verdict = ("SPLIT VALID (central~q_t, tail~0)" if central / qt > 0.95
                   else "split still leaking" if central / qt > 0.5
                   else "*** central misses bulk -- window too narrow ***")
        print(f"   => {verdict}")
