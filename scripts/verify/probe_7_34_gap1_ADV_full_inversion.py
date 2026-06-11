#!/usr/bin/env python3
"""
ADVERSARIAL (decisive): compute q_t THREE ways and compare, to see whether the
contour-shift on the SMALL window eps_strip actually reproduces q_t or whether
the verticals/middle-zone spoil it.

 (M0) q_t exact = (1/2pi) Re int_{-pi}^{pi} Phi(s) e^{-ist} ds   (fine grid)
 (M1) saddle prediction: (2 pi v_x)^{-1/2} e^{-beta^2/2}
 (M2) contour-shifted central piece on |s|<=eps_strip ALONG s+i s0, PLUS the two
      vertical ends, PLUS the real-axis tail |s|>eps_strip.   Should equal (M0).
 We report:
   - central(shifted) contribution
   - vertical-ends contribution  (the thing the draft claims is O(e^{-cn}))
   - real tail |s|in[eps_strip,pi] contribution
   - their sum vs q_t exact, and vs the saddle prediction M1.
The KEY diagnostic: is |vertical ends| << |central|?  And is |real tail| << central?
If the verticals or the tail are NOT negligible at the chosen eps_strip, the
draft's split is INVALID at that window size.
"""
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 45


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


def integ_line(lP, t, a, b, npts, y):
    """integral of Phi(s+iy) e^{-i t (s+iy)} ds over s in [a,b], complex."""
    ss = np.linspace(a, b, npts)
    vals = []
    for s in ss:
        z = mp.mpf(s) + 1j * mp.mpf(y)
        vals.append(complex(mp.e ** (lP(z) - 1j * t * z)))
    vals = np.array(vals)
    return np.trapezoid(vals, ss)


def integ_vertical(lP, t, sedge, y0, npts):
    """integral over the vertical segment s=sedge, y from 0 to y0, ds = i dy."""
    ys = np.linspace(0, y0, npts)
    vals = []
    for y in ys:
        z = mp.mpf(sedge) + 1j * mp.mpf(y)
        vals.append(complex(mp.e ** (lP(z) - 1j * t * z)))
    vals = np.array(vals)
    # ds = i dy
    return 1j * np.trapezoid(vals, ys)


if __name__ == "__main__":
    p = 0.4; D = 0.9 * Dc(p)
    A = D * (1 - D) / (1 - 2 * D)
    r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
    eps_strip = r_lb / (8 * A)
    print("=" * 92)
    print(f"FULL inversion 3 ways.  p={p} D=0.9Dc={D:.5f}  eps_strip={eps_strip:.5f}")
    print("=" * 92)
    rng = np.random.default_rng(21)
    for n in (200, 400, 800):
        t = math.floor(n * D)
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
        # M0: exact, fine grid (dense near 0)
        w = 1.0 / math.sqrt(n * 0.05)
        sf = np.linspace(-15 * w, 15 * w, 2001)
        st = np.concatenate([np.linspace(-math.pi, -15 * w, 600, endpoint=False),
                             np.linspace(15 * w, math.pi, 600)])
        sg = np.sort(np.concatenate([sf, st]))
        Phi = np.array([complex(mp.e ** lP(s)) for s in sg])
        qt_exact = float(np.real(np.trapezoid(Phi * np.exp(-1j * sg * t), sg)) / (2 * math.pi))
        # M1 saddle
        qt_saddle = 1.0 / math.sqrt(2 * math.pi * v) * math.exp(-beta * beta / 2)
        # M2: central shifted + verticals + real tail
        central = integ_line(lP, t, -eps_strip, eps_strip, 801, s0) / (2 * math.pi)
        vL = integ_vertical(lP, t, -eps_strip, s0, 201) / (2 * math.pi)  # left edge, 0->s0
        vR = integ_vertical(lP, t, eps_strip, s0, 201) / (2 * math.pi)
        # the rectangle: real axis [-eps,eps] = shifted line + verticals (sign bookkeeping)
        # closed contour: real[-eps,eps] - vert_left(0->s0) + shifted(-eps->eps at y=s0) + vert_right(s0->0)
        # => real_central = shifted_central + (vR_contribution with proper orientation)
        # we instead DIRECTLY compute real-axis central and compare to shifted+verticals.
        real_central = integ_line(lP, t, -eps_strip, eps_strip, 801, 0.0) / (2 * math.pi)
        # real tail
        tail = (integ_line(lP, t, eps_strip, math.pi, 1500, 0.0)
                + integ_line(lP, t, -math.pi, -eps_strip, 1500, 0.0)) / (2 * math.pi)
        print(f"\n n={n} t={t} beta={beta:+.3f} s0={s0:+.5f} v={v:.3f}")
        print(f"   q_t exact (M0)            = {qt_exact:.6e}")
        print(f"   saddle pred (M1)          = {qt_saddle:.6e}   (ratio M1/M0 = {qt_saddle/qt_exact:.4f})")
        print(f"   real central |s|<=eps     = {real_central.real:+.6e} (+{real_central.imag:+.2e} i)")
        print(f"   real tail eps<|s|<=pi     = {tail.real:+.6e} (+{tail.imag:+.2e} i)  "
              f"[|tail|/|central| = {abs(tail)/abs(real_central):.3e}]")
        print(f"   shifted central (y=s0)    = {central.real:+.6e} (+{central.imag:+.2e} i)")
        print(f"   vertical end R (0->s0)    = {vR.real:+.6e} (+{vR.imag:+.2e} i)")
        print(f"   vertical end L (0->s0)    = {vL.real:+.6e} (+{vL.imag:+.2e} i)")
        # Cauchy check: real_central = shifted_central + (vR - vL) with orientation
        cauchy_resid = real_central - (central + (vR - vL))
        print(f"   Cauchy residual real_central-(shifted+vR-vL) = {abs(cauchy_resid):.2e}")
        # decisive: are verticals negligible vs central?
        vmag = max(abs(vR), abs(vL))
        print(f"   *** |vertical|/|central| = {vmag/abs(real_central):.3e}  "
              f"=> {'NEGLIGIBLE' if vmag < 0.05*abs(real_central) else '*** NOT NEGLIGIBLE -- split invalid at eps_strip ***'}")
