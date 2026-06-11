#!/usr/bin/env python3
"""
ADVERSARIAL (4-deep): the vertical-end segments, DONE PROPERLY with sign.

The draft's Step 4(a) writes the vertical travel as
   integral_0^{|s0|} |Phi(+-eps+iy)| e^{ty} dy  <=  |s0| 2 e^{-cn} e^{C3 n|s0|} e^{t|s0|}
and then claims t|s0| = O(1).  We show:
  (A) t|s0| = Theta(sqrt n), NOT O(1)  -- the draft's "= O(|beta|)" is wrong by sqrt n.
  (B) BUT the actual vertical contribution has the FAVORABLE sign: the real
      integrand on the vertical segment is |Phi(+-eps+iy)| * e^{ty} with y going
      from 0 to s0 = (mu-t)/v.  When beta>0, s0<0, so e^{ty}<1 (decaying).
      Even worst-case over sign, the question is whether
        -log|Phi(eps)| (= Theta(n)) beats the GROWTH from the y-travel.
  (C) The DECISIVE quantity is whether the *whole* vertical integrand magnitude
      |Phi(eps+iy)| e^{ty} stays << the central peak ~ 1/sqrt(v_x) at q_t scale.
      We measure |Phi(eps+iy)| e^{ty} DIRECTLY for real words and compare to the
      central mass q_t.

We compute everything on real shell words via the O(n) recursion (exact).
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


def make_recursion(x, p):
    n = len(x)
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
    return logS


def logPhi_factory(x, p, D):
    logS = make_recursion(x, p)
    n = len(x)

    def lP(s):
        return n * mp.log(C_ratio(p, D, s)) + logS(eta_closed(p, D, s))
    return lP


def cumulants(lP, h=mp.mpf('1e-6')):
    f0 = lP(0); fp = lP(h); fm = lP(-h)
    d1 = (fp - fm) / (2 * h)
    d2 = (fp - 2 * f0 + fm) / (h * h)
    return float((d1 / 1j).real), float((-d2).real)


if __name__ == "__main__":
    p = 0.4; D = 0.9 * Dc(p)
    rng = np.random.default_rng(5)
    print("=" * 80)
    print("(A)+(B)+(C): vertical-end real magnitude vs the draft's worst-case bound")
    print(f"p={p} D=0.9Dc={D:.5f}")
    print("=" * 80)
    for n in (200, 400, 800):
        eps = 0.013  # the fixed strip window from the budget probe
        t = math.floor(n * D)
        # build a shell word
        while True:
            x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
            sw = int(np.sum(x[1:] != x[:-1]))
            if abs(sw - (n - 1) * p) <= 2 * math.sqrt(n * p * (1 - p)):
                break
        x = [int(b) for b in x]
        lP = logPhi_factory(x, p, D)
        mu, v = cumulants(lP)
        beta = (t - mu) / math.sqrt(v)
        s0 = (mu - t) / v
        # the CENTRAL mass scale ~ 1/sqrt(2 pi v) (the peak of q_t integrand contribution)
        central_scale = 1.0 / math.sqrt(2 * math.pi * v)
        # |Phi(eps)| (the tail-edge value)
        absPhi_eps = float(abs(mp.e ** lP(eps)))
        # the FULL vertical integrand sup over y in [0, s0]:  |Phi(eps+iy)| e^{t y}
        # scan y between 0 and s0
        ys = np.linspace(0, s0, 25)
        vert_vals = []
        for y in ys:
            val = float(abs(mp.e ** lP(eps + 1j * y))) * math.exp(t * float(y))
            vert_vals.append(val)
        vert_sup = max(vert_vals)
        # draft's worst-case bound factor:  2 e^{-cn} e^{C3 n|s0|} e^{t|s0|}
        # we don't have c,C3 cleanly; instead compare the REAL vert_sup to central_scale
        # and report t|s0| to expose the sqrt-n growth.
        t_s0 = t * abs(s0)
        print(f" n={n:>4} t={t:>4} beta={beta:+.3f} s0={s0:+.5f} | "
              f"t|s0|={t_s0:.3f} (draft claims O(1); sqrt(n)={math.sqrt(n):.1f})")
        print(f"        |Phi(eps)|={absPhi_eps:.3e}  central_scale 1/sqrt(2piv)={central_scale:.3e}")
        print(f"        REAL vert integrand sup_y |Phi(eps+iy)|e^(ty) = {vert_sup:.3e}  "
              f"(ratio to central = {vert_sup/central_scale:.3e})")
        # vertical CONTRIBUTION to I (length |s0| times sup), vs central mass q_t~central_scale
        vert_contrib = abs(s0) * vert_sup
        print(f"        vertical contribution ~ |s0|*sup = {vert_contrib:.3e}  "
              f"vs |I|~central_scale={central_scale:.3e}  "
              f"=> {'NEGLIGIBLE' if vert_contrib < 0.01*central_scale else '*** NOT negligible ***'}")
        print()
    print("CONCLUSION on t|s0|:  t|s0| GROWS ~ sqrt(n) (NOT O(1) as the draft writes).")
    print("Decisive: whether |Phi(eps+iy)| e^{ty} stays exponentially below central_scale.")
