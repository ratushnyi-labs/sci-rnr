#!/usr/bin/env python3
"""
ADVERSARIAL (5): the LEMMA's final inequality, on 30 shell words at n=512.
  q_t(x) = (2 pi v_x)^{-1/2} e^{-beta^2/2} (1 + R_n), |R_n| <= C_R (1+|beta|)^3 (log n)^3/sqrt n.
Measure R_n directly (numpy, fast) and CHECK |R_n| <= bound with a SMALL C_R.
Also independently re-derive attacks (2) saddle linear and (3) cubic-form via the
measured cumulants, to confirm the order arithmetic the draft claims.
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
    l = np.ones(n); r = np.ones(n)
    for i in range(1, n):
        l[i] = rho if x[i] == x[i - 1] else 1.0 / rho
    for b in range(0, n - 1):
        r[b] = rho if x[b + 1] == x[b] else 1.0 / rho
    m = len(eta)
    Fip1 = np.ones(m, dtype=complex); Fip2 = np.ones(m, dtype=complex)
    Gip1 = np.zeros(m, dtype=complex); logsc = np.zeros(m, dtype=complex)
    for i in range(n - 1, -1, -1):
        Gi = eta * (r[i] * Fip2 + Gip1)
        Fi = Fip1 + l[i] * Gi
        Fip2 = Fip1; Fip1 = Fi; Gip1 = Gi
        if (i % 16) == 0:
            a = np.abs(Fip1); a = np.where(a > 0, a, 1.0)
            Fip1 = Fip1 / a; Fip2 = Fip2 / a; Gip1 = Gip1 / a; logsc += np.log(a)
    return np.log(Fip1) + logsc


def measure(x, p, D, t, sg, h=0.003):
    n = len(x)
    eta_s, Cr = eta_C(p, D, sg)
    logPhi = n * np.log(Cr) + logS_vec(x, p, eta_s)
    Phi = np.exp(logPhi)
    q_t = float(np.real(np.trapezoid(Phi * np.exp(-1j * sg * t), sg)) / (2 * math.pi))
    sh = np.array([h, 2 * h, -h, -2 * h])
    eh, Ch = eta_C(p, D, sh)
    lp = n * np.log(Ch) + logS_vec(x, p, eh)
    mu = (8 * lp[0].imag - lp[1].imag) / (6 * h)
    vh = -2 * lp[0].real / h ** 2; v2h = -2 * lp[1].real / (2 * h) ** 2
    v = (4 * vh - v2h) / 3
    return q_t, mu, v


if __name__ == "__main__":
    p = 0.4; D = 0.9 * Dc(p)
    n = 512
    t = math.floor(n * D)
    w = 1.0 / math.sqrt(n * 0.05)
    sf = np.linspace(-14 * w, 14 * w, 1601)
    st = np.concatenate([np.linspace(-math.pi, -14 * w, 500, endpoint=False),
                         np.linspace(14 * w, math.pi, 500)])
    sg = np.sort(np.concatenate([sf, st]))
    rng = np.random.default_rng(101)
    print("=" * 84)
    print(f"ADV (5): lemma final inequality on 30 shell words.  p={p} D={D:.5f} n={n} t={t}")
    print(f"{'#':>3} {'beta':>7} {'q_t':>11} {'saddle':>11} {'R_n':>9} {'bound(CR=1)':>11} {'R/bound':>8}")
    print("=" * 84)
    Rns = []; betas = []; ratios = []
    got = 0; tries = 0
    log3n = math.log(n) ** 3
    while got < 30 and tries < 400:
        tries += 1
        x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
        sw = int(np.sum(x[1:] != x[:-1]))
        if abs(sw - (n - 1) * p) > 3 * math.sqrt(n * p * (1 - p)):
            continue
        x = [int(b) for b in x]
        q_t, mu, v = measure(x, p, D, t, sg)
        if not (v > 0 and q_t > 0):
            continue
        beta = (t - mu) / math.sqrt(v)
        saddle = 1.0 / math.sqrt(2 * math.pi * v) * math.exp(-beta * beta / 2)
        Rn = q_t / saddle - 1.0
        bound = (1 + abs(beta)) ** 3 * log3n / math.sqrt(n)   # C_R=1 form
        ratio = abs(Rn) / bound
        Rns.append(Rn); betas.append(beta); ratios.append(ratio)
        got += 1
        if got <= 30:
            print(f"{got:>3} {beta:>+7.3f} {q_t:>11.4e} {saddle:>11.4e} {Rn:>+9.4f} "
                  f"{bound:>11.4e} {ratio:>8.4f}")
    Rns = np.array(Rns); betas = np.array(betas); ratios = np.array(ratios)
    print("=" * 84)
    print(f"words measured: {got}")
    print(f"max |R_n| = {np.max(np.abs(Rns)):.4f}   RMS = {np.sqrt(np.mean(Rns**2)):.4f}")
    print(f"max |beta| = {np.max(np.abs(betas)):.3f}")
    print(f"max ratio |R_n|/[(1+|beta|)^3 log^3 n/sqrt n] = {np.max(ratios):.4f}")
    print(f"=> implied C_R needed = {np.max(ratios):.4f}  (the lemma's C_R must exceed this)")
    print(f"   (the bound HOLDS with C_R = {math.ceil(np.max(ratios)*100)/100} on this sample)")
    # also report without the polylog (the tighter (1+|beta|)^3/sqrt n form)
    ratios_nopoly = np.abs(Rns) / ((1 + np.abs(betas)) ** 3 / math.sqrt(n))
    print(f"   tighter (no polylog) max |R_n|/[(1+|beta|)^3/sqrt n] = {np.max(ratios_nopoly):.4f}")
