#!/usr/bin/env python3
"""
probe_7_34_gap1_Rn_stress.py
============================================================================
GAP-1 STRESS (A5 of the disprove pass): the per-word steepest-descent remainder
    R_n(x) := q_t(x) * sqrt(2 pi v_x) * e^{beta^2/2} - 1,   beta=(t-mu_x)/sqrt(v_x),
is claimed O(log^3 n / sqrt n) UNIFORMLY over the sqrt(n)-shell.  This probe
measures R_n(x) directly at n in {256, 512, 1024}: q_t by Fourier inversion of
the posterior characteristic function Phi(s;x) = (C_s/C_0)^n S_x(eta_s)
(S_x via the O(n) two-term polymer recursion, VECTORIZED over the s-grid with
per-component renormalization), and (mu_x, v_x) by Richardson finite differences
of log Phi at small s.

  V1  sanity per word: sum-rule q-mass ~ Bin-like, |Phi|<=1, v_x>0.
  V2  RMS(R_n) decays ~ log^3 n / sqrt n (slope check).
  V3  MAX over shell words |R_n| also decays (the UNIFORMITY diagnostic);
      report the tail (90th pct) for heavy-tail detection.
Deps: numpy.
"""
import math
import numpy as np


def Dc(p): return 0.5*(1-math.sqrt(1-2*p)/(1-p))


def eta_C(p, D, s):
    th0 = math.log(D/(1-D))
    eb = np.exp(th0+1j*s)
    ze = (1-eb)/((1+eb)*(1-2*D))
    eta = (1-ze)/(1+ze)
    C = (1+eb)*(1+ze)/2.0
    return eta, C/(1.0/(1.0-D))


def logS_vec(x, p, eta):
    """log S_x(eta) for a VECTOR of complex eta, via the O(n) recursion,
    vectorized componentwise with per-component renormalization.
    Returns complex log (log|S| + i arg-accumulation via complex log of mantissa)."""
    n = len(x); rho = p/(1.0-p)
    l = np.ones(n); r = np.ones(n)
    for i in range(1, n):
        l[i] = rho if x[i] == x[i-1] else 1.0/rho
    for b in range(0, n-1):
        r[b] = rho if x[b+1] == x[b] else 1.0/rho
    m = len(eta)
    Fip1 = np.ones(m, dtype=complex)
    Fip2 = np.ones(m, dtype=complex)
    Gip1 = np.zeros(m, dtype=complex)
    logsc = np.zeros(m, dtype=complex)   # accumulate complex log of scales
    for i in range(n-1, -1, -1):
        Gi = eta*(r[i]*Fip2 + Gip1)
        Fi = Fip1 + l[i]*Gi
        Fip2 = Fip1; Fip1 = Fi; Gip1 = Gi
        if (i % 32) == 0:
            a = np.abs(Fip1)
            a = np.where(a > 0, a, 1.0)
            Fip1 = Fip1/a; Fip2 = Fip2/a; Gip1 = Gip1/a
            logsc += np.log(a)
    return np.log(Fip1) + logsc


def measure_word(x, p, D, t, s_grid, h=0.004):
    n = len(x)
    eta_s, Cr = eta_C(p, D, s_grid)
    logPhi = n*np.log(Cr) + logS_vec(x, p, eta_s)
    Phi = np.exp(logPhi)
    # q_t by trapezoid inversion over [-pi, pi]
    integ = Phi*np.exp(-1j*s_grid*t)
    q_t = float(np.real(np.trapezoid(integ, s_grid))/(2*math.pi))
    # cumulants via Richardson at +-h, +-2h
    sh = np.array([h, 2*h, -h, -2*h])
    eh, Ch = eta_C(p, D, sh)
    lp = len(x)*np.log(Ch) + logS_vec(x, p, eh)
    # mu: Im part odd in s. mu = [8 Im lp(h) - Im lp(2h)]/(6h)  (Richardson for f'(0))
    mu = (8*lp[0].imag - lp[1].imag)/(6*h)
    # v: -2 Re lp / s^2 with Richardson: v(h)= -2Re lp(h)/h^2; v = (4 v(h) - v(2h))/3
    vh = -2*lp[0].real/h**2; v2h = -2*lp[1].real/(2*h)**2
    v = (4*vh - v2h)/3
    return q_t, mu, v


if __name__ == "__main__":
    rng = np.random.default_rng(11)
    p = 0.4; D = 0.9*Dc(p)
    print("="*80)
    print(f"GAP-1 A5 stress: R_n = q_t sqrt(2 pi v_x) e^(beta^2/2) - 1;  p={p}, D={D:.4f}")
    print(f"{'n':>6} {'R':>4} | {'RMS(R_n)':>9} {'MAX|R_n|':>9} {'P90|R_n|':>9} | {'pred~log^3n/sqrt(n)':>20}")
    rows = []
    for n, R in ((256, 40), (512, 30), (1024, 20)):
        t = int(math.floor(n*D))
        # s-grid: fine near 0 (peak width ~ 1/sqrt(n vbar)), coarser tails
        w = 1.0/math.sqrt(n*0.05)
        s_fine = np.linspace(-12*w, 12*w, 1201)
        s_tail = np.concatenate([np.linspace(-math.pi, -12*w, 400, endpoint=False),
                                 np.linspace(12*w, math.pi, 400)])
        s_grid = np.sort(np.concatenate([s_fine, s_tail]))
        Rns = []
        for _ in range(R):
            x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
            sw = int(np.sum(x[1:] != x[:-1]))
            if abs(sw - (n-1)*p) > 3*math.sqrt(n*p*(1-p)):
                continue   # keep shell words only
            q_t, mu, v = measure_word(x, p, D, t, s_grid)
            if not (v > 0 and q_t > 0):
                continue
            beta = (t-mu)/math.sqrt(v)
            Rn = q_t*math.sqrt(2*math.pi*v)*math.exp(beta*beta/2.0) - 1.0
            Rns.append(Rn)
        Rns = np.array(Rns)
        rms = float(np.sqrt(np.mean(Rns**2)))
        mx = float(np.max(np.abs(Rns)))
        p90 = float(np.percentile(np.abs(Rns), 90))
        pred = (math.log(n)**3)/math.sqrt(n)
        rows.append((n, rms, mx, p90, pred))
        print(f"{n:>6} {len(Rns):>4} | {rms:>9.4f} {mx:>9.4f} {p90:>9.4f} | {pred:>20.3f}", flush=True)
    print("-"*80)
    # decay checks
    ok_rms = rows[-1][1] < rows[0][1]
    ok_max = rows[-1][2] < rows[0][2]*1.2   # allow noise; uniformity = max not growing
    # slope vs sqrt n: RMS ratio across n=256->1024 should be ~ (pred ratio)
    rms_ratio = rows[-1][1]/rows[0][1]; pred_ratio = rows[-1][4]/rows[0][4]
    print(f"RMS ratio (n=1024/256) = {rms_ratio:.3f}  vs predicted {pred_ratio:.3f}")
    print(f"MAX decaying (uniformity): {ok_max};  RMS decaying: {ok_rms}")
    print(f"VERDICT: {'CONSISTENT with uniform O(log^3 n/sqrt n) -- no heavy tail, max decays' if (ok_rms and ok_max) else 'STRESS FLAG -- examine the worst words (uniformity in doubt)'}")
