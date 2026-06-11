#!/usr/bin/env python3
"""
ADVERSARIAL verification of the GAP-1 contour-shift analysis (Remark 7.34b).
Independent re-derivation, NOT trusting the committed/draft text.

Attack vectors:
 (1) complex-strip analyticity: is |eta(z)| <= r_KP genuinely maintained on the
     shifted contour, and is the "4A(D)|z|" complex-slope constant correct/usable?
 (2) saddle: residual linear term order arithmetic.
 (3) cubic-error integral: re-do independently.
 (4) vertical segments: sign/absolute bound; does e^{C2 n s0^2} eat e^{-cn}?
 (5) lemma final inequality on shell words at n=512.

Uses /Users/para/.venvs/rnr/bin/python (numpy + mpmath).
"""
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 50


def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))


# ---------------------------------------------------------------------------
# The dual closed forms.  s is the spectral variable; the posterior char fn is
#   Phi(s;x) = (C_s/C_0)^n S_x(eta_s),  with eta_s, C_s/C_0 below.
# The committed forms (from tex line 12690):
#   eta_s = D(1-D)(w-1)/((1-D)^2 - D^2 w),   w = e^{is}
#   C_s/C_0 = ((1-D)^2 - D^2 w)/(1-2D)
# We use COMPLEX s = sr + i*si throughout (the shifted contour has si = s0).
# ---------------------------------------------------------------------------
def eta_closed(p, D, s):
    """eta_s with COMPLEX s, using the committed closed form directly."""
    w = mp.e ** (1j * s)
    num = D * (1 - D) * (w - 1)
    den = (1 - D) ** 2 - D ** 2 * w
    return num / den


def C_ratio_closed(p, D, s):
    w = mp.e ** (1j * s)
    return ((1 - D) ** 2 - D ** 2 * w) / (1 - 2 * D)


def A_of_D(D):
    return D * (1 - D) / (1 - 2 * D)


def r_lb(p):
    # KP analytic lower bound used in the draft (Step 5 DAG)
    return p ** 2 / (12 * math.e * (1 - p) ** 2)


# ===========================================================================
# ATTACK (1): the complex slope.  The draft claims |eta_z| <= 4 A(D) |z| on the
# strip Sigma.  Test whether this holds, and whether a usable bound keeping
# |eta_z| <= r_KP is achievable on the contour we ACTUALLY shift to.
# ===========================================================================
def attack1_complex_slope(p, D, deltas, simaxes, verbose=True):
    A = A_of_D(D)
    rlb = r_lb(p)
    print("=" * 78)
    print(f"ATTACK (1) complex slope & strip containment:  p={p}, D={D:.5f}")
    print(f"  A(D)={A:.5f}   r_lb(p)={rlb:.6f}")
    worst_ratio = 0.0
    worst_at = None
    for delta in deltas:
        for simax in simaxes:
            # scan the rectangle |sr|<=delta, |si|<=simax (the strip we use)
            grid_sr = np.linspace(-delta, delta, 41)
            grid_si = np.linspace(-simax, simax, 41)
            maxeta = 0.0
            max_ratio = 0.0
            for sr in grid_sr:
                for si in grid_si:
                    z = mp.mpf(sr) + 1j * mp.mpf(si)
                    if abs(z) == 0:
                        continue
                    e = abs(eta_closed(p, D, z))
                    maxeta = max(maxeta, float(e))
                    ratio = float(e / abs(z))  # |eta_z|/|z|: the realized slope
                    if ratio > max_ratio:
                        max_ratio = ratio
            if max_ratio > worst_ratio:
                worst_ratio = max_ratio
                worst_at = (delta, simax)
            # the draft's claim: maxeta <= 4 A |z|_edge, and <= r_KP if window small
            z_edge = math.hypot(delta, simax)
            bound_4A = 4 * A * z_edge
            inside_rkp = maxeta < rlb
            if verbose:
                print(f"  delta={delta:.4f} simax={simax:.4f} | "
                      f"max|eta|={maxeta:.5f}  realized_slope={max_ratio:.4f}  "
                      f"4A={4*A:.4f}  4A|z_edge|={bound_4A:.5f}  "
                      f"|eta|<r_lb? {inside_rkp}")
    print(f"  --> WORST realized slope |eta_z|/|z| over scanned strips: {worst_ratio:.4f} "
          f"at (delta,simax)={worst_at};  compare 4A={4*A:.4f}, A={A:.4f}")
    # Verdict on the slope constant
    if worst_ratio <= A * 1.001:
        print(f"  SLOPE VERDICT: realized slope <= A(D) (the REAL-axis slope) -- "
              f"the '4A' bound is LOOSE but valid; real slope ~ A even off-axis.")
    elif worst_ratio <= 4 * A * 1.001:
        print(f"  SLOPE VERDICT: realized slope in (A, 4A] -- 4A bound holds, A alone would FAIL.")
    else:
        print(f"  SLOPE VERDICT: *** realized slope EXCEEDS 4A *** -- the draft's "
              f"complex-slope constant is INVALID.")
    return worst_ratio


# ===========================================================================
# Independent eval of log Phi(s;x) via the O(n) polymer recursion (complex s).
# Cross-check the closed-form eta against the recursion's S_x.
# ===========================================================================
def logS_mp(x, p, eta):
    n = len(x)
    rho = mp.mpf(p) / (1 - p)
    l = [mp.mpf(1)] * n
    r = [mp.mpf(1)] * n
    for i in range(1, n):
        l[i] = rho if x[i] == x[i - 1] else 1 / rho
    for b in range(0, n - 1):
        r[b] = rho if x[b + 1] == x[b] else 1 / rho
    Fip1 = mp.mpf(1)
    Fip2 = mp.mpf(1)
    Gip1 = mp.mpf(0)
    logsc = mp.mpf(0)
    for i in range(n - 1, -1, -1):
        Gi = eta * (r[i] * Fip2 + Gip1)
        Fi = Fip1 + l[i] * Gi
        Fip2 = Fip1
        Fip1 = Fi
        Gip1 = Gi
        a = abs(Fip1)
        if a > 0:
            Fip1 /= a
            Fip2 /= a
            Gip1 /= a
            logsc += mp.log(a)
    return mp.log(Fip1) + logsc


def logPhi(x, p, D, s):
    eta = eta_closed(p, D, s)
    Cr = C_ratio_closed(p, D, s)
    return len(x) * mp.log(Cr) + logS_mp(x, p, eta)


def cumulants(x, p, D, h=mp.mpf('1e-6')):
    f0 = logPhi(x, p, D, 0)
    fp = logPhi(x, p, D, h)
    fm = logPhi(x, p, D, -h)
    d1 = (fp - fm) / (2 * h)            # = i mu
    d2 = (fp - 2 * f0 + fm) / (h * h)   # = -v
    fpp = logPhi(x, p, D, 2 * h)
    fmm = logPhi(x, p, D, -2 * h)
    d3 = (fpp - 2 * fp + 2 * fm - fmm) / (2 * h ** 3)  # = -i v3 (third cumulant up to i)
    mu = float((d1 / 1j).real)
    v = float((-d2).real)
    c3 = float(abs(d3))   # |Phi'''(0)|
    return mu, v, c3


# ===========================================================================
# ATTACK (2): saddle residual linear-term order arithmetic.
#   Lambda'(0) = -1/2 psi'''(xi) s0^2,  |Lambda'(0)| <= 1/2 C3 n s0^2 = O(beta^2/sqrt n).
#   Relative contribution to I: Lambda'(0) * Gaussian width (1/sqrt v_x)
#     = O(beta^2 C3 / (vbar^{3/2} sqrt n)).
# ===========================================================================
def attack2_saddle_linear(x, p, D):
    n = len(x)
    mu, v, c3 = cumulants(x, p, D)
    t = math.floor(n * D)
    beta = (t - mu) / math.sqrt(v)
    s0 = (mu - t) / v
    # measured residual linear term after shift
    h = mp.mpf('1e-6')
    g = lambda s: logPhi(x, p, D, s + 1j * s0) - 1j * t * (s + 1j * s0)
    dpsi = (g(h) - g(-h)) / (2 * h)
    meas = float(abs(dpsi))
    pred = c3 * s0 * s0 / 2  # = 1/2 |Phi'''| s0^2
    # relative contribution: Lambda'(0) * (1/sqrt v)  (the Gaussian width)
    rel = meas / math.sqrt(v)
    rel_pred_form = beta * beta * c3 / (n) / (v / n) ** 1.5 / math.sqrt(n)
    return dict(n=n, mu=mu, v=v, c3=c3, beta=beta, s0=s0,
                lin_meas=meas, lin_pred=pred, rel=rel,
                rel_form_pred=rel_pred_form)


# ===========================================================================
# ATTACK (3): cubic error integral, re-done INDEPENDENTLY (analytic, not
# trusting probe_7_34_gap1_cubic_remainder_form.py).
#   rho_1 = (1/Z) int |E(s)| e^{|E|} e^{-v s^2/2} ds,   E = (1/6) C3 n s^3 (worst).
# Claim:  |rho_1| <= (2 C3 n / 6) * sqrt(v)/sqrt(2pi) * int |s|^3 e^{-v s^2/2} ds
#                  = (2 C3 n /6)*(2/v^2)*sqrt(v)/sqrt(2pi) = O(C3/(vbar^{3/2} sqrt n)).
# Re-derive the closed form and compare to direct numeric integral.
# ===========================================================================
def attack3_cubic_integral(C3, vbar, n, beta=0.0, efactor=2.0):
    v = vbar * n
    C3n = C3 * n
    s0 = abs(beta) / math.sqrt(v)
    # direct numeric integral of |E| e^{|E|} e^{-v s^2/2}, E=(1/6)C3n(|s|+|s0|)^3
    s = np.linspace(-8 / math.sqrt(v), 8 / math.sqrt(v), 8001)
    E = (1 / 6.0) * C3n * (np.abs(s) + s0) ** 3
    integ = np.abs(E) * np.exp(np.minimum(np.abs(E), 30)) * np.exp(-v * s * s / 2)
    Z = np.trapezoid(np.exp(-v * s * s / 2), s)
    rho_direct = np.trapezoid(integ, s) / Z
    # closed-form analytic majorant (beta=0, e^{|E|}<=efactor):
    #   rho <= efactor * (C3n/6) * E[|s|^3]_{N(0,1/v)} = efactor*(C3n/6)*(2/sqrt(pi))*(2/v)^{3/2}/...
    # E|s|^3 for N(0,sigma^2), sigma^2=1/v: = 2 sqrt(2/pi) sigma^3 = 2 sqrt(2/pi) v^{-3/2}
    Eabs3 = 2 * math.sqrt(2 / math.pi) * v ** (-1.5)
    rho_closed = efactor * (C3n / 6.0) * Eabs3
    return dict(n=n, v=v, beta=beta, rho_direct=float(rho_direct),
                rho_closed=rho_closed,
                rho_direct_times_sqrtn=float(rho_direct) * math.sqrt(n),
                rho_closed_times_sqrtn=rho_closed * math.sqrt(n))


# ===========================================================================
# ATTACK (4): vertical segments.  Does e^{C3 n |s0|} * e^{t|s0|} eat e^{-cn}?
#   n |s0|^2 = O(1) (s0 ~ beta/sqrt(vbar n)), but the draft writes the vertical
#   travel as e^{C3 n |s0|} (NOT n s0^2!).  C3 n |s0| = C3 n beta/sqrt(vbar n)
#     = C3 sqrt(n) beta/sqrt(vbar) -- this GROWS like sqrt(n).  Check whether the
#   zone-2 domination -vbar n eps^2/4 (Theta(n)) beats it.
# ===========================================================================
def attack4_vertical(p, D, C3, vbar, eps, n, beta):
    A = A_of_D(D)
    v = vbar * n
    s0 = abs(beta) / math.sqrt(v)
    # vertical exponent budget at edge |s|=eps:
    zone2 = -0.25 * vbar * n * eps ** 2     # -log|Phi(eps)| >= vbar n eps^2/4
    travel_C3 = C3 * n * s0                  # C3 n |s0|  (the "dangerous" term)
    travel_t = math.floor(n * D) * s0        # t|s0| ~ O(beta)
    total = zone2 + travel_C3 + travel_t
    return dict(n=n, beta=beta, eps=eps, zone2=zone2,
                travel_C3=travel_C3, travel_t=travel_t, total=total,
                closes=(total < 0))


if __name__ == "__main__":
    print("#" * 78)
    print("# GAP-1 ADVERSARIAL CONTOUR VERIFICATION")
    print("#" * 78)

    # ---- ATTACK (1): complex slope, several p,D ----
    for p in (0.4, 0.25, 0.1, 0.05):
        D = 0.9 * Dc(p)
        attack1_complex_slope(p, D,
                              deltas=[0.02, 0.05, 0.10],
                              simaxes=[0.02, 0.05],
                              verbose=(p == 0.4))
        print()

    # also test near the boundary D = D_c (where the strip can be tightest)
    print("### near-boundary D = 0.99 D_c (slope stress) ###")
    for p in (0.4, 0.1):
        D = 0.99 * Dc(p)
        attack1_complex_slope(p, D, deltas=[0.05, 0.10], simaxes=[0.05], verbose=True)
        print()
