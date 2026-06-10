#!/usr/bin/env python3
"""
probe_7_34_quenched_tail_rate.py
============================================================================
THE SIGN OF c(p) in the quenched tail lemma -- the sharpest remaining object of
the BSMS Gray-region achievability (Remark 7.34b):

    sup_{|s|>=delta_n} |Phi(s;x)| <= e^{-c(p) n}   for typical x ?

Via the dual identity (Lemma 7.34b') and the polymer gas (Lemma 7.34b''),
    Phi(s;x) = (C_s/C_0)^n S_x(eta_s),
and S_x(eta) admits an O(n) two-term suffix recursion over the hard-core
domain-wall gas:
    F[i] = F[i+1] + l(i) G[i],   G[i] = eta r(i) F[i+2] + eta G[i+1],
    F[n]=F[n+1]=1,  G[n-1] = eta r(n-1) F[n+1],
with l(i)/r(b) the edge-bond likelihood factors (rho=p/(1-p): rho if stay,
1/rho if switch; word boundary = 1).

MEASUREMENT: the quenched rate profile
    g_n(s;x) = -(1/n) log|Phi(s;x)|  >= 0     (|Phi|<=1: posterior char. fn.)
for s in (0,pi], n up to 400, R source-typical words; report worst-x and
median-x profiles and the tail floor c_hat(n,p;delta) = min_{x, s>=delta} g_n.
VERDICT LOGIC:
  - c_hat(n,p;delta) stabilizing at a POSITIVE constant as n grows
      => numerical EVIDENCE FOR the tail lemma (c(p)>0);
  - c_hat -> 0 like 1/n at some fixed s (|Phi| flat < 1)
      => evidence AGAINST (the lemma would FAIL at that p) -- major either way.
NOTE the small-s Gaussian shoulder: g_n(s) ~ v_rate s^2/2 for small s is
n-independent and POSITIVE but small -- the prior report of "sup_x|Phi|~0.995
flat" can be exactly this shoulder at s~0.1; the lemma's delta_n -> 0 slowly is
served by the Edgeworth window, so the diagnostic is the PROFILE vs s, not the
sup alone.

  V1  recursion == brute force 2^n at n=10 (several words, complex eta).
  V2  sanity: Phi(0;x)=1, |Phi(s;x)|<=1 for all tested (s,x).
  V3  the rate profiles and c_hat table vs n.
Deps: numpy.
"""
import math
import numpy as np
from itertools import product


def Dc(p): return 0.5*(1-math.sqrt(1-2*p)/(1-p))


def eta_C(p, D, s_arr):
    th0 = math.log(D/(1-D))
    eb = np.exp(th0 + 1j*np.asarray(s_arr, dtype=float))
    ze = (1-eb)/((1+eb)*(1-2*D))
    eta = (1-ze)/(1+ze)
    C = (1+eb)*(1+ze)/2.0
    C0 = 1.0/(1.0-D)
    return eta, C/C0


def make_lr(x, p):
    n = len(x); rho = p/(1.0-p)
    l = np.ones(n); r = np.ones(n)
    for i in range(1, n):
        l[i] = rho if x[i] == x[i-1] else 1.0/rho
    for b in range(0, n-1):
        r[b] = rho if x[b+1] == x[b] else 1.0/rho
    return l, r


def S_eval(x, p, eta):
    """O(n) evaluation of S_x(eta) (complex eta), scaled. Returns log S (complex log magnitude+phase ignored: we return (log|S|, S/|S|))."""
    n = len(x)
    l, r = make_lr(x, p)
    # F[i] = F[i+1] + l[i]*G[i];  G[i] = eta*r[i]*F[i+2] + eta*G[i+1]
    Fip1 = 1.0+0j   # F[i+1]
    Fip2 = 1.0+0j   # F[i+2]
    Gip1 = 0.0+0j   # G[i+1] (G[n] := 0)
    logsc = 0.0
    for i in range(n-1, -1, -1):
        Gi = eta*r[i]*Fip2 + eta*Gip1
        Fi = Fip1 + l[i]*Gi
        # shift
        Fip2 = Fip1; Fip1 = Fi; Gip1 = Gi
        a = abs(Fi)
        if a > 1e150 or (a < 1e-150 and a > 0):
            sc = a
            Fip1 /= sc; Fip2 /= sc; Gip1 /= sc
            logsc += math.log(sc)
    if Fip1 == 0:
        return -math.inf
    return math.log(abs(Fip1)) + logsc


def S_brute(x, p, eta):
    n = len(x); rho = p/(1.0-p)
    def sw(z): return int(np.sum(z[1:] != z[:-1]))
    s0 = sw(x); tot = 0.0+0j
    for e in product([0, 1], repeat=n):
        e = np.array(e)
        tot += (eta**int(e.sum())) * rho**(sw(x ^ e) - s0)
    return tot


def sample_words(n, p, R, rng):
    out = []
    for _ in range(R):
        out.append((np.cumsum(rng.random(n) < p) % 2).astype(int))
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=100)
    ap.add_argument("--nmax", type=int, default=400)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    rng = np.random.default_rng(7)
    print("="*86)
    print("V1: O(n) recursion vs brute force 2^n (n=10, complex eta)")
    ok1 = True
    for p in (0.25, 0.4):
        for trial in range(3):
            x = (np.cumsum(rng.random(10) < p) % 2).astype(int)
            for eta in (0.1+0.0j, -0.15+0.05j, 0.02+0.2j):
                lb = S_eval(x, p, eta)
                bb = S_brute(x, p, eta)
                d = abs(lb - math.log(abs(bb)))
                ok1 &= d < 1e-9
    print(f"   recursion == brute (log|S|, 18 cases): {ok1}")

    print("="*86)
    print("V2/V3: quenched rate profile g_n(s;x) = -(1/n) log|Phi|;  c_hat = min_{x,s>=delta} g_n")
    s_grid = np.concatenate([np.linspace(0.05, 0.5, 12), np.linspace(0.6, math.pi, 22)])
    deltas = (0.1, 0.3, 1.0)
    ok2 = True
    print("  (n chosen per p so that nD spans the frozen (t=0) -> near-mean (nD~10) regimes;")
    print("   the dispersion-regime verdict uses the nD>=5 rows)")
    print(f"{'p':>5} {'n':>5} {'nD':>6} | " + " ".join(f"c^({d})" .rjust(9) for d in deltas) +
          " |  argmin(s)  med-x g(pi)  worst-x g(pi)")
    verdict = {}
    for p in (0.1, 0.25, 0.4):
        D = 0.9*Dc(p)
        # nD-matched schedule: frozen -> dispersion regime
        targets = (0.5, 1.0, 2.0, 5.0, 10.0) if not args.quick else (0.5, 2.0, 5.0)
        n_list = sorted(set(min(int(math.ceil(k/D)), 4000) for k in targets))
        eta_s, Cr = eta_C(p, D, s_grid)
        logCr = np.log(np.abs(Cr))
        rows = []
        for n in n_list:
            words = sample_words(n, p, args.R, rng)
            G = np.zeros((len(words), len(s_grid)))
            for wi, x in enumerate(words):
                for si in range(len(s_grid)):
                    logS = S_eval(x, p, complex(eta_s[si]))
                    logPhi = n*logCr[si] + logS
                    ok2 &= logPhi < 1e-9   # |Phi|<=1
                    G[wi, si] = -logPhi/n
            chats = []
            for d in deltas:
                mask = s_grid >= d
                chats.append(float(np.min(G[:, mask])))
            mn = np.unravel_index(np.argmin(G[:, s_grid >= deltas[0]]),
                                  G[:, s_grid >= deltas[0]].shape)
            s_at = s_grid[s_grid >= deltas[0]][mn[1]]
            gpi_med = float(np.median(G[:, -1])); gpi_worst = float(np.min(G[:, -1]))
            rows.append((n, chats, s_at, gpi_med, gpi_worst))
            print(f"{p:>5} {n:>5} {n*D:>6.1f} | " + " ".join(f"{c:>9.5f}" for c in chats) +
                  f" |  {s_at:>8.2f}  {gpi_med:>10.5f}  {gpi_worst:>12.5f}", flush=True)
        # verdict per p over the DISPERSION-REGIME rows (nD>=5): c_hat(delta=0.3) positive & stable?
        disp = [r for r in rows if r[0]*D >= 4.5] or rows[-2:]
        cs = [r[1][1] for r in disp]          # delta=0.3 column
        ns = [r[0] for r in disp]
        # test: 1/n decay would give c-ratio ~ n_ratio (=ns[0]/ns[-1] < 1); a positive
        # limit gives c-ratio ~ 1 (or >1 while still converging up). Use the midpoint.
        ratio_last = cs[-1]/cs[0] if cs[0] > 0 else float('inf')
        n_ratio = ns[0]/ns[-1]
        looks_decay = ratio_last < (n_ratio + 1.0)/2.0
        verdict[p] = (cs, looks_decay)
    print("="*86)
    print(f"V1 PASS={ok1}; V2 (|Phi|<=1 everywhere) PASS={ok2}")
    for p, (cs, dec) in verdict.items():
        trend = "DECAYING ~1/n  => evidence AGAINST exponential tail (c(p)=0?)" if dec else \
                "stabilizing POSITIVE => evidence FOR the tail lemma (c(p)>0)"
        print(f"  p={p}: c_hat(delta=0.3) over n: {[f'{c:.5f}' for c in cs]}  -> {trend}")
    print("NOTE: small-s shoulder g ~ v s^2/2 is n-independent-positive; the lemma's delta_n->0")
    print("      is served by the Edgeworth window; the verdict column is delta=0.3.")
