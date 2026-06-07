#!/usr/bin/env python3
"""
probe_7_34_route_b_crux.py
============================================================================
THE CRUX TEST for Route B (P2c) of the BSMS Gray-region RD-dispersion
ACHIEVABILITY (Remark 7.34b, sole open gap).

Route B target: Var_{x~BSMS}(log2 q_t(X^n)) = O(1) uniformly in n, where
q_t(x) = tilted boundary local mass at t=floor(nD). With E2-E4 this gives
Var(G_n)/n -> V_lossless, hence V_op <= V_lossless (achievability) modulo
the deterministic (1/2)log n prefactor.

The proof route is the AVERAGED-DIRICHLET Markov-Poincare inequality
    Var(f) <= (1/gap) * E[ sum_i Var(f | X_{-i}) ]            (*)
with f = log2 q_t, gap = O(1) [1D Gibbs / Martinelli], and the CRUX (P2c) is
    E[ sum_i Var(log2 q_t | X_{-i}) ] = O(1)  uniformly in n.

This probe measures, IN THE NEAR-MEAN regime (t=floor(nD)>=8), via the exact
O(n^2) Walsh/Krawtchouk transfer-matrix local mass:

 (1) The AVERAGED-DIRICHLET (conditional-variance) energy of the path law,
     EXACT for the 1D Gibbs measure: conditioned on X_{-i}=x_{-i}, the single
     coordinate X_i is a 2-point law with flip-prob r_i (computed exactly from
     the BSMS nearest-neighbour conditional P(X_i | X_{i-1},X_{i+1})). Then
        Var(f | X_{-i}) = r_i (1-r_i) (D_i f)^2,   D_i f = f(x) - f(x^{flip i}),
     and the energy is  Edir = E_x[ sum_i r_i(1-r_i) (D_i f)^2 ].
     [This is the EXACT object in (*); it is <= (1/4) E[sum (D_i f)^2], the
      hard-flip Efron-Stein energy -- same order, but the conditional-variance
      form is what the Poincare inequality uses.]

 (2) The per-coordinate influence DISTRIBUTION |D_i log2 q_t|: quantiles
     (median, 90%, 99%, max) and the fraction of (x,i) pairs with |D_i|>c.
     CRUX claim: |D_i| = O(n^{-1/2}) on TYPICAL configs (median ~ 1/sqrt n),
     with Theta(1) influence only on a SMALL fraction of (x,i) pairs.
     If the energy is O(1) it must be that the typical |D_i| ~ n^{-1/2} (n of
     them, each^2 ~ 1/n => sum O(1)) and the heavy tail has vanishing weight.

 (3) Direct Var_x(log2 q_t) over the source shell vs n (the target quantity).

VERDICT LOGIC:
  - Edir BOUNDED (O(1), slope~0) AND Var(log2 q_t) bounded => (*) is CONSISTENT:
    the averaged-Dirichlet Poincare CAN deliver the O(1) bound (modulo proving
    gap=O(1) and the energy bound rigorously).
  - If Edir = Theta(n): Route B FAILS the same way as the worst-case route.

Run:  /tmp/rnr_venv/bin/python scripts/verify/probe_7_34_route_b_crux.py [p]
Deps: numpy
"""
import sys, math
import numpy as np
from math import comb


def D_c(p):
    return 0.5 * (1.0 - math.sqrt(1.0 - 2.0 * p) / (1.0 - p))


def V_lossless(p):
    return p * (1 - p) * (math.log2((1 - p) / p)) ** 2


def H_coeffs_f(x, p):
    n = len(x); half = 0.5; x0 = int(x[0])

    def mul(c, s):
        m = len(c); o = [0.0] * (m + 1); o[0] = c[0]
        for k in range(1, m):
            o[k] = c[k] + s * c[k - 1]
        o[m] = s * c[m - 1]; return o

    v = [mul([half], 1 if x0 == 0 else -1), mul([half], 1 if x0 == 1 else -1)]
    P = p; Q = 1 - p
    for i in range(1, n):
        xi = int(x[i]); vb, va = v[0], v[1]; L = len(vb)
        a0 = [Q * vb[k] + P * va[k] for k in range(L)]
        a1 = [P * vb[k] + Q * va[k] for k in range(L)]
        v = [mul(a0, 1 if xi == 0 else -1), mul(a1, 1 if xi == 1 else -1)]
    Qc = [v[0][k] + v[1][k] for k in range(len(v[0]))]
    while len(Qc) < n + 1:
        Qc.append(0.0)
    return np.array(Qc[:n + 1])


_KC = {}
def kraw(n):
    if n in _KC:
        return _KC[n]
    K = np.zeros((n + 1, n + 1))
    for j in range(n + 1):
        for k in range(n + 1):
            s = 0; lo = max(0, k - (n - j)); hi = min(k, j)
            for l in range(lo, hi + 1):
                s += ((-1) ** l) * comb(j, l) * comb(n - j, k - l)
            K[k, j] = s
    _KC[n] = K
    return K


def logqt(x, p, D, Km, th0, t):
    n = len(x); H = H_coeffs_f(x, p)
    inv = 1.0 / (1 - 2 * D); ip = inv ** np.arange(n + 1)
    pis = (Km * ip[None, :]) @ H / (2.0 ** n)
    w = pis * np.exp(th0 * np.arange(n + 1))
    M = w.sum()
    if M <= 0 or pis[t] <= 0:
        return None
    qt = pis[t] * math.exp(th0 * t) / M
    if qt <= 0:
        return None
    return math.log2(qt)


def sample(n, p, rng):
    st = (rng.random(n) < p).astype(np.int8); st[0] = rng.random() < 0.5
    return np.cumsum(st) % 2


def flip_prob(x, i, p):
    """EXACT BSMS conditional flip-probability r_i = P(X_i = 1-x_i | X_{-i}).
    BSMS: P(X) propto prod transitions, transition prob p (switch) or 1-p (stay),
    first symbol uniform. For interior i: the two neighbours (x_{i-1},x_{i+1}) and
    the two transitions (i-1->i, i->i+1) determine the 2-point law of X_i.
    For an endpoint, only one neighbour/transition. Returns flip-prob in [0,1]."""
    n = len(x)
    def edge(a, b):  # weight of transition a->b
        return p if a != b else (1 - p)
    xi = int(x[i]); xf = 1 - xi
    if i == 0:
        # uniform * transition to x[1]
        w_cur = 0.5 * edge(xi, int(x[1])) if n > 1 else 0.5
        w_flip = 0.5 * edge(xf, int(x[1])) if n > 1 else 0.5
    elif i == n - 1:
        w_cur = edge(int(x[i - 1]), xi)
        w_flip = edge(int(x[i - 1]), xf)
    else:
        w_cur = edge(int(x[i - 1]), xi) * edge(xi, int(x[i + 1]))
        w_flip = edge(int(x[i - 1]), xf) * edge(xf, int(x[i + 1]))
    return w_flip / (w_cur + w_flip)


def main():
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
    D = 0.9 * D_c(p)
    th0 = math.log(D / (1 - D))
    print(f"# BSMS p={p}  D=0.9 D_c={D:.5f}  V_lossless={V_lossless(p):.5f}  theta0={th0:.3f}")
    print(f"# ROUTE B CRUX: averaged-Dirichlet (conditional-variance) energy of log2 q_t,")
    print(f"#   Edir = E[sum_i r_i(1-r_i)(D_i log2 q_t)^2], r_i = exact BSMS flip-prob.")
    print(f"#   Poincare: Var(log2 q_t) <= (1/gap) Edir. Want Edir=O(1) & |D_i|~n^-1/2 typ.")
    print("#" * 100)
    print(f"{'n':>4} {'t':>3} {'nD':>6} {'kE':>3} | {'Edir':>8} {'Edir/n':>8} {'Ehard':>8} | "
          f"{'Var(lq)':>8} | {'|Di| med':>8} {'90%':>7} {'99%':>7} {'max':>7} "
          f"{'sqrt(n)*med':>11} {'P(|Di|>.1)':>10}")
    rows = []
    for n in (80, 100, 120, 140, 160, 200):
        t = int(math.floor(n * D))
        if t < 8:
            print(f"{n:>4} t={t}<8 skip"); continue
        Km = kraw(n)
        rng = np.random.default_rng(13000 + n)
        # direct variance over many samples
        lqs = []
        for _ in range(1000):
            x = sample(n, p, rng); v = logqt(x, p, D, Km, th0, t)
            if v is not None:
                lqs.append(v)
        Vlq = float(np.var(lqs, ddof=1))
        # influence + conditional-variance energy over fewer samples (n flips each)
        Edir_list = []; Ehard_list = []
        allDi = []
        nrep = 80 if n <= 140 else 50
        for _ in range(nrep):
            x = sample(n, p, rng); lq0 = logqt(x, p, D, Km, th0, t)
            if lq0 is None:
                continue
            Di = np.zeros(n); rr = np.zeros(n); ok = True
            for i in range(n):
                xf = x.copy(); xf[i] ^= 1; lf = logqt(xf, p, D, Km, th0, t)
                if lf is None:
                    ok = False; break
                Di[i] = lq0 - lf
                rr[i] = flip_prob(x, i, p)
            if not ok:
                continue
            cv = rr * (1 - rr) * Di ** 2          # conditional variance per coord
            Edir_list.append(float(np.sum(cv)))
            Ehard_list.append(float(np.sum(Di ** 2)))
            allDi.append(np.abs(Di))
        Edir = float(np.mean(Edir_list)); Ehard = float(np.mean(Ehard_list))
        allDi = np.concatenate(allDi)
        med = float(np.median(allDi)); q90 = float(np.quantile(allDi, 0.90))
        q99 = float(np.quantile(allDi, 0.99)); mx = float(allDi.max())
        pheavy = float(np.mean(allDi > 0.1))
        rows.append((n, t, Edir, Ehard, Vlq, med, q90, q99, mx, pheavy))
        print(f"{n:>4} {t:>3} {n*D:>6.2f} {len(Edir_list):>3} | {Edir:>8.4f} {Edir/n:>8.5f} "
              f"{Ehard:>8.4f} | {Vlq:>8.4f} | {med:>8.4f} {q90:>7.4f} {q99:>7.4f} {mx:>7.4f} "
              f"{math.sqrt(n)*med:>11.4f} {pheavy:>10.4f}", flush=True)
    if len(rows) >= 3:
        ns = np.array([r[0] for r in rows], float)
        Edir = np.array([r[2] for r in rows]); Vlq = np.array([r[4] for r in rows])
        med = np.array([r[5] for r in rows]); pheavy = np.array([r[9] for r in rows])
        sE = np.polyfit(ns, Edir, 1)[0]; sV = np.polyfit(ns, Vlq, 1)[0]
        # is sqrt(n)*median bounded (=> median ~ c/sqrt(n))?
        snmed = math.sqrt(1.0) * (np.sqrt(ns) * med)
        print("-" * 100)
        print(f"Edir slope vs n = {sE:+.5f}  (O(1) if ~0; values {[round(r[2],3) for r in rows]})")
        print(f"Edir BOUNDED (=> averaged-Dirichlet Poincare can give Var=O(1)): "
              f"{'YES' if abs(sE) < 0.02 else 'NO -> Theta(n)'}")
        print(f"Var(log2 q_t) slope vs n = {sV:+.6f}  (bounded if ~0)")
        print(f"sqrt(n)*median|D_i| = {[round(float(v),3) for v in (np.sqrt(ns)*med)]}  "
              f"(BOUNDED => typical influence ~ c/sqrt(n))")
        print(f"P(|D_i|>0.1) = {[round(float(v),4) for v in pheavy]}  "
              f"(vanishing => heavy influence on o(1) set)")
        print()
        print("READING:")
        print("  (P2c) requires the averaged-Dirichlet energy Edir = E[sum_i Var(log2 q_t|X_-i)]")
        print("  to be O(1). The conditional-variance form r_i(1-r_i)(D_i)^2 is the EXACT object the")
        print("  Markov-Poincare inequality (*) integrates. If Edir is bounded AND the typical |D_i|")
        print("  scales as c/sqrt(n) (n of them, each^2 ~ 1/n => sum O(1)) with heavy influence only")
        print("  on an o(1)-probability set, the mechanism in (P2c) is numerically supported: the")
        print("  Theta(1) per-coordinate influence at atypical (frozen-radius/boundary-straddling)")
        print("  configs is tolerated by the AVERAGED (L2) energy because those configs have")
        print("  vanishing probability. This is exactly the gap between the averaged-Dirichlet")
        print("  Poincare and the worst-case McDiarmid (which sees the Theta(1) sup-influence).")


if __name__ == "__main__":
    main()
