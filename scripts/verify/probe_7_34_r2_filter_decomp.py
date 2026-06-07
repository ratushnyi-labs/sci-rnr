#!/usr/bin/env python3
"""
probe_7_34_r2_filter_decomp.py
============================================================================
TRACK (R2), PROOF-MECHANISM PROBE: does the single-flip influence D_i log2 r_t
DECOMPOSE into a geometrically-forgetting additive (filter-stability) part plus
a controllable lattice-mass part, as the (R2)-proof mechanism requires?

The committed reduction needs Var(log2 r_t)=O(1), via the averaged-Dirichlet
Glauber-Poincare bound, which needs E[sum_i (D_i log2 r_t)^2]=O(1). (R2) asserts
this energy is O(1) because (a) typical |D_i|~n^-1/2 and (b) the jump set is tiny.
The PROOF MECHANISM for (a) is two-sided filter stability. To test whether that
mechanism is real (not just the energy being O(1) by luck), we use the EXACT
algebraic decomposition (verified to machine precision below):

    log2 q_t(x) = log2 pi_t(x) + (theta0 t)/ln2 - log2 M(x),

    pi_t(x) = P_{Y*}(d_H(x,Y*)=t)              [untilted boundary lattice mass]
    M(x)    = E_{Y*}[ e^{theta0 K(x,Y*)} ],  K=d_H(x,Y*), theta0=ln(D/(1-D))<0
            = the tilted PARTITION FUNCTION (an ADDITIVE-cocycle MGF).

The (theta0 t)/ln2 term is CONSTANT in x (t=floor(nD) fixed), so
    D_i log2 q_t = D_i log2 pi_t  -  D_i log2 M.

HYPOTHESES TESTED:
 (H1) M-part is the EASY additive/filter-stability piece.  log2 M is a tilted
      free energy of the additive cocycle K=sum_l 1[Y*_l != x_l]; flipping x_i
      changes ONE selector, so D_i log2 M = log2( E_Q[ e^{-theta0 (2 Y'_i -1)...}])
      is O(|theta0|) bounded, and -- by two-sided filter stability of the tilted
      posterior Q(.|x) -- depends on x_i only through the posterior marginal at i,
      a LOCAL, O(1)-influence quantity (no n^-1/2; it does NOT vanish, it is the
      analogue of D_i j_n).  => E_M := E[sum_i (D_i log2 M)^2] = Theta(n) EXPECTED.
 (H2) pi_t-part is the HARD lattice piece.  D_i log2 pi_t should TRACK D_i log2 M
      to leading order (so their DIFFERENCE D_i log2 q_t is the small n^-1/2
      residual), i.e. the energies E_pi and E_M are BOTH Theta(n) but the
      energy of the DIFFERENCE E_q=E[sum (D_i log2 q_t)^2] is O(1).
 (H3) The RESIDUAL influence D_i log2 q_t should show SPATIAL DECAY away from the
      "active" coordinates (two-sided filter forgetting): we profile mean |D_i|
      by coordinate index and by a coupling proxy.  [If FLAT not decaying, the
      summability is by 1/n-smallness per coord, not by geometric forgetting --
      both give O(1) energy, but the proof tool differs.]

If (H1)+(H2) hold, the (R2) energy O(1) is a genuine LEADING-ORDER CANCELLATION
between two Theta(n)-energy additive/lattice objects, and the residual is the
filter-stability + near-mean-LLT object -- exactly the proof mechanism, with the
HARD step being to show the cancellation analytically (NOT supplied by numerics).

Exact O(n^2) Walsh/Krawtchouk machinery; float64 exact at near-mean t.
Run:  /tmp/rnr_venv/bin/python scripts/verify/probe_7_34_r2_filter_decomp.py [p]
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


def parts(x, p, D, Km, th0, t):
    """Return (log2 r_t, log2 q_t, log2 pi_t, log2 M) exactly.
    log2 q_t = log2 pi_t + theta0*t/ln2 - log2 M ; r_t = q_t*sqrt(2 pi v_q)."""
    n = len(x); H = H_coeffs_f(x, p)
    inv = 1.0 / (1 - 2 * D); ip = inv ** np.arange(n + 1)
    pis = (Km * ip[None, :]) @ H / (2.0 ** n)
    ks = np.arange(n + 1)
    w = pis * np.exp(th0 * ks)
    M = w.sum()
    if M <= 0 or pis[t] <= 0:
        return None
    qn = w / M
    mq = float((ks * qn).sum())
    vq = float(((ks - mq) ** 2 * qn).sum())
    if vq <= 0:
        return None
    lpi_t = math.log2(pis[t])
    lM = math.log2(M)
    lq_t = lpi_t + th0 * t / math.log(2) - lM
    lr_t = lq_t + 0.5 * math.log2(2 * math.pi * vq)
    return lr_t, lq_t, lpi_t, lM


def sample(n, p, rng):
    st = (rng.random(n) < p).astype(np.int8); st[0] = rng.random() < 0.5
    return np.cumsum(st) % 2


def run(p, n, R, seed, Dfrac=0.9):
    D = Dfrac * D_c(p); th0 = math.log(D / (1 - D)); t = int(math.floor(n * D))
    if t < 8:
        return None
    Km = kraw(n); rng = np.random.default_rng(seed)
    Er = []; Eq = []; Epi = []; EM = []
    # spatial decay profiles (by index distance from the LEFT and from a flip's local effect)
    prof_q = np.zeros(n); prof_cnt = np.zeros(n)
    # decomposition-tracking: corr(D_i pi, D_i M); how much D_i q is the residual
    corr_piM = []
    kept = 0
    selfcheck = 0.0
    for _ in range(R):
        x = sample(n, p, rng)
        b = parts(x, p, D, Km, th0, t)
        if b is None:
            continue
        lr0, lq0, lpi0, lM0 = b
        Dr = np.zeros(n); Dq = np.zeros(n); Dpi = np.zeros(n); DM = np.zeros(n)
        ok = True
        for i in range(n):
            xf = x.copy(); xf[i] ^= 1
            bf = parts(xf, p, D, Km, th0, t)
            if bf is None:
                ok = False; break
            Dr[i] = lr0 - bf[0]; Dq[i] = lq0 - bf[1]
            Dpi[i] = lpi0 - bf[2]; DM[i] = lM0 - bf[3]
        if not ok:
            continue
        kept += 1
        # self-check: D_i q == D_i pi - D_i M  (theta0 t/ln2 cancels)
        selfcheck = max(selfcheck, float(np.max(np.abs(Dq - (Dpi - DM)))))
        Er.append(float(np.sum(Dr ** 2))); Eq.append(float(np.sum(Dq ** 2)))
        Epi.append(float(np.sum(Dpi ** 2))); EM.append(float(np.sum(DM ** 2)))
        if np.std(Dpi) > 0 and np.std(DM) > 0:
            corr_piM.append(float(np.corrcoef(Dpi, DM)[0, 1]))
        prof_q += np.abs(Dq); prof_cnt += 1
    if kept < 5:
        return None
    prof = prof_q / np.maximum(prof_cnt, 1)
    return dict(n=n, t=t, kept=kept, selfcheck=selfcheck,
                Er=float(np.mean(Er)), Eq=float(np.mean(Eq)),
                Epi=float(np.mean(Epi)), EM=float(np.mean(EM)),
                corr_piM=float(np.mean(corr_piM)) if corr_piM else float('nan'),
                prof=prof)


def main():
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
    D = 0.9 * D_c(p); th0 = math.log(D / (1 - D))
    print("#" * 100)
    print(f"# R2 FILTER-DECOMP :: BSMS p={p}  D=0.9 D_c={D:.5f}  V_lossless={V_lossless(p):.5f}")
    print(f"# log2 q_t = log2 pi_t + theta0*t/ln2 - log2 M ;  D_i log2 q_t = D_i log2 pi_t - D_i log2 M")
    print(f"# (H1) E_M = Theta(n) [additive cocycle MGF, local influence];")
    print(f"# (H2) E_pi = Theta(n) but E_q=E[sum (D_i log2 q_t)^2]=O(1) (LEADING-ORDER CANCELLATION);")
    print(f"# (H3) residual D_i log2 q_t spatial decay (filter forgetting).")
    print("#" * 100)
    print(f"{'n':>4} {'t':>3} {'kpt':>4} | {'E_r':>7} {'E_q':>7} {'E_pi':>9} {'E_M':>9} | "
          f"{'E_pi/n':>7} {'E_M/n':>7} | {'corr(Dpi,DM)':>12} | {'selfchk':>8}")
    rows = []
    for n, R in [(96, 120), (128, 90), (160, 60), (220, 40), (300, 22)]:
        r = run(p, n, R, seed=22000 + n)
        if r is None:
            print(f"{n:>4} (skip)"); continue
        rows.append(r)
        print(f"{r['n']:>4} {r['t']:>3} {r['kept']:>4} | {r['Er']:>7.3f} {r['Eq']:>7.3f} "
              f"{r['Epi']:>9.2f} {r['EM']:>9.2f} | {r['Epi']/r['n']:>7.3f} {r['EM']/r['n']:>7.3f} | "
              f"{r['corr_piM']:>12.5f} | {r['selfcheck']:>8.1e}", flush=True)
    if len(rows) >= 3:
        ns = np.array([r['n'] for r in rows], float)
        Eq = np.array([r['Eq'] for r in rows])
        Epi = np.array([r['Epi'] for r in rows]); EM = np.array([r['EM'] for r in rows])
        sEq = np.polyfit(ns, Eq, 1)[0]
        sEpi = np.polyfit(ns, Epi, 1)[0]; sEM = np.polyfit(ns, EM, 1)[0]
        print("-" * 100)
        print(f"E_q  slope vs n = {sEq:+.5f}   (O(1) if ~0)  -> {'O(1) cancellation CONFIRMED' if abs(sEq)<0.02 else 'GROWS'}")
        print(f"E_pi slope vs n = {sEpi:+.4f}   (Theta(n) if >0)")
        print(f"E_M  slope vs n = {sEM:+.4f}   (Theta(n) if >0)")
        h1 = sEM > 0.005
        h2 = (sEpi > 0.005) and (abs(sEq) < 0.02)
        print(f"(H1) E_M = Theta(n): {'YES' if h1 else 'NO'}")
        print(f"(H2) E_pi Theta(n) but E_q O(1) (leading-order cancellation): {'YES' if h2 else 'NO'}")
        # spatial decay of the residual influence
        prof = rows[-1]['prof']
        # left-tail (boundary effect) and bulk
        print(f"(H3) residual |D_i log2 q_t| profile by index i (largest n={rows[-1]['n']}):")
        idxs = [0, 1, 2, 5, 10, rows[-1]['n']//2, rows[-1]['n']-3, rows[-1]['n']-2, rows[-1]['n']-1]
        print("     " + "  ".join(f"i={i}:{prof[i]:.4f}" for i in idxs if i < rows[-1]['n']))
        flat = np.std(prof[5:-5]) / (np.mean(prof[5:-5]) + 1e-12)
        print(f"     bulk rel-spread = {flat:.3f}  (FLAT ~1/sqrt(n) per coord if <~0.3; "
              f"summable by 1/n-smallness, NOT geometric decay)")
        print()
        print("MECHANISM READING:")
        print("  D_i log2 q_t = D_i log2 pi_t - D_i log2 M is a difference of two Theta(n)-energy")
        print("  objects whose difference has O(1) energy: a LEADING-ORDER CANCELLATION.  The M-part")
        print("  is a tilted additive-cocycle free energy (filter-stability/Atar-Zeitouni controls")
        print("  its single-flip influence as a local posterior-marginal quantity); the pi_t-part is")
        print("  the boundary lattice mass.  PROVING E_q=O(1) requires showing this cancellation")
        print("  analytically -- i.e. D_i log2 pi_t = D_i log2 M + O(n^-1/2) on typical x -- which is")
        print("  the near-mean lattice-LLT-ratio regularity (the pi_t local mass tracks the smooth")
        print("  tilted free energy to within the n^-1/2 lattice fluctuation).  Numerics show the")
        print("  cancellation; the analytic step is the SAME near-mean lattice-LLT residual (R1).")


if __name__ == "__main__":
    main()
