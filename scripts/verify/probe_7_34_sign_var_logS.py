#!/usr/bin/env python3
"""
probe_7_34_sign_var_logS.py
===========================================================================
PROBE THE SIGN of the open BSMS Hamming RD-dispersion ACHIEVABILITY residual.

Binary symmetric Markov source (BSMS), switch prob p, Hamming distortion.
In Gray's region 0 < D <= D_c(p), D_c(p) = 0.5*(1 - sqrt(1-2p)/(1-p)), the
optimal output Y* is the exact BSC(D) deconvolution of the source:
    X = Y* XOR Z,   Z ~ iid Bern(D),
so in the WALSH domain  Phat_Y(w) = Phat_X(w) (1-2D)^{-|w|}.

For a source word x^n (sampled annealed, x^n ~ BSMS):
    G_n(x)  = -log2 P_{Y*}(B_{nD}(x^n)),  B_t = {y : d_H(x,y) <= t}, t=floor(nD)
    i_n(x)  = -log2 P_X(x^n)
    j_n(x)  = i_n(x) - n*h(D),   h = binary entropy
    log2 S_n(x) = j_n(x) - G_n(x)   (so G_n = j_n - log2 S_n)

THE QUESTION (annealed): does  Var(log2 S_n)/n -> 0  or -> positive const c>0 ?
  -> 0           : V_op = V_lossless = p(1-p) log2^2((1-p)/p) numerically CERTIFIED
                   (achievability closes at the conjectured value).
  -> c > 0       : V_op = V_lossless + c +/- cross-term  => V_op != V_lossless,
                   known converse V_op >= V_lossless would be LOOSE.

EXACT METHOD (O(n^2) per sample, no 2^n blowup):
    P_{Y*}(B_t(x)) = (1/2^n) sum_{j=0}^n (1-2D)^{-j} Gt(j) H_j(x)
  Gt(j) = sum_{i=0}^{t} K_i(j;n)            partial Krawtchouk sum (depends on j only)
  K_i(j;n) = sum_l (-1)^l C(j,l) C(n-j, i-l)  Krawtchouk polynomial
  H_j(x)   = z^j coefficient of  Q_x(z) = sum_{x'} P_X(x') (1+z)^{n-dH(x,x')} (1-z)^{dH(x,x')}
           computed as a 2x2 transfer-matrix product over x' (order-1 Markov):
             Q_x(z) = [1/2, 1/2] (prod_i T_i(z)) [1,1]^T,
             (T_i(z))_{a,b} = P(X'_i=b|X'_{i-1}=a) * ((1+z) if b==x_i else (1-z)).
           For i=1 the transition factor is the stationary pi=[1/2,1/2].

NUMERICS: (1-2D)^{-j} grows and Gt(j),H_j(x) oscillate => catastrophic
cancellation at large n. We use mpmath at high precision and DETECT precision
loss per sample (cancellation magnitude vs working precision); n is capped to
the reliable range automatically and reported honestly.

Dependencies: numpy, mpmath  (e.g.  python -m venv venv && venv/bin/pip install numpy mpmath)
Run:    python scripts/verify/probe_7_34_sign_var_logS.py          # full (~20 min, 7 procs)
Quick:  python scripts/verify/probe_7_34_sign_var_logS.py --quick  # smoke + validation only
"""
import math
import time
import argparse
import numpy as np
import mpmath as mp

# --------------------------------------------------------------------------
# basic scalars
# --------------------------------------------------------------------------
def h2(x):
    if x <= 0 or x >= 1:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)

def D_c(p):
    return 0.5 * (1.0 - math.sqrt(1.0 - 2.0 * p) / (1.0 - p))

def V_lossless(p):
    return p * (1 - p) * (math.log2((1 - p) / p)) ** 2

# --------------------------------------------------------------------------
# stationary BSMS log-prob (exact, cheap)
# --------------------------------------------------------------------------
def i_n_bits(x, p):
    """-log2 P_X(x^n), stationary BSMS switch prob p. X_1 uniform."""
    n = len(x)
    flips = int(np.sum(x[1:] != x[:-1]))
    # P = 0.5 * p^flips * (1-p)^(n-1-flips)
    logP = math.log2(0.5) + flips * math.log2(p) + (n - 1 - flips) * math.log2(1 - p)
    return -logP

# --------------------------------------------------------------------------
# H_j(x): coefficients of Q_x(z) via 2x2 transfer matrix in mpmath polynomials
# --------------------------------------------------------------------------
def _mul_deg1(coeffs, sign):
    """
    Multiply a Python-list polynomial (coeffs low->high, mpf entries) by
    (1 + sign*z), sign in {+1,-1}.  (poly*(1+sz))_k = c_k + sign*c_{k-1}.
    O(deg) shift-and-add — NOT a general convolution.
    NB: plain Python lists, NOT numpy object arrays (the latter trigger a
    pathological mpmath scalar*ndarray convert path that is ~50x slower).
    """
    m = len(coeffs)
    out = [None] * (m + 1)
    out[0] = coeffs[0]
    for k in range(1, m):
        out[k] = coeffs[k] + sign * coeffs[k - 1]
    out[m] = sign * coeffs[m - 1]
    return out

def H_coeffs(x, p):
    """
    Return Python list [H_0,...,H_n] (mpmath mpf) = coefficients of
        Q_x(z) = sum_{x'} P_X(x') prod_i ((1+z) if x'_i==x_i else (1-z)).
    Transfer-matrix product over x' (order-1 Markov, X'_1 uniform):
      v[b] = polynomial (in z) accumulating the partial product for X'_i=b.
      step i: v_new[b] = (sum_a T(a->b) v[a]) * ((1+z) if b==x_i else (1-z)).
    Each factor is degree-1 => O(n) per step, O(n^2) total.
    """
    n = len(x)
    half = mp.mpf('0.5')
    x0 = int(x[0])
    # v[b] = poly for X'_1 = b, weighted by pi=1/2 and the (1+-z) factor.
    v = [_mul_deg1([half], 1 if x0 == 0 else -1),
         _mul_deg1([half], 1 if x0 == 1 else -1)]
    mp_p = mp.mpf(p)
    mp_q = mp.mpf(1) - mp_p
    for i in range(1, n):
        xi = int(x[i])
        vb, va = v[0], v[1]   # current state polys
        # acc_b = q*v[b] + p*v[1-b], elementwise (lists, same length here)
        L = len(vb)
        acc0 = [mp_q * vb[k] + mp_p * va[k] for k in range(L)]   # b=0
        acc1 = [mp_p * vb[k] + mp_q * va[k] for k in range(L)]   # b=1
        v = [_mul_deg1(acc0, 1 if xi == 0 else -1),
             _mul_deg1(acc1, 1 if xi == 1 else -1)]
    Q = [v[0][k] + v[1][k] for k in range(len(v[0]))]
    while len(Q) < n + 1:
        Q.append(mp.mpf(0))
    return Q[:n + 1]

# --------------------------------------------------------------------------
# Krawtchouk partial sums Gt(j) = sum_{i=0}^t K_i(j;n), exact integers
# --------------------------------------------------------------------------
def krawtchouk_partial(n, t, jmax=None):
    """
    Return dict-free list Gt[j] for j=0..n (Python ints, exact).
    K_i(j;n) = sum_{l=0}^{i} (-1)^l C(j,l) C(n-j, i-l).
    Gt(j) = sum_{i=0}^{t} K_i(j;n).
    """
    if jmax is None:
        jmax = n
    from math import comb
    Gt = [0] * (n + 1)
    for j in range(jmax + 1):
        s = 0
        for i in range(0, t + 1):
            Ki = 0
            lo = max(0, i - (n - j))
            hi = min(i, j)
            for l in range(lo, hi + 1):
                Ki += ((-1) ** l) * comb(j, l) * comb(n - j, i - l)
            s += Ki
        Gt[j] = s
    return Gt

# --------------------------------------------------------------------------
# ball probability via the formula (mpmath)
# --------------------------------------------------------------------------
def ball_prob_formula(x, p, D, t, Gt, dps):
    """
    P_{Y*}(B_t(x)) = (1/2^n) sum_j (1-2D)^{-j} Gt(j) H_j(x).
    Returns (prob, cancellation_ratio): ratio = sum|term| / |sum| (precision diag).
    """
    n = len(x)
    H = H_coeffs(x, p)
    inv = mp.mpf(1) / (mp.mpf(1) - 2 * mp.mpf(D))   # (1-2D)^{-1}
    total = mp.mpf(0)
    abssum = mp.mpf(0)
    invj = mp.mpf(1)
    for j in range(n + 1):
        term = invj * Gt[j] * H[j]
        total += term
        abssum += abs(term)
        invj *= inv
    prob = total / (mp.mpf(2) ** n)
    cancel = (abssum / abs(total)) if total != 0 else mp.inf
    return prob, cancel

# --------------------------------------------------------------------------
# BRUTE FORCE (n<=16) — three independent ball-prob computations
# --------------------------------------------------------------------------
def brute_ball_probs(x, p, D, t):
    """
    Returns (P_formula_free, P_walsh, P_direct):
      method (ii) walsh: Phat_Y(w)=Phat_X(w)(1-2D)^{-|w|}, inverse transform to
                  P_{Y*}(y) for all 2^n y, sum over ball.
      method (iii) direct: P_{Y*}(y) = (1/2^n) sum_w Phat_Y(w)(-1)^{<w,y>} (== ii,
                  but computed with explicit per-y Walsh sum to cross-check).
    Also returns min_y P_{Y*}(y) (deconvolution validity).
    Uses float64 (n<=16 fine; cross-checked against mpmath formula separately).
    """
    n = len(x)
    N = 1 << n
    # enumerate all y as bit arrays
    ys = np.arange(N, dtype=np.int64)
    bits = ((ys[:, None] >> np.arange(n)[::-1]) & 1).astype(np.int8)  # N x n, MSB first
    # P_X(y) for all y (stationary BSMS)
    flips = np.sum(bits[:, 1:] != bits[:, :-1], axis=1)
    PX = 0.5 * (p ** flips) * ((1 - p) ** (n - 1 - flips))
    # Walsh transform: Phat_X(w) = sum_y P_X(y) (-1)^{<w,y>}, w ranges all 2^n
    # We need Phat over all w; do via Hadamard (fast) using sign matrix in chunks.
    # For n<=16, N<=65536: N x N float is 4e9 -> too big. Use FWHT.
    def fwht(a):
        a = a.astype(np.float64).copy()
        h = 1
        while h < len(a):
            for i in range(0, len(a), h * 2):
                xx = a[i:i + h].copy()
                yy = a[i + h:i + 2 * h].copy()
                a[i:i + h] = xx + yy
                a[i + h:i + 2 * h] = xx - yy
            h *= 2
        return a
    # FWHT computes sum_y a[y] (-1)^{<w,y>} with natural (Sylvester) ordering;
    # index ordering matches integer bit dot product. PX indexed by integer y.
    PhatX = fwht(PX)                       # = sum_y P_X(y)(-1)^{<w,y>}, indexed by w int
    wt = np.array([bin(w).count("1") for w in range(N)])
    PhatY = PhatX * ((1 - 2 * D) ** (-wt.astype(np.float64)))
    # inverse Walsh: P_Y(y) = (1/N) sum_w PhatY(w)(-1)^{<w,y>} = (1/N) FWHT(PhatY)
    PY = fwht(PhatY) / N                   # method (ii)
    # method (iii): direct per-y for the words in the ball only (cross-check a few)
    x_int = 0
    for b in x:
        x_int = (x_int << 1) | int(b)
    dist = np.array([bin(int(yy) ^ x_int).count("1") for yy in ys])
    inball = dist <= t
    P_walsh = PY[inball].sum()
    # method (iii) direct: for each ball y, recompute P_Y(y) = sum_w PhatY[w](-1)^{<w,y>}/N
    ball_ys = ys[inball]
    P_direct = 0.0
    allw = np.arange(N)
    wbits = ((allw[:, None] >> np.arange(n)[::-1]) & 1).astype(np.int8)
    for yint in ball_ys:
        ybits = ((yint >> np.arange(n)[::-1]) & 1).astype(np.int8)
        dot = (wbits * ybits).sum(axis=1) & 1
        s = np.sum(PhatY * np.where(dot == 0, 1.0, -1.0)) / N
        P_direct += s
    return P_walsh, P_direct, PY.min()

# --------------------------------------------------------------------------
# VALIDATION
# --------------------------------------------------------------------------
def validation(verbose=True):
    print("=" * 78)
    print("VALIDATION (mandatory): formula vs brute force, n<=16, several x and (p,D)")
    print("=" * 78)
    ok = True
    rng = np.random.default_rng(12345)
    cases = [(0.25, None), (0.4, None), (0.1, None)]
    for p, _ in cases:
        Dc = D_c(p)
        D = 0.9 * Dc
        for n in (8, 10, 12, 14, 16):
            t = int(math.floor(n * D))
            # need t>=0; at small n*D, t may be 0 -> ball is just {x}. still valid test.
            Gt = krawtchouk_partial(n, t)
            # a few random x ~ BSMS plus deterministic ones
            xs = []
            for _ in range(3):
                xs.append(sample_one_bsms(n, p, rng))
            xs.append(np.zeros(n, dtype=np.int8))
            xs.append(np.array([i % 2 for i in range(n)], dtype=np.int8))
            for x in xs:
                with mp.workdps(60):
                    Pf, cancel = ball_prob_formula(x, p, D, t, Gt, 60)
                    Pf = float(Pf)
                P_walsh, P_direct, miny = brute_ball_probs(x, p, D, t)
                e1 = abs(Pf - P_walsh)
                e2 = abs(Pf - P_direct)
                e3 = abs(P_walsh - P_direct)
                good = (e1 < 1e-9 and e2 < 1e-9 and e3 < 1e-9 and miny > -1e-12)
                ok = ok and good
                if not good or (n == 16):
                    print(f"  p={p:.2f} D={D:.4f} n={n:2d} t={t} | "
                          f"Pf={Pf:.3e} Pw={P_walsh:.3e} Pd={P_direct:.3e} | "
                          f"|f-w|={e1:.1e} |f-d|={e2:.1e} |w-d|={e3:.1e} "
                          f"min_yPY={miny:.2e} cancel~{float(cancel):.1e} "
                          f"{'OK' if good else 'FAIL<<<'}")
    # Var(j_n)/n -> V_lossless check (deferred to main run, but a quick small-n note)
    print("-" * 78)
    print(f"  (deconvolution validity min_y P_Y(y) >= 0 confirmed in-Gray above)")
    print(f"VALIDATION: {'PASS' if ok else 'FAIL'}")
    print("=" * 78)
    return ok

def sample_one_bsms(n, p, rng):
    steps = (rng.random(n) < p).astype(np.int8)
    steps[0] = (rng.random() < 0.5)
    return np.cumsum(steps) % 2

# --------------------------------------------------------------------------
# MAIN annealed measurement (per-sample worker, optionally parallel)
# --------------------------------------------------------------------------
def _sample_chunk(args):
    """
    Worker: compute (j_n, G_n, log2S) for `R` samples ~ BSMS(p), returning
    raw arrays + (dropped, worst_cancel).  Precision tracked per sample.
    `Gt` (Krawtchouk partial sums) passed in (built once in parent).
    """
    p, D, n, t, Gt, R, dps, keepdig, child_seed = args
    rng = np.random.default_rng(child_seed)
    nh_D = n * h2(D)
    jn, Gn, lS = [], [], []
    dropped = 0
    worst_cancel = 0.0
    for _ in range(R):
        x = sample_one_bsms(n, p, rng)
        i_n = i_n_bits(x, p)
        j_n = i_n - nh_D
        with mp.workdps(dps):
            P, cancel = ball_prob_formula(x, p, D, t, Gt, dps)
            cancel_f = float(cancel)
            if cancel_f > worst_cancel:
                worst_cancel = cancel_f
            lost = math.log10(cancel_f) if cancel_f > 1 else 0.0
            if (dps - lost) < keepdig or P <= 0:
                dropped += 1
                continue
            G_n = float(-mp.log(P, 2))
        jn.append(j_n); Gn.append(G_n); lS.append(j_n - G_n)
    return (np.array(jn), np.array(Gn), np.array(lS), dropped, worst_cancel)

def measure(p, D, n, R, dps, keepdig, base_seed, procs=1, pool=None):
    """
    Sample R words x^n ~ BSMS(p) (annealed). Compute j_n, G_n, log2S per sample,
    drop samples whose cancellation eats the precision budget. Parallel over
    `procs` worker chunks. Return stats dict (or None if too few kept).
    """
    t = int(math.floor(n * D))
    Gt = krawtchouk_partial(n, t)
    # split R into `procs` chunks, distinct child seeds
    nchunks = max(1, procs)
    base = R // nchunks
    sizes = [base + (1 if i < R - base * nchunks else 0) for i in range(nchunks)]
    ss = np.random.SeedSequence(base_seed + n * 7919 + int(1000 * p))
    children = ss.spawn(nchunks)
    jobs = [(p, D, n, t, Gt, sizes[i], dps, keepdig, children[i])
            for i in range(nchunks) if sizes[i] > 0]
    if pool is not None and len(jobs) > 1:
        outs = pool.map(_sample_chunk, jobs)
    else:
        outs = [_sample_chunk(j) for j in jobs]
    jn = np.concatenate([o[0] for o in outs]) if outs else np.array([])
    Gn = np.concatenate([o[1] for o in outs]) if outs else np.array([])
    lS = np.concatenate([o[2] for o in outs]) if outs else np.array([])
    dropped = sum(o[3] for o in outs)
    worst_cancel = max((o[4] for o in outs), default=0.0)
    kept = len(jn)
    if kept < 30:
        return None
    def vstats(a):
        v = a.var(ddof=1)
        return v, v * math.sqrt(2.0 / (kept - 1))  # SE(var) ~ var*sqrt(2/(m-1)), normal approx
    Vj, seVj = vstats(jn)
    VG, seVG = vstats(Gn)
    VlS, seVlS = vstats(lS)
    cov = np.cov(jn, lS, ddof=1)[0, 1]
    # SE of covariance (delta-method-ish): sd of product deviations / sqrt(m)
    prod = (jn - jn.mean()) * (lS - lS.mean())
    se_cov = np.std(prod, ddof=1) / math.sqrt(kept)
    return dict(n=n, t=t, kept=kept, dropped=dropped, worst_cancel=worst_cancel,
                Vj_n=Vj / n, seVj_n=seVj / n,
                VG_n=VG / n, seVG_n=seVG / n,
                VlS_n=VlS / n, seVlS_n=seVlS / n,
                cov_n=cov / n, secov_n=se_cov / n,
                meanG_n=Gn.mean() / n, meanj_n=jn.mean() / n, meanlS_n=lS.mean() / n)

def run_one_p(p, nR_schedule, dps, keepdig, seed, pool, procs):
    """Run all (n,R) cells for a single p; print the table; return (D,Vloss,rows)."""
    Dc = D_c(p)
    D = 0.9 * Dc
    Vloss = V_lossless(p)
    print(f"\n### p={p}  D=0.9*D_c={D:.6f}  (D_c={Dc:.6f})  "
          f"V_lossless={Vloss:.5f} bits^2  (1-2D={1-2*D:.4f}) ###")
    print(f"{'n':>5} {'t':>3} {'nD':>7} {'kept':>6} {'drop':>5} "
          f"{'Var(j)/n':>14} {'Var(G)/n':>14} {'Var(log2S)/n':>16} {'Cov/n':>10} "
          f"{'wcancel':>9}")
    rows = []
    for n, R in nR_schedule:
        res = measure(p, D, n, R, dps, keepdig, seed, procs=procs, pool=pool)
        if res is None:
            print(f"{n:>5}  -- insufficient kept samples (precision ceiling) --")
            continue
        rows.append(res)
        print(f"{res['n']:>5} {res['t']:>3} {n*D:>7.2f} {res['kept']:>6} "
              f"{res['dropped']:>5} "
              f"{res['Vj_n']:>8.4f}+-{res['seVj_n']:.3f} "
              f"{res['VG_n']:>8.4f}+-{res['seVG_n']:.3f} "
              f"{res['VlS_n']:>9.5f}+-{res['seVlS_n']:.4f} "
              f"{res['cov_n']:>+9.4f} {res['worst_cancel']:>9.1e}", flush=True)
    if rows:
        # identity check Var(G) = Var(j) - 2 Cov(j,log2S) + Var(log2S)
        print("  identity check  Var(G)/n  vs  Var(j)/n - 2Cov/n + Var(log2S)/n :")
        for r in rows:
            rhs = r['Vj_n'] - 2 * r['cov_n'] + r['VlS_n']
            print(f"    n={r['n']:>4}: {r['VG_n']:.5f}  vs  {rhs:.5f}  "
                  f"(diff {abs(r['VG_n']-rhs):.1e})")
    return D, Vloss, rows

def verdict_for_p(p, D, Vloss, rows):
    if not rows:
        print(f"  p={p}: no reliable rows."); return
    # trust rows past the integer-floor artifact: t = floor(nD) >= 4 (nD>=~5)
    good = [r for r in rows if r['t'] >= 4]
    tag = "t>=4 (past integer ball-floor artifact, nD>=~5)"
    if len(good) < 3:
        good = [r for r in rows if r['t'] >= 2]
        tag = "t>=2 (partial floor escape)"
    if len(good) < 3:
        good = rows; tag = "ALL rows (could not escape floor)"
    Vj = [r['Vj_n'] for r in good]
    VlS = [r['VlS_n'] for r in good]
    seVlS = [r['seVlS_n'] for r in good]
    ns = [r['n'] for r in good]
    print(f"\n  p={p}  D={D:.6f}  V_lossless={Vloss:.4f}   [{tag}]")
    print(f"    n             : {ns}")
    print(f"    nD            : {[f'{n*D:.1f}' for n in ns]}")
    print(f"    Var(j)/n      : {[f'{v:.4f}' for v in Vj]}   (target -> {Vloss:.4f})")
    print(f"    Var(log2S)/n  : {[f'{v:.5f}' for v in VlS]}")
    print(f"    SE            : {[f'{s:.5f}' for s in seVlS]}")
    if len(ns) >= 3:
        # fit Var(log2S)/n ~ a/n + b ; if dispersion residual -> 0, b ~ 0 and decreasing
        invn = np.array([1.0 / n for n in ns])
        yv = np.array(VlS)
        wts = 1.0 / np.array(seVlS) ** 2
        A = np.vstack([invn, np.ones_like(invn)]).T
        Aw = A * np.sqrt(wts)[:, None]; yw = yv * np.sqrt(wts)
        coef, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
        a, b = coef
        # crude SE of intercept from weighted LSQ covariance
        cov = np.linalg.inv(Aw.T @ Aw)
        se_b = math.sqrt(cov[1, 1])
        print(f"    weighted fit Var(log2S)/n ~ a/n + b :  a={a:+.4f},  "
              f"b(extrap n->inf)={b:+.5f} +- {se_b:.5f}")
        # also: is the SEQUENCE itself decreasing toward 0?
        last, lastse = VlS[-1], seVlS[-1]
        decreasing = VlS[-1] < VlS[0]
        if b < 2 * se_b and (last < 4 * lastse + 0.01):
            print(f"    => Var(log2S)/n trends to 0 (extrap b~0) "
                  f"=> ACHIEVABILITY-CONSISTENT: V_op = V_lossless certified-numerically.")
        elif b > 3 * se_b and b > 0.02:
            print(f"    => Var(log2S)/n trends to POSITIVE b={b:.4f}+-{se_b:.4f} "
                  f"=> CONVERSE LOOSE: V_op = V_lossless + c (c>0).")
        else:
            print(f"    => INCONCLUSIVE in reliable range "
                  f"(extrap b={b:+.4f}+-{se_b:.4f}; "
                  f"last={last:.4f}+-{lastse:.4f}; "
                  f"{'decreasing' if decreasing else 'NOT monotone-decreasing'}).")

# --------------------------------------------------------------------------
def build_schedule(p, nmax, R_small, R_large):
    """
    Per-p n-grid + sample budget. D_c(0.25)~0.029 vs D_c(0.4)~0.127, so to reach
    nD ~ 5-10 we need much larger n at p=0.25. Bigger R at small n (cheap), smaller
    R at large n (expensive). Returns list of (n, R).
    """
    Dc = D_c(p)
    D = 0.9 * Dc
    # candidate n's, geometric-ish, capped at nmax
    cands = [16, 24, 32, 48, 64, 96, 128, 160, 200, 240, 280, 320, 400, 480, 560]
    sched = []
    for n in cands:
        if n > nmax:
            break
        # cost ~ n^2 per sample; scale R down as n grows
        R = max(R_small if n <= 64 else int(R_large * (96.0 / n) ** 1.3), 400)
        R = min(R, R_small)
        sched.append((n, R))
    return sched

def main():
    import multiprocessing as mpc
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=6000, help="samples at small n")
    ap.add_argument("--Rlarge", type=int, default=6000, help="base for large-n schedule")
    ap.add_argument("--dps", type=int, default=60)
    ap.add_argument("--keepdig", type=int, default=15)
    ap.add_argument("--seed", type=int, default=20260606)
    ap.add_argument("--nmax", type=int, default=320)
    ap.add_argument("--procs", type=int, default=max(1, (mpc.cpu_count() or 2) - 1))
    ap.add_argument("--skip-validation", action="store_true")
    ap.add_argument("--quick", action="store_true", help="small R, few n, for smoke")
    args = ap.parse_args()

    t0 = time.time()
    if not args.skip_validation:
        ok = validation()
        if not ok:
            print("\n*** VALIDATION FAILED — aborting before scale-up. ***")
            return

    print("=" * 118)
    print(f"MAIN: annealed Var/n trends. mpmath dps={args.dps}, keep if retained "
          f"digits >= {args.keepdig}. procs={args.procs} seed={args.seed}")
    print("=" * 118)

    p_list = [0.25, 0.4]
    pool = None
    if args.procs > 1:
        ctx = mpc.get_context("fork")
        pool = ctx.Pool(args.procs)
    try:
        all_rows = {}
        for p in p_list:
            if args.quick:
                sched = [(16, 400), (24, 400), (32, 300)]
            else:
                sched = build_schedule(p, args.nmax, args.R, args.Rlarge)
            D, Vloss, rows = run_one_p(p, sched, args.dps, args.keepdig,
                                       args.seed, pool, args.procs)
            all_rows[p] = (D, Vloss, rows)
    finally:
        if pool is not None:
            pool.close(); pool.join()

    print("\n" + "=" * 118)
    print("VERDICT")
    print("=" * 118)
    for p, (D, Vloss, rows) in all_rows.items():
        verdict_for_p(p, D, Vloss, rows)
    print("\n  HONEST CAVEATS:")
    print("   - Integer ball-floor t=floor(nD): at small nD, t is frozen (0,1,2) and")
    print("     Var(log2S)/n SAWTOOTHS at each radius jump; trust only nD>=~5 (t>=4) rows.")
    print("   - Catastrophic cancellation in (1-2D)^{-j}Gt(j)H_j(x): mild in-Gray")
    print("     (worst_cancel column); samples below the digit budget are DROPPED.")
    print("   - Finite-sample SE on each variance ~ var*sqrt(2/(m-1)); shown as +-.")
    print(f"   - Elapsed {time.time()-t0:.0f}s.")

if __name__ == "__main__":
    main()
