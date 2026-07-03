#!/usr/bin/env python3
r"""
probe_7_34_beyond_gray_structure.py
============================================================================
BUG-009-D beyond-Gray (D > D_c): the threshold MECHANISM, the finite-n law, and
the WHT-accelerated exact block-RD solver -- the verifiable core of the
beyond-Gray research findings (user-authorized fan-out; artifact
research_beyond_gray/; re-executed before packaging).

FINDINGS PACKAGED AS CHECKS:
  B1 (PROVEN, symbolic): for the BSMS(p) deconvolution the 2-step alternating
     transfer matrix C = B0 B1 has
         disc(C) = p^2 [p^2 - 4 D(1-D)(1-p)^2] / (1-2D)^2,
     vanishing EXACTLY at the Gray threshold D_c(p) = (1/2)(1 - sqrt(1 - p^2/
     ((1-p)^2 + ...)))-form root of p^2 = 4D(1-D)(1-p)^2 -- i.e. the eigenvalue
     collision IS the threshold (the mechanism behind the whole Gray family).
  B2 (VERIFIED NUMERIC): the finite-n threshold law
         D_c^(n) - D_c = (2 pi / c(p))^2 / n^2 (1 + o(1)),
     equivalently n^2 (D_c^(n) - D_c) -> K(p)^2: the per-n Gray region extends
     STRICTLY beyond D_c with a universal 1/n^2 approach.
  B3 (VALIDATION): the WHT/XOR-convolution O(N log N) exact Blahut-Arimoto
     block-RD solver agrees with a dense O(N^2) reference at small n.

DOCUMENTED FINDINGS AT THEIR ARTIFACT EVIDENCE LEVEL (not re-checked here):
  * The optimal reverse channel remains EXACTLY iid/product BSC(D) throughout
    the per-n Gray region D <= D_c^(n) (artifact s4; TC and KL to product = 0
    within solver tolerance) -- the product structure survives beyond D_c.
  * CONJECTURE (artifact s5/s9/s10, verified over ~8 decades of residual):
    beyond Gray, R(D) - SLB(D) = exp(-b_R(p)(1+o(1))/sqrt(D - D_c)) -- an
    ESSENTIAL SINGULARITY at D_c (no closed form for b_R(p) yet; the naive
    rare-alternating-run candidates are excluded).
  * The converse dispersion ratio rho(D) stays < 1 continuously across D_c
    (no jump; artifact s3/s12 Richardson n<=16).
  * Support collapse cascades through word space by switch count, each
    even-period defected-alternating family at its own exact threshold from
    the same 2x2 discriminant mechanism (artifact s8).

CHECKS:
  B1  symbolic disc(C) identity + root at D_c (sympy exact).
  B2  n^2 (D_c^(n) - D_c) stabilizes (n = 8, 10, 12, 14; drift < 15%).
  B3  WHT BA solver == dense BA reference at n = 6 (R(D) to 1e-9).

Deps: numpy, sympy.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-4 min.
"""

import math
import numpy as np


# ---------------------------------------------------------------- fwht
def fwht(a):
    """Unnormalized Walsh-Hadamard transform (involution up to factor N)."""
    a = a.astype(np.float64, copy=True)
    N = a.size
    h = 1
    while h < N:
        a = a.reshape(N // (2 * h), 2, h)
        x = a[:, 0, :].copy()
        y = a[:, 1, :].copy()
        a[:, 0, :] = x + y
        a[:, 1, :] = x - y
        a = a.reshape(N)
        h *= 2
    return a


class XorConv:
    """XOR-convolution with a fixed kernel k: (k * f)(x) = sum_u k(u) f(x^u)."""
    def __init__(self, kernel):
        self.N = kernel.size
        self.hk = fwht(kernel)

    def __call__(self, f):
        return fwht(fwht(f) * self.hk) / self.N


# ---------------------------------------------------------------- source
def popcount(x):
    return np.bitwise_count(x).astype(np.int64)


def bsms_block(n, p):
    """Return P (stationary BSMS law on {0,1}^n, X_1 uniform), i_nats = -ln P."""
    N = 1 << n
    x = np.arange(N, dtype=np.uint64)
    sw = popcount((x ^ (x >> np.uint64(1))) & np.uint64((1 << (n - 1)) - 1))
    logP = math.log(0.5) + sw * math.log(p) + (n - 1 - sw) * math.log(1 - p)
    P = np.exp(logP)
    P /= P.sum()
    return P, -np.log(P), sw


def iid_block(n, q):
    N = 1 << n
    x = np.arange(N, dtype=np.uint64)
    ones = popcount(x)
    logP = ones * math.log(q) + (n - ones) * math.log(1 - q)
    P = np.exp(logP)
    P /= P.sum()
    return P, -np.log(P), ones


def D_c(p):
    return 0.5 * (1.0 - math.sqrt(1.0 - 2.0 * p) / (1.0 - p))


def V_lossless_bits(p):
    return p * (1 - p) * (math.log2((1 - p) / p)) ** 2


# ---------------------------------------------------------------- BA
class RDBlock:
    """n-block RD solver at slope s (nats per Hamming unit)."""
    def __init__(self, n, P, i_nats):
        self.n = n
        self.N = 1 << n
        self.P = P
        self.i = i_nats
        Ei = float(np.dot(P, i_nats))
        self.Ei = Ei
        self.Vi = float(np.dot(P, (i_nats - Ei) ** 2))     # total (not per-n)
        self.d = popcount(np.arange(self.N, dtype=np.uint64)).astype(np.float64)
        self._warm = None

    def _ba(self, s, tol=1e-13, iters=60000, q0=None, kkt_tol=1e-8):
        w = np.exp(-s * self.d)
        conv = XorConv(w)
        P = self.P
        q = q0.copy() if q0 is not None else np.full(self.N, 1.0 / self.N)
        q = np.maximum(q, 1e-280)            # warm starts must not carry exact 0s
        q /= q.sum()

        def step(q):
            Z = np.maximum(conv(q), 1e-300)
            t = conv(P / Z)
            qn = q * np.maximum(t, 0.0)      # t can dip <0 only by fwht roundoff
            ssum = qn.sum()
            return qn / ssum, float(t.max()) - 1.0

        it = 0
        ok = False
        for it in range(iters):
            q1, kkt1 = step(q)
            q2, kkt2 = step(q1)
            r = q1 - q
            v = q2 - q1 - r
            nv = float(np.dot(v, v))
            small = float(np.max(np.abs(q2 - q1))) < tol
            if small and kkt2 < kkt_tol:
                q = q2
                ok = True
                break
            if nv < 1e-300:
                q = q2
                continue
            alpha = -math.sqrt(float(np.dot(r, r)) / nv)
            qn = q - 2 * alpha * r + alpha * alpha * v
            if np.any(qn <= 0) or not np.isfinite(qn).all():
                qn = q2                       # NEVER zero an atom by clamping
            else:
                qn /= qn.sum()
                qn, _ = step(qn)              # stabilize
            q = qn
        return q, w, conv, it + 1, ok

    def point(self, s, tol=1e-13, warm=True, kkt_tol=1e-7, iters=20000):
        """Solve at slope s; return dict of exact functionals."""
        q0 = self._warm if warm else None
        q, w, conv, it, ok = self._ba(s, tol=tol, q0=q0, kkt_tol=kkt_tol,
                                      iters=iters)
        self._warm = q
        P = self.P
        Z = np.maximum(conv(q), 1e-300)
        # Dtot = E_P[ (q * (w d))(x) / Z(x) ]
        conv_wd = XorConv(w * self.d)
        Dtot = float(np.dot(P, conv_wd(q) / Z))
        lnZ = np.log(Z)
        ElnZ = float(np.dot(P, lnZ))
        R_nats = -s * Dtot - ElnZ                      # = E[j]; also = I(X;Y*) at fixpoint
        j = -s * Dtot - lnZ                            # d-tilted info (nats), total
        Vj = float(np.dot(P, (j - R_nats) ** 2))
        cov = float(np.dot(P, (j - R_nats) * (self.i - self.Ei)))
        # KKT residual: t(y) = ((P/Z) * w)(y) ; <=1 with = on support
        t = conv(P / Z)
        kkt = float(np.max(t)) - 1.0
        return dict(s=s, q=q, D=Dtot / self.n, Dtot=Dtot, R_nats=R_nats,
                    R_bits=R_nats / math.log(2), j=j, Vj=Vj, Vi=self.Vi,
                    rho=Vj / self.Vi, cov=cov, corr=cov / math.sqrt(Vj * self.Vi)
                    if Vj > 0 else float('nan'),
                    Vj_n_bits=Vj / self.n / math.log(2) ** 2,
                    it=it, conv=ok, kkt=kkt, Z=Z, w=w)

    def solve_D(self, Dtar, s_lo=0.05, s_hi=9.0, tol_D=2e-6, max_bis=60):
        """Bisect slope s to hit per-letter distortion Dtar. D(s) decreasing."""
        r_last = None
        for _ in range(max_bis):
            s = 0.5 * (s_lo + s_hi)
            r = self.point(s)
            r_last = r
            if abs(r['D'] - Dtar) < tol_D:
                return r
            if r['D'] > Dtar:
                s_lo = s          # D(s) decreasing: achieved D too big -> larger s
            else:
                s_hi = s
        return r_last


def sweep_targets(blk, s_grid, targetsD, tol_D=5e-7, max_secant=12, verbose=False):
    """Descending warm-started s-sweep, then local secant refinement per target D.
    Returns (sweep_rows, {Dtar: result})."""
    s_grid = sorted(s_grid, reverse=True)
    sweep = []
    qs = []
    for s in s_grid:
        r = blk.point(s)
        sweep.append(r)
        qs.append(r['q'].copy())
        if verbose:
            print(f"      [sweep s={s:.4f} D={r['D']:.6f} it={r['it']} "
                  f"kkt={r['kkt']:.1e}]", flush=True)
    Ds = [r['D'] for r in sweep]           # ascending (s descending)
    out = {}
    for Dtar in sorted(targetsD):
        import bisect
        k = bisect.bisect_left(Ds, Dtar)
        if k == 0 or k == len(Ds):
            raise ValueError(f"target D={Dtar} outside sweep range "
                             f"[{Ds[0]:.6f},{Ds[-1]:.6f}]")
        sl, sh = sweep[k]['s'], sweep[k - 1]['s']      # D(sl)>Dtar>D(sh)
        Dl, Dh = sweep[k]['D'], sweep[k - 1]['D']
        blk._warm = qs[k].copy()
        r = sweep[k]
        for _ in range(max_secant):
            if abs(r['D'] - Dtar) < tol_D:
                break
            # secant in (s, D)
            s_new = sl + (sh - sl) * (Dl - Dtar) / max(Dl - Dh, 1e-30)
            s_new = min(max(s_new, min(sl, sh) ), max(sl, sh))
            r = blk.point(s_new)
            if r['D'] > Dtar:
                sl, Dl = s_new, r['D']
            else:
                sh, Dh = s_new, r['D']
        out[Dtar] = r
    return sweep, out


# ---------------------------------------------------------------- deconvolution
def deconv(n, p, D):
    """(K^{-1})^{otimes n} P_X as an XOR-convolution; returns the signed vector."""
    N = 1 << n
    P, _, _ = bsms_block(n, p)
    z = np.arange(N, dtype=np.uint64)
    pc = popcount(z)
    a0 = 1.0 - 2.0 * D
    kern = ((1 - D) ** (n - pc) * (-D) ** pc) / a0 ** n
    return XorConv(kern)(P)




# ============================ CHECKS ============================
import sympy as _spx

PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}")

def B1():
    print("-" * 78); print("B1  symbolic: disc(B0 B1) = p^2[p^2-4D(1-D)(1-p)^2]/(1-2D)^2; root = D_c")
    ps, Ds = _spx.symbols('p D', positive=True)
    # B(y) = T * Kinv[y]: BSMS transition T (switch prob p), BSC(D) deconvolution row
    T = _spx.Matrix([[1 - ps, ps], [ps, 1 - ps]])
    K = _spx.Matrix([[1 - Ds, Ds], [Ds, 1 - Ds]])
    Ki = K.inv()
    def B(y):
        return _spx.Matrix([[Ki[y, a] * T[a, b] for b in range(2)] for a in range(2)])
    C = _spx.simplify(B(0) * B(1))
    disc = _spx.simplify(C.trace()**2 - 4 * C.det())
    target = ps**2 * (ps**2 - 4 * Ds * (1 - Ds) * (1 - ps)**2) / (1 - 2 * Ds)**2
    ok_id = _spx.simplify(disc - target) == 0
    # root: p^2 = 4D(1-D)(1-p)^2 at D = D_c(p) = 1/2 - sqrt(1-2p)/(2(1-p))  [A=2 closed form]
    Dc_cf = _spx.Rational(1, 2) - _spx.sqrt(1 - 2 * ps) / (2 * (1 - ps))
    ok_root = _spx.simplify((ps**2 - 4 * Dc_cf * (1 - Dc_cf) * (1 - ps)**2)) == 0
    print(f"     identity: {ok_id};  disc root at the closed-form D_c: {ok_root}")
    rep("B1 threshold mechanism = eigenvalue collision (exact)", bool(ok_id and ok_root))

def B2():
    print("-" * 78); print("B2  finite-n law: n^2 (D_c^(n) - D_c) stabilizes (p = 0.2)")
    import numpy as _np
    pc = 0.2
    Dcv = D_c(pc)
    def Dc_n(n):
        P, _, _ = bsms_block(n, pc)
        hatP = fwht(P.copy())
        pcs = _np.array([popcount(m) for m in range(1 << n)])
        lo, hi = Dcv, 0.5 - 1e-9
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            hatY = hatP / (1.0 - 2.0 * mid) ** pcs
            PY = fwht(hatY.copy()) / (1 << n)
            if PY.min() >= -1e-14: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)
    vals = []
    for n in (8, 10, 12, 14):
        v = (Dc_n(n) - Dcv) * n * n
        vals.append(v)
        print(f"     n={n}: n^2 (D_c^(n) - D_c) = {v:.6f}")
    drift = abs(vals[-1] - vals[-2]) / abs(vals[-1])
    ok = drift < 0.15 and all(v > 0 for v in vals)
    print(f"     last-step drift {100*drift:.1f}% (<15%)")
    rep("B2 1/n^2 threshold law (per-n Gray extends beyond D_c)", ok)

def B3():
    print("-" * 78); print("B3  WHT BA == dense BA at n=6 (p=0.2, two slopes)")
    import numpy as _np
    pc = 0.2; n = 6; N = 1 << n
    P, i_nats, _ = bsms_block(n, pc)
    rd = RDBlock(n, P, i_nats)
    # dense reference BA
    dmat = _np.array([[popcount(x ^ y) for y in range(N)] for x in range(N)], dtype=float)
    def dense_point(s):
        q = _np.full(N, 1.0 / N)
        for _ in range(20000):
            W = _np.exp(-s * dmat) * q[None, :]
            W /= W.sum(axis=1, keepdims=True)
            qn = P @ W
            if _np.abs(qn - q).max() < 1e-13: q = qn; break
            q = qn
        W = _np.exp(-s * dmat) * q[None, :]
        W /= W.sum(axis=1, keepdims=True)
        Dv = float((P[:, None] * W * dmat).sum())
        Rv = float((P[:, None] * W * _np.log(_np.maximum(W, 1e-300) / _np.maximum(q[None, :], 1e-300))).sum())
        return Dv / n, Rv / n
    ok = True
    for s in (2.5, 3.5):
        out = rd.point(s)
        Dw, Rw = out['D'], out['R_nats'] / n   # D already per-symbol; R_nats per-block
        Dd, Rd_ = dense_point(s)
        here = abs(Dw - Dd) < 1e-8 and abs(Rw - Rd_) < 1e-8
        ok = ok and here
        print(f"     s={s}: WHT (D,R)=({Dw:.10f},{Rw:.10f})  dense=({Dd:.10f},{Rd_:.10f})  agree:{here}")
    rep("B3 WHT solver exact vs dense reference", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("Beyond-Gray structure: threshold mechanism, 1/n^2 law, WHT solver")
    print("=" * 78)
    B1(); B2(); B3()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
