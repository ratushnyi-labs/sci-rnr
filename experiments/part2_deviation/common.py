"""Shared machinery for the Part II deviation-hierarchy simulation studies.

Sources: i.i.d. over a finite alphabet and finite-state Markov chains, with
closed-form entropy rate h, varentropy rate V (fundamental-matrix CLT variance),
scaled cumulant generating function Lambda(t) (base 2), Legendre rate function
E(R) = sup_t [tR - Lambda(t)], max-cycle-mean endpoint, and exponentially tilted
importance sampling for deep-tail overflow probabilities (unbiased; the tilted
kernel for Markov chains uses the Perron right eigenvector of the tilted transfer
matrix).

All laws under test are those of Part II (random access + deviation hierarchy),
Theorems 7.31-7.33, 7.35, 7.36; the exact statements are mirrored from
scripts/verify/thm_7_3*.py.
"""
import csv
import math
import os

import numpy as np

LN2 = math.log(2.0)
EULER = 0.5772156649015329  # Euler-Mascheroni

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# Fixed-order categorical palette (Okabe-Ito; colorblind-safe).
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]


# ----------------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------------
def report(tag, ok, msg, results=None):
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}: {msg}"
    print(line, flush=True)
    if results is not None:
        results.append((tag, bool(ok), msg))
    return bool(ok)


def save_csv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    print(f"  [csv] {path}", flush=True)


def style_axes(ax):
    """Recessive grid/axes for diagnostic figures."""
    ax.grid(True, alpha=0.25, linewidth=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save_fig(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(fig)
    print(f"  [png] {path}", flush=True)


# ----------------------------------------------------------------------------
# i.i.d. source
# ----------------------------------------------------------------------------
class IIDSource:
    """i.i.d. source with per-symbol code cost l(x) = -log2 M(x) (model M,
    default matched M = P). h below is the cross-entropy rate (= entropy rate
    when matched); V the variance of the per-symbol cost."""

    def __init__(self, p, model=None):
        self.p = np.asarray(p, dtype=float)
        self.A = len(self.p)
        self.m = self.p if model is None else np.asarray(model, dtype=float)
        self.l = -np.log2(self.m)
        self.h = float(self.p @ self.l)
        self.V = float(self.p @ (self.l - self.h) ** 2)

    def sample_counts(self, N, R, rng):
        return rng.multinomial(N, self.p, size=R)

    def lengths_from_counts(self, counts):
        return counts @ self.l

    def sample_lengths(self, N, R, rng):
        return self.lengths_from_counts(self.sample_counts(N, R, rng))

    def Lambda(self, t):
        """Base-2 scaled CGF: Lambda(t) = log2 sum_x P(x) M(x)^{-t}."""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        val = np.log2((self.p[None, :] * self.m[None, :] ** (-t[:, None])).sum(axis=1))
        return val if val.size > 1 else float(val[0])

    def legendre_E(self, R, t_hi=60.0, n_grid=120001):
        """E(R) = sup_{t>=0} [tR - Lambda(t)]  (base 2)."""
        ts = np.linspace(0.0, t_hi, n_grid)
        return float(np.max(ts * R - self.Lambda(ts)))

    def tilt_root(self, R, t_hi=60.0, n_grid=120001):
        """argmax_t [tR - Lambda(t)]  (= t* solving Lambda'(t*)=R on the interior)."""
        ts = np.linspace(0.0, t_hi, n_grid)
        return float(ts[np.argmax(ts * R - self.Lambda(ts))])

    def is_tail_log2(self, N, thresh, R_samp, rng, t=None):
        """Importance-sampled log2 P(sum_i l_i >= thresh) via the exponential
        tilt q_t(x) = p(x) m(x)^{-t} / Z_t.  Unbiased; returns (log2 p_hat, hits,
        rel_se) where rel_se is the relative standard error of p_hat."""
        if t is None:
            t = self.tilt_root(thresh / N)
        w = self.p * self.m ** (-t)
        Z = w.sum()
        q = w / Z
        counts = rng.multinomial(N, q, size=R_samp)
        L = counts @ self.l
        log2w = N * np.log2(Z) - t * L  # per-sample log2 importance weight
        sel = L >= thresh - 1e-12
        if not sel.any():
            return -np.inf, 0, np.inf
        mx = float(log2w[sel].max())
        vals = np.zeros(R_samp)
        vals[sel] = np.exp((log2w[sel] - mx) * LN2)
        mean = vals.mean()
        se = vals.std(ddof=1) / math.sqrt(R_samp)
        return mx + math.log2(mean), int(sel.sum()), float(se / mean)


# ----------------------------------------------------------------------------
# finite-state Markov source
# ----------------------------------------------------------------------------
def markov_stationary(P):
    P = np.asarray(P, dtype=float)
    n = len(P)
    pi = np.full(n, 1.0 / n)
    for _ in range(200000):
        nx = pi @ P
        if np.max(np.abs(nx - pi)) < 1e-15:
            return nx
        pi = nx
    return pi


class MarkovSource:
    """Finite-state irreducible Markov chain P with per-step code cost
    g(x,x') = -log2 M(x'|x) on edges P(x,x')>0 (model M, default matched).
    h = long-run mean cost (cross-entropy rate); V = long-run variance of the
    additive cost functional (fundamental-matrix Markov CLT variance, the same
    varentropy constant as Theorems 7.28/7.29/7.35/7.36)."""

    def __init__(self, P, model=None):
        self.P = np.asarray(P, dtype=float)
        self.n = len(self.P)
        self.M = self.P if model is None else np.asarray(model, dtype=float)
        self.pi = markov_stationary(self.P)
        with np.errstate(divide="ignore"):
            g = np.where(self.P > 0, -np.log2(np.where(self.M > 0, self.M, 1.0)), 0.0)
        if np.any((self.P > 0) & (self.M <= 0)):
            raise ValueError("model must be positive on every source edge")
        self.g = g
        self.h, self.V = self._h_V()
        self.cdf = np.cumsum(self.P, axis=1)

    def _h_V(self):
        P, pi, g, n = self.P, self.pi, self.g, self.n
        f = (P * g).sum(axis=1)  # conditional mean cost per state
        h = float(pi @ f)
        Z = np.linalg.inv(np.eye(n) - P + np.outer(np.ones(n), pi))
        d = np.where(P > 0, g - h, 0.0)
        fc = f - h
        term = float(pi @ (P * d ** 2).sum(axis=1))
        u = Z @ fc
        cross = 2.0 * float(pi @ (P * (d * u[None, :])).sum(axis=1))
        return h, term + cross

    # -- path simulation -----------------------------------------------------
    def sample_lengths(self, N, R, rng, kernel=None, init=None, return_paths=False):
        """R paths of N steps; returns total cost per path (under self.g).
        kernel/init override the transition law (used by importance sampling).
        If return_paths, also returns (x0, xN)."""
        K = self.P if kernel is None else np.asarray(kernel, dtype=float)
        cdf = np.cumsum(K, axis=1)
        p0 = self.pi if init is None else np.asarray(init, dtype=float)
        cur = rng.choice(self.n, size=R, p=p0)
        x0 = cur.copy()
        L = np.zeros(R)
        for _ in range(N):
            u = rng.random(R)
            nxt = np.minimum((u[:, None] > cdf[cur]).sum(axis=1), self.n - 1)
            L += self.g[cur, nxt]
            cur = nxt
        if return_paths:
            return L, x0, cur
        return L

    def sample_cost_paths(self, N, R, rng):
        """R x N array of per-step costs (for partial-sum / buffer studies)."""
        cur = rng.choice(self.n, size=R, p=self.pi)
        out = np.empty((R, N))
        for t in range(N):
            u = rng.random(R)
            nxt = np.minimum((u[:, None] > self.cdf[cur]).sum(axis=1), self.n - 1)
            out[:, t] = self.g[cur, nxt]
            cur = nxt
        return out

    # -- tilted transfer matrix / CGF ----------------------------------------
    def tilted_matrix(self, t):
        return np.where(self.P > 0, self.P * 2.0 ** (t * self.g), 0.0)

    def perron(self, t):
        """(rho, r): Perron root and positive right eigenvector of A_t."""
        A = self.tilted_matrix(t)
        w, vr = np.linalg.eig(A)
        i = int(np.argmax(w.real))
        rho = float(w[i].real)
        r = vr[:, i].real
        r = r * np.sign(r[np.argmax(np.abs(r))])
        r = np.maximum(r, 1e-300)
        return rho, r

    def Lambda(self, t):
        """Base-2 scaled CGF: log2 rho(A_t), A_t(x,x') = P(x,x') M(x'|x)^{-t}."""
        if np.ndim(t) == 0:
            return math.log2(self.perron(float(t))[0])
        return np.array([math.log2(self.perron(float(tt))[0]) for tt in np.asarray(t)])

    def legendre_E(self, R, t_hi=40.0, n_grid=4001):
        ts = np.linspace(0.0, t_hi, n_grid)
        return float(np.max(ts * R - self.Lambda(ts)))

    def tilt_root(self, R, t_hi=40.0, n_grid=4001):
        ts = np.linspace(0.0, t_hi, n_grid)
        return float(ts[np.argmax(ts * R - self.Lambda(ts))])

    def lambda_prime(self, t, dt=1e-4):
        return (self.Lambda(t + dt) - self.Lambda(t - dt)) / (2 * dt)

    # -- endpoint machinery ---------------------------------------------------
    def max_cycle_mean(self):
        """Max over simple cycles (edges with P>0) of the mean cost -- the LDP
        upper endpoint r_+ = lim_{t->oo} Lambda'(t).  Exhaustive for small n."""
        from itertools import permutations

        n = len(self.P)
        best = -np.inf
        for k in range(1, n + 1):
            for cyc in permutations(range(n), k):
                if cyc[0] != min(cyc):  # canonical rotation only
                    continue
                edges = list(zip(cyc, cyc[1:] + cyc[:1]))
                if all(self.P[i, j] > 0 for i, j in edges):
                    mean = sum(self.g[i, j] for i, j in edges) / k
                    best = max(best, mean)
        return best

    def max_path_total(self, N):
        """DP: maximum total cost of any positive-probability N-step path."""
        dp = np.where(self.pi > 0, 0.0, -np.inf)
        for _ in range(N):
            nxt = np.full(self.n, -np.inf)
            for j in range(self.n):
                cand = dp + np.where(self.P[:, j] > 0, self.g[:, j], -np.inf)
                nxt[j] = cand.max()
            dp = nxt
        return float(dp.max())

    def edge_max(self):
        """Max single-step cost over positive-probability edges (ell_max)."""
        return float(np.max(np.where(self.P > 0, self.g, -np.inf)))

    # -- tilted-kernel importance sampling ------------------------------------
    def is_tail_log2(self, N, thresh, R_samp, rng, t=None):
        """Importance-sampled log2 P(L_N >= thresh) via the tilted kernel
        Q_t(x'|x) = A_t(x,x') r(x') / (rho r(x)).  Per-path log2 weight:
        N log2 rho - t L + log2 r(x0) - log2 r(xN)  (initial dist kept = pi).
        Returns (log2 p_hat, hits, rel_se)."""
        if t is None:
            t = self.tilt_root(thresh / N)
        rho, r = self.perron(t)
        Q = self.tilted_matrix(t) * r[None, :] / (rho * r[:, None])
        Q = Q / Q.sum(axis=1, keepdims=True)  # numerical renorm
        L, x0, xN = self.sample_lengths(N, R_samp, rng, kernel=Q, return_paths=True)
        log2w = N * math.log2(rho) - t * L + np.log2(r[x0]) - np.log2(r[xN])
        sel = L >= thresh - 1e-12
        if not sel.any():
            return -np.inf, 0, np.inf
        mx = float(log2w[sel].max())
        vals = np.zeros(R_samp)
        vals[sel] = np.exp((log2w[sel] - mx) * LN2)
        mean = vals.mean()
        se = vals.std(ddof=1) / math.sqrt(R_samp)
        return mx + math.log2(mean), int(sel.sum()), float(se / mean)


# ----------------------------------------------------------------------------
# small statistical helpers
# ----------------------------------------------------------------------------
def gumbel_norm(m):
    """(a_m, b_m) Gumbel norming constants: a_m = sqrt(2 ln m),
    b_m = a_m - (ln ln m + ln 4 pi)/(2 a_m)."""
    a = math.sqrt(2 * math.log(m))
    b = a - (math.log(math.log(m)) + math.log(4 * math.pi)) / (2 * a)
    return a, b


def fit_affine_with_offset(N_arr, y_arr, offset_fn=None):
    """Least-squares fit y = slope*N + c (+ offset_fn(N) imposed).  Returns slope."""
    N_arr = np.asarray(N_arr, dtype=float)
    y = np.asarray(y_arr, dtype=float)
    if offset_fn is not None:
        y = y - np.array([offset_fn(N) for N in N_arr])
    X = np.vstack([N_arr, np.ones_like(N_arr)]).T
    slope, _c = np.linalg.lstsq(X, y, rcond=None)[0]
    return float(slope)
