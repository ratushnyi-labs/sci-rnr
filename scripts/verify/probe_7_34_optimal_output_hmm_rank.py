#!/usr/bin/env python3
"""
PROBE the 7.34 achievability crux: is the RD-OPTIMAL OUTPUT process for the BSMS a FINITE-STATE
HMM?  If yes, the optimal n-letter d-tilted info j_n is a log-transfer-matrix product (asymptotically
additive over finitely many hidden states) => Markov Berry-Esseen (Zhou-Tan-Motani template) =>
V_op = V_conv < V_lossless, refuting the candidate.  If the optimal output needs UNBOUNDED hidden
state, the residual stands.

Test: the HANKEL matrix H[u,v] = P(prefix=u, suffix=v) of a stationary process has rank = the number
of predictive (causal) states; a finite-state HMM has FINITE Hankel rank.  We form the block-measure
Hankel by splitting the n-letter law into two contiguous halves (u = first n/2 symbols, v = last n/2)
and take its numerical (SVD) rank.

Baseline/validation: the SOURCE (order-1 Markov) must give Hankel rank 2.  We then compare the
OPTIMAL OUTPUT q* (from the exact n-letter Blahut-Arimoto at distortion D) -- small stable rank =>
finite-state HMM (=> V_op=V_conv likely); rank growing with n => not finite-state.
"""
import math
import numpy as np

def popcount_table(N):
    return np.array([bin(i).count("1") for i in range(N)], dtype=np.int64)

def bsms_P(n, p):
    N = 1 << n
    x = np.arange(N, dtype=np.int64)
    shifted = (x ^ (x >> 1)) & ((1 << n) - 1)
    mask = (1 << (n-1)) - 1
    sw = np.array([bin(int(v & mask)).count("1") for v in shifted], dtype=np.int64)
    stays = (n-1) - sw
    logP = math.log2(0.5) + stays*math.log2(1-p) + sw*math.log2(p)
    P = np.exp(logP*math.log(2)); P /= P.sum()
    return P

def blahut_arimoto(P, Dm, s, iters=500, tol=1e-13):
    N = Dm.shape[0]
    W = np.exp(-s * Dm.astype(np.float64))
    q = np.full(N, 1.0/N)
    for _ in range(iters):
        Z = np.maximum(W @ q, 1e-300)
        qn = q * (W.T @ (P / Z)); qn /= qn.sum()
        if np.max(np.abs(qn - q)) < tol:
            q = qn; break
        q = qn
    return q

def hankel_rank(prob_vec, n, rel_tol=1e-9):
    """Split n-letter law into first n//2 / last n//2 symbols; SVD rank of [P(u,v)]."""
    h = n // 2
    M = prob_vec.reshape(1 << (n-h), 1 << h)   # rows=high (n-h) bits, cols=low h bits
    sv = np.linalg.svd(M, compute_uv=False)
    if sv[0] <= 0: return 0, sv
    rank = int(np.sum(sv > rel_tol * sv[0]))
    return rank, sv

if __name__ == "__main__":
    p = 0.1; s = 2.0   # slope -> D ~ 0.05 (Gray region)
    print("="*78)
    print(f"BSMS p={p}, BA slope s={s} (D~0.05). Hankel rank = # predictive states (finite => HMM).")
    print(f"{'n':>3} {'src rank':>9} {'src top-sv':>26} {'q* rank':>8} {'q* top singular values':>30}")
    for n in (8, 10, 12):
        N = 1 << n
        P = bsms_P(n, p)
        pc = popcount_table(N); x = np.arange(N, dtype=np.int64)
        Dm = pc[(x[:, None] ^ x[None, :])]
        q = blahut_arimoto(P, Dm, s)
        rs, svs = hankel_rank(P, n)
        rq, svq = hankel_rank(q, n)
        srcsv = " ".join(f"{v:.2e}" for v in svs[:4])
        qsv = " ".join(f"{v:.2e}" for v in svq[:6])
        print(f"{n:>3} {rs:>9} {srcsv:>26} {rq:>8}  {qsv}")
    print("="*78)
    print("INTERPRETATION:")
    print("  - source Hankel rank should be 2 (order-1 Markov) -- validates the probe.")
    print("  - q* rank small & STABLE across n (e.g. 2-8) => optimal output is a FINITE-STATE HMM")
    print("    => j_n asymptotically additive (log-transfer-matrix) => V_op=V_conv (candidate refuted).")
    print("  - q* rank GROWING with n / many comparable singular values => not finite-state => residual.")
    print("  NOTE: finite-block BA has boundary effects; this is SUGGESTIVE evidence, not a proof.")
