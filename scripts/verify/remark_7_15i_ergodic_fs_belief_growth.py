#!/usr/bin/env python3
"""
remark_7_15i_ergodic_fs_belief_growth.py
======================================================================
Regime (II) of the §7.15 belief-growth trichotomy (Remark 7.15h) is NON-TRIVIAL:
a STATIONARY ERGODIC FILTER-STABLE hidden process can have EXPONENTIAL exact
belief-growth.  Filter stability bounds the EPS-distinct beliefs to poly (the
7.15c truncation), but does NOT bound the EXACT-distinct beliefs, which proliferate
as 2^N on the fractal Blackwell attractor.  Hence the exact belief-DP (cost
Theta(N*B(N)), 7.15h) is blocked in (II) just as in (III); ergodic forgetting does
NOT force exact-poly.  (Existence of ergodic-FS fractal-Blackwell processes:
Jurgens-Crutchfield, Chaos 31:083114 (2021), arXiv:2102.10487.)

Witness: a generic 2-hidden-state, binary-output, NON-UNIFILAR HMM.
  hidden transition A (row-stochastic, strictly positive => irreducible+aperiodic
  => ergodic, geometric forgetting); emission B[s,y]=P(Y=y|S=s) (non-unifilar:
  the symbol does not pin the next hidden state).
  Forward belief b=P(S=1|observed past).  Bayes+transition update b->F(b,y).

  V1 ERGODIC: A strictly positive (Doeblin) => unique stationary pi, mixing.
  V2 FILTER-STABLE: two OPPOSITE initial beliefs b=0,1 fed a COMMON random word
     contract in |.|: |b_t^{(0)}-b_t^{(1)}| -> 0 geometrically, rate rho<1.
  V3 EXPONENTIAL EXACT belief-growth: # distinct beliefs reachable by length-N
     words = 2^{N+1}-1 (doubles each step; non-unifilar => no exact collapse).
  V4 CONTRAST with regime I (renewal/unifilar): a UNIFILar HMM (emission pins the
     next state) has the SAME ergodicity+forgetting but POLY (here O(1)) exact
     belief-growth -- so it is exact-poly. FS is orthogonal to exact-belief-growth.

Deps: numpy.
"""
import numpy as np


def belief_update(b, y, A, B):
    """b=P(S=1|past) scalar. Posterior given Y=y, then push through transition A.
    prior over states p=[1-b, b]; posterior ~ p[s]*B[s,y]; predict s'=sum_s post[s]*A[s,:]."""
    p = np.array([1.0 - b, b])
    post = p * B[:, y]
    Z = post.sum()
    if Z <= 0:
        return b  # unreachable y; leave unchanged
    post = post / Z
    pred = post @ A          # next-state distribution
    return float(pred[1])    # P(S'=1)


def word_prob_and_belief(word, b0, A, B, pi):
    """P(word) and final belief, starting from stationary pi (b0=pi[1])."""
    b = b0
    logp = 0.0
    p = np.array([1.0 - b, b])
    for y in word:
        post = p * B[:, y]
        Z = post.sum()
        logp += np.log(max(Z, 1e-300))
        b = float((post / max(Z, 1e-300)) @ A @ np.array([0.0, 1.0]))
        p = (post / max(Z, 1e-300)) @ A
    return logp, b


def stationary(A):
    w, V = np.linalg.eig(A.T)
    i = np.argmin(np.abs(w - 1.0))
    pi = np.real(V[:, i]); pi = pi / pi.sum()
    return pi


def distinct_beliefs(N, A, B, pi, tol=1e-9):
    """# distinct beliefs over all 2^N words of length N (rounded), starting from pi."""
    b0 = pi[1]
    beliefs = set()
    # BFS over words; track (belief) reachable. Use exact recursion over the tree.
    cur = {round(b0, 10)}
    total = set(cur)
    for t in range(N):
        nxt = set()
        for b in cur:
            for y in (0, 1):
                # only count y reachable (positive prob) from belief b
                p = np.array([1.0 - b, b])
                if (p * B[:, y]).sum() > 1e-12:
                    nb = round(belief_update(b, y, A, B), 10)
                    nxt.add(nb)
        total |= nxt
        cur = nxt
    return len(total)


def reachable_count_by_depth(N, A, B, pi):
    b0 = pi[1]
    counts = []
    cur = {round(b0, 10)}
    seen = set(cur)
    for t in range(N):
        nxt = set()
        for b in cur:
            for y in (0, 1):
                p = np.array([1.0 - b, b])
                if (p * B[:, y]).sum() > 1e-12:
                    nxt.add(round(belief_update(b, y, A, B), 10))
        seen |= nxt
        counts.append(len(seen))
        cur = nxt
    return counts


if __name__ == "__main__":
    print("=" * 76)
    print("Remark 7.15i: ergodic FILTER-STABLE process with EXPONENTIAL exact belief-growth")
    print("=" * 76)

    # --- NON-UNIFILAR ergodic 2-state HMM ---
    A = np.array([[0.7, 0.3],
                  [0.4, 0.6]])          # row-stochastic, strictly positive => ergodic
    B = np.array([[0.8, 0.2],
                  [0.35, 0.65]])        # emission P(Y|S); non-unifilar (both states emit both symbols)
    pi = stationary(A)
    print(f"\n[NON-UNIFILAR HMM]  pi={pi.round(4)}")

    v1 = (A > 0).all()
    print(f"  V1 ERGODIC (A strictly positive, Doeblin): {v1}")

    # V2 filter stability: opposite initial beliefs, common random word
    rng = np.random.default_rng(0)
    diffs = []
    for trial in range(200):
        # generate a word from the stationary process
        s = rng.choice(2, p=pi); word = []
        for _ in range(40):
            y = rng.choice(2, p=B[s]); word.append(y); s = rng.choice(2, p=A[s])
        b0, b1 = 0.0, 1.0
        seq = []
        for y in word:
            b0 = belief_update(b0, y, A, B); b1 = belief_update(b1, y, A, B)
            seq.append(abs(b0 - b1))
        diffs.append(seq)
    diffs = np.array(diffs)
    med = np.median(diffs, axis=0)
    # contraction rate from the geometric decay (median |b0-b1|)
    valid = med[med > 1e-12]
    rho = float(np.exp(np.mean(np.diff(np.log(valid[:15]))))) if len(valid) > 3 else float('nan')
    v2 = (med[15] < 1e-3)
    print(f"  V2 FILTER-STABLE: |b_t^(0)-b_t^(1)| median at t=5,10,15 = "
          f"{med[5]:.2e},{med[10]:.2e},{med[15]:.2e}; contraction rate rho~{rho:.3f}  => {v2}")

    # V3 exponential exact belief-growth
    counts = reachable_count_by_depth(12, A, B, pi)
    expected = [2 ** (t + 2) - 1 for t in range(12)]   # cumulative incl. root: 1+2+4+...
    # cumulative distinct = 1 + 2 + 4 + ... + 2^(t+1) = 2^(t+2)-1
    v3 = counts == expected
    print(f"  V3 EXACT belief-growth B(N) (cumulative distinct), N=1..12:")
    print(f"      got      = {counts}")
    print(f"      2^(N+2)-1= {expected}   match={v3}  (=> EXPONENTIAL, 2^N proliferation)")

    # --- UNIFILAR contrast (regime I): emission pins next state => O(1) beliefs ---
    Au = np.array([[0.7, 0.3], [0.4, 0.6]])
    Bu = np.array([[1.0, 0.0], [0.0, 1.0]])   # unifilar: Y reveals S exactly
    piu = stationary(Au)
    cu = reachable_count_by_depth(12, Au, Bu, piu)
    v4 = max(cu) <= 4
    print(f"\n[UNIFILAR contrast]  belief counts N=1..12: {cu}  (bounded O(1)) => {v4}")
    print(f"  V4 same ergodicity+forgetting but POLY exact-belief-growth => exact-poly (regime I).")

    print("\n" + "=" * 76)
    allok = v1 and v2 and v3 and v4
    print(f"RESULT: {'ALL PASS' if allok else 'CHECK'} -- a STATIONARY ERGODIC FILTER-STABLE process")
    print("  (rho<1 forgetting) can have EXPONENTIAL exact belief-growth B(N)=2^N (non-unifilar),")
    print("  so the exact belief-DP is blocked: ergodic forgetting does NOT force exact-poly.")
    print("  Filter stability bounds EPS-distinct beliefs (7.15c approx), not EXACT-distinct.")
    print("  => Regime II of the 7.15h trichotomy is non-trivial; its exact complexity is OPEN.")
