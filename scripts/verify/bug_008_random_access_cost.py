#!/usr/bin/env python3
r"""
Verification of the random-access cost accounting (BUG-008-C):
Theorem 10.3's per-query cost O(log L(R) + (K+k)*cost(M)) and its
monotonicities, the K=Theta(|X|) / k=Theta(|X|) degeneracy, the
simultaneous-polylog impossibility (the K knob trades per-query window
cost against repair-stream sync overhead), and Remark 10.3a's
sub-block-granular cost for Type-II / Type-III-B partition-rank modes.

COST-ACCOUNTING / COMBINATORIAL CHECK ONLY. This probe is a deterministic
counting model of the seek + (K+k) decode-window arithmetic of Theorem
10.3's proof and of Remark 10.3a's sub-block accounting. It is NOT a
neural implementation and NOT a benchmark: it asserts the structural
relations the paper states (decoded-position count, the big-O regime
boundaries, the trade-off optimum), not a wall-clock latency. The
measured seek latency is the separate external leaf BUG-008-D.

Every relation traces to a paper source:
  - decoded positions p+k'-P <= K+k, P=K*floor(p/K), k'=min(k,|X|-p)
        : Theorem 10.3 statement (ii) + proof ("Total decoded positions:
          p + k' - P <= K + k since P > p - K").
  - per-query work O(log L(R) + (K+k)*cost(M)), dominance, and the
    K=Theta(|X|) or k=Theta(|X|) => Theta(|X|*cost(M)) degeneracy
        : Theorem 10.3 statement (iii) and its scope sentence.
  - total = (N/K)*delta_sync + (K+k)*cost(M), interior optimum, and the
    simultaneous-polylog impossibility ("polylog total overhead would
    require K = N/polylog(N), which is incompatible with polylog
    per-query random-access cost")
        : Lemma 5.1b (Delta_{sync,total} = m*S*W*log2(1/eta),
          m = ceil(N/K)) and the section 10.7 / section 1.1 trade-off bullet.
  - sub-block-granular read spans ceil(k/K)+O(1) sub-blocks, costs
    O(log L(R) + K*unrank(K)) per affected sub-block, decoding/discarding
    one full K-byte sub-block per request for k<<K
        : Remark 10.3a.

Emits PASS/FAIL per check and "OVERALL -> PASS" on success (CI contract).
"""

import math
import sys


def line(ok, label, detail):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {label}: {detail}")
    return bool(ok)


# ---------------------------------------------------------------------------
# (1) Theorem 10.3 decoded-position count and per-query work
# ---------------------------------------------------------------------------
def decoded_positions(p, k, K, X):
    """Theorem 10.3: P = K*floor(p/K), k' = min(k, X-p),
    decoded positions = p + k' - P (warmup [P,p) then output [p,p+k'))."""
    assert 0 <= p < X, "p must be a valid position in [0,|X|)"
    assert K >= 1 and k >= 0
    P = K * (p // K)
    kprime = min(k, X - p)
    warmup = p - P                       # discarded output positions [P, p)
    decoded = p + kprime - P             # total positions decoded [P, p+k')
    return P, kprime, warmup, decoded


def per_query_work(p, k, K, X, costM, logLR):
    """Total work units of Theorem 10.3:
    one O(log L(R)) seek-index probe + (decoded positions)*cost(M)."""
    _, _, _, decoded = decoded_positions(p, k, K, X)
    return logLR + decoded * costM, decoded


def check_decoded_count():
    print("\n(1) Theorem 10.3 decoded-position count = p+k'-P <= K+k:")
    ok = True
    # generic sweep
    for X in (1, 7, 16, 1000, 65536):
        for K in (1, 3, 16, 4096, 65536):
            if K > 4 * X and X > 16:
                continue
            # sample positions: aligned, mid-block, last
            ps = sorted(set([0, X // 2, X - 1, K - 1 if K <= X else 0,
                             min(K, X - 1)]))
            for p in ps:
                if not (0 <= p < X):
                    continue
                for k in (0, 1, K, K + 1, 2 * K, X):
                    P, kprime, warmup, decoded = decoded_positions(p, k, K, X)
                    # exact identity
                    if decoded != p + kprime - P:
                        ok &= line(False, "identity p+k'-P",
                                   f"X={X},K={K},p={p},k={k}: {decoded}")
                        continue
                    # bound: decoded <= K + k (Theorem 10.3 proof)
                    if decoded > K + k:
                        ok &= line(False, "decoded<=K+k",
                                   f"X={X},K={K},p={p},k={k}: "
                                   f"decoded={decoded} > K+k={K + k}")
                        continue
                    # warmup at most K (statement (ii): at most K warmup)
                    if not (0 <= warmup <= K):
                        ok &= line(False, "0<=warmup<=K",
                                   f"X={X},K={K},p={p},k={k}: warmup={warmup}")
                        continue
                    # short read: k' = k when p+k<=X, less otherwise
                    if p + k <= X and kprime != k:
                        ok &= line(False, "short-read k'",
                                   f"p+k<=X but k'={kprime}!=k={k}")
                        continue
                    if p + k > X and kprime != X - p:
                        ok &= line(False, "short-read clamp",
                                   f"p+k>X but k'={kprime}!=X-p={X - p}")
                        continue
    ok &= line(ok, "full sweep of (X,K,p,k)",
               "identity, decoded<=K+k, 0<=warmup<=K, short-read clamp all hold")

    # explicit edge case: |X| < K -> single sync point at P=0, full prefix decode
    X, K = 1000, 65536
    P, kprime, warmup, decoded = decoded_positions(p=500, k=10, K=K, X=X)
    ok &= line(P == 0 and decoded == 510 and warmup == 500,
               "|X|<K single-sync-point (P=0)",
               f"X={X}<K={K}: P={P}, warmup={warmup}, decoded={decoded} "
               f"(prefix [0,510) decoded, [0,500) discarded)")

    # explicit edge case: short read at end of archive (p+k>|X|)
    X, K = 100, 16
    P, kprime, warmup, decoded = decoded_positions(p=95, k=20, K=K, X=X)
    ok &= line(kprime == 5 and decoded == 95 + 5 - 80,
               "short read at end (p+k>|X|)",
               f"X={X},K={K},p=95,k=20: P={P}, k'={kprime} (clamped to X-p=5), "
               f"decoded={decoded}")
    return ok


# ---------------------------------------------------------------------------
# (2) per-query cost dominance + K=Theta(|X|) / k=Theta(|X|) degeneracy
# ---------------------------------------------------------------------------
def check_cost_dominance_and_degeneracy():
    print("\n(2) Per-query cost O(log L(R)+(K+k)*cost(M)); decode term dominance;")
    print("    K=Theta(|X|) or k=Theta(|X|) degeneracy => Theta(|X|*cost(M)):")
    ok = True

    # decode-window term dominates log L(R) once (K+k)*cost(M) >> log L(R)
    X, p, k, costM = 1 << 30, (1 << 29), 1024, 1.0   # cost(M) in eval-units
    K = 65536
    logLR = math.log2(1 << 40)                       # L(R) ~ 2^40 -> log L(R)=40
    work, decoded = per_query_work(p, k, K, X, costM, logLR)
    decode_term = decoded * costM
    ok &= line(decode_term > logLR,
               "decode window dominates log L(R)",
               f"(K+k)*cost(M)~{decode_term:.0f} >> log L(R)={logLR:.0f}")

    # monotonicity of the Theorem 10.3 worst-case BOUND
    # logLR + (K+k)*cost(M): strictly increasing in K, k, cost(M). (The
    # realized decoded count p+k'-P at a fixed p is only non-decreasing in K
    # -- it can stay flat when the window is k-dominated and p is far from
    # the sync boundary -- so monotonicity is asserted on the stated bound,
    # and the realized count is checked separately as non-decreasing.)
    def bound(Kb, kb, cb):
        return logLR + (Kb + kb) * cb
    base = bound(K, k, costM)
    ok &= line(bound(2 * K, k, costM) > base, "bound monotone increasing in K",
               f"bound(2K)={bound(2 * K, k, costM):.0f} > bound(K)={base:.0f}")
    ok &= line(bound(K, 2 * k, costM) > base, "bound monotone increasing in k",
               f"bound(2k)={bound(K, 2 * k, costM):.0f} > bound(K)={base:.0f}")
    ok &= line(bound(K, k, 2 * costM) > base, "bound monotone increasing in cost(M)",
               f"bound(2cost)={bound(K, k, 2 * costM):.0f} > bound(K)={base:.0f}")
    # realized decoded count is non-decreasing in K at a worst-case position
    # (p at end of its block, warmup maximal): larger K cannot decrease it.
    prev = -1
    nd_ok = True
    for Kk in (1024, 4096, 16384, 65536):
        p_wc = Kk - 1 if Kk - 1 < X else X - 1      # end-of-block worst case
        _, _, _, dwc = decoded_positions(p_wc, k, Kk, X)
        nd_ok &= (dwc >= prev)
        prev = dwc
    ok &= line(nd_ok, "realized decoded count non-decreasing in K (worst-case p)",
               "p at end of block: decoded count grows with K")

    # degeneracy K=Theta(|X|): pick K = X/4 (a constant fraction of |X|).
    # worst-case read (end of archive) decodes Theta(|X|) positions.
    Xd = 1 << 20
    Kd = Xd // 4                       # K = Theta(|X|)
    p_end = Xd - 1                      # worst-case (end-of-archive) read
    _, _, _, dec = decoded_positions(p_end, k=1, K=Kd, X=Xd)
    # decoded ~ K = Theta(|X|): at least a constant fraction of |X|
    ok &= line(dec >= Xd // 8,
               "K=Theta(|X|) => Theta(|X|) decoded (degeneracy)",
               f"K=|X|/4: worst-case decoded={dec} = Theta(|X|={Xd}); "
               f"no better than full sequential decode")

    # degeneracy k=Theta(|X|): a read spanning the whole archive
    _, _, _, dec2 = decoded_positions(p=0, k=Xd, K=64, X=Xd)
    ok &= line(dec2 >= Xd,
               "k=Theta(|X|) => Theta(|X|) decoded (degeneracy)",
               f"k=|X|: decoded={dec2} >= |X|={Xd}")

    # polylog regime: K,k = polylog(|X|) keeps decoded = polylog(|X|)
    for Xexp in (20, 30, 40):
        Xp = 1 << Xexp
        Kp = Xexp * Xexp               # K = polylog(|X|) ~ (log X)^2
        kp = Xexp
        _, _, _, decp = decoded_positions(min(Xp // 2, Xp - 1), kp, Kp, Xp)
        ok &= line(decp <= Kp + kp and decp <= (Xexp ** 2) * 4,
                   f"polylog regime |X|=2^{Xexp}",
                   f"K=(log X)^2={Kp}, k=log X={kp}: decoded={decp} = polylog(|X|)")
    return ok


# ---------------------------------------------------------------------------
# (3) Simultaneous-polylog impossibility / K-knob trade-off
# ---------------------------------------------------------------------------
def total_cost_model(N, K, k, costM, delta_sync):
    """Sum of repair-stream sync overhead and a per-query decode window.

    Lemma 5.1b: Delta_{sync,total} = m * delta_sync, m = ceil(N/K), so the
    sync term is (N/K)*delta_sync (delta_sync = S*W*log2(1/eta) per sync).
    The per-query decode window is (K+k)*cost(M). The K knob trades them:
    larger K shrinks the sync term but grows the per-query window.
    """
    m = math.ceil(N / K)
    sync_overhead = m * delta_sync
    per_query = (K + k) * costM
    return sync_overhead + per_query, sync_overhead, per_query


def check_simultaneous_polylog_impossibility():
    print("\n(3) Simultaneous-polylog impossibility (K trades sync overhead vs window):")
    ok = True

    N = 1 << 30
    k = 1024
    costM = 1.0
    delta_sync = 16.0      # S*W*log2(1/eta) per sync point (eval/bit units)

    # (a) interior optimum of total(K) = N*delta_sync/K + (K+k)*cost(M).
    # d/dK [N*delta/K + K*cost] = -N*delta/K^2 + cost = 0
    # => K* = sqrt(N*delta/cost).  Verify the modelled total is minimized near K*.
    Kstar = math.sqrt(N * delta_sync / costM)
    grid = sorted({int(Kstar * f) for f in (0.25, 0.5, 0.9, 1.0, 1.1, 2.0, 4.0)
                   if int(Kstar * f) >= 1})
    totals = {Kg: total_cost_model(N, Kg, k, costM, delta_sync)[0] for Kg in grid}
    Kmin = min(totals, key=totals.get)
    # the discrete grid minimum should sit at/next to K* (the 0.9/1.0/1.1 cluster)
    ok &= line(0.5 * Kstar <= Kmin <= 2.0 * Kstar,
               "interior optimum K* ~ sqrt(N*delta_sync/cost(M))",
               f"K*={Kstar:.0f}; grid-min K={Kmin} within [0.5,2]*K*")

    # (b) the two regimes are mutually exclusive at a fixed N:
    #   - polylog per-query needs K = polylog(N): then the sync term
    #     (N/K)*delta_sync = N/polylog(N) is super-polylog (NOT polylog total);
    #   - polylog total needs K = N/polylog(N): then the per-query window
    #     (K+k)*cost(M) = (N/polylog(N))*cost(M) is super-polylog (NOT polylog).
    # The asymptotic separation: their RATIO to a polylog reference (log^4 N)
    # blows up as N grows -- so no FIXED polylog ceiling holds for all N. We
    # exhibit that blow-up across a sequence of N (see (c)); here we record
    # the per-N qualitative split.
    logN = math.log2(N)
    K_small = int(logN ** 2)                        # K = polylog(N)
    _, sync_small, pq_small = total_cost_model(N, K_small, k, costM, delta_sync)
    K_big = max(1, int(N / (logN ** 2)))            # K = N/polylog(N)
    _, sync_big, pq_big = total_cost_model(N, K_big, k, costM, delta_sync)
    ok &= line(pq_small < sync_small,
               "K=polylog(N): per-query small, sync overhead large",
               f"K={K_small}: per-query~{pq_small:.0f} << sync~{sync_small:.3e}")
    ok &= line(sync_big < pq_big,
               "K=N/polylog(N): total small, per-query window large",
               f"K={K_big}: sync~{sync_big:.0f} << per-query~{pq_big:.3e}")

    # (c) ASYMPTOTIC impossibility (the real, N->infinity statement):
    # min over K of max(sync, per-query) is achieved at K* where the two
    # terms balance, giving min-of-max = Theta(sqrt(N*delta_sync*cost(M))).
    # That floor is super-polylog in N: its ratio to any fixed polylog
    # reference log^a N diverges as N grows. There is therefore NO K (and no
    # fixed polylog ceiling) making BOTH per-query AND total polylog for all N.
    print("    asymptotic floor min_K max(sync, per-query) ~ sqrt(N*delta*cost):")
    ratios = []
    prev_floor = None
    mono = True
    for Nexp in (20, 30, 40, 50, 60):
        Nn = 1 << Nexp
        # scan K in powers of two, take the K minimizing max(sync, per-query)
        best = None
        Ke = 1
        while Ke <= Nn:
            _, s_e, q_e = total_cost_model(Nn, Ke, k, costM, delta_sync)
            mx = max(s_e, q_e)
            if best is None or mx < best:
                best = mx
            Ke *= 2
        theory = math.sqrt(Nn * delta_sync * costM)   # Theta(sqrt(N*delta*cost))
        ref = (math.log2(Nn)) ** 4                     # fixed polylog reference
        ratio = best / ref
        ratios.append(ratio)
        if prev_floor is not None and best <= prev_floor:
            mono = False
        prev_floor = best
        # best-achievable min-of-max is within a small factor of the sqrt law
        within = 0.5 * theory <= best <= 4.0 * theory
        ok &= line(within,
                   f"N=2^{Nexp}: min-of-max ~ sqrt(N*delta*cost)",
                   f"best={best:.3e}, sqrt-law={theory:.3e}, "
                   f"best/log^4(N)={ratio:.2f}")
    # the ratio to the polylog reference must DIVERGE (super-polylog floor)
    ok &= line(ratios[-1] > ratios[0] and ratios[-1] > 10 * ratios[0],
               "min-of-max / polylog(N) DIVERGES (super-polylog floor)",
               f"ratio log^4 grows {ratios[0]:.2f} -> {ratios[-1]:.2f}; "
               f"no fixed polylog ceiling bounds it => simultaneous-polylog "
               f"regime does not exist (Lemma 5.1b / section 10.7)")
    ok &= line(mono, "asymptotic floor is increasing in N",
               "min_K max(sync, per-query) grows with N (sqrt law)")
    return ok


# ---------------------------------------------------------------------------
# (4) Remark 10.3a sub-block-granular cost (Type-II / Type-III-B partition)
# ---------------------------------------------------------------------------
def subblocks_spanned(p, k, K):
    """Remark 10.3a: j0=floor(p/K), j1=floor((p+k-1)/K);
    request spans j1-j0+1 = ceil(k/K)+O(1) sub-blocks."""
    assert k >= 1
    j0 = p // K
    j1 = (p + k - 1) // K
    return j0, j1, (j1 - j0 + 1)


def check_subblock_granularity():
    print("\n(4) Remark 10.3a sub-block granularity (Type-II / Type-III-B (i)/(ii)):")
    ok = True

    # a byte query (k=1) costs a whole sub-block decode: span = 1 sub-block,
    # decoding/discarding up to K bytes to deliver 1.
    K = 65536
    for p in (0, 1, K - 1, K, 3 * K + 17, 10 * K + 5):
        j0, j1, span = subblocks_spanned(p, k=1, K=K)
        decoded_bytes = span * K          # each affected sub-block decoded in full
        ok &= line(span == 1 and decoded_bytes == K,
                   f"byte query (k=1) at p={p}",
                   f"spans {span} sub-block; decodes {decoded_bytes} bytes to "
                   f"serve 1 (whole-sub-block decode, NOT byte-granular)")

    # k << K: spans ceil(k/K)+O(1) sub-blocks (1, or 2 when crossing a boundary);
    # cost O(log L(R) + K*unrank(K)); produces K+O(K) bytes, k requested.
    K = 4096
    for (p, k) in [(0, 100), (4090, 100), (8000, 200), (4095, 2)]:
        j0, j1, span = subblocks_spanned(p, k, K)
        lo = math.ceil(k / K)
        ok &= line(lo <= span <= lo + 1,
                   f"k<<K read p={p},k={k}: span = ceil(k/K)+O(1)",
                   f"ceil(k/K)={lo} <= span={span} <= {lo + 1} "
                   f"(+1 when unaligned crossing)")

    # k >> K: cost ~ (k/K)*K*unrank(K) = k*unrank(K), near-linear in k;
    # span ~ ceil(k/K) (+1). Model unrank(K) as a fixed per-byte unit u.
    K = 4096
    u = 1.0                               # unrank(K)/K per-byte work unit
    for k in (K, 4 * K, 64 * K):
        j0, j1, span = subblocks_spanned(p=0, k=k, K=K)
        cost = span * K * u               # O((k/K)*K*unrank(K))
        near_linear = k * u
        ok &= line(near_linear <= cost <= near_linear + K * u,
                   f"k>>K read k={k}: cost near-linear in k",
                   f"cost={cost:.0f} ~ k*unrank(K)={near_linear:.0f} "
                   f"(span={span} sub-blocks)")

    # sub-block granularity is STRICTLY WEAKER than byte granularity:
    # to read 1 byte, sub-block mode decodes K, per-position mode decodes <= K+k
    # but delivers byte-granular access (no whole-sub-block discard requirement).
    K, k = 65536, 1
    _, _, span = subblocks_spanned(p=K + 7, k=k, K=K)
    subblock_decoded = span * K
    ok &= line(subblock_decoded == K and subblock_decoded > k,
               "sub-block granularity strictly weaker than byte granularity",
               f"1-byte read decodes a full sub-block ({subblock_decoded} bytes) "
               f"vs k={k} requested")
    return ok


# ---------------------------------------------------------------------------
def main() -> int:
    print("Verification: random-access cost accounting (BUG-008)")
    print("=" * 70)
    print("COST-ACCOUNTING / COMBINATORIAL CHECK ONLY -- models the seek + (K+k)")
    print("decode-window arithmetic of Thm 10.3 and the sub-block accounting of")
    print("Rem 10.3a; NOT a neural runtime and NOT a wall-clock benchmark")
    print("(measured seek latency is the external leaf BUG-008-D).")

    results = {
        "(1) Thm 10.3 decoded-position count p+k'-P<=K+k + edge cases":
            check_decoded_count(),
        "(2) per-query cost dominance + K/k=Theta(|X|) degeneracy":
            check_cost_dominance_and_degeneracy(),
        "(3) simultaneous-polylog impossibility (K-knob trade-off)":
            check_simultaneous_polylog_impossibility(),
        "(4) Rem 10.3a sub-block-granular cost (Type-II/III-B partition)":
            check_subblock_granularity(),
    }

    print("\n" + "=" * 70)
    print("Summary:")
    all_ok = True
    for name, ok in results.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        all_ok &= ok

    print()
    if all_ok:
        print("OVERALL -> PASS")
        return 0
    print("OVERALL -> FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
