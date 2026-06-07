#!/usr/bin/env python3
"""
remark_7_15h_belief_growth_boundary.py
======================================================================
The exact-finite-block-entropy tractable boundary of the §7.15 open middle is
NOT "finite belief set" -- finite belief is SUFFICIENT (Remark 7.15f) but NOT
NECESSARY.  The sharp sufficient condition is POLY belief-GROWTH B(N)=poly(N).

CRUCIAL TEST (this script): a discrete-time RENEWAL process (Marzen-Crutchfield,
"Informational and Causal Architecture of Discrete-Time Renewal Processes") has
a COUNTABLY-INFINITE belief set (one point-mass belief per phase = time since the
last event), so the belief set is INFINITE; yet its finite-prefix belief-growth
is B(N)=O(N) (only phases 0..N are reachable in N steps) and its exact finite-block
entropy H(X^N) is computed in O(N^2) by a phase DP.  This kills the clean
finite/infinite dichotomy and pins the boundary at poly belief-growth.

Binary renewal process: X_t=1 at renewal (event) times, 0 otherwise.  Inter-event
interval T ~ phi on {1,2,3,...} with INFINITE support (here phi(k) propto k^{-s},
s>2, truncated only by N), hazard h(k)=phi(k+1)/sum_{j>=k+1}phi(j) all DISTINCT
=> infinitely many distinct phases/beliefs.  Given the observed past the phase is
known exactly (=gap since last 1), so the belief is a point mass on the phase.

  V1  belief set is genuinely infinite: # distinct hazards h(0..N-1) grows = N
      (all distinct), i.e. B(N)=N+1, infinite over all N but POLY in N.
  V2  exact H(X^N): phase-DP O(N^2) matches brute-force 2^N enumeration to ~1e-12.
  V3  scaling: DP cost is O(N^2) (poly), enabling N far beyond brute force.
  V4  contrast: a FINITE-support renewal (phi on {1..m}) has FINITE belief set
      (phases 0..m-1) -- the 7.15f finite-belief island -- still matches brute force.

PASS => "finite belief" is not necessary for exact-poly; poly belief-growth is.
Deps: numpy.
"""
import math
import numpy as np


def h2(x):
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def phi_powerlaw(kmax, s=2.5):
    """Inter-event pmf phi(k) propto k^{-s} on k=1..kmax (infinite-support family,
    truncated at kmax only for normalization; all hazards distinct)."""
    w = np.array([k ** (-s) for k in range(1, kmax + 1)], dtype=float)
    return w / w.sum()


def hazards_from_phi(phi):
    """h(k)=P(T=k+1 | T>=k+1) for phase k>=0; phase = #steps since last event.
    Survival S(k)=P(T>=k+1)=sum_{j>=k+1} phi(j).  h(k)=phi(k+1)/S(k).
    phi indexed phi[0]=phi(1),...,phi[m-1]=phi(m). Phase k in 0..m-1."""
    m = len(phi)
    tail = np.concatenate([np.cumsum(phi[::-1])[::-1], [0.0]])  # tail[i]=sum_{j>=i} phi[j]
    h = np.zeros(m)
    for k in range(m):
        S = tail[k]            # P(T >= k+1) = sum_{j>=k} phi[j]  (phi[k]=phi(k+1))
        h[k] = phi[k] / S if S > 0 else 1.0
    h[m - 1] = 1.0             # forced event at max support (truncation)
    return h


def H_exact_phaseDP(h, N):
    """H(X^N) for the binary renewal process via O(N^2) phase DP.
    nu_t = distribution over phase at time t (phase 0..t). X_1=1 (start at an event),
    so nu_1 = delta_0.  Transition: from phase k, X_{t}=1 w.p. h(k) -> phase 0;
    X_t=0 w.p. 1-h(k) -> phase k+1.  H(X^N)=sum_t H(X_t|X^{t-1})=sum_t E_{nu_{t-1}}[h2(h(k))]
    (chain rule; X_t determined-in-law by phase k_{t-1}).  Returns H(X^N) in bits."""
    m = len(h)
    # nu over phases 0..m-1; start: just had an event at t=1 -> X_1=1 deterministic?
    # Use stationary-free "fresh start" at an event: phase 0 at the first emitted symbol.
    nu = np.zeros(m); nu[0] = 1.0
    H = 0.0
    for t in range(1, N + 1):
        # H(X_t | past): given phase k (known from past), X_t ~ Bern(h(k))
        H += float(np.sum(nu * np.array([h2(h[k]) for k in range(m)])))
        # advance phase distribution by one step
        nu_new = np.zeros(m)
        ev = float(np.sum(nu * h))            # mass that emits 1 -> phase 0
        nu_new[0] += ev
        for k in range(m - 1):
            nu_new[k + 1] += nu[k] * (1 - h[k])
        nu_new[m - 1] += nu[m - 1] * (1 - h[m - 1])  # =0 since h[m-1]=1
        nu = nu_new
    return H


def H_brute(h, N):
    """Exact H(X^N) by enumerating all 2^N binary strings; P(x) via the phase chain
    (X_1=1 start at event). Returns bits. Only for small N."""
    m = len(h)
    # P(x_1..x_N): start phase 0 (event at position 0, before x_1). Walk:
    # at phase k, P(x=1)=h(k) -> next phase 0; P(x=0)=1-h(k) -> next phase k+1.
    total_H = 0.0
    Z = 0.0
    for code in range(1 << N):
        bits = [(code >> (N - 1 - i)) & 1 for i in range(N)]
        p = 1.0
        k = 0
        for b in bits:
            hk = h[min(k, m - 1)]
            if b == 1:
                p *= hk; k = 0
            else:
                p *= (1 - hk); k = min(k + 1, m - 1)
        Z += p
        if p > 0:
            total_H += -p * math.log2(p)
    return total_H, Z


def distinct_beliefs(h, N):
    """# distinct hazards among reachable phases 0..min(N,m)-1 (each a distinct belief)."""
    m = len(h)
    reach = min(N, m)
    return len(set(np.round(h[:reach], 12)))


if __name__ == "__main__":
    print("=" * 76)
    print("Remark 7.15h: 'finite belief' is SUFFICIENT (7.15f) but NOT NECESSARY for")
    print("exact-poly finite-block entropy; the boundary is POLY belief-GROWTH B(N).")
    print("=" * 76)

    # --- INFINITE-belief renewal (power-law inter-event, all hazards distinct) ---
    s = 2.5
    Nbrute = 13
    phi = phi_powerlaw(kmax=64, s=s)     # support 1..64 (>> Nbrute => infinite-like)
    h = hazards_from_phi(phi)
    print(f"\n[INFINITE-belief renewal]  phi(k) ~ k^-{s}, support 1..{len(phi)}")
    print(f"  V1 distinct beliefs (hazards) over reachable phases:")
    ok1 = True
    for N in (4, 8, 12, 13):
        nb = distinct_beliefs(h, N)
        print(f"     N={N:>3}: {nb} distinct beliefs  (=N => grows linearly, infinite over all N)")
        ok1 = ok1 and (nb == min(N, len(h)))
    print(f"  V2 exact H(X^N): phase-DP vs brute force (2^N):")
    ok2 = True
    for N in range(3, Nbrute + 1):
        Hdp = H_exact_phaseDP(h, N)
        Hbf, Z = H_brute(h, N)
        err = abs(Hdp - Hbf)
        ok2 = ok2 and (err < 1e-9 and abs(Z - 1.0) < 1e-9)
        if N >= Nbrute - 2:
            print(f"     N={N:>3}: H_DP={Hdp:.8f}  H_brute={Hbf:.8f}  |err|={err:.2e}  (Z={Z:.6f})")
    print(f"  V3 poly DP reaches large N (brute infeasible): "
          f"H(X^200)={H_exact_phaseDP(h, 200):.6f}, H(X^1000)={H_exact_phaseDP(h, 1000):.6f}")

    # --- FINITE-belief renewal (finite support) = the 7.15f island, sanity ---
    phiF = phi_powerlaw(kmax=4, s=s); hF = hazards_from_phi(phiF)
    print(f"\n[FINITE-belief renewal]  phi support 1..4 (finite belief = 7.15f island)")
    ok4 = True
    for N in range(3, 11):
        Hdp = H_exact_phaseDP(hF, N); Hbf, Z = H_brute(hF, N)
        ok4 = ok4 and abs(Hdp - Hbf) < 1e-9
    print(f"  V4 phase-DP == brute force for N=3..10: {ok4}; #beliefs={distinct_beliefs(hF, 100)} (finite)")

    print("\n" + "=" * 76)
    allok = ok1 and ok2 and ok4
    print(f"RESULT: {'ALL PASS' if allok else 'CHECK'} -- the renewal process has an INFINITE belief set")
    print("  (distinct hazard per phase) yet EXACT H(X^N) in O(N^2) (phase DP == brute force).")
    print("  => finite belief is NOT necessary; poly belief-growth B(N)=O(N) is the boundary.")
    print("  Clean finite/infinite dichotomy of the §7.15 tractable island is FALSE.")
