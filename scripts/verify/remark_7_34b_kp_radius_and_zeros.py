#!/usr/bin/env python3
"""
remark_7_34b_kp_radius_and_zeros.py
============================================================================
NEGATIVE-BRANCH DOCUMENTATION for the dual-identity route (Lemmas 7.34b'/b''):
the proposed CLUSTER-EXPANSION CLOSURE of the Gray-region achievability is
REFUTED.  Records the three refutation anchors so the route is not re-attempted:

  N1  the Kotecky-Preiss convergence radius r_KP(p) of the polymer gas is
      4-10x SMALLER than the contour maximum max_s |eta_s| at D=D_c, for EVERY
      p in (0,1/2).  [The single-flip activity threshold (p/(1-p))^2 quoted in
      the breakthrough memo is NOT the cluster radius.]
  N2  S_x(eta) has genuine complex zeros (Fisher / Lee-Yang at negative
      fugacity) lying ON or NEAR the eta-contour at s ~ +-pi for a TYPICAL
      fraction of source words => log S_x(eta_s) is singular on the contour
      => no contour-wide cluster expansion. (The zeros sit at the benign
      anti-periodic endpoint where |Phi(pi;x)| -> 0 -- a span-1 LLT signature
      -- so they break the METHOD, not the lattice LLT itself.)
  N3  |C(theta0+is)/C(theta0)| > 1 on (0,pi] (peaks at s=pi): the deterministic
      factor GROWS; all decay of |Phi| must come from S_x -- the Theta(n)-vs-
      Theta(n) cancellation again (the committed 7.34b cancellation obstruction,
      now in dual form).
KP criterion used (1-D hard-core polymers, activities bounded by
B^2 u^L, B=((1-p)/p)^2, u=|eta| e^mu):  sum over polymers touching a site of
|z_I| e^{mu |I|} <= mu, i.e.  B(2u/(1-u) + u/(1-u)^2) <= mu  for some mu>0.
"""
import math
import numpy as np

def Dc(p): return 0.5*(1-math.sqrt(1-2*p)/(1-p))

def kp_feasible(eta_abs, p):
    B = ((1-p)/p)**2
    for mu in np.linspace(1e-4, 6.0, 6000):
        u = eta_abs*math.exp(mu)
        if u >= 1: continue
        if B*(2*u/(1-u) + u/(1-u)**2) <= mu:
            return True
    return False

def r_cluster(p):
    lo, hi = 1e-6, 0.999
    for _ in range(50):
        mid = 0.5*(lo+hi)
        if kp_feasible(mid, p): lo = mid
        else: hi = mid
    return lo

def eta_contour(p, D, s):
    th0 = math.log(D/(1-D)); eb = np.exp(th0+1j*s)
    ze = (1-eb)/((1+eb)*(1-2*D)); return (1-ze)/(1+ze)

def C_ratio(p, D, s):
    th0 = math.log(D/(1-D))
    eb = np.exp(th0+1j*s); eb0 = math.exp(th0)
    ze = (1-eb)/((1+eb)*(1-2*D)); ze0 = (1-eb0)/((1+eb0)*(1-2*D))
    C = (1+eb)*(1+ze)/2; C0 = (1+eb0)*(1+ze0)/2
    return np.abs(C/C0)

def block_dp_poly(x, n, ratio):
    memo=[None]*(n+2); memo[n]=np.array([1.0])
    for i in range(n-1,-1,-1):
        poly=memo[i+1].copy()
        for b in range(i,n):
            L=b-i+1; w=1.0
            if i>=1: w*=ratio if (x[i]==x[i-1]) else 1.0/ratio
            if b+1<=n-1: w*=ratio if (x[b+1]==x[b]) else 1.0/ratio
            nxt=memo[b+2] if b+2<=n else np.array([1.0])
            term=np.concatenate([np.zeros(L), w*nxt])
            m=max(len(poly),len(term)); pp=np.zeros(m); pp[:len(poly)]+=poly; pp[:len(term)]+=term; poly=pp
        memo[i]=poly
    return memo[0]

if __name__ == "__main__":
    print("="*84)
    print("Dual-route NEGATIVE anchors: KP radius vs contour; Fisher zeros; |C-ratio|>1")
    print("="*84)
    # N1: KP radius vs contour max
    print(f"  N1: {'p':>5} {'r_KP':>9} {'max|eta|@Dc':>12} {'ratio':>6}")
    n1 = True
    for p in (0.05, 0.1, 0.25, 0.4, 0.49):
        rkp = r_cluster(p)
        s = np.linspace(-math.pi, math.pi, 4001)
        me = float(np.max(np.abs(eta_contour(p, Dc(p), s))))
        n1 &= me > rkp     # the refutation: contour EXCEEDS the KP disc
        print(f"      {p:>5} {rkp:>9.5f} {me:>12.5f} {me/rkp:>6.1f}x")
    print(f"      => contour exceeds the KP disc at every p: {n1} (closure via contour-wide KP: REFUTED)")
    # N2: Fisher zeros near the contour (typical words)
    rng = np.random.default_rng(11); n2_hit = 0; R = 200; n = 24; p = 0.4
    D = 0.9*Dc(p); ratio = p/(1-p)
    s = np.linspace(-math.pi, math.pi, 1001); et = eta_contour(p, D, s)
    minvals = []
    for _ in range(R):
        bits = (np.cumsum(rng.random(n) < p) % 2).astype(int)
        c = block_dp_poly(bits, n, ratio)
        Sx = np.polyval(c[::-1], et)
        minvals.append(float(np.min(np.abs(Sx))))
        r = np.roots(c[::-1])
        if np.min(np.abs(r[:, None]-et[None, :])) < 0.03: n2_hit += 1
    n2 = (np.median(minvals) < 0.01) and (n2_hit > R//10)
    print(f"  N2: p={p}, n={n}: median min_s|S_x(eta_s)| = {np.median(minvals):.1e}; "
          f"{n2_hit}/{R} words have a zero within 0.03 of the contour => log S_x singular: {n2}")
    # N3: |C-ratio| > 1 on (0,pi]
    s = np.linspace(0.05, math.pi, 500)
    n3 = True
    for p in (0.25, 0.4):
        cr = C_ratio(p, 0.9*Dc(p), s)
        n3 &= bool(np.max(cr) > 1.0)
        print(f"  N3: p={p}: max_s |C_s/C_0| = {np.max(cr):.4f} at s={s[np.argmax(cr)]:.2f} "
              f"(>1: deterministic factor GROWS; decay must come from S_x)")
    allok = n1 and n2 and n3
    print("="*84)
    print(f"RESULT: {'ALL ANCHORS CONFIRMED' if allok else 'CHECK'} -- the contour-wide cluster-expansion")
    print("  closure is REFUTED (do not re-attempt); the residual is the Theta(n)-vs-Theta(n)")
    print("  cancellation (now: C^n growth vs S_x Fisher zeros), i.e. the same open quenched")
    print("  shell-uniform lattice tail bound, in the (elliptic) dual representation.")
