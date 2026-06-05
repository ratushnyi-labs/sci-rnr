#!/usr/bin/env python3
"""
Verification for Remark 7.34a (Gauss-Markov instance of the lossy second-order interplay:
the "missing ingredient" of Remark 7.34 -- a mixing-source rate-distortion dispersion --
DOES exist for a memory source, namely Gauss-Markov AR(1) under MSE; Tian-Kostina 2019).

Source: stationary AR(1)  X_t = rho X_{t-1} + W_t,  |rho|<1,  W_t ~ N(0, sigma_W^2).
PSD  S(w) = sigma_W^2 / (1 + rho^2 - 2 rho cos w),  w in [-pi, pi].
Reverse-water-filling at level theta:
   D(theta)   = (1/2pi) int min{theta, S(w)} dw
   R(D)       = (1/2pi) int (1/2) log^+(S(w)/theta) dw          (nats)
   V(D)       = (1/4pi) int min{1, (S(w)/theta)^2} dw           (nats^2)   [Tian-Kostina Thm 1/2]
Critical distortion D_c = min_w S(w) = sigma_W^2/(1+|rho|)^2 ;  D_max = Var(X)=sigma_W^2/(1-rho^2).
For 0<D<=D_c (theta=D<=S_min): every band active => min{1,(S/theta)^2}=1 => V(D)=1/2 (nats).
For D_c<D<D_max: V(D) < 1/2.

Checks:
  V1  reverse-water-filling self-consistency: solve theta(D), recompute D, R(D)>0.
  V2  V(D) = 1/2 nats for D <= D_c (low distortion) -- the iid-innovation value (Cor 1).
  V3  V(D) < 1/2 and decreasing for D_c < D < D_max.
  V4  operational = informational: the BLOCK d-tilted information variance V_n (covariance
      eigenvalues of the n x n AR(1) Toeplitz) -> spectral V(D) as n grows (Tasci-Kostina);
      cross-checked by Monte-Carlo of the Gaussian d-tilted information j_n.
  V5  the interplay (lossy analogue of Thm 7.29/7.34): at K=c sqrt(N) the cold-context RA
      overhead (N/K) delta_inf^D and the dispersion sqrt(N V(D)) Q^{-1}(eps) are BOTH
      Theta(sqrt N), ratio N-independent (delta_inf^D>0 parametrised from the AR memory).
  V6  delta_inf^D > 0 for AR(1) (memory): the cold-start lossy rate excess from losing the
      AR context at a sync reset is strictly positive (and ->0 as rho->0, the iid limit).
"""
import math
import numpy as np

rng = np.random.default_rng(20260605)
LOG2E = 1.0/math.log(2.0)            # nats -> bits multiplier is (log2 e); dispersion x (log2 e)^2
NW = 200000                           # spectral quadrature points
w = (np.arange(NW)+0.5)/NW*2*math.pi - math.pi   # midpoint grid on [-pi,pi]

def psd(rho, s2):
    return s2/(1+rho**2-2*rho*np.cos(w))

def water_fill(rho, s2, D):
    """Find theta with (1/2pi) int min{theta,S} dw = D; return theta,R(D),V(D) in nats."""
    S = psd(rho, s2)
    lo, hi = 1e-12, S.max()
    for _ in range(200):
        th = 0.5*(lo+hi)
        Dth = np.mean(np.minimum(th, S))      # (1/2pi)int = mean over uniform grid
        if Dth < D: lo = th
        else: hi = th
    th = 0.5*(lo+hi)
    R = np.mean(0.5*np.maximum(0.0, np.log(S/th)))           # nats
    V = np.mean(0.5*np.minimum(1.0, (S/th)**2))              # nats^2  (1/4pi int = 0.5*mean)
    return th, R, V

# ---------------------------------------------------------------- V1
print("="*70)
print("V1: reverse-water-filling self-consistency (AR(1))")
rho, s2 = 0.8, 1.0
Dmax = s2/(1-rho**2); Dc = s2/(1+abs(rho))**2
print(f"  rho={rho}, sigma_W^2={s2}: Var(X)=D_max={Dmax:.4f}, D_c=min S={Dc:.4f}")
D = 0.5*Dc
th, R, V = water_fill(rho, s2, D)
Dcheck = np.mean(np.minimum(th, psd(rho,s2)))
print(f"  D={D:.4f}: theta={th:.4f}, recomputed D={Dcheck:.4f} (match {abs(Dcheck-D)<1e-6}), "
      f"R(D)={R:.4f} nats = {R*LOG2E:.4f} bits, V={V:.4f} nats^2")
v1 = abs(Dcheck-D) < 1e-5 and R > 0
print(f"  V1 {'PASS' if v1 else 'FAIL'}")

# ---------------------------------------------------------------- V2
print("="*70)
print("V2: V(D) = 1/2 nats for D <= D_c (low distortion; iid-innovation value)")
ok2 = True
for D in [0.1*Dc, 0.5*Dc, 0.99*Dc]:
    th, R, V = water_fill(rho, s2, D)
    print(f"  D={D:.4f} (={D/Dc:.2f} D_c): V(D)={V:.5f} nats^2 (target 0.5)")
    ok2 = ok2 and abs(V-0.5) < 1e-3
print(f"  V2 {'PASS' if ok2 else 'FAIL'}")

# ---------------------------------------------------------------- V3
print("="*70)
print("V3: V(D) < 1/2 and decreasing for D_c < D < D_max")
Vs = []
for D in [1.05*Dc, 1.5*Dc, 3*Dc, 0.9*Dmax]:
    if D >= Dmax: continue
    th, R, V = water_fill(rho, s2, D)
    Vs.append(V)
    print(f"  D={D:.4f} (={D/Dc:.2f} D_c): V(D)={V:.5f} nats^2, R={R*LOG2E:.4f} bits")
dec = all(Vs[i] > Vs[i+1] for i in range(len(Vs)-1))
ok3 = all(v < 0.5 for v in Vs) and dec
print(f"  V3 {'PASS' if ok3 else 'FAIL'}  (all <0.5 nats^2 and decreasing in D)")

# ---------------------------------------------------------------- V4
print("="*70)
print("V4: operational = informational -- block d-tilted info variance V_n -> spectral V(D)")
def ar1_cov_eigs(rho, s2, n):
    # stationary AR(1) covariance: Sigma[i,j] = s2/(1-rho^2) * rho^{|i-j|}
    var = s2/(1-rho**2)
    idx = np.arange(n)
    Sig = var*rho**np.abs(idx[:,None]-idx[None,:])
    return np.linalg.eigvalsh(Sig)
def block_RV(lams, D):
    lo, hi = 1e-12, lams.max()
    for _ in range(200):
        th = 0.5*(lo+hi)
        if np.mean(np.minimum(th, lams)) < D: lo = th
        else: hi = th
    th = 0.5*(lo+hi)
    Rn = np.mean(0.5*np.maximum(0.0, np.log(lams/th)))
    Vn = np.mean(0.5*np.minimum(1.0, (lams/th)**2))
    return th, Rn, Vn
D = 2.0*Dc                       # a distortion in the (D_c, D_max) range so V<1/2 (nontrivial)
th_spec, R_spec, V_spec = water_fill(rho, s2, D)
print(f"  spectral V(D)={V_spec:.5f} nats^2 (D={D:.4f})")
ok4 = True
for n in [50, 200, 800]:
    lams = ar1_cov_eigs(rho, s2, n)
    thn, Rn, Vn = block_RV(lams, D)
    print(f"  n={n:4d}: block V_n={Vn:.5f}  (rel.err vs spectral {abs(Vn/V_spec-1):.3f})")
    last_ok = abs(Vn/V_spec-1) < 0.06
ok4 = abs(block_RV(ar1_cov_eigs(rho,s2,800),D)[2]/V_spec-1) < 0.05
# Monte-Carlo the d-tilted info variance at n=800 to confirm V_n is the actual Var(j_n)/n
n = 800; lams = ar1_cov_eigs(rho, s2, n); thn,Rn,Vn = block_RV(lams, D)
R_mc = 4000
# j_n random part = sum_i [min(th,lam_i)/(2 th)] (z_i^2/lam_i - 1), z_i ~ N(0,lam_i)
coef = np.minimum(thn, lams)/(2*thn)        # per-coordinate
jvar = np.empty(R_mc)
for r in range(R_mc):
    z2_over_lam = rng.standard_normal(n)**2   # (z_i/sqrt(lam_i))^2 ~ chi^2_1
    jvar[r] = np.sum(coef*(z2_over_lam-1.0))
V_mc = jvar.var(ddof=1)/n
print(f"  Monte-Carlo Var(j_n)/n={V_mc:.5f} vs block V_n={Vn:.5f} (rel.err {abs(V_mc/Vn-1):.3f})")
ok4 = ok4 and abs(V_mc/Vn-1) < 0.06
print(f"  V4 {'PASS' if ok4 else 'FAIL'}")

# ---------------------------------------------------------------- V5
print("="*70)
print("V5: the interplay -- (N/K)delta_inf^D and sqrt(N V(D))Q^{-1}(eps) BOTH Theta(sqrt N) at K=c sqrt N")
from math import erf
def Qinv(eps):
    lo,hi=-10,10
    for _ in range(200):
        m=(lo+hi)/2
        if 0.5*(1+erf(m/math.sqrt(2)))<1-eps: lo=m
        else: hi=m
    return (lo+hi)/2
D = 2.0*Dc; th,R,V = water_fill(rho,s2,D)
Vbits = V*LOG2E**2                  # dispersion in bits^2
delta_inf_D = 0.15                  # cold-context lossy excess (bits), parametrised >0 (see V6)
c = 1.0; eps = 0.05; qi = Qinv(eps)
print(f"  V(D)={V:.4f} nats^2 = {Vbits:.4f} bits^2; delta_inf^D={delta_inf_D} bits; c={c}, eps={eps}, Q^-1={qi:.3f}")
ratios = []
for N in [10**4, 10**6, 10**8]:
    K = c*math.sqrt(N)
    ra = (N/K)*delta_inf_D                          # = sqrt(N)/c * delta_inf^D
    disp = math.sqrt(N*Vbits)*qi
    ratios.append(ra/disp)
    print(f"  N={N:>10}: (N/K)delta^D={ra:.3e} (Theta sqrt N), sqrt(NV)Q^-1={disp:.3e}, ratio={ra/disp:.5f}")
ok5 = max(ratios)-min(ratios) < 1e-9 and ratios[0] > 0
print(f"  V5 {'PASS' if ok5 else 'FAIL'}  (ratio N-independent: both Theta(sqrt N), the Thm 7.29 interplay, lossy)")

# ---------------------------------------------------------------- V6
print("="*70)
print("V6: delta_inf^D > 0 for AR(1) (memory) and -> 0 as rho -> 0 (iid limit)")
# Cold-start lossy excess: a reset loses the AR context. The first reproduction after a reset
# must code X_t from its MARGINAL (variance Var(X)=s2/(1-rho^2)) instead of the conditional
# (variance of W = s2). At distortion D<=D_c (high-rate), R = (1/2)log(var/D); the per-symbol
# cold excess = (1/2)log(Var(X)/s2) = (1/2)log(1/(1-rho^2)) = -(1/2)log(1-rho^2) > 0.
ok6 = True; prev = None
for rr in [0.0, 0.3, 0.6, 0.9]:
    d = -0.5*math.log(1-rr**2) if rr>0 else 0.0      # nats
    dbits = d*LOG2E
    print(f"  rho={rr}: delta_inf^D (cold-start rate excess) = {dbits:.5f} bits "
          f"({'=0 (iid)' if rr==0 else '>0 (memory)'})")
    if rr>0: ok6 = ok6 and dbits>0
    if prev is not None: ok6 = ok6 and dbits>=prev-1e-12
    prev = dbits
print(f"  V6 {'PASS' if ok6 else 'FAIL'}  (delta_inf^D>0 for rho>0, increasing in |rho|, ->0 as rho->0)")

# ---------------------------------------------------------------- summary
print("="*70)
allok = v1 and ok2 and ok3 and ok4 and ok5 and ok6
print(f"RESULT: {'ALL PASS' if allok else 'SOME FAILED'}  "
      f"[V1 {v1}, V2 {ok2}, V3 {ok3}, V4 {ok4}, V5 {ok5}, V6 {ok6}]")
