#!/usr/bin/env python3
"""
Verification for Theorem 7.36 (Functional CLT for the RNR repair process; the
streaming-buffer high-water-mark).

Setup: repair lengths ell_i = -log2 P(x_i | context). Centered partial sum
D_t = sum_{i<=t} (ell_i - h) is the occupancy of a streaming encoder buffer that
emits ell_i bits per symbol while the channel drains h bits/symbol. h = entropy rate,
V = varentropy rate (long-run Var of the information density; same V as Thm 7.28/7.29/7.35).

Claim (functional CLT / Donsker for exp-mixing sequences, Merlevede-Peligrad-Rio):
    W_N(t) := D_{floor(Nt)} / sqrt(N V)  =>  standard Brownian motion W(t) on D[0,1].
Consequences (continuous mapping):
    (1) max_t D_t / sqrt(NV)  =>  sup_{[0,1]} W  =_d  |N(0,1)|   (half-normal; reflection),
        E = sqrt(2/pi) ~ 0.79788  -- the NET-SURPLUS high-water-mark (signed buffer).
    (2) max_t |D_t| / sqrt(NV) =>  sup_{[0,1]} |W|  -- two-sided buffer range (over+under).
    (3) D_N / sqrt(NV) = W_N(1) => N(0,1)  -- recovers the Thm 7.28 marginal CLT.

Checks:
  V1  finite-dim distributions of W_N: mean~0, Var(W_N(t))~t, Cov(W_N(s),W_N(t))~min(s,t).
  V2  buffer high-water-mark max_t D_t / sqrt(NV) -> half-normal: mean sqrt(2/pi),
      and tail P(>a)=2(1-Phi(a)) at a quantile.
  V3  two-sided max |D_t|/sqrt(NV) -> sup|W| (reference mean by fine BM simulation).
  V4  Markov source: same half-normal limit with the long-run (covariance) V.
  V5  distributional convergence (not just mean): KS-style quantile match of the
      rescaled high-water-mark to half-normal across growing N.
  V6  same-V unification: V here = varentropy = the Thm 7.29 dispersion constant; the
      three buffers (archive dispersion, RA read buffer 7.35, streaming buffer 7.36) share V.
"""
import math
import numpy as np

rng = np.random.default_rng(20260605)

def H_iid(p):
    p = np.asarray(p); return float(-(p*np.log2(p)).sum())
def V_iid(p):
    p = np.asarray(p); h = H_iid(p); return float((p*(np.log2(p))**2).sum() - h*h)

def Phi(x):  # standard normal CDF
    return 0.5*(1+math.erf(x/math.sqrt(2)))

HALF_NORMAL_MEAN = math.sqrt(2/math.pi)   # E|N(0,1)| ~ 0.79788

def iid_paths_maxes(p, N, R):
    """R realizations of length N; return centered partial-sum running maxes and |.|."""
    p = np.asarray(p); h = H_iid(p)
    surpr = -np.log2(p); cump = np.cumsum(p); cump[-1] = 1.0
    maxpos = np.empty(R); maxabs = np.empty(R); endval = np.empty(R)
    # also collect W_N at t=0.3,0.5,0.7 for V1 (indices)
    i3, i5, i7 = int(0.3*N), int(0.5*N), int(0.7*N)
    w3 = np.empty(R); w5 = np.empty(R); w7 = np.empty(R)
    for r in range(R):
        u = rng.random(N)
        idx = np.clip(np.searchsorted(cump, u, side='right'), 0, len(p)-1)
        D = np.cumsum(surpr[idx] - h)         # centered partial sums D_t
        maxpos[r] = D.max(); maxabs[r] = np.abs(D).max(); endval[r] = D[-1]
        w3[r] = D[i3-1]; w5[r] = D[i5-1]; w7[r] = D[i7-1]
    return maxpos, maxabs, endval, (w3, w5, w7)

# reference: sup and sup|.| of standard BM on [0,1] by fine simulation (vectorized, batched)
def bm_sup_refs(steps=2000, R=20000, batch=4000):
    dt = 1.0/steps
    sup1s = []; supabss = []
    for b0 in range(0, R, batch):
        nb = min(batch, R-b0)
        incr = rng.normal(0, math.sqrt(dt), size=(nb, steps))
        w = np.cumsum(incr, axis=1)
        sup1s.append(np.maximum(0.0, w.max(axis=1)))
        supabss.append(np.abs(w).max(axis=1))
    return float(np.concatenate(sup1s).mean()), float(np.concatenate(supabss).mean())

# ---------------------------------------------------------------- V1
print("="*70)
print("V1: finite-dim distributions of W_N(t) ~ Brownian motion")
p = [0.5, 0.25, 0.15, 0.10]
h = H_iid(p); V = V_iid(p)
N = 3000; R = 12000
maxpos, maxabs, endval, (w3, w5, w7) = iid_paths_maxes(p, N, R)
s = math.sqrt(N*V)
W3 = w3/s; W5 = w5/s; W7 = w7/s; W1 = endval/s
print(f"  h={h:.5f}, V={V:.5f}, N={N}")
print(f"  Var W_N(0.3)={W3.var():.4f} (exp 0.3); Var W_N(0.5)={W5.var():.4f} (exp 0.5); "
      f"Var W_N(1)={W1.var():.4f} (exp 1.0)")
cov37 = float(np.cov(W3, W7)[0,1])
print(f"  mean W_N(0.5)={W5.mean():+.4f} (exp 0); Cov(W_N(0.3),W_N(0.7))={cov37:.4f} (exp 0.3=min)")
v1 = (abs(W3.var()-0.3)<0.02 and abs(W5.var()-0.5)<0.03 and abs(W1.var()-1.0)<0.05
      and abs(W5.mean())<0.03 and abs(cov37-0.3)<0.03)
print(f"  V1 {'PASS' if v1 else 'FAIL'}")

# ---------------------------------------------------------------- V2
print("="*70)
print("V2: buffer high-water-mark max_t D_t / sqrt(NV) -> half-normal |N(0,1)|")
M = maxpos/s
emp_mean = M.mean()
print(f"  E[max_t D_t/sqrt(NV)]={emp_mean:.4f}  vs half-normal mean sqrt(2/pi)={HALF_NORMAL_MEAN:.4f}")
# tail check: P(M>a)=2(1-Phi(a)); test at a=1.0
a = 1.0; emp_tail = float((M>a).mean()); pred_tail = 2*(1-Phi(a))
print(f"  P(M>{a})={emp_tail:.4f}  vs 2(1-Phi({a}))={pred_tail:.4f}")
v2 = abs(emp_mean-HALF_NORMAL_MEAN)<0.02 and abs(emp_tail-pred_tail)<0.03
print(f"  V2 {'PASS' if v2 else 'FAIL'}")

# ---------------------------------------------------------------- V3
print("="*70)
print("V3: two-sided max |D_t|/sqrt(NV) -> sup_{[0,1]}|W| (reference by fine BM sim)")
sup1_ref, supabs_ref = bm_sup_refs()
Mabs = maxabs/s
print(f"  E[max|D_t|/sqrt(NV)]={Mabs.mean():.4f}  vs BM-sim E[sup|W|]={supabs_ref:.4f}")
print(f"  (also one-sided E[sup W^+] BM-sim={sup1_ref:.4f} vs sqrt(2/pi)={HALF_NORMAL_MEAN:.4f})")
v3 = abs(Mabs.mean()-supabs_ref)<0.03 and abs(sup1_ref-HALF_NORMAL_MEAN)<0.02
print(f"  V3 {'PASS' if v3 else 'FAIL'}")

# ---------------------------------------------------------------- V4  Markov
print("="*70)
print("V4: Markov source -- same half-normal high-water-mark with long-run V")

def markov_stationary(P):
    P = np.asarray(P); n=len(P); pi=np.full(n,1.0/n)
    for _ in range(5000):
        nx = pi@P
        if np.max(np.abs(nx-pi))<1e-15: pi=nx; break
        pi=nx
    return pi

def markov_h_V(P, pi):
    P=np.asarray(P); pi=np.asarray(pi); n=len(P)
    with np.errstate(divide='ignore'):
        g = np.where(P>0,-np.log2(P),0.0)
    f = (P*g).sum(axis=1); h=float(pi@f)
    A = np.eye(n)-P+np.outer(np.ones(n),pi); Z=np.linalg.inv(A)
    d=g-h; fc=f-h
    term=float(pi@(P*d**2).sum(axis=1)); u=Z@fc
    cross=float(pi@(P*(d*u[None,:])).sum(axis=1))*2.0
    return h, term+cross

def markov_buffer_maxes(P, pi, N, R):
    P=np.asarray(P); n=len(P)
    with np.errstate(divide='ignore'):
        g=np.where(P>0,-np.log2(P),0.0)
    cdf=np.cumsum(P,axis=1)
    hM,_=markov_h_V(P,pi)
    cur=rng.choice(n,size=R,p=pi)
    D=np.zeros(R); runmax=np.full(R,-1e18)
    for _ in range(N):
        u=rng.random(R)
        nxt=np.minimum((u[:,None]>cdf[cur]).sum(axis=1),n-1)
        D += g[cur,nxt]-hM
        runmax=np.maximum(runmax,D)
        cur=nxt
    return runmax

P=[[0.7,0.2,0.1],[0.1,0.6,0.3],[0.2,0.3,0.5]]
pi=markov_stationary(P); hM,VM=markov_h_V(P,pi)
N=3000; R=8000
rm = markov_buffer_maxes(P,pi,N,R)
Mm = rm/math.sqrt(N*VM)
print(f"  Markov h={hM:.5f}, V(long-run)={VM:.5f}")
print(f"  E[max_t D_t/sqrt(NV)]={Mm.mean():.4f}  vs sqrt(2/pi)={HALF_NORMAL_MEAN:.4f}")
v4 = abs(Mm.mean()-HALF_NORMAL_MEAN)<0.03
print(f"  V4 {'PASS' if v4 else 'FAIL'}")

# ---------------------------------------------------------------- V5
print("="*70)
print("V5: distributional convergence to half-normal (quantile match, growing N)")
# half-normal quantiles: F^{-1}(q) = Phi^{-1}((1+q)/2)
def halfnormal_quantile(q):
    from math import sqrt
    # invert Phi via bisection
    target=(1+q)/2; lo,hi=0.0,10.0
    for _ in range(100):
        mid=(lo+hi)/2
        if Phi(mid)<target: lo=mid
        else: hi=mid
    return (lo+hi)/2
qs=[0.5,0.9,0.99]
ok5=True
for Ntest in [1000, 6000]:
    mp,_,_,_ = iid_paths_maxes(p, Ntest, 8000)
    Mt = mp/math.sqrt(Ntest*V)
    line=f"  N={Ntest:5d}: "
    for q in qs:
        emp=np.quantile(Mt,q); pred=halfnormal_quantile(q)
        line+=f"q{q}={emp:.3f}/{pred:.3f} "
        if abs(emp-pred)>0.08: ok5=False
    print(line)
print(f"  V5 {'PASS' if ok5 else 'FAIL'}  (empirical/half-normal quantiles match, sharper at large N)")

# ---------------------------------------------------------------- V6
print("="*70)
print("V6: same-V unification -- V here = varentropy (7.29 dispersion = 7.35 EVT constant)")
# iid: V = Var(surprisal); confirm the buffer rescaling uses exactly this V.
emp_var_surpr = None
u=rng.random(200000); idx=np.clip(np.searchsorted(np.cumsum(p),u,side='right'),0,len(p)-1)
sv=(-np.log2(np.asarray(p)))[idx]
emp_var_surpr=float(sv.var(ddof=1))
print(f"  iid: V(analytic)={V:.5f}  empirical Var(surprisal)={emp_var_surpr:.5f}  ratio={emp_var_surpr/V:.4f}")
print(f"  Markov: V(long-run, fundamental matrix)={VM:.5f} (same constant feeding 7.29/7.35/7.36)")
v6 = abs(emp_var_surpr/V-1)<0.03
print(f"  V6 {'PASS' if v6 else 'FAIL'}")

# ---------------------------------------------------------------- V7  sub-block drift
print("="*70)
print("V7: SUB-BLOCK cold-start resets -> drifted buffer sup(W+beta s), NOT half-normal")
# (codex round-1 fix: parts (a),(b) are the IDEAL stream; the practical sub-block coder resets
# context every K, paying the cold-start excess delta_inf per block = a deterministic drift.)
# 2-state SYMMETRIC Markov: marginal pi=(0.5,0.5), so cold-start (first symbol of each block
# coded from pi instead of the conditional) pays exactly delta_inf = H(pi)-h = 1-h per block.
p_stay = 0.9
Pm = np.array([[p_stay, 1-p_stay],[1-p_stay, p_stay]])
pim = markov_stationary(Pm)            # (0.5,0.5)
hm, Vm = markov_h_V(Pm, pim)
H_marg = -float((pim*np.log2(pim)).sum())
delta_inf = H_marg - hm                # cold-start excess per block (order-1: one symbol)
print(f"  2-state p_stay={p_stay}: h={hm:.4f}, V={Vm:.4f}, H(pi)={H_marg:.4f}, "
      f"delta_inf=H(pi)-h={delta_inf:.4f}")

def subblock_buffer_max(P, pi, K, N, R):
    """Sub-block coder: first symbol of each length-K block coded from the marginal pi
    (cold), the rest from the conditional. Returns max_t D_t (D_t=sum(ell_i - h))."""
    P=np.asarray(P); n=len(P)
    with np.errstate(divide='ignore'):
        g=np.where(P>0,-np.log2(P),0.0)
    smarg = -np.log2(pi)               # cold marginal surprisal
    cdf=np.cumsum(P,axis=1); hh=markov_h_V(P,pi)[0]
    cur=rng.choice(n,size=R,p=pi)
    D=np.zeros(R); runmax=np.full(R,-1e18)
    for t in range(N):
        u=rng.random(R)
        nxt=np.minimum((u[:,None]>cdf[cur]).sum(axis=1),n-1)
        if t % K == 0:                 # cold start of a block: code nxt from marginal
            ell = smarg[nxt]
        else:
            ell = g[cur,nxt]
        D += ell - hh
        runmax=np.maximum(runmax,D)
        cur=nxt
    return runmax

# reference: E[sup_{[0,1]} (W(s)+beta s)] by fine Brownian-with-drift simulation
def sup_drift_ref(beta, steps=4000, R=30000, batch=6000):
    dt=1.0/steps; ts=(np.arange(1,steps+1))*dt; sups=[]
    for b in range(0,R,batch):
        nb=min(batch,R-b)
        w=np.cumsum(rng.normal(0,math.sqrt(dt),size=(nb,steps)),axis=1)+beta*ts[None,:]
        sups.append(np.maximum(0.0, w.max(axis=1)))
    return float(np.concatenate(sups).mean())

c = 1.0
N = 10000; K = int(round(c*math.sqrt(N))); R = 6000
beta = delta_inf/(c*math.sqrt(Vm))
rm = subblock_buffer_max(Pm, pim, K, N, R)
emp = (rm/math.sqrt(N*Vm)).mean()
ref_drift = sup_drift_ref(beta)
print(f"  K=c*sqrt(N)={K}, beta=delta_inf/(c sqrt(V))={beta:.4f}")
print(f"  E[max_t D_t/sqrt(NV)] (sub-block)={emp:.4f}")
print(f"  vs drifted-sup E[sup(W+beta s)]={ref_drift:.4f}  (half-normal sqrt(2/pi)={HALF_NORMAL_MEAN:.4f})")
# the sub-block high-water-mark matches the DRIFTED sup, and is clearly ABOVE the half-normal
v7 = abs(emp-ref_drift) < 0.05 and emp > HALF_NORMAL_MEAN + 0.15
print(f"  V7 {'PASS' if v7 else 'FAIL'}  (drifted sup, not half-normal: drift beta matters)")

# ---------------------------------------------------------------- V8  physical Lindley queue
print("="*70)
print("V8: PHYSICAL empty-queue high-water-mark (Lindley draw-up) -> sup|W|, NOT half-normal")
# (codex round-2 fix: the SIGNED net-surplus max_t D_t -> half-normal (V2); a PHYSICAL queue
# with no negative occupancy has high-water-mark = max draw-up max_{j<=t}(D_t-D_j) = Lindley
# max_t Q_t, which by Levy's theorem =_d sup_{[0,1]}|W| (theta-function, mean sqrt(pi/2)~1.2533),
# distinct from the half-normal.)
SUP_ABS_MEAN = math.sqrt(math.pi/2)            # E[sup_{[0,1]}|W|] = sqrt(pi/2)
N = 4000; R = 12000
hq = H_iid(p); Vq = V_iid(p)
surpr = -np.log2(np.asarray(p)); cump = np.cumsum(np.asarray(p)); cump[-1] = 1.0
qmax = np.empty(R); dmax = np.empty(R)
for r in range(R):
    u = rng.random(N)
    idx = np.clip(np.searchsorted(cump, u, side='right'), 0, len(p)-1)
    xi = surpr[idx] - hq                        # centered increments (mean 0, ideal drain h)
    D = np.cumsum(xi)
    # Lindley empty-start queue: Q_t = max(0, Q_{t-1}+xi_t); max_t Q_t = max draw-up
    Q = 0.0; qm = 0.0
    # vectorized running cumulative-min draw-up: max_t (D_t - min_{s<=t} D_s)
    runmin = np.minimum.accumulate(np.concatenate(([0.0], D)))[1:]   # min_{s<=t} D_s (incl 0)
    drawup = D - runmin
    qmax[r] = drawup.max()
    dmax[r] = D.max()
emp_q = (qmax/math.sqrt(N*Vq)).mean()
emp_d = (dmax/math.sqrt(N*Vq)).mean()
print(f"  physical-queue (draw-up) E[max Q/sqrt(NV)]={emp_q:.4f}  vs sup|W| mean sqrt(pi/2)={SUP_ABS_MEAN:.4f}")
print(f"  net-surplus     E[max D/sqrt(NV)]={emp_d:.4f}  vs half-normal sqrt(2/pi)={HALF_NORMAL_MEAN:.4f}")
# draw-up -> sup|W| (theta, ~1.25, from below at finite N); clearly ABOVE the half-normal net-surplus
v8 = (emp_q > emp_d + 0.15) and abs(emp_q - SUP_ABS_MEAN) < 0.06 and abs(emp_d - HALF_NORMAL_MEAN) < 0.03
print(f"  V8 {'PASS' if v8 else 'FAIL'}  (physical queue = sup|W| theta-law, NOT the net-surplus half-normal)")

# ---------------------------------------------------------------- summary
print("="*70)
allok = v1 and v2 and v3 and v4 and ok5 and v6 and v7 and v8
print(f"RESULT: {'ALL PASS' if allok else 'SOME FAILED'}  "
      f"[V1 {v1}, V2 {v2}, V3 {v3}, V4 {v4}, V5 {ok5}, V6 {v6}, V7 {v7}, V8 {v8}]")
