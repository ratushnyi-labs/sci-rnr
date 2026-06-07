#!/usr/bin/env python3
"""
probe_7_34_R1_avg_vs_worst_influence.py
============================================================================
TRACK (R1), discriminator A: AVERAGED vs WORST-CASE per-coordinate single-flip
influence of log2 r_t, in the GENUINE NEAR-MEAN regime (t=floor(nD) large).

  avg sumE_i      = E_x[ sum_i (D_i log2 r_t)^2 ]   (averaged Dirichlet energy)
  max|D_i|        = sup over words & coords of |D_i log2 r_t|
  sum(worst_i)^2  = sum_i (max_x |D_i|)^2           (McDiarmid worst-case budget)
  #big            = avg #coords per word with |D_i|>0.3

FINDING (p=0.4): the averaged energy is O(1) and DECREASING (0.49->0.20 as
n:60->100) -- the ANNEALED Efron-Stein/Glauber-Poincare object. The Theta(1)
per-coordinate jumps are CONFINED to the atypical frozen-radius (t small, n<=60)
words (#big ~0.7); once t enters the near-mean zone (n>=80, t>=9) #big=0.00 and
max|D_i| drops to ~0.24, consistent with every single-flip increment being
Theta(n^{-1/2}). So the worst-case single-flip influence is NOT Theta(1) in the
genuine regime -- but it is SPREAD over Theta(n) coordinates (participation ratio
~0.5; see probe_7_34_R1_chaining_fast.py), so the chaining BUDGET over the shell
still grows (see probe_7_34_R1_chaining_growth.py for the decisive sqrt(n) growth).

Run: /tmp/rnr_venv/bin/python scripts/verify/probe_7_34_R1_avg_vs_worst_influence.py [p]
"""
import math, numpy as np, mpmath as mp
from math import comb
def D_c(p): return 0.5*(1.0-math.sqrt(1.0-2.0*p)/(1.0-p))
def _mul(c,s):
    m=len(c); o=[None]*(m+1); o[0]=c[0]
    for k in range(1,m): o[k]=c[k]+s*c[k-1]
    o[m]=s*c[m-1]; return o
def H_coeffs(x,p):
    n=len(x); half=mp.mpf('0.5'); x0=int(x[0])
    v=[_mul([half],1 if x0==0 else -1),_mul([half],1 if x0==1 else -1)]
    P=mp.mpf(p); Qq=1-P
    for i in range(1,n):
        xi=int(x[i]); vb,va=v[0],v[1]; L=len(vb)
        a0=[Qq*vb[k]+P*va[k] for k in range(L)]; a1=[P*vb[k]+Qq*va[k] for k in range(L)]
        v=[_mul(a0,1 if xi==0 else -1),_mul(a1,1 if xi==1 else -1)]
    Q=[v[0][k]+v[1][k] for k in range(len(v[0]))]
    while len(Q)<n+1: Q.append(mp.mpf(0))
    return Q[:n+1]
_KC={}
def kraw(n):
    if n in _KC: return _KC[n]
    K=[[0]*(n+1) for _ in range(n+1)]
    for j in range(n+1):
        for k in range(n+1):
            s=0; lo=max(0,k-(n-j)); hi=min(k,j)
            for l in range(lo,hi+1): s+=((-1)**l)*comb(j,l)*comb(n-j,k-l)
            K[k][j]=s
    _KC[n]=K; return K
def pik(x,p,D,Km):
    n=len(x); H=H_coeffs(x,p); inv=mp.mpf(1)/(1-2*mp.mpf(D))
    ip=[mp.mpf(1)]*(n+1)
    for j in range(1,n+1): ip[j]=ip[j-1]*inv
    out=[mp.mpf(0)]*(n+1)
    for k in range(n+1):
        s=mp.mpf(0)
        for j in range(n+1): s+=ip[j]*Km[k][j]*H[j]
        out[k]=s/(mp.mpf(2)**n)
    return out
def lr(x,p,D,Km,dps=50):
    n=len(x); t=int(math.floor(n*D)); th0=float(mp.log(mp.mpf(D)/(1-mp.mpf(D))))
    with mp.workdps(dps):
        pis=pik(x,p,D,Km); M=sum(pis[k]*mp.e**(mp.mpf(th0)*k) for k in range(n+1))
        if M<=0: return None
        q=[float(pis[k]*mp.e**(mp.mpf(th0)*k)/M) for k in range(n+1)]
    if q[t]<=0: return None
    mq=sum(k*q[k] for k in range(n+1)); vq=sum((k-mq)**2*q[k] for k in range(n+1))
    if vq<=0: return None
    rt=q[t]*math.sqrt(2*math.pi*vq)
    if rt<=0: return None
    return math.log2(rt),mq
def samp(n,p,rng):
    st=(rng.random(n)<p).astype(np.int8); st[0]=(rng.random()<0.5); return np.cumsum(st)%2
# decisive: AVERAGED per-coord energy E[D_i^2] vs WORST per-coord (max_x |D_i|), summed
p=0.4; D=0.9*D_c(p)
print(f"p={p} D={D:.5f}  AVERAGED vs WORST-CASE per-coordinate flip energy of log2 r_t")
print(f"{'n':>4} {'t':>3}{'kept':>5} | {'avg sumE_i':>11} {'max|D_i| overall':>16} {'sum(worst_i)^2':>15} {'#big(|D_i|>0.3)/word':>20}")
for n in (40,60,80,100):
    t=int(math.floor(n*D))
    if t<4: continue
    Km=kraw(n); rng=np.random.default_rng(3000+n)
    R=120 if n<=60 else (70 if n==80 else 40)
    avgE=np.zeros(n); worst=np.zeros(n); kept=0; gmax=0.0; bigcnt=[]
    while kept<R:
        x=samp(n,p,rng); b=lr(x,p,D,Km)
        if b is None: continue
        l0,mq=b
        if abs(mq-n*D)>1.0+0.1*n*D: continue
        Di=np.zeros(n); ok=True
        for i in range(n):
            xf=x.copy(); xf[i]^=1; bf=lr(xf,p,D,Km)
            if bf is None: ok=False; break
            Di[i]=l0-bf[0]
        if not ok: continue
        kept+=1; a=np.abs(Di); avgE+=Di**2
        worst=np.maximum(worst,a); gmax=max(gmax,a.max())
        bigcnt.append(int((a>0.3).sum()))
    avgE/=kept
    print(f"{n:>4} {t:>3}{kept:>5} | {avgE.sum():>11.4f} {gmax:>16.3f} {float((worst**2).sum()):>15.3f} {np.mean(bigcnt):>20.2f}", flush=True)
