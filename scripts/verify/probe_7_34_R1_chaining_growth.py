#!/usr/bin/env python3
"""
probe_7_34_R1_chaining_growth.py
============================================================================
TRACK (R1), DECISIVE test: does the CONTINUITY/CHAINING proof technique bridge
the ANNEALED boundedness of r_t to the QUENCHED (shell-uniform) boundedness (R1)?

For two near-mean shell schedules x, x' we compare:
  TRUE gap     = |log2 r_t(x) - log2 r_t(x')|     (the actual quantity; annealed
                 concentration predicts O(1) since both r_t ~ 0.94)
  chain bound  = sum_{i in diff} |single-flip increment|  along a path x->x'
                 (the triangle-inequality bound a continuity/McDiarmid/chaining
                 argument can deliver)

FINDING (p=0.4, n=40,60,80): TRUE gap mean ~0.10-0.12, max ~0.3-0.6 = O(1)
[boundedness HOLDS], but mean Hamming distance between shell words is Theta(n)
(~n/2) and the chaining bound grows like Theta(sqrt n): chain/sqrt(n) ~ 0.19-0.20
CONSTANT across n. => the triangle-inequality continuity argument over the shell
yields only O(sqrt n), i.e. Var = O(n), NOT the uniform O(1) that (R1) requires.
The single-flip increments are non-localized (spread over Theta(n) coords) and do
not telescope, so chaining over-counts by sqrt(n). CONCLUSION: route (a)
continuity/chaining provably does NOT upgrade annealed -> quenched; the quenched
(R1) needs a genuinely new (sequential/inhomogeneous, density-floor-free) lattice
LLT (route (b)), NOT obtainable from annealed + (R2) + continuity.

Run: /tmp/rnr_venv/bin/python scripts/verify/probe_7_34_R1_chaining_growth.py
"""
# DECISIVE: chaining bound over the shell. Two near-mean shell schedules differ in Theta(n)
# coords. The TRIANGLE-INEQUALITY chaining bound sum_{i in diff}|incr| -- does it grow like
# sqrt(n) (=> continuity FAILS to give O(1)) or stay O(1) (=> continuity closes it)?
# We measure: (a) typical Hamming distance between two shell words, (b) the actual
# |log r_t(x)-log r_t(x')| (the TRUE gap, which annealed concentration says is O(1)),
# (c) the chaining UPPER BOUND = sum of single-flip increments along a path x->x'.
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
    n=len(x); H=H_coeffs(x,p); inv=mp.mpf(1)/(1-2*mp.mpf(D)); ip=[mp.mpf(1)]*(n+1)
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
p=0.4; D=0.9*D_c(p)
print(f"p={p} D={D:.5f}: TRUE gap |lr(x)-lr(x')| vs CHAINING bound (path of single flips), shell words")
print(f"{'n':>4} | {'meanHamm':>9} | {'TRUE gap mean':>13} {'TRUE gap max':>12} | {'chain bound mean':>16} {'chain/sqrt(n)':>13}")
for n in (40,60,80):
    t=int(math.floor(n*D))
    Km=kraw(n); rng=np.random.default_rng(7000+n)
    R=60 if n<=60 else 35
    def nm(x):
        b=lr(x,p,D,Km); 
        return b if (b is not None and abs(b[1]-n*D)<1.0+0.1*n*D) else None
    hamms=[]; truegaps=[]; chains=[]
    cnt=0; tries=0
    while cnt<R and tries<R*40:
        tries+=1
        x=samp(n,p,rng); bx=nm(x)
        if bx is None: continue
        xp=samp(n,p,rng); bxp=nm(xp)
        if bxp is None: continue
        diff=np.where(x!=xp)[0]
        if len(diff)<1: continue
        # chaining path: flip diff coords one at a time, sum |single-flip increments|
        cur=x.copy(); lcur=bx[0]; cb=0.0; ok=True
        for i in diff:
            cur[i]^=1; bc=lr(cur,p,D,Km)
            if bc is None: ok=False; break
            cb+=abs(bc[0]-lcur); lcur=bc[0]
        if not ok: continue
        cnt+=1; hamms.append(len(diff)); truegaps.append(abs(bx[0]-bxp[0])); chains.append(cb)
    if cnt<5:
        print(f"{n:>4} | few"); continue
    h=np.array(hamms); tg=np.array(truegaps); ch=np.array(chains)
    print(f"{n:>4} | {h.mean():>9.1f} | {tg.mean():>13.4f} {tg.max():>12.4f} | {ch.mean():>16.4f} {ch.mean()/math.sqrt(n):>13.4f}", flush=True)
print()
print("KEY: TRUE gap is O(1) (annealed concentration). If chaining bound GROWS like sqrt(n)")
print("(chain/sqrt(n) ~ const), the triangle-inequality continuity argument over the shell")
print("gives only O(sqrt n), NOT O(1): continuity does NOT bridge annealed->quenched.")
