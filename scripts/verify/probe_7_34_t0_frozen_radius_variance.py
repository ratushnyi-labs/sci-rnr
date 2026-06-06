#!/usr/bin/env python3
"""
At t=0 the ball is the singleton {x}, so
  log2 S_n(x) = j_n(x) - G_n(x) = log2[ P_{Y*}(x)/P_X(x) ] - n h(D).
We can compute P_{Y*}(x) EXACTLY in O(n^2) via the transfer-matrix Walsh identity
  P_{Y*}(x) = 2^{-n} sum_j (1-2D)^{-j} H_j(x),   H_j(x) = [z^j] Q_x(z),
  Q_x(z) = sum_{x'} P_X(x') prod_i ((1+z) if x'_i==x_i else (1-z)).
(This is the t=0 case of the paper's formula: K_{<=0}(j)=K_0(j)=1.)

GOAL: directly measure Var(log2 S_n) at t=0 to LARGE n (mpmath, O(n^2)), to see
whether Var(log2 S_n) is o(n) or Theta(n) IN THIS REGIME. This separates:
  - if Var/n -> 0 here too: the Efron-Stein energy E_F=Theta(n) is just a SLACK
    upper bound (influences cancel in the variance) => Efron-Stein cannot close it.
  - if Var/n -> const>0 here: then AT t=0 the residual is genuinely Theta(n), and
    the o(n) only emerges for t->infty (a regime brute force can't reach) -- meaning
    the t=0 Efron-Stein reading is not even probing the right asymptotic.
Either way it adjudicates the Efron-Stein route.
"""
import math, numpy as np, mpmath as mp

def h2(x):
    if x<=0 or x>=1: return 0.0
    return -x*math.log2(x)-(1-x)*math.log2(1-x)
def D_c(p): return 0.5*(1.0-math.sqrt(1.0-2.0*p)/(1.0-p))
def V_lossless(p): return p*(1-p)*(math.log2((1-p)/p))**2

def _mul_deg1(coeffs, sign):
    m=len(coeffs); out=[None]*(m+1); out[0]=coeffs[0]
    for k in range(1,m): out[k]=coeffs[k]+sign*coeffs[k-1]
    out[m]=sign*coeffs[m-1]; return out

def H_coeffs(x,p):
    n=len(x); half=mp.mpf('0.5'); x0=int(x[0])
    v=[_mul_deg1([half],1 if x0==0 else -1),_mul_deg1([half],1 if x0==1 else -1)]
    mp_p=mp.mpf(p); mp_q=mp.mpf(1)-mp_p
    for i in range(1,n):
        xi=int(x[i]); vb,va=v[0],v[1]; Ln=len(vb)
        acc0=[mp_q*vb[k]+mp_p*va[k] for k in range(Ln)]
        acc1=[mp_p*vb[k]+mp_q*va[k] for k in range(Ln)]
        v=[_mul_deg1(acc0,1 if xi==0 else -1),_mul_deg1(acc1,1 if xi==1 else -1)]
    Q=[v[0][k]+v[1][k] for k in range(len(v[0]))]
    while len(Q)<n+1: Q.append(mp.mpf(0))
    return Q[:n+1]

def py_star_point(x,p,D,dps=80):
    """P_{Y*}(x) exact via t=0 formula, mpmath."""
    with mp.workdps(dps):
        H=H_coeffs(x,p); n=len(x)
        inv=mp.mpf(1)/(mp.mpf(1)-2*mp.mpf(D)); total=mp.mpf(0); abss=mp.mpf(0); invj=mp.mpf(1)
        for j in range(n+1):
            term=invj*H[j]; total+=term; abss+=abs(term); invj*=inv
        prob=total/(mp.mpf(2)**n)
        cancel=(abss/abs(total)) if total!=0 else mp.inf
        return prob, float(cancel)

def i_n_bits(x,p):
    n=len(x); flips=int(np.sum(x[1:]!=x[:-1]))
    logP=math.log2(0.5)+flips*math.log2(p)+(n-1-flips)*math.log2(1-p)
    return -logP

def sample_bsms(n,p,rng):
    steps=(rng.random(n)<p).astype(np.int8); steps[0]=rng.random()<0.5
    return np.cumsum(steps)%2

def measure_t0(p,n,R,seed,Dfrac=0.9,dps=80,keepdig=15):
    D=Dfrac*D_c(p); nhD=n*h2(D)
    rng=np.random.default_rng(seed)
    lS=[]; jn=[]; dropped=0; wc=0.0
    for _ in range(R):
        x=sample_bsms(n,p,rng)
        prob,cancel=py_star_point(x,p,D,dps)
        wc=max(wc,cancel)
        lost=math.log10(cancel) if cancel>1 else 0.0
        if (dps-lost)<keepdig or prob<=0:
            dropped+=1; continue
        G=float(-mp.log(prob,2))                 # G_n = -log2 P_{Y*}(x)
        i_n=i_n_bits(x,p); j=i_n-nhD; S=j-G
        lS.append(S); jn.append(j)
    lS=np.array(lS); jn=np.array(jn); kept=len(lS)
    if kept<30: return None
    return dict(n=n,kept=kept,dropped=dropped,wc=wc,
        VlS=lS.var(ddof=1), VlS_se=lS.var(ddof=1)*math.sqrt(2/(kept-1)),
        Vj=jn.var(ddof=1),
        meanlS=lS.mean())

if __name__=="__main__":
    print("DIRECT Var(log2 S_n) at t=0 (singleton ball), exact O(n^2), to large n.")
    print("Adjudicates whether Efron-Stein E_F=Theta(n) is slack (true Var=o(n)) or tight.\n")
    for p in (0.25,0.4):
        D=0.9*D_c(p)
        print(f"### p={p}  D={D:.4f}  V_lossless={V_lossless(p):.4f}  (t=0 forced: nD<1) ###")
        print(f"{'n':>5} {'kept':>5} {'drop':>4} | {'Var(j)/n':>10} | {'Var(log2S)':>12} {'Var(log2S)/n':>14} | {'wc':>8}")
        rows=[]
        for n in (16,24,32,48,64,96,128,160,200,256,320):
            R = 4000 if n<=64 else (2000 if n<=128 else 1000)
            r=measure_t0(p,n,R,seed=4242+n)
            if r is None:
                print(f"{n:>5}  -- insufficient (precision) --"); continue
            rows.append(r)
            print(f"{n:>5} {r['kept']:>5} {r['dropped']:>4} | "
                  f"{r['Vj']/n:>10.4f} | {r['VlS']:>12.5f} {r['VlS']/n:>14.6f} | {r['wc']:>8.1e}", flush=True)
        if len(rows)>=4:
            ns=np.array([r['n'] for r in rows]); y=np.array([r['VlS']/r['n'] for r in rows])
            se=np.array([r['VlS_se']/r['n'] for r in rows])
            invn=1.0/ns; A=np.vstack([invn,np.ones_like(invn)]).T
            w=1/se**2; Aw=A*np.sqrt(w)[:,None]; yw=y*np.sqrt(w)
            coef,*_=np.linalg.lstsq(Aw,yw,rcond=None); a,b=coef
            cov=np.linalg.inv(Aw.T@Aw); seb=math.sqrt(cov[1,1])
            print(f"  fit Var(log2S)/n ~ a/n + b :  a={a:+.4f}  b(intercept)={b:+.6f} +- {seb:.6f}")
            # also fit raw Var(log2S) ~ c + d*n to see if it's bounded or linear
            yr=np.array([r['VlS'] for r in rows])
            ser=np.array([r['VlS_se'] for r in rows]); wr=1/ser**2
            Ar=np.vstack([ns,np.ones_like(ns)]).T; Arw=Ar*np.sqrt(wr)[:,None]; yrw=yr*np.sqrt(wr)
            cr,*_=np.linalg.lstsq(Arw,yrw,rcond=None)
            print(f"  fit Var(log2S) ~ d*n + c :  d(slope)={cr[0]:+.6f}  c={cr[1]:+.4f}")
            print(f"  => raw Var(log2S) {'GROWS ~ linearly (slope>0)' if cr[0]>3e-3 else 'is ~BOUNDED (slope~0)'}")
        print()
