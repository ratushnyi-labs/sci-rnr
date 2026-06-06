#!/usr/bin/env python3
"""
Faster influence probe + the KEY decomposition test.

For each schedule x and each coordinate i:
   D_i j_n  = j_n(x) - j_n(x flip i)      [local, O(1), depends on 3 bits around i]
   D_i G_n  = G_n(x) - G_n(x flip i)      [global shell-slab object]
   D_i F    = D_i j_n - D_i G_n           [F = log2 S_n]

Efron-Stein / bounded-difference energy:
   E_j  = E[ sum_i (D_i j_n)^2 ]   (expect Theta(n): additive)
   E_G  = E[ sum_i (D_i G_n)^2 ]   (expect Theta(n): tracks j_n)
   E_F  = E[ sum_i (D_i F)^2 ]     (THE question: o(n)? then Var(log2 S_n)=o(n))

Also: CORRELATION between D_i j_n and D_i G_n. If D_i G_n ~ D_i j_n (leading cancellation),
the discrepancy D_i F is the "second-order" residual. We report:
   - E_F / n   (should -> 0 if Efron-Stein closes; const if not)
   - per-coordinate mean |D_i F| profile vs distance from chain "switch" structure.

Vectorized ball: enumerate all 2^n y once, compute P_Y(y) once, and for the flip
recompute ball mass by adjusting the boundary shells (gain shell_{t+1} with y_i!=x_i,
lose shell_t with y_i=x_i) -- exact, O(2^n) per (x) total not per (x,i).
"""
import math
import numpy as np

def h2(x):
    if x<=0 or x>=1: return 0.0
    return -x*math.log2(x)-(1-x)*math.log2(1-x)
def D_c(p): return 0.5*(1.0-math.sqrt(1.0-2.0*p)/(1.0-p))
def V_lossless(p): return p*(1-p)*(math.log2((1-p)/p))**2

def fwht(a):
    a=a.astype(np.float64).copy(); h=1
    while h<len(a):
        for i in range(0,len(a),h*2):
            xx=a[i:i+h].copy(); yy=a[i+h:i+2*h].copy()
            a[i:i+h]=xx+yy; a[i+h:i+2*h]=xx-yy
        h*=2
    return a

def setup(n,p,D):
    N=1<<n
    ys=np.arange(N,dtype=np.int64)
    bits=((ys[:,None]>>np.arange(n)[::-1])&1).astype(np.int8)  # N x n
    flips=np.sum(bits[:,1:]!=bits[:,:-1],axis=1)
    PX=0.5*(p**flips)*((1-p)**(n-1-flips))
    PhatX=fwht(PX)
    wt=np.array([bin(w).count("1") for w in range(N)])
    PhatY=PhatX*((1-2*D)**(-wt.astype(np.float64)))
    PY=fwht(PhatY)/N
    return ys, bits, PY

def analyze(p,n,R,seed,Dfrac=0.9):
    D=Dfrac*D_c(p); t=int(math.floor(n*D)); nhD=n*h2(D)
    ys,ybits,PY=setup(n,p,D)
    N=1<<n
    rng=np.random.default_rng(seed)
    L=math.log2((1-p)/p)
    Ej=[]; EG=[]; EF=[]; corr=[]
    for _ in range(R):
        steps=(rng.random(n)<p).astype(np.int8); steps[0]=rng.random()<0.5
        x=np.cumsum(steps)%2
        # distance of every y to x
        dist=np.sum(ybits!=x[None,:],axis=1)   # N
        inball = dist<=t
        ballmass = PY[inball].sum()
        G0 = -math.log2(ballmass) if ballmass>0 else math.inf
        # i_n
        f0 = int(np.sum(x[1:]!=x[:-1]))
        i0 = -(math.log2(0.5)+f0*math.log2(p)+(n-1-f0)*math.log2(1-p))
        j0 = i0 - nhD; S0 = j0 - G0
        # precompute shell membership masks
        shell_t   = (dist==t)
        shell_tp1 = (dist==t+1)
        dj=np.zeros(n); dG=np.zeros(n); dF=np.zeros(n)
        for i in range(n):
            # D_i G: ball(x') = ballmass + mass(shell_{t+1} & y_i!=x_i) - mass(shell_t & y_i=x_i)
            yi = ybits[:,i]
            gain = PY[shell_tp1 & (yi!=x[i])].sum()
            loss = PY[shell_t   & (yi==x[i])].sum()
            ballmass_f = ballmass + gain - loss
            Gf = -math.log2(ballmass_f) if ballmass_f>0 else math.inf
            # D_i j: flip x_i changes flip count
            xf = x.copy(); xf[i]^=1
            ff = int(np.sum(xf[1:]!=xf[:-1]))
            i_f = -(math.log2(0.5)+ff*math.log2(p)+(n-1-ff)*math.log2(1-p))
            jf = i_f - nhD
            dj[i]=j0-jf; dG[i]=j0 - Gf - (j0-G0) if False else (G0-Gf)
            dF[i]=(j0-jf)-(G0-Gf)
        Ej.append(np.sum(dj**2)); EG.append(np.sum(dG**2)); EF.append(np.sum(dF**2))
        # correlation of the influence vectors
        if np.std(dj)>0 and np.std(dG)>0:
            corr.append(np.corrcoef(dj,dG)[0,1])
    Ej=np.array(Ej);EG=np.array(EG);EF=np.array(EF)
    return dict(D=D,t=t,
        Ej=Ej.mean(),Ej_se=Ej.std()/math.sqrt(R),
        EG=EG.mean(),EG_se=EG.std()/math.sqrt(R),
        EF=EF.mean(),EF_se=EF.std()/math.sqrt(R),
        corr=np.mean(corr) if corr else float('nan'))

if __name__=="__main__":
    print("KEY TEST: Efron-Stein influence energy E[sum_i (D_i F)^2], F=log2 S_n.")
    print("Var(log2 S_n) <~ C * E_F (geometric-mixing bounded-diff). Want E_F = o(n).\n")
    for p in (0.25,0.4):
        print(f"### p={p}  D=0.9 D_c={0.9*D_c(p):.4f}  V_lossless={V_lossless(p):.4f} ###")
        print(f"{'n':>4} {'t':>3} {'nD':>6} | {'E_j':>14} {'E_j/n':>7} | {'E_G':>14} | {'E_F':>14} {'E_F/n':>7} | corr(dj,dG)")
        for n in (12,14,16,18,20):
            R = 3000 if n<=14 else (1200 if n==16 else (500 if n==18 else 200))
            r=analyze(p,n,R,seed=999+n)
            print(f"{n:>4} {r['t']:>3} {n*r['D']:>6.2f} | "
                  f"{r['Ej']:>8.3f}+-{r['Ej_se']:.2f} {r['Ej']/n:>7.3f} | "
                  f"{r['EG']:>8.3f}+-{r['EG_se']:.2f} | "
                  f"{r['EF']:>8.4f}+-{r['EF_se']:.3f} {r['EF']/n:>7.4f} | {r['corr']:>+.4f}", flush=True)
        print()
