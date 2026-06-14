#!/usr/bin/env python3
"""
probe_7_34p_lossy_devhier.py
============================================================================
Remark 7.34p -- the LOSSY deviation hierarchy on the Gray region: §7.34 (x) §7.29-7.36.

The all-order dual identity j_m = i_m - m c (Remark 7.34n') holds for EVERY block
length m (the dual identity is per-length; the per-symbol constant
c = lam* D log2(e) + log2(gamma) is m-INDEPENDENT). The shift is therefore a
DETERMINISTIC LINEAR DRIFT, so the CENTERED d-tilted-information process equals
the centered surprisal process PATHWISE:
    j_m - m R(D) = i_m - m h    (every m, on the Gray region; R(D)=h-c).

Consequently the ENTIRE §7 deviation hierarchy -- proved there for the lossless
surprisal/archive process -- transfers verbatim to the lossy rate-distortion
setting on the Gray region, via i_m -> j_m = i_m - m c:
  * Functional CLT / streaming buffer (Thm 7.36): centred process pathwise
    identical => same Brownian limit, same V=V_lossless, same buffer laws.
  * Large-deviation overflow exponent (Thm 7.31): E_lossy(R) = E_lossless(R+c).
  * Moderate deviations / RA fade (Thm 7.33): same V, shifted mean.
  * Extreme-value worst-case query (Thm 7.35): max block d-tilted info = max
    surprisal - K c => Gumbel law and read-buffer M_m = K R(D) + sqrt(2KV ln(N/K))
    + ... with h -> R(D).

CHECKS:
  P1  j_m = i_m - m c for each m (spread ~0); per-symbol c is m-independent.
  P2  centered process identity: {j_m - m R(D)} = {i_m - m h} per word (pathwise).
  P3  consequence sample -- the overflow/tail shift: the empirical distribution of
      j_m equals that of i_m shifted by -m c (so any tail/LDP/EVT functional of
      j_m equals that of i_m at the shifted argument).

Deps: numpy.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools, math
import numpy as np
PASS=True
def rep(name, ok):
    global PASS; PASS=PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok
def stat2(a,b): s=a+b; return np.array([b/s, a/s])
def Kinv(D): a0=1-2*D; return (1/a0)*np.array([[1-D,-D],[-D,1-D]])
def Dc(a,b): return 0.5*(1-math.sqrt(1-4*a*b/(2-a-b)**2))
def block(a,b,D,m):
    pi=stat2(a,b); T=np.array([[1-a,a],[b,1-b]]); Ki=Kinv(D)
    idx=list(itertools.product([0,1],repeat=m)); PX=np.zeros([2]*m)
    for w in idx: PX[w]=pi[w[0]]*np.prod([T[w[t-1],w[t]] for t in range(1,m)])
    Pt=PX
    for ax in range(m): Pt=np.tensordot(Ki,Pt,axes=([1],[ax])); Pt=np.moveaxis(Pt,0,ax)
    nu=D/(1-D); lam=math.log((1-D)/D); rows=[]
    for wx in idx:
        if PX[wx]<=0: continue
        M=sum(Pt[wy]*nu**sum(1 for t in range(m) if wx[t]!=wy[t]) for wy in idx)
        jm=(-lam*m*D-math.log(M))/math.log(2); im=-math.log2(PX[wx])
        rows.append((wx, im, jm, PX[wx]))
    return rows

def P1():
    print("-"*78); print("P1  j_m = i_m - m c for each m; per-symbol c m-independent")
    a=b=0.2; D=0.6*Dc(a,b); cs=[]; ok=True
    for m in range(2,9):
        r=block(a,b,D,m); d=[jm-im for _,im,jm,_ in r]
        spread=max(d)-min(d); cs.append(-np.mean(d)/m)
        ok=ok and spread<1e-10
    print(f"     spreads ~{spread:.1e}; per-symbol c across m=2..8 std={np.std(cs):.1e} (c={cs[0]:.6f})")
    return rep("P1 j_m=i_m-mc, c m-independent (deterministic linear drift)", ok and np.std(cs)<1e-9)

def P2():
    print("-"*78); print("P2  centered process pathwise identical: j_m - m R(D) == i_m - m h")
    a=b=0.2; D=0.6*Dc(a,b)
    # h = entropy rate; R(D)=h-c
    pi=stat2(a,b); T=np.array([[1-a,a],[b,1-b]])
    h=-sum(pi[i]*T[i,j]*math.log2(T[i,j]) for i in range(2) for j in range(2) if T[i,j]>0)
    ok=True; worst=0
    for m in (4,6,8):
        r=block(a,b,D,m); c=-np.mean([jm-im for _,im,jm,_ in r])/m; R=h-c
        for _,im,jm,_ in r:
            worst=max(worst, abs((jm-m*R)-(im-m*h)))
    print(f"     max| (j_m-mR(D)) - (i_m-mh) | over words = {worst:.2e}  (h={h:.4f})")
    return rep("P2 centered j-process = centered i-process pathwise (=> FCLT 7.36 transfers)", worst<1e-9)

def P3():
    print("-"*78); print("P3  distribution shift: dist(j_m) = dist(i_m) - m c (=> LDP/MDP/EVT transfer)")
    a=b=0.2; D=0.6*Dc(a,b); m=8; r=block(a,b,D,m)
    c=-np.mean([jm-im for _,im,jm,_ in r])/m
    # compare sorted (j_m) vs sorted(i_m)-mc weighted by P; and the MAX (EVT) and a tail (LDP)
    ims=np.array([im for _,im,jm,_ in r]); jms=np.array([jm for _,im,jm,_ in r])
    P=np.array([p for *_,p in r]); P/=P.sum()
    shift_match=np.max(np.abs(jms-(ims-m*c)))
    max_shift=abs(jms.max()-(ims.max()-m*c))   # EVT: max j = max i - mc
    print(f"     max|j_m-(i_m-mc)|={shift_match:.2e}; EVT max-shift residual={max_shift:.2e}")
    print(f"     => overflow P(j_m>t)=P(i_m>t+mc): lossy LDP exponent E(R)=E_lossless(R+c)")
    return rep("P3 dist(j_m)=dist(i_m)-mc: tail/LDP/EVT functionals transfer", shift_match<1e-9)

if __name__=="__main__":
    print("="*78); print("Remark 7.34p -- lossy deviation hierarchy on Gray (7.34 x 7.29-7.36)")
    print("="*78)
    P1(); P2(); P3()
    print("="*78); print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
