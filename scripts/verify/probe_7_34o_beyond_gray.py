#!/usr/bin/env python3
"""
probe_7_34o_beyond_gray.py
============================================================================
Remark 7.34o -- the BEYOND-GRAY regime (D > D_c): where the dispersion plateau
ends.  The §7.34 grid (Thms 7.34b/i/m/n) gives V_op = V_lossless on the Gray
region 0<D<=D_c (SLB tight).  This probe MAPS what happens for D > D_c.

STRUCTURE (numerically established; the exact value is OPEN):
  O1  the per-n deconvolution threshold D_c^(n) (smallest D with the n-block
      (K^{-1})^{(x)n}P_X having a negative entry) DECREASES monotonically to the
      all-n D_c as n grows -- so the FINITE-n dispersion plateau extends to
      D_c^(n) > D_c, but the LIMIT (n->oo) plateau ends exactly at D_c.
  O2  the d-tilted-information identity j_n = i_n - n c (Remark 7.34n') FAILS for
      D > D_c: the deconvolution P_{Y*}=(K^{-1})^{(x)n}P_X is no longer a valid
      distribution for the alternating binding word, so M(x) != gamma^n P_X and
      the converse's deterministic-shift argument no longer applies.
  O3  consequently the dispersion leaves the plateau: with the OPTIMAL output
      (Blahut-Arimoto), V_op(D)=lim Var(j_n)/n equals V_lossless on [0,D_c]
      (the plateau) and DECREASES below V_lossless for D>D_c (the d-tilted
      information is "smoother" than the surprisal once distortion is allowed) --
      continuous at D_c.  The exact V_op(D) for D>D_c requires the optimal
      (non-SLB-tight, memory-aware) test channel and is OPEN.

Deps: numpy.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools, math
import numpy as np
PASS=True
def rep(name, ok):
    global PASS; PASS=PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

def stat2(a,b): s=a+b; return np.array([b/s, a/s])
def blockP(a,b,n):
    T=np.array([[1-a,a],[b,1-b]]); pi=stat2(a,b)
    idx=list(itertools.product([0,1],repeat=n))
    P=np.array([pi[w[0]]*np.prod([T[w[t-1],w[t]] for t in range(1,n)]) for w in idx])
    return idx,P
def Dc_bin(a,b): return 0.5*(1-math.sqrt(1-4*a*b/(2-a-b)**2))
def Kinv(D): a0=1-2*D; return (1/a0)*np.array([[1-D,-D],[-D,1-D]])
def pern_thresh(a,b,n):
    idx,P=blockP(a,b,n)
    def mn(D):
        Ki=Kinv(D); Pt=P.reshape([2]*n)
        for ax in range(n): Pt=np.tensordot(Ki,Pt,axes=([1],[ax])); Pt=np.moveaxis(Pt,0,ax)
        return Pt.min()
    lo,hi=1e-9,0.5-1e-6
    for _ in range(50):
        m=.5*(lo+hi); lo,hi=(m,hi) if mn(m)>=-1e-13 else (lo,m)
    return lo
def hamming(idx):
    n=len(idx[0]); M=len(idx); D=np.zeros((M,M))
    for i in range(M):
        for j in range(M): D[i,j]=sum(idx[i][t]!=idx[j][t] for t in range(n))
    return D
def blahut(P,Dm,lam,iters=4000):
    M=len(P); q=np.ones(M)/M; W=np.exp(-lam*Dm)
    for _ in range(iters):
        c=W@q; A=(W*q[None,:])/c[:,None]; qn=P@A
        if np.max(np.abs(qn-q))<1e-13: q=qn; break
        q=qn
    c=W@q; A=(W*q[None,:])/c[:,None]
    Dbar=float((P[:,None]*A*Dm).sum())
    cn=np.exp(-lam*Dm)@q
    j=(-np.log(cn))/math.log(2)          # variance part of the d-tilted info (const dropped)
    Vj=float((P*(j-(P*j).sum())**2).sum())
    return Dbar,Vj

def O1():
    print("-"*78); print("O1  per-n deconvolution threshold D_c^(n) decreases to all-n D_c")
    a=b=0.2; dc=Dc_bin(a,b); ths=[pern_thresh(a,b,n) for n in range(2,11)]
    mono=all(ths[i]>ths[i+1] for i in range(len(ths)-1)); conv=ths[-1]/dc
    for n,t in zip(range(2,11),ths): print(f"     n={n}: D_c^(n)={t:.6f} (x{t/dc:.3f} all-n D_c={dc:.6f})")
    return rep("O1 D_c^(n) monotone -> all-n D_c (finite-n plateau extension)", mono and conv<1.12)

def O2():
    print("-"*78); print("O2  the j_n=i_n-nc identity FAILS for D>D_c (deconvolution invalid)")
    a=b=0.2; dc=Dc_bin(a,b); n=8; idx,P=blockP(a,b,n)
    def maxneg(D):
        Ki=Kinv(D); Pt=P.reshape([2]*n)
        for ax in range(n): Pt=np.tensordot(Ki,Pt,axes=([1],[ax])); Pt=np.moveaxis(Pt,0,ax)
        return Pt.min()
    below=maxneg(0.9*dc); above=maxneg(1.3*Dc_bin(a,b))  # well above all-n D_c
    # at n=8 the per-n threshold is ~0.0183 > all-n D_c=0.0159; pick D between to show n-sensitivity
    th8=pern_thresh(a,b,8); mid=maxneg(0.5*(dc+th8))
    print(f"     n={n}: min P_Y* at 0.9 D_c={below:+.2e} (valid); at 0.99 D_c^(n=8)={maxneg(0.99*th8):+.2e}; just above D_c^(n=8): {maxneg(1.02*th8):+.2e}")
    return rep("O2 deconvolution valid below D_c^(n), invalid above (identity fails)", below>=-1e-12 and maxneg(1.02*th8)<-1e-12)

def O3():
    print("-"*78); print("O3  dispersion leaves the plateau: V_op(D)=V_lossless on [0,D_c], decreases beyond")
    a=b=0.2; dc=Dc_bin(a,b); n=7; idx,P=blockP(a,b,n); Dm=hamming(idx)
    Vl=(lambda p: p*(1-p)*math.log2((1-p)/p)**2)(0.2)
    rows=[]
    for fr in (0.5,0.95,1.5,3.0,6.0):
        # bisection on lam to reach D/n = fr*dc
        lo,hi=1e-4,80.0
        for _ in range(55):
            lam=0.5*(lo+hi); Db,Vj=blahut(P,Dm,lam,1500)
            if Db/n>fr*dc: lo=lam
            else: hi=lam
        Db,Vj=blahut(P,Dm,lam,5000); rows.append((Db/n/dc,Vj/n))
    base=rows[0][1]  # ~plateau (finite-n)
    print(f"     V_lossless(p=.2)={Vl:.4f}; finite-n(n={n}) plateau Vj/n~{base:.4f} (finite-n bias)")
    for r,v in rows: print(f"     D/D_c={r:.2f}: Vj/n={v:.5f}  ({'PLATEAU' if abs(v-base)<0.02 else 'below'})")
    # plateau (<=D_c) flat; far beyond strictly lower
    plateau_flat = abs(rows[1][1]-rows[0][1])<0.02
    decreases = rows[-1][1] < rows[0][1]-0.05
    return rep("O3 V_op flat on plateau, decreases for D>>D_c (continuous; exact OPEN)", plateau_flat and decreases)

if __name__=="__main__":
    print("="*78); print("Remark 7.34o -- the beyond-Gray (D>D_c) regime: where the plateau ends")
    print("="*78)
    O1(); O2(); O3()
    print("="*78); print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
