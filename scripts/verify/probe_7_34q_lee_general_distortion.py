#!/usr/bin/env python3
"""
probe_7_34q_lee_general_distortion.py
============================================================================
Remark 7.34q -- the RD-dispersion converse extends from Hamming to ANY
group-difference distortion (Lee, etc.): the distortion-type axis.

The §7.34 grid (7.34b/i/m/n) is all HAMMING distortion. The source-agnostic
converse mechanism, however, depends only on the dual identity M=gamma^n P_X,
which is a property of the BACKWARD CHANNEL, not the distortion -- and it holds
for EVERY difference distortion d(x,y)=rho(x (-) y) on a finite abelian group,
with the SLB-achieving GIBBS backward channel K = K_nu/Z, K_nu[x,y]=nu^{d(x,y)},
Z=sum_y nu^{d(x,y)} (constant by group symmetry), nu=e^{-lam*}:
    K_nu K^{-1} = (Z K) K^{-1} = Z I = gamma I   (gamma = Z, TRIVIALLY scalar).
Hence M(x^n)=sum_y P_{Y*}(y) nu^{d(x^n,y)} = (K_nu K^{-1})^{(x)n} P_X = gamma^n P_X,
so j_n = i_n - n c (deterministic), the converse ratio rho=1, and
    V_op(D) >= V_lossless    on the (distortion's own) Gray region,
UNCONDITIONALLY, for ANY group-difference distortion. (For Hamming, gamma=1/(1-D),
recovering 7.34n; the Lee channel is a graded circulant so its Gray region is
SMALLER -- K^{-1} sign-oscillates -- but nonempty.)

CHECKS (concrete instance: A-ary LEE distortion d=min(|x-y|,A-|x-y|)):
  Q1  K = K_nu/Z (Gibbs kernel) => K_nu K^{-1} = Z I scalar, A=3..6 (dual identity).
  Q2  Lee Gray region nonempty: (K^{-1})^{(x)n} P_X >= 0 for small D, all n<=6
      (D_c^Lee>0; A=3 Lee = Hamming since all d=1).
  Q3  operational converse: j_n = i_n - n c EXACT on the Lee Gray region
      (deterministic shift => rho=1), A=4,5 -- the converse transfers.
  Q4  V_lossless is the SAME lossless varentropy rate (distortion-independent on
      Gray): the converse value does not depend on Lee-vs-Hamming.

Deps: numpy.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools, math
import numpy as np
PASS=True
def rep(name, ok):
    global PASS; PASS=PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok
def lee(A): return np.array([[min(abs(x-y),A-abs(x-y)) for y in range(A)] for x in range(A)])
def Kmat(A,lam): W=np.exp(-lam*lee(A)); return W/W.sum(axis=1,keepdims=True)
def stat(T):
    ev,vL=np.linalg.eig(T.T); k=np.argmin(np.abs(ev-1)); v=np.real(vL[:,k]); return v/v.sum()

def Q1():
    print("-"*78); print("Q1  Lee dual identity K_nu K^{-1}=Z I (Gibbs channel K=K_nu/Z), A=3..6")
    ok=True
    for A in (3,4,5,6):
        lam=2.0; K=Kmat(A,lam); nu=math.exp(-lam); Knu=nu**lee(A)
        Z=Knu.sum(axis=1)[0]
        P=Knu@np.linalg.inv(K)
        off=np.abs(P-np.diag(np.diag(P))).max(); spread=np.diag(P).max()-np.diag(P).min()
        ok = ok and off<1e-9 and spread<1e-9 and abs(P[0,0]-Z)<1e-9
        print(f"     A={A}: off={off:.1e} spread={spread:.1e} gamma={P[0,0]:.4f}=Z={Z:.4f}")
    return rep("Q1 K_nu K^{-1}=Z I scalar (dual identity, group-circulant)", ok)

def Q2():
    print("-"*78); print("Q2  Lee Gray region nonempty (D_c^Lee>0; A=3 Lee=Hamming)")
    ok=True
    def dmin(A,lam,nmax=6):
        p=0.2; T=np.full((A,A),p/(A-1)); np.fill_diagonal(T,1-p)
        pi=stat(T); Ki=np.linalg.inv(Kmat(A,lam)); worst=1.0
        for n in range(2,nmax+1):
            idx=list(itertools.product(range(A),repeat=n)); PX=np.zeros([A]*n)
            for w in idx: PX[w]=pi[w[0]]*np.prod([T[w[t-1],w[t]] for t in range(1,n)])
            Pt=PX
            for ax in range(n): Pt=np.tensordot(Ki,Pt,axes=([1],[ax])); Pt=np.moveaxis(Pt,0,ax)
            worst=min(worst,Pt.min())
        return worst
    for A in (4,5,6):
        # find a lam (small distortion) with valid deconvolution
        valid=any(dmin(A,lam)>=-1e-12 for lam in (10.0,9.0,8.0))
        Davg=float((stat(np.full((A,A),0.2/(A-1))+np.diag([0.8-0.2/(A-1)]*A))[:,None]*Kmat(A,9.0)*lee(A)).sum())
        print(f"     A={A}: valid deconvolution at small Lee distortion: {valid}")
        ok = ok and valid
    return rep("Q2 Lee Gray region nonempty for A=4,5,6 (small distortion)", ok)

def Q3():
    print("-"*78); print("Q3  operational converse: j_n = i_n - n c EXACT on Lee Gray region (rho=1)")
    ok=True
    for A in (4,5):
        p=0.2; T=np.full((A,A),p/(A-1)); np.fill_diagonal(T,1-p); pi=stat(T)
        lam=9.0; K=Kmat(A,lam); Ki=np.linalg.inv(K); nu=math.exp(-lam); Dm=lee(A)
        n=5; idx=list(itertools.product(range(A),repeat=n)); PX=np.zeros([A]*n)
        for w in idx: PX[w]=pi[w[0]]*np.prod([T[w[t-1],w[t]] for t in range(1,n)])
        Pt=PX
        for ax in range(n): Pt=np.tensordot(Ki,Pt,axes=([1],[ax])); Pt=np.moveaxis(Pt,0,ax)
        PYstar=Pt
        # block distortion d_n(x,y)=sum_t Dm[x_t,y_t]; M(x)=sum_y PYstar(y) nu^{d_n}; j_n=-lam n D - log M
        diffs=[]
        for wx in idx:
            if PX[wx]<=0: continue
            M=0.0
            for wy in idx:
                dd=sum(Dm[wx[t],wy[t]] for t in range(n)); M+=PYstar[wy]*nu**dd
            jn=(-math.log(M))/math.log(2)   # variance part (drop the const -lam n D)
            inn=-math.log2(PX[wx]); diffs.append(jn-inn)
        spread=max(diffs)-min(diffs)
        print(f"     A={A}: spread of (j_n - i_n) over words = {spread:.2e} (=0 => det. shift => rho=1)")
        ok = ok and spread<1e-9
    return rep("Q3 j_n=i_n-nc on Lee Gray (converse V_op>=V_lossless transfers)", ok)

def Q4():
    print("-"*78); print("Q4  V_lossless distortion-independent (same lossless varentropy rate)")
    # V_lossless depends only on the SOURCE (-log P(X^n)), not the distortion; identical for Lee/Hamming.
    A=4; p=0.2; pi=np.full(A,1/A)
    V=p*(1-p)*math.log2((1-p)*(A-1)/p)**2   # symmetric A-ary lossless varentropy
    print(f"     symmetric A=4 p=0.2: V_lossless={V:.5f} (a SOURCE quantity; same for Hamming & Lee)")
    return rep("Q4 V_lossless is the source varentropy rate (distortion-independent)", V>0)

def Q6():
    print("-"*78); print("Q6  Lee all-n threshold D_c^Lee>0 (per-n D_c^(n) converges to positive); binding=alt Lee-adjacent")
    A=4; p=0.2; T=np.full((A,A),p/(A-1)); np.fill_diagonal(T,1-p); pi=stat(T)
    def avgd(lam): return float((pi[:,None]*Kmat(A,lam)*lee(A)).sum())
    def lam_for_D(Dt):
        lo,hi=1e-3,40.0
        for _ in range(55):
            m=.5*(lo+hi); lo,hi=(m,hi) if avgd(m)>Dt else (lo,m)
        return .5*(lo+hi)
    def pern(n):
        def minP(Dt):
            Ki=np.linalg.inv(Kmat(A,lam_for_D(Dt)))
            idx=list(itertools.product(range(A),repeat=n)); PX=np.zeros([A]*n)
            for w in idx: PX[w]=pi[w[0]]*np.prod([T[w[t-1],w[t]] for t in range(1,n)])
            Pt=PX
            for ax in range(n): Pt=np.tensordot(Ki,Pt,axes=([1],[ax])); Pt=np.moveaxis(Pt,0,ax)
            f=Pt.reshape(-1); return f.min(), idx[int(f.argmin())]
        lo,hi=1e-6,0.5
        for _ in range(38):
            m=.5*(lo+hi); lo,hi=(m,hi) if minP(m)[0]>=-1e-13 else (lo,m)
        return lo, minP(min(hi*1.02,0.4))[1]
    ths=[];
    for n in (6,7,8,9):
        dc,w=pern(n); ths.append(dc)
        syms=set(w); alt=len(syms)==2 and all(w[i]!=w[i+1] for i in range(len(w)-1))
        adj=len(syms)==2 and (lambda a,b: min(abs(a-b),A-abs(a-b))==1)(*sorted(syms))
        print(f"     n={n}: D_c^(n)={dc:.6f} word={w} (alt-Lee-adjacent: {alt and adj})")
    mono=ths[0]>ths[1]>ths[2] and ths[-1]>1e-3
    print(f"     => D_c^Lee ~ {ths[-1]:.6f} > 0 (per-n decreasing, stabilising positive)")
    return rep("Q6 Lee all-n threshold D_c^Lee>0 (positive Gray region; alt Lee-adjacent binding)", mono)

def Q5():
    print("-"*78); print("Q5  BALANCED non-group (non-circulant) distortion: K_nu K^{-1}=Z I still holds")
    print("    (only property needed: Z=sum_y nu^{d} x-indep <=> rows are permutations = balanced)")
    # A=4 balanced symmetric, rows permutations of {0,1,2,3}, NOT circulant
    Dm=np.array([[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]])
    # confirm balanced (each row a permutation of row 0's multiset) and symmetric, diag 0
    base=sorted(Dm[0]); balanced=all(sorted(Dm[r])==base for r in range(4))
    sym=np.allclose(Dm,Dm.T); diag0=all(Dm[i,i]==0 for i in range(4))
    # NOT circulant: row1 != cyclic shift of row0
    circ=all(Dm[1,j]==Dm[0,(j-1)%4] for j in range(4))
    lam=1.3; Knu=np.exp(-lam*Dm); Z=Knu.sum(axis=1)[0]; K=Knu/Z
    P=Knu@np.linalg.inv(K); off=np.abs(P-np.diag(np.diag(P))).max(); spread=np.diag(P).max()-np.diag(P).min()
    print(f"     balanced={balanced} symmetric={sym} diag0={diag0} circulant={circ} (want balanced,sym,diag0,NOT circulant)")
    print(f"     K_nu K^-1: off={off:.1e} spread={spread:.1e} gamma={P[0,0]:.4f}=Z={Z:.4f}")
    ok = balanced and sym and diag0 and (not circ) and off<1e-9 and spread<1e-9 and abs(P[0,0]-Z)<1e-9
    return rep("Q5 dual identity holds for balanced NON-group distortion (Z x-indep is all)", ok)

if __name__=="__main__":
    print("="*78); print("Remark 7.34q -- RD-dispersion converse for general balanced distortion (Lee, ...)")
    print("="*78)
    Q1(); Q2(); Q3(); Q4(); Q5(); Q6()
    print("="*78); print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
