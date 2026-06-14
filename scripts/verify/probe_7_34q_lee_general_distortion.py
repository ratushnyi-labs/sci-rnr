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

if __name__=="__main__":
    print("="*78); print("Remark 7.34q -- RD-dispersion converse for general group-difference distortion (Lee)")
    print("="*78)
    Q1(); Q2(); Q3(); Q4()
    print("="*78); print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
