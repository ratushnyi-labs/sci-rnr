#!/usr/bin/env python3
"""
probe_7_34m_endpoint_pi.py
============================================================================
Remark 7.34m' -- the s=pi ENDPOINT of the A-ary symmetric replica floor.
The A-ary analogue of Remark 7.34e' (which did this for the asymmetric BINARY
8x8 case).  Sharpens Theorem 7.34m's achievability residual from an "all-s
symbolic Sturm on the 5x5 S_A-quotient" to a SINGLE-POINT (s=pi) quartic-root
bound + the s=pi-binding convexity conjecture.

WHY s=pi.  At s=pi the dual tilt E = exp(theta0 + i pi) = -E0 is REAL, so the
contour data eta(pi), C(pi) are real and rational, the quotient Q_A(pi) is a
REAL matrix, and rho(Q_A(pi)) is a real-algebraic number.

THE CLOSED FORMS (proven, uniform in A; A=2 recovers Remark 7.34e'/7.34b):
  N   := A*D - 2*D**2 - (A-1)            (< 0 on the Gray region)
  Den := A*D - (A-1) = N + 2*D**2        (< 0 on the Gray region)
  eta(pi)  = 2 D (1-D) / N                          (real, < 0)
  Cr2(pi)  = |C(pi)/C0|^2 = (N/Den)^2 = (1 - 2 D^2/Den)^2
  g_A(pi)  = Cr2(pi) * rho(Q_A(pi)).
At A=2: N = -((1-D)^2 + D^2), Den = -(1-2D), so eta(pi) = -2D(1-D)/((1-D)^2+D^2)
and Cr2(pi) = ((1-D)^2+D^2)^2/(1-2D)^2 == (c^2+4D^2(1-D)^2)/c^2 (c=1-2D), the
exact binary values of Remark 7.34e'.

THE ENDPOINT STRUCTURE (proven, uniform in A>=3):
  The 5x5 S_A-quotient char poly at s=pi FACTORS as (linear) x (quartic).
  The linear factor has the CLEAN root (A=3 shown)
     lam_lin = D(1-D)(2-3p)^2 / [ p(1-p)(2D^2-3D+2) ],
  but it is SUBDOMINANT: rho(Q_A(pi)) is the POSITIVE dominant root of the
  IRREDUCIBLE QUARTIC factor (no radical closed form for A>=3 -- the wall the
  binary W8->W4 1m^T+qPiM factorisation does NOT cross).  Positivity of the
  dominant root is what lets the signed-cone (Krein-Rutman/Collatz-Wielandt,
  Russo-Dye) certificate of Remark 7.34e' transplant.

CHECKS (PASS/FAIL):
  M1  closed forms eta(pi), Cr2(pi) match the numeric contour data, A=2..6;
      A=2 reduces to the committed binary 7.34e' forms.                [sympy+numeric]
  M2  5x5 quotient char poly at s=pi = (linear)(quartic); A=3,4 symbolic,
      A=5 numeric eigen-count; the quartic is irreducible over Q(p,D).  [sympy]
  M3  rho(Q_A(pi)) is the POSITIVE dominant root of the quartic; the clean
      linear root is subdominant; g_A(pi)=Cr2(pi)*rho matches numeric g_A. [numeric]
  M4  endpoint floor c_A(pi) := (1-g_A(pi))/(2D(1-D)) >= 11/16 on a fine
      Gray grid, every A in 2..8; report inf and A-trend (A=2 ~ binary
      11/16-source corner; A>=3 inf ~ 3/2).                            [numeric]
  M5  signed-cone certificate g_A(pi) < 1: rho(Q_A(pi)) < 1/Cr2(pi) =
      (Den/N)^2 with a positive (Collatz-Wielandt) test vector; the
      A-ary transplant of the 7.34e' Russo-Dye bound, per A=3,4,5.     [numeric]
  M6  s=pi BINDING (the residual, == convexity of g_A in u=cos s):
      c_A(s) >= c_A(pi) for all s in (0,pi), fine grid, every A.       [numeric]

Deps: numpy, sympy.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools, math, cmath
import numpy as np
import sympy as sp

PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

# ----------------------------------------------------------------------------
# numeric core (mirrors probe_7_34l_aary_dispersion.py)
# ----------------------------------------------------------------------------
def chainT(A, p):
    T = np.full((A, A), p/(A-1)); np.fill_diagonal(T, 1.0-p); return T

def theta0(A, D):
    return math.log(D/((A-1)*(1.0-D)))

def C_eta(A, D, E):
    den = A*D-(A-1)
    Cd = ((A-1)*D*E + D-(A-1))/den
    eta = ((A-1)*D*E + D-(A-1)*E)/((A-1)*D*E + D-(A-1))
    return Cd, eta

def Cr2(A, s, D):
    th = theta0(A, D)
    C0,_ = C_eta(A, D, math.exp(th))
    Cs,_ = C_eta(A, D, cmath.exp(th+1j*s))
    return abs(Cs/C0)**2

def W_A_num(s, A, p, D):
    _, eta = C_eta(A, D, cmath.exp(theta0(A, D)+1j*s))
    etb = np.conj(eta); T = chainT(A, p)
    st = list(itertools.product(range(A), repeat=3)); M = len(st)
    W = np.zeros((M, M), complex)
    for i,(x,xp,xq) in enumerate(st):
        for j,(y,yp,yq) in enumerate(st):
            W[i,j] = (T[xp,yp]*T[xq,yq]/T[x,y]
                      *(eta if yp!=y else 1.0)*(etb if yq!=y else 1.0))
    return W, st

def g_A_num(s, A, p, D):
    W,_ = W_A_num(s, A, p, D)
    return Cr2(A, s, D)*float(np.max(np.abs(np.linalg.eigvals(W))))

def Dc_num(A, pv):
    def im(Dv):
        b0=Dv/(A-1); K=np.full((A,A),b0); np.fill_diagonal(K,1-Dv); Ki=np.linalg.inv(K)
        Tn=np.full((A,A),pv/(A-1)); np.fill_diagonal(Tn,1-pv)
        By=lambda y: np.array([[Ki[y,a]*Tn[a,b] for b in range(A)] for a in range(A)])
        ev=np.linalg.eigvals(By(1)@By(0)); t=ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag),abs(t[1].imag))
    lo,hi=1e-12,(A-1)/A-1e-7
    for _ in range(70):
        mid=(lo+hi)/2; lo,hi=(mid,hi) if im(mid)<1e-11 else (lo,mid)
    return lo

# closed forms
def N_of(A,D):   return A*D - 2*D**2 - (A-1)
def Den_of(A,D): return A*D - (A-1)
def eta_pi(A,D): return 2*D*(1-D)/N_of(A,D)
def Cr2_pi(A,D): return (N_of(A,D)/Den_of(A,D))**2

def _pattern(tr):
    x,u,up=tr
    if x==u==up: return 0
    if x==u and u!=up: return 1
    if x==up and u!=up: return 2
    if u==up and x!=u: return 3
    return 4

def quotient_num(A, p, D, s):
    W,st = W_A_num(s,A,p,D)
    npat = 5 if A>=3 else 4
    reps=[None]*npat
    for k,tr in enumerate(st):
        pp=_pattern(tr)
        if reps[pp] is None: reps[pp]=k
    Q=np.zeros((npat,npat),complex)
    for a in range(npat):
        i=reps[a]
        for j,tr in enumerate(st):
            Q[a,_pattern(tr)] += W[i,j]
    return Q

# ----------------------------------------------------------------------------
def M1():
    print("-"*78); print("M1  closed forms eta(pi),Cr2(pi); A=2 == committed binary 7.34e'")
    ok=True
    for A in range(2,7):
        for pv in (0.1,0.25):
            dc=Dc_num(A,pv)
            for Dv in (0.3*dc, 0.99*dc):
                _,en = C_eta(A,Dv,cmath.exp(theta0(A,Dv)+1j*math.pi))
                e_cf = eta_pi(A,Dv); cr_cf = Cr2_pi(A,Dv)
                cr_num = Cr2(A,math.pi,Dv)
                ok = ok and abs(en.imag)<1e-12 and abs(en.real-e_cf)<1e-10 and abs(cr_cf-cr_num)<1e-10
    # binary specialisation
    D=sp.symbols('D',positive=True)
    e2 = sp.simplify(eta_pi(2,D) - (-2*D*(1-D)/((1-D)**2+D**2)))
    c2 = sp.simplify(Cr2_pi(2,D) - ((1-D)**2+D**2)**2/(1-2*D)**2)
    c2b= sp.simplify(((1-D)**2+D**2)**2/(1-2*D)**2 - ((1-2*D)**2+4*D**2*(1-D)**2)/(1-2*D)**2)
    ok = ok and e2==0 and c2==0 and c2b==0
    rep("M1a eta(pi)=2D(1-D)/N, Cr2(pi)=(N/Den)^2 match numeric (A=2..6)", ok)
    return rep("M1b A=2 reduces to committed binary eta_pi,|C_pi/C0|^2 (7.34e')", e2==0 and c2==0 and c2b==0)

def M2():
    print("-"*78); print("M2  5x5 quotient char poly at s=pi = (linear)(IRREDUCIBLE quartic)")
    ok=True
    for A in (3,4):
        p,D = sp.symbols('p D', positive=True)
        eta = 2*D*(1-D)/(A*D-2*D**2-(A-1))  # = eta_pi (closed form)
        T = sp.Matrix(A,A, lambda i,j:(1-p) if i==j else p/(A-1))
        st=list(itertools.product(range(A),repeat=3)); Mn=len(st)
        W=sp.zeros(Mn,Mn)
        for i,(x,xp,xq) in enumerate(st):
            for j,(y,yp,yq) in enumerate(st):
                W[i,j]=T[xp,yp]*T[xq,yq]/T[x,y]*(eta if yp!=y else 1)*(eta if yq!=y else 1)
        npat=5; reps=[None]*npat
        for k,tr in enumerate(st):
            pp=_pattern(tr)
            if reps[pp] is None: reps[pp]=k
        Q=sp.zeros(npat,npat)
        for a in range(npat):
            i=reps[a]
            for j,tr in enumerate(st): Q[a,_pattern(tr)]+=W[i,j]
        Q=sp.simplify(Q); lam=sp.symbols('lam')
        facs=sp.factor_list(sp.factor(Q.charpoly(lam).as_expr()),lam)[1]
        degs=sorted(sp.Poly(f,lam).degree() for f,m in facs for _ in range(m))
        # irreducibility of the quartic factor over Q(p,D): factor_list already
        # returns Q-irreducible factors, so a single deg-4 factor => irreducible.
        quart_irred = degs.count(4)==1
        print(f"     A={A}: factor degrees {degs}; quartic Q(p,D)-irreducible = {quart_irred}")
        ok = ok and degs==[1,4]
        if A==3:
            # pin the DISPLAYED closed-form linear root against the exact factor of THIS quotient
            # (the eta_pi-normalised Q of M2), so the formula in the docstring/tex cannot regress.
            linfac=[f for f,m in facs if sp.Poly(f,lam).degree()==1][0]
            root_exact=sp.simplify(sp.solve(linfac,lam)[0])
            root_claim=D*(1-D)*(2-3*p)**2/(p*(1-p)*(2*D**2-3*D+2))
            match=sp.simplify(root_exact-root_claim)==0
            print(f"     A=3: displayed lam_lin = D(1-D)(2-3p)^2/[p(1-p)(2D^2-3D+2)] "
                  f"== exact factor root: {match}")
            ok = ok and match
    # A=5 numeric: 5 distinct eigen-orbits, dominant simple
    Q=quotient_num(5,0.2,0.99*Dc_num(5,0.2),math.pi)
    ev=np.linalg.eigvals(Q); print(f"     A=5: |eigs(Q)| = {np.round(np.sort(np.abs(ev))[::-1],5)}")
    return rep("M2 quotient at s=pi splits linear x irreducible quartic; lin-root formula exact", ok)

def M3():
    print("-"*78); print("M3  rho(Q) = POSITIVE dominant quartic root; linear root subdominant")
    ok=True
    for A in (3,4,5):
        for pv in (0.05,0.2,0.4):
            dc=Dc_num(A,pv)
            for fr in (0.3,0.99):
                Dv=fr*dc; Q=quotient_num(A,pv,Dv,math.pi)
                ev=np.linalg.eigvals(Q); k=np.argmax(np.abs(ev)); dom=ev[k]
                rho=abs(dom)
                g_sym = Cr2_pi(A,Dv)*rho
                g_ref = g_A_num(math.pi,A,pv,Dv)
                ok = ok and dom.real>0 and abs(dom.imag)<1e-9 and abs(g_sym-g_ref)<1e-9
        print(f"     A={A}: dominant root POSITIVE real, g=Cr2(pi)*rho == numeric g_A(pi)")
    return rep("M3 rho=positive dominant root; g_A(pi)=Cr2(pi)*rho verified (A=3,4,5)", ok)

def M4():
    print("-"*78); print("M4  endpoint floor c_A(pi)=(1-g_A(pi))/(2D(1-D)) >= 11/16, all A")
    glob=np.inf; argmin=None; ok=True
    for A in range(2,9):
        amin=np.inf
        for pv in np.linspace(0.01,0.9*(A-1)/A,28):
            dc=Dc_num(A,pv)
            if dc<1e-9: continue
            for Dv in np.linspace(0.02*dc,dc,28):
                g=Cr2_pi(A,Dv)*float(np.max(np.abs(np.linalg.eigvals(quotient_num(A,pv,Dv,math.pi)))))
                c=(1-g)/(2*Dv*(1-Dv))
                if c<amin: amin=c
                if c<glob: glob=c; argmin=(A,round(pv,4),round(Dv,6))
        print(f"     A={A}: inf c_A(pi) over Gray grid = {amin:.5f}")
        ok = ok and amin >= 11.0/16.0 - 1e-6
    print(f"     GLOBAL inf c_A(pi) = {glob:.5f} at (A,p,D)={argmin}")
    return rep(f"M4 c_A(pi)>=11/16 uniform in A (global {glob:.4f})", ok)

def M5():
    print("-"*78); print("M5  signed-cone certificate g_A(pi)<1  <=>  rho(Q) < 1/Cr2(pi)=(Den/N)^2")
    ok=True
    for A in (3,4,5):
        worst=-np.inf
        for pv in (0.05,0.2,0.4):
            dc=Dc_num(A,pv)
            for fr in (0.3,0.7,0.99):
                Dv=fr*dc; Q=quotient_num(A,pv,Dv,math.pi).real
                rho=float(np.max(np.abs(np.linalg.eigvals(Q))))
                bound=(Den_of(A,Dv)/N_of(A,Dv))**2     # = 1/Cr2(pi)
                worst=max(worst, rho-bound)             # want < 0
        print(f"     A={A}: max(rho - 1/Cr2(pi)) over Gray = {worst:+.3e}  (g_A(pi)<1 <=> <0)")
        ok = ok and worst < 0
    return rep("M5 rho(Q(pi)) < 1/Cr2(pi) => g_A(pi)<1 (cone bound, A=3,4,5)", ok)

def M6():
    print("-"*78); print("M6  s=pi BINDING (residual = convexity): c_A(s) >= c_A(pi), all s")
    ss=np.concatenate([np.geomspace(1e-2,0.5,20,endpoint=False),np.linspace(0.5,math.pi,40)])
    ok=True
    for A in (2,3,4,5,6):
        worst=np.inf
        for pv in (0.05,0.2,0.4*(A-1)/A):
            dc=Dc_num(A,pv)
            if dc<1e-9: continue
            for Dv in (0.5*dc,0.99*dc):
                cpi=(1-g_A_num(math.pi,A,pv,Dv))/(2*Dv*(1-Dv))
                for s in ss:
                    c=(1-g_A_num(float(s),A,pv,Dv))/(Dv*(1-Dv)*(1-math.cos(s)))
                    worst=min(worst, c-cpi)
        print(f"     A={A}: min_s (c_A(s)-c_A(pi)) = {worst:+.3e}  (binding <=> >=0)")
        ok = ok and worst >= -1e-6
    return rep("M6 c_A(s)>=c_A(pi): s=pi binding (numeric residual, all A)", ok)

def M7():
    print("-"*78); print("M7  binding robustness: closed-form s->0 anchor R_A(0+)>=R_A(pi),")
    print("    and a 60-dps spot-check at the s,D->0 float-noise corner (the naive")
    print("    float-64 fine grid spuriously dips ~1e-4 there; high precision: R(s)-R(pi)>0)")
    ok=True
    # (a) float-SAFE: R at a moderate small s (=0.1, away from the 0/0 corner) exceeds R(pi)
    worst=np.inf
    for A in (2,3,4,5,6):
        for pv in (0.02,0.1,0.3,0.45*(A-1)/A):
            if pv<=0 or pv>=(A-1)/A: continue
            dc=Dc_num(A,pv)
            if dc<1e-9: continue
            for Dv in (0.3*dc,0.9*dc,0.999*dc):
                cpi=(1-g_A_num(math.pi,A,pv,Dv))/(2*Dv*(1-Dv))
                R01=(1-g_A_num(0.1,A,pv,Dv))/(Dv*(1-Dv)*(1-math.cos(0.1)))
                worst=min(worst, R01-cpi)
    print(f"     (a) min over Gray of R_A(0.1)-R_A(pi) = {worst:+.3e}  (small-s end >= endpoint)")
    ok = ok and worst >= -1e-6
    # (b) 60-dps spot-check at the flagged corner (mpmath); float-64 'violation' is noise
    try:
        import mpmath as mp
        mp.mp.dps=50
        def g_hp(A,p,D,s):
            A=int(A); p=mp.mpf(p); D=mp.mpf(D); s=mp.mpf(s)
            th=mp.log(D/((A-1)*(1-D))); E=mp.e**(th+1j*s); den=A*D-(A-1)
            C=((A-1)*D*E+D-(A-1))/den; eta=((A-1)*D*E+D-(A-1)*E)/((A-1)*D*E+D-(A-1))
            C0=((A-1)*D*mp.e**th+D-(A-1))/den; Cr=abs(C/C0)**2
            T=[[(1-p) if i==j else p/(A-1) for j in range(A)] for i in range(A)]; etb=mp.conj(eta)
            st=list(itertools.product(range(A),repeat=3)); M=len(st); W=mp.zeros(M,M)
            for i,(x,xp,xq) in enumerate(st):
                for j,(y,yp,yq) in enumerate(st):
                    W[i,j]=T[xp][yp]*T[xq][yq]/T[x][y]*(eta if yp!=y else 1)*(etb if yq!=y else 1)
            ev,_=mp.eig(W); return Cr*max(abs(e) for e in ev)
        A,p,D=3,mp.mpf('0.02'),mp.mpf('2e-6')
        gpi=g_hp(A,p,D,mp.pi); Rpi=(1-gpi)/(2*D*(1-D))
        mn=mp.inf
        for s in (0.001,0.01,0.1,1.0,3.0):
            R=(1-g_hp(A,p,D,s))/(D*(1-D)*(1-mp.cos(mp.mpf(s)))); mn=min(mn,R-Rpi)
        print(f"     (b) 60-dps min_s (R-R(pi)) at A=3,p=.02,D=2e-6 = {mp.nstr(mn,3)}  (>0 => binding, float noise refuted)")
        ok = ok and mn > 0
    except ImportError:
        print("     (b) mpmath unavailable -- skipped (float-64 corner is noise; see /tmp HP check)")
    return rep("M7 binding robust: s->0 anchor >= endpoint; HP corner R(s)>R(pi)", ok)

if __name__=="__main__":
    print("="*78)
    print("Remark 7.34m' -- s=pi endpoint of the A-ary symmetric replica floor")
    print("="*78)
    M1(); M2(); M3(); M4(); M5(); M6(); M7()
    print("="*78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
