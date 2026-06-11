import math, numpy as np, mpmath as mp
mp.mp.dps = 40

def Dc(p): return 0.5*(1-math.sqrt(1-2*p)/(1-p))

def eta_C_mp(p, D, z):
    # z complex; beta = theta0 + i z  (z = s + i*y, the spectral variable on the strip)
    th0 = mp.log(D/(1-D))
    eb = mp.e**(th0 + 1j*z)
    ze = (1-eb)/((1+eb)*(1-2*D))
    eta = (1-ze)/(1+ze)
    C = (1+eb)*(1+ze)/2
    C0 = 1/(1-D)
    return eta, C/C0

def logS_mp(x, p, eta):
    # O(n) polymer/transfer recursion for S_x(eta), mpmath complex
    n=len(x); rho=mp.mpf(p)/(1-p)
    l=[mp.mpf(1)]*n; r=[mp.mpf(1)]*n
    for i in range(1,n):
        l[i]= rho if x[i]==x[i-1] else 1/rho
    for b in range(0,n-1):
        r[b]= rho if x[b+1]==x[b] else 1/rho
    Fip1=mp.mpf(1); Fip2=mp.mpf(1); Gip1=mp.mpf(0)
    logsc=mp.mpf(0)
    for i in range(n-1,-1,-1):
        Gi=eta*(r[i]*Fip2+Gip1)
        Fi=Fip1+l[i]*Gi
        Fip2=Fip1; Fip1=Fi; Gip1=Gi
        a=abs(Fip1)
        if a>0:
            Fip1/=a; Fip2/=a; Gip1/=a; logsc+=mp.log(a)
    return mp.log(Fip1)+logsc

def logPhi(x,p,D,z):
    eta,Cr = eta_C_mp(p,D,z)
    return len(x)*mp.log(Cr) + logS_mp(x,p,eta)

# pick a shell word
rng=np.random.default_rng(3)
p=0.4; D=0.9*Dc(p); n=200
while True:
    x=(np.cumsum(rng.random(n)<p)%2).astype(int)
    sw=int(np.sum(x[1:]!=x[:-1]))
    if abs(sw-(n-1)*p)<=2*math.sqrt(n*p*(1-p)): break
x=[int(b) for b in x]

# cumulants by complex-step / finite diff at s=0 of log Phi(theta0+i s) ... here z plays role of s
h=mp.mpf('1e-6')
f0=logPhi(x,p,D,0)
fp=logPhi(x,p,D,h); fm=logPhi(x,p,D,-h)
# mu_x = -d/ds Im? log Phi'(0) = i * mu (since Phi=E e^{isK}), so log Phi(s)= i mu s - v s^2/2+...
# log Phi'(0) = i mu ; log Phi''(0) = -v
d1=(fp-fm)/(2*h)        # = i mu
d2=(fp-2*f0+fm)/(h*h)   # = -v
mu=float((d1/1j).real); v=float((-d2).real)
print(f"n={n}  mu_x={mu:.4f}  v_x={v:.4f}  v/n={v/n:.4f}  (vbar pred ~ {D*(1-D)*(1-((1-2*p)**2/(p*p*(1-p)**2))*D*(1-D)/(1-2*D)**2):.4f})")

t=math.floor(n*D)
beta=(t-mu)/math.sqrt(v)
s0=(mu-t)/v   # committed z0 = i*s0 shift; s0 = -beta/sqrt(v)
print(f"t={t}  beta={beta:.4f}  s0=(mu-t)/v={s0:.5f}  |s0|=|beta|/sqrt(v)={abs(beta)/math.sqrt(v):.5f}")

# RISK (ii): does the saddle shift s -> s + i s0 kill the linear term?
# Define psi(s) = log Phi(s + i s0) - i t (s + i s0).  Its derivative at s=0 should ~0.
def psi_deriv0():
    g = lambda s: logPhi(x,p,D, s + 1j*s0) - 1j*t*(s+1j*s0)
    return (g(h)-g(-h))/(2*h)
dpsi=psi_deriv0()
print(f"RISK(ii) linear term after shift: |psi'(0)| = {abs(dpsi):.3e}  "
      f"(raw |logPhi'(0)-i t| = {abs(d1-1j*t):.3e}; ratio killed = {abs(dpsi)/abs(d1-1j*t):.2e})")
# the residual linear term should be O(C3 n s0^2) = O(third deriv * s0^2)
d3=( logPhi(x,p,D,2*h)-2*logPhi(x,p,D,h)+2*logPhi(x,p,D,-h)-logPhi(x,p,D,-2*h) )/(2*h**3)
print(f"  third cumulant log Phi'''(0): {complex(d3)}  |.|/n = {abs(d3)/n:.4f}  (C3 ~ this/n)")
print(f"  predicted residual linear ~ |d3|*s0^2/2 = {abs(d3)*s0*s0/2:.3e}")

# RISK (i): strip analyticity. Check |eta(z)| stays < r_KP on a vertical segment Im up to |s0|max ~ A ln n/sqrt n region.
# Here just verify |eta| at the shifted contour endpoints is small/inside.
rho=p/(1-p); B=((1-p)/p)**2
r_lb=p**2/(12*math.e*(1-p)**2)
# scan complex z = s + i*s0 for s in [-delta, delta]
A=D*(1-D)/(1-2*D)
delta=0.05
print(f"RISK(i) strip: r_lb(KP analytic lower)={r_lb:.5f}, A(D)={A:.4f}, |s0|={abs(s0):.5f}")
maxeta=0.0
for s in np.linspace(-delta,delta,21):
    z=s+1j*s0
    eta,_=eta_C_mp(p,D,z)
    maxeta=max(maxeta,float(abs(eta)))
print(f"  max|eta| on shifted window |s|<={delta}, Im=s0: {maxeta:.5f}  < r_lb? {maxeta<r_lb}  "
      f"(4A|z| bound at edge: {4*A*math.hypot(delta,abs(s0)):.5f})")
