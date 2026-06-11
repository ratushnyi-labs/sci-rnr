import math, numpy as np
# Verify the cubic-remainder integral that yields the (1+|beta|)^3/sqrt(n) form.
# On the shifted contour the integrand is e^{-v s^2/2} * e^{E(s)}, with
#   E(s) = (1/6) Phi'''(saddle) (i s)^3 + linear_residual*(i s) + higher,
# linear_residual = O(C3 n s0^2), |Phi'''| <= C3 n.  Using |e^E - 1| <= |E| e^{|E|}:
#   contribution = integral |E| e^{|E|} e^{-v s^2/2} ds.
# The cubic term integrates: int |s|^3 e^{-v s^2/2} ds = 2/v^2; times C3 n /6 = (C3 n)/(3 v^2).
# With v ~ vbar n: (C3 n)/(3 (vbar n)^2) = C3/(3 vbar^2 n).  Times the 1/sqrt(2 pi v) normalization
# the *relative* remainder R_n ~ C3/(3 vbar^2) * (1/n) * sqrt(v) ... plus the |beta|^3 from the
# linear residual on the vertical end.  Track the (1+|beta|)^3 scaling numerically by sweeping beta.
def remainder_model(C3n, v, s0, smax):
    # numerically integrate |E| e^{|E|} e^{-v s^2/2} over [-smax,smax] on the shifted line,
    # E = (1/6) C3n (|s|+|s0|)^3 (worst-case modulus bound on the cubic+residual)
    s=np.linspace(-smax,smax,4000)
    E=(1/6.0)*C3n*(np.abs(s)+abs(s0))**3
    integ=np.abs(E)*np.exp(np.minimum(np.abs(E),20))*np.exp(-v*s*s/2)
    num=np.trapezoid(integ,s)
    den=np.trapezoid(np.exp(-v*s*s/2),s)   # the Gaussian normalizer
    return num/den

# scale test: fix vbar, let n grow; beta fixed -> s0 ~ beta/sqrt(v) ~ beta/sqrt(vbar n)
vbar=0.09; C3=0.07
print(f"{'n':>7} {'beta':>5} | {'R_model':>10} {'R*sqrt(n)/(1+|b|)^3':>20}")
for beta in (0.0, 1.0, 3.0):
    for n in (256,1024,4096,16384):
        v=vbar*n; C3n=C3*n; s0=abs(beta)/math.sqrt(v); smax=6/math.sqrt(v)
        R=remainder_model(C3n,v,s0,smax)
        print(f"{n:>7} {beta:>5.1f} | {R:>10.5f} {R*math.sqrt(n)/(1+abs(beta))**3:>20.5f}")
    print()
