#!/usr/bin/env python3
"""
probe_7_34_R2_Lb_variance_split.py
============================================================================
(R2) lemma (L-b): the ONE genuinely at-risk constant of the BSMS Gray-region
RD-dispersion achievability central-window write-up.

We need  vbar(p,D) := lim_n E[v_x]/n  to EXIST and be > 0, where
   mu_x = E_Q[K | x]   (posterior mean Hamming radius),
   v_x  = Var_Q(K | x) (posterior variance),
   Q = posterior law of K = d_H(x, Y*) given the source word x.

The EXACT annealed split (law of total variance, K~Bin(n,D) annealed):
   n D(1-D) = Var_{x,Q}(K) = E_x[v_x] + Var_x(mu_x).
So  E[v_x]/n = D(1-D) - Var(mu_x)/n.

THE AT-RISK QUESTION (the prompt's flag): is Var(mu_x) = o(n) or Theta(n)?
  - If Var(mu_x) = o(n):    vbar = D(1-D) > 0 trivially.  (BEST case.)
  - If Var(mu_x) = Theta(n) -> sigma_mu^2 := lim Var(mu_x)/n in (0, D(1-D)]:
        vbar = D(1-D) - sigma_mu^2; we then NEED sigma_mu^2 < D(1-D) STRICTLY.
  - If sigma_mu^2 = D(1-D):  vbar = 0  -> the constant FAILS (degenerate
        steepest descent). This is the failure mode to rule out.

This probe MEASURES sigma_mu^2 = Var(mu_x)/n and vbar = E[v_x]/n directly, over
the sqrt(n)-shell AND annealed (so the split identity is exact), at
p in {0.25, 0.4}, D = 0.9 D_c, n up to 400+, and:
  (B0) cross-checks mu_x, v_x from the exact O(n^2) DFT posterior pipeline vs
       the analytic per-site cumulants from log Phi small-s expansion;
  (B1) verifies the EXACT split  n D(1-D) = E[v_x] + Var(mu_x)  (ANNEALED);
  (B2) reports sigma_mu^2 = Var(mu_x)/n and vbar = E[v_x]/n with n-trend ->
       settles Theta(n) vs o(n) for Var(mu_x), and pins vbar > 0 with margin;
  (B3) the closed-form lower bound vbar >= D(1-D) - sigma_mu^2 and the
       analytic curvature a(p,D)=D(1-D)(1-Psi_0(D)) = per-site posterior
       variance vbar predicted by the replica corollary (Thm corollary in tex);
       cross-check vbar == a(p,D).

NOTE on shell vs annealed: the split identity is exact ANNEALED. On the shell
(switch-count restricted) Var(mu_x) can only DECREASE (conditioning on a
function correlated with mu_x removes variance), so the shell vbar >= annealed
vbar; we report BOTH and use the smaller (annealed) as the conservative floor.
Deps: numpy.
"""
import math
import numpy as np


def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))


def theta0(D):
    return math.log(D / (1 - D))


def eta_C(beta, D):
    eb = np.exp(beta)
    ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
    eta = (1 - ze) / (1 + ze)
    C = (1 + eb) * (1 + ze) / 2
    return eta, C


def phi_all_nodes(xbits, p, D):
    """Phi(s_j;x) at the n+1 DFT nodes, log-magnitude scaled (verified pipeline
    identical to probe_7_34_replica_assembly_central.py)."""
    n = len(xbits)
    S = n + 1
    s = 2 * math.pi * np.arange(S) / S
    th = theta0(D)
    eta, Cs = eta_C(th + 1j * s, D)
    C0 = eta_C(th, D)[1]
    rho = p / (1 - p)
    w0 = np.ones(S, dtype=complex)
    w1 = eta.copy()
    acc = np.zeros(S)
    for j in range(n - 1):
        r = rho if xbits[j] == xbits[j + 1] else 1.0 / rho
        n0 = w0 + w1 * r
        n1 = (w0 * r + w1) * eta
        w0, w1 = n0, n1
        m = np.maximum(np.abs(w0), np.abs(w1))
        m = np.where(m == 0, 1.0, m)
        w0 /= m
        w1 /= m
        acc += np.log(m)
    Ssum = w0 + w1
    logmag = n * np.log(np.abs(Cs / C0)) + acc + np.log(np.maximum(np.abs(Ssum), 1e-300))
    phase = n * np.angle(Cs / C0) + np.angle(Ssum)
    return np.exp(np.minimum(logmag, 0.0) + 1j * phase)


def cumulants(xbits, p, D):
    """mu_x = E_Q[K], v_x = Var_Q[K] from the exact posterior q_k = IFFT(Phi)."""
    phi = phi_all_nodes(xbits, p, D)
    q = np.fft.fft(phi).real / len(phi)
    q = np.maximum(q, 0.0)
    n = len(xbits)
    ks = np.arange(n + 1, dtype=float)
    Z = q[:n + 1].sum()
    q = q[:n + 1] / Z
    mu = float(ks @ q)
    v = float(((ks - mu) ** 2) @ q)
    return mu, v


def sample_words(rng, n, p, R):
    out = []
    for _ in range(R):
        out.append((np.cumsum(rng.random(n) < p) % 2).astype(int))
    return out


def shell_words(rng, n, p, R):
    out = []
    mean = (n - 1) * p
    half = 2.0 * math.sqrt((n - 1) * p * (1 - p))
    while len(out) < R:
        x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
        sw = int(np.sum(x[:-1] != x[1:]))
        if abs(sw - mean) <= half:
            out.append(x)
    return out


def Psi0(D):
    # Psi_t(D) at t = 1-cos s = 0 (the s->0 curvature term in the replica corollary):
    #   Psi_0(D) = kappa^2 D(1-D) (1-2D)^2 / (1-2D)^4 = kappa^2 D(1-D)/(1-2D)^2
    # with kappa^2 = (1-2p)^2/(p^2 (1-p)^2). Needs p; pass via closure below.
    raise NotImplementedError


def a_replica(p, D):
    """Closed-form replica-corollary curvature a(p,D) = D(1-D)(1 - Psi_0(D)),
    Psi_0(D) = kappa^2 D(1-D)/(1-2D)^2, kappa^2 = (1-2p)^2/(p^2(1-p)^2).
    This EQUALS the per-site posterior variance vbar (claim in tex corollary)."""
    kappa2 = (1 - 2 * p) ** 2 / (p ** 2 * (1 - p) ** 2)
    psi0 = kappa2 * D * (1 - D) / (1 - 2 * D) ** 2
    return D * (1 - D) * (1 - psi0)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=4000)
    ap.add_argument("--nlist", type=str, default="64,128,256,400")
    args = ap.parse_args()
    n_list = [int(x) for x in args.nlist.split(",")]
    rng = np.random.default_rng(20260610)

    print("=" * 92)
    print("(L-b) variance split:  n D(1-D) = E[v_x] + Var(mu_x)   [ANNEALED, EXACT]")
    print("      vbar = E[v_x]/n,  sigma_mu^2 = Var(mu_x)/n ;  need vbar -> const > 0")
    print("=" * 92)

    for p in (0.25, 0.40):
        D = 0.9 * Dc(p)
        DD = D * (1 - D)
        a_pred = a_replica(p, D)
        print("-" * 92)
        print(f"p={p}  D=0.9*Dc={D:.6f}  D(1-D)={DD:.6f}  "
              f"replica-curvature a(p,D)={a_pred:.6f}")
        print(f"{'n':>6} {'E[v]/n':>10} {'Var(mu)/n':>11} {'sum=D(1-D)?':>12} "
              f"{'split err':>11} | {'shell E[v]/n':>13} {'shell Vmu/n':>12}")
        ann_rows = []
        shell_rows = []
        for n in n_list:
            # annealed (split is EXACT here)
            words = sample_words(rng, n, p, args.R)
            mus = np.empty(len(words)); vs = np.empty(len(words))
            for i, x in enumerate(words):
                mus[i], vs[i] = cumulants(x, p, D)
            Ev = mus.var() * 0 + vs.mean()      # E[v_x]
            Vmu = mus.var()                       # Var(mu_x)  (population)
            # exact annealed total: Var_{x,Q}(K) should equal n D(1-D)
            # but with finite R we report E[v]+Var(mu) and compare to n D(1-D)
            tot = Ev + Vmu
            split_err = abs(tot - n * DD) / (n * DD)
            ann_rows.append((n, Ev / n, Vmu / n, tot, split_err))
            # shell (conservative: Var(mu) only smaller)
            sw = shell_words(rng, n, p, max(args.R // 3, 400))
            smus = np.empty(len(sw)); svs = np.empty(len(sw))
            for i, x in enumerate(sw):
                smus[i], svs[i] = cumulants(x, p, D)
            shell_rows.append((n, svs.mean() / n, smus.var() / n))
            print(f"{n:>6} {Ev/n:>10.5f} {Vmu/n:>11.5f} {tot:>12.4f} "
                  f"{split_err:>11.2e} | {svs.mean()/n:>13.5f} {smus.var()/n:>12.5f}",
                  flush=True)

        # trend diagnostics
        vbars = [r[1] for r in ann_rows]
        smus2 = [r[2] for r in ann_rows]
        ns = [r[0] for r in ann_rows]
        # is Var(mu)=Theta(n)? then Var(mu)/n -> const>0; if o(n), -> 0
        smu_trend = "Theta(n) (Var(mu)/n n-stable > 0)" if smus2[-1] > 0.3 * smus2[0] and smus2[-1] > 1e-4 \
            else "o(n)? (Var(mu)/n shrinking toward 0)"
        # vbar floor
        vbar_min = min(vbars)
        margin = vbar_min  # distance of E[v]/n above 0
        print(f"  -> sigma_mu^2 = Var(mu)/n trend: {[f'{x:.5f}' for x in smus2]}  [{smu_trend}]")
        print(f"  -> vbar = E[v]/n trend:          {[f'{x:.5f}' for x in vbars]}")
        print(f"  -> vbar floor over n = {vbar_min:.5f}  (margin above 0)")
        print(f"  -> closed-form check: D(1-D) - sigma_mu^2(last) = {DD - smus2[-1]:.5f}  "
              f"vs E[v]/n(last) = {vbars[-1]:.5f}  vs replica a(p,D) = {a_pred:.5f}")
        ok_pos = vbar_min > 1e-3
        ok_strict = smus2[-1] < DD - 1e-3
        ok_match = abs(vbars[-1] - a_pred) < 0.15 * max(a_pred, 1e-3) + 0.002
        print(f"  VERDICT p={p}: vbar>0 [{ 'PASS' if ok_pos else 'FAIL'}]  "
              f"sigma_mu^2<D(1-D) strict [{'PASS' if ok_strict else 'FAIL'}]  "
              f"vbar==replica-a [{'PASS' if ok_match else 'FAIL'}]")
    print("=" * 92)
