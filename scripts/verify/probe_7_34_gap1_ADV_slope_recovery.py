#!/usr/bin/env python3
"""
ADVERSARIAL (1c): is the narrow-window pathology purely an artifact of the LOOSE
4A complex-slope constant?  The draft sets y_max = (1/2) r_KP/(4A) and requires
eps <= r_KP/(8A) = y_max.  But attack(1) showed the REAL complex slope is ~A
(not 4A) on the relevant tiny strip (|Im z| <= |s0| = O(1/sqrt n)).  With the true
slope, the strip-containment |eta_z|<=r_KP permits eps up to ~ r_KP/A - |s0| ~ delta_KP.
Then crossover n ~ 1/(vbar delta_KP^2) is far smaller.

We DIRECTLY measure the largest eps such that |eta(s+i s0)| <= r_KP for s in [-eps,eps]
at the actual shell |s0|=O(1/sqrt n), and report the resulting #peak-widths and crossover.
"""
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 40


def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))


def eta_closed(p, D, s):
    w = mp.e ** (1j * s)
    return D * (1 - D) * (w - 1) / ((1 - D) ** 2 - D ** 2 * w)


if __name__ == "__main__":
    print("=" * 90)
    print("Max eps with |eta(s+i s0)|<=r_KP on the ACTUAL shifted strip (|s0|=beta/sqrt(vbar n))")
    print("=" * 90)
    for p in (0.4, 0.25, 0.1):
        D = 0.9 * Dc(p)
        A = D * (1 - D) / (1 - 2 * D)
        r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
        kappa2 = (1 - 2 * p) ** 2 / (p ** 2 * (1 - p) ** 2)
        vbar = D * (1 - D) * (1 - kappa2 * D * (1 - D) / (1 - 2 * D) ** 2)
        eps_draft = r_lb / (8 * A)
        delta_KP = r_lb / A
        print(f"\np={p} D=0.9Dc={D:.5f}: A={A:.4f} r_KP>={r_lb:.5f} vbar={vbar:.4f}")
        print(f"  draft eps=r_KP/(8A)={eps_draft:.5f}   delta_KP=r_KP/A={delta_KP:.5f}")
        for n in (10**3, 10**4, 10**5):
            beta = 2.0  # typical shell
            s0 = beta / math.sqrt(vbar * n)
            # find max eps: scan s up, check |eta(s+i s0)| <= r_lb
            eps_max = 0.0
            for s in np.linspace(0, delta_KP * 1.2, 600):
                e = float(abs(eta_closed(p, D, s + 1j * s0)))
                if e <= r_lb:
                    eps_max = s
                else:
                    break
            pw = 1.0 / math.sqrt(vbar * n)
            nwidths_draft = eps_draft / pw
            nwidths_true = eps_max / pw
            cross_draft = 36 / (vbar * eps_draft ** 2)
            cross_true = 36 / (vbar * eps_max ** 2) if eps_max > 0 else float('inf')
            print(f"   n={n:>7}: |s0|={s0:.5f}  eps_max(true slope)={eps_max:.5f}  "
                  f"(vs draft {eps_draft:.5f})  #widths: true={nwidths_true:.1f} draft={nwidths_draft:.1f}")
        print(f"   crossover n (6 widths): TRUE-slope ~{36/(vbar*eps_max**2):.2e}  "
              f"DRAFT ~{36/(vbar*eps_draft**2):.2e}  (improvement ~{(eps_max/eps_draft)**2:.0f}x)")
    print("\nFINDING: with the TRUE complex slope (~A, attack 1), eps reaches ~delta_KP,")
    print("shrinking the window-vs-peak crossover by ~64x.  The draft's narrow window is")
    print("an artifact of the loose 4A constant; the architecture is repairable to a")
    print("reasonable n0, but the SUBMITTED draft (eps=r_KP/8A) has the huge-n0 pathology.")
