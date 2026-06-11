#!/usr/bin/env python3
"""
ADVERSARIAL (2)+(3): saddle residual-linear and cubic-error order arithmetic,
re-derived INDEPENDENTLY (analytic, symbolic-ish), confirming/refuting the draft.

(2) The s0-shift residual linear term:
    Lambda'(0) = -1/2 psi'''(xi) s0^2,  |.| <= 1/2 C3 n s0^2.
    s0 = -beta/sqrt(v_x), v_x ~ vbar n => s0^2 = beta^2/(vbar n).
    => |Lambda'(0)| <= 1/2 C3 n beta^2/(vbar n) = C3 beta^2/(2 vbar) = O(beta^2)  [n-INDEPENDENT!]
    Its CONTRIBUTION to the Gaussian integral: relative shift of the peak,
    exponent gain Lambda'(0)^2/(2 v_x') = O(beta^4)/(vbar n) = O(beta^4/n) -> 0.
    DRAFT claims |Lambda'(0)| = O(beta^2/sqrt n).  *** CHECK: is it /sqrt n or O(1)? ***
    Recompute: 1/2 C3 n s0^2 with s0^2 = beta^2/(vbar n):
       = 1/2 C3 n * beta^2/(vbar n) = C3 beta^2/(2 vbar)  -- this is O(beta^2), NOT O(beta^2/sqrt n)!
    The draft's "O(beta^2/sqrt n)" REQUIRES an extra 1/sqrt n that is NOT present
    in 1/2 C3 n s0^2.  Investigate: is the draft using s0 = O(1/n) (not 1/sqrt n)?
    -> s0 = (mu-t)/v.  On the shell, mu-t = O(sqrt(v)) = O(sqrt n) (beta=O(1)), v=O(n),
       so s0 = O(sqrt n)/O(n) = O(1/sqrt n).  So s0^2 = O(1/n), n s0^2 = O(1). CONFIRMS
       |Lambda'(0)| = O(1), NOT O(1/sqrt n).  The draft's extra 1/sqrt n is an ERROR,
       but it is HARMLESS: the contribution to log q_t is Lambda'(0)^2/(2v) = O(1)/O(n)
       = O(1/n), subdominant either way.

(3) Cubic remainder rho_1 = O(C3/(vbar^{3/2} sqrt n)).  Re-derive:
    rho_1 ~ (C3 n/6) * E|s|^3 with s~N(0,1/v_x), E|s|^3 = 2 sqrt(2/pi) v_x^{-3/2}.
    = (C3 n/6)(2 sqrt(2/pi))(vbar n)^{-3/2} = C3/(3 sqrt(pi/2)) * vbar^{-3/2} * n^{-1/2}.
    CONFIRMS O(1/sqrt n), polylog-free.  The |beta|^3 enters via the shifted-saddle
    cubic at the displaced center, |s+s0|^3, giving (1+|beta|)^3.
"""
import math


def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))


def vbar_closed(p, D):
    kappa2 = (1 - 2 * p) ** 2 / (p ** 2 * (1 - p) ** 2)
    Psi0 = kappa2 * D * (1 - D) / (1 - 2 * D) ** 2
    return D * (1 - D) * (1 - Psi0)


if __name__ == "__main__":
    print("=" * 80)
    print("(2) residual linear term |Lambda'(0)| = 1/2 |psi'''| s0^2  order arithmetic")
    print("=" * 80)
    p = 0.4; D = 0.9 * Dc(p); vbar = vbar_closed(p, D)
    C3 = 0.10
    print(f"p={p} D=0.9Dc={D:.5f} vbar={vbar:.4f} (taking C3~{C3})")
    print(f"{'n':>8} {'beta':>5} {'s0^2':>12} {'n s0^2':>10} {'|Lam(0)|=C3 n s0^2/2':>22} "
          f"{'draft O(b^2/sqrt n)':>20} {'contrib to logq=Lam^2/2v':>26}")
    for n in (256, 1024, 4096, 16384, 65536):
        for beta in (1.0, 3.0):
            v = vbar * n
            s0sq = beta ** 2 / v
            lam = 0.5 * C3 * n * s0sq           # = C3 beta^2/(2 vbar) -- n-INDEPENDENT
            draft = beta ** 2 / math.sqrt(n)    # the draft's claimed scaling (x const)
            contrib = lam ** 2 / (2 * v)        # gain to log q_t
            print(f"{n:>8} {beta:>5.1f} {s0sq:>12.3e} {n*s0sq:>10.4f} {lam:>22.5f} "
                  f"{draft:>20.5f} {contrib:>26.3e}")
    print("\n  FINDING: |Lambda'(0)| = C3 beta^2/(2 vbar) is n-INDEPENDENT (O(beta^2)),")
    print("  NOT O(beta^2/sqrt n) as the draft states.  HARMLESS: its contribution to")
    print("  log q_t is Lambda'(0)^2/(2 v_x) = O(beta^4/n) -> 0 regardless.")

    print("\n" + "=" * 80)
    print("(3) cubic remainder rho_1 order arithmetic (independent)")
    print("=" * 80)
    print(f"{'n':>8} {'rho_1 closed':>14} {'rho_1*sqrt n':>14}")
    for n in (256, 1024, 4096, 16384, 65536):
        v = vbar * n
        Eabs3 = 2 * math.sqrt(2 / math.pi) * v ** (-1.5)
        rho = 2.0 * (C3 * n / 6.0) * Eabs3   # e^{|E|}<=2 factor
        print(f"{n:>8} {rho:>14.6e} {rho*math.sqrt(n):>14.6f}")
    print("\n  FINDING: rho_1 ~ const * n^{-1/2}, polylog-FREE.  rho_1*sqrt n -> const")
    print(f"  (= 2 C3/(3 sqrt(pi/2) vbar^{{3/2}}) = {2*C3/(3*math.sqrt(math.pi/2)*vbar**1.5):.4f}).")
    print("  The (1+|beta|)^3 enters from the displaced cubic |s+s0|^3; the (log n)^3 in")
    print("  the LEMMA is a LOOSE majorant (matches the draft's own honesty note).")
