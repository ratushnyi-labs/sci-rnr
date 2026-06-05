#!/usr/bin/env python3
"""
Verification for the gamma_SGP sharpening (symbol-count smallest-grammar inapproximability):
8569/8568 (Charikar et al. 2005, via Berman-Karpinski 145/144) can be sharpened by substituting
the stronger Chlebik-Chlebikova bounded-degree Vertex-Cover hardness into the SAME Charikar
reduction (optimum grammar size OPT_SGP(sigma_G) = 15|V| + 3|E| + tau(G), Charikar et al. Thm 1).

For a d-regular hard VC family with inapproximability factor 1+delta and tau >= |V|/2,
using |E| = d|V|/2:
   gamma_SGP >= (15|V| + 3(d|V|/2) + (1+delta)(|V|/2)) / (15|V| + 3(d|V|/2) + (|V|/2))
             = 1 + delta/(31 + 3d).
Checks:
  V1 the gap formula gamma = 1 + delta/(31+3d) (symbolic, from OPT=15|V|+3|E|+tau, |E|=d|V|/2, tau=|V|/2).
  V2 3-regular CC VC (100/99, delta=1/99): gamma = 3961/3960 (DIRECTLY valid: 3-regular subset of
     Charikar's max-degree-3 + |E|=1.5|V|>=|V| input).
  V3 4-regular CC VC (53/52, delta=1/52): gamma = 2237/2236 (CONDITIONAL on Charikar's exact-size
     lemma extending to simple 4-regular graphs).
  V4 both EXCEED the baseline 8569/8568 (i.e. 1/3960 > 1/8568 and 1/2236 > 1/8568).
  V5 reconciliation: the baseline 8569/8568 = 1 + delta/(2*59.5) with delta=1/144 corresponds to a
     DENSER instance (c_E=|E|/|V|=4.75) than the CC sparse families; the formula is consistent.
"""
from fractions import Fraction as F

def gamma_d_regular(delta, d):
    # gamma - 1 = (delta/2) / (15 + 3d/2 + 1/2) = delta/(31+3d)
    return 1 + delta/F(31+3*d)

print("="*68)
print("V1: gap formula gamma = 1 + delta/(31+3d) from OPT=15|V|+3|E|+tau, d-regular")
# symbolic check: (15 + 3d/2 + (1+delta)/2) / (15 + 3d/2 + 1/2) - 1 == delta/(31+3d)
for d in (3,4,5):
    for delta in (F(1,99), F(1,52), F(1,144)):
        num = F(15) + 3*F(d,2) + (1+delta)/2
        den = F(15) + 3*F(d,2) + F(1,2)
        lhs = num/den
        rhs = gamma_d_regular(delta, d)
        assert lhs == rhs, (d, delta, lhs, rhs)
print("  formula verified symbolically for d in {3,4,5}, delta in {1/99,1/52,1/144}: PASS")
v1 = True

print("="*68)
print("V2: 3-regular CC VC (100/99): gamma_SGP (directly valid)")
g3 = gamma_d_regular(F(1,99), 3)
print(f"  gamma = 1 + (1/99)/40 = {g3} = {float(g3):.7f}")
v2 = (g3 == F(3961,3960))
print(f"  == 3961/3960 ? {v2}")

print("="*68)
print("V3: 4-regular CC VC (53/52): gamma_SGP (conditional on degree-4 lemma)")
g4 = gamma_d_regular(F(1,52), 4)
print(f"  gamma = 1 + (1/52)/43 = {g4} = {float(g4):.7f}")
v3 = (g4 == F(2237,2236))
print(f"  == 2237/2236 ? {v3}")

print("="*68)
print("V4: both exceed the baseline 8569/8568")
base = F(8569,8568)
print(f"  8569/8568 = {float(base):.7f}; 3961/3960 = {float(g3):.7f}; 2237/2236 = {float(g4):.7f}")
v4 = (g3 > base) and (g4 > base)
print(f"  3961/3960 > 8569/8568 and 2237/2236 > 8569/8568 ? {v4}")

print("="*68)
print("V5: reconciliation via the paper's EXACT extreme-point method rho=(15+3cE+(1+d)cT)/(15+3cE+cT)")
def rho_extreme(delta, cE, cT):   # the paper's adversarial-family ratio at (|E|/|V|=cE, k/|V|=cT)
    return (15 + 3*cE + (1+delta)*cT) / (15 + 3*cE + cT)
# baseline (paper lines 5854-5855): Berman-Karpinski 145/144, max-degree-3 |E|>=|V| extreme:
#   c_E = 3/2 (max edges), c_T = 1/3 (k>=|V|/3 from |E|>=|V|, each vertex covers <=3 edges).
base_chk = rho_extreme(F(1,144), F(3,2), F(1,3))
print(f"  baseline (delta=1/144, c_E=3/2, c_T=1/3) = {base_chk} = {float(base_chk):.7f}  == 8569/8568 ? {base_chk==base}")
# improvement: Chlebik-Chlebikova 3-regular VC 100/99; 3-regular forces c_E=3/2 and k>=|E|/3=|V|/2 => c_T=1/2.
imp3 = rho_extreme(F(1,99), F(3,2), F(1,2))
print(f"  CC 3-regular (delta=1/99, c_E=3/2, c_T=1/2) = {imp3} = {float(imp3):.7f}  == 3961/3960 ? {imp3==F(3961,3960)}")
# 4-regular CC 53/52: c_E=2, k>=|E|/4=|V|/2 => c_T=1/2 (conditional on Charikar's lemma at degree 4).
imp4 = rho_extreme(F(1,52), F(2), F(1,2))
print(f"  CC 4-regular (delta=1/52, c_E=2,   c_T=1/2) = {imp4} = {float(imp4):.7f}  == 2237/2236 ? {imp4==F(2237,2236)}")
v5 = (base_chk==base) and (imp3==F(3961,3960)) and (imp4==F(2237,2236)) and (imp3>base) and (imp4>base)
print(f"  baseline reproduces 8569/8568, CC families give larger valid bounds (3961/3960, 2237/2236): {v5}")

print("="*68)
allok = v1 and v2 and v3 and v4 and v5
print(f"RESULT: {'ALL PASS' if allok else 'SOME FAILED'}  [V1 {v1}, V2 {v2}, V3 {v3}, V4 {v4}, V5 {v5}]")
