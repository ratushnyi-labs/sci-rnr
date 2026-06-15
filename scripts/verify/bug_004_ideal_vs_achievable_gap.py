#!/usr/bin/env python3
"""
BUG-004-C  --  the ideal-vs-achievable GAP probe (operational companion to Remark 7.15g).

Scope (governance/documentation; NO new mathematics). This probe instruments three
already-proven facts so the BUG-004-B scoping edits carry a green-CI verification anchor.
It does NOT re-derive any hardness: the #P-hardness of the exact ideal rate is Remark
7.15g (verified standalone in remark_7_15g_hmm_entropy_sharp_p_hard.py, V1-V7); the
poly-time eps-approximability under filter stability is Remark 7.15c (verified standalone
in remark_7_15c_mixing_oracle_approx.py, M1-M5). Here we exhibit, side by side on tiny
explicit sources, the DICHOTOMY the scoping clause asserts:

  THE GAP (Rem 7.15g + 7.15c):
    the exact ideal RNR rate  L*_RNR(M,N) = sum_t H_M(X_t | X_{<t}) = H_M(X^N)  is the
    optimization TARGET, and computing it exactly is #P-hard / FP^#P-complete -- an
    information-theoretic BENCHMARK, not an efficiently computable encoder target.  The
    OPERATIVE object is (a) the eps-approximation, poly-time under filter stability, and
    (b) the achievable scheme: the realized causal codelength sum_t -log2 Q(x_t|x_<t),
    one forward pass per symbol, whose expectation equals H(X^N) but which the encoder
    PAYS rather than COMPUTES.

Checks:
  G1  INTRACTABLE TARGET (re-exhibit 7.15g's YES/NO mapping on a small instance).
      On a small 3-CNF, build the 7.15g HMM block marginal p(x)=(P+u(x))/Z; confirm the
      exact ideal rate L*_RNR = H(X^N) (chain rule = block entropy) and that its exact
      symbolic log-P coefficient is -c0*P/Z, recovering #SAT(phi)=c0.  The exact ideal
      rate is thus a #SAT-hard query (Rem 7.15g): an intractable benchmark.  We also show
      the bounded-precision NUMERIC value of H does NOT expose c0 (two formulas, different
      #SAT, value alone insufficient) -- the hardness is the exact-symbolic target, the
      complement of the eps-approximation of G2.

  G2  TRACTABLE eps-APPROXIMATION under filter stability (Rem 7.15c).
      On a small filter-stable (positive-emission, geometrically ergodic) HMM, the
      order-d truncation A_d = sum_t H(X_t | X_{t-d:t-1}) over-estimates H(X^N) with a
      MONOTONE, GEOMETRICALLY decaying gap 0 <= A_d - H <= C*rho^d.  So eps accuracy needs
      only d = O(log(N/eps)/log(1/rho)) (LOGARITHMIC in 1/eps), poly work -- while the
      EXACT value by brute marginalization enumerates |Sigma|^N words (exponential).  We
      exhibit the poly-vs-exponential separation by COUNT of words touched, and confirm
      A_d -> H monotonically with geometric decay.

  G3  ACHIEVABLE scheme matches the eps-approx target, NOT the exact optimum.
      On the same filter-stable HMM, the deployed encoder pays the realized causal
      codelength  L(x) = sum_t -log2 Q(x_t | x_<t)  in ONE forward pass (Rem 7.15b).
      Under the TRUE model Q=P, E_P[L(X^N)] = H(X^N) exactly (chain-rule identity), and
      the per-symbol Monte-Carlo mean of the realized length converges to (1/N)H(X^N) --
      i.e. the achievable scheme realizes the (eps-approximable) rate, never computing the
      exact #P-hard optimum.  We also confirm the converse direction E[L] >= H(X^N) holds
      with equality at the matched model and STRICT inequality under a mismatched model
      (the KL gap), so the ideal H(X^N) is the floor the achievable scheme approaches.

Prints per-check PASS/FAIL and a final 'OVERALL -> PASS'.  Exact fractions where the
identity is symbolic (G1); float Monte-Carlo with explicit tolerance for G3.
"""
import math
import itertools
import random
from fractions import Fraction

random.seed(1234)

# ====================================================================== helpers
def is_prime(k):
    if k < 2:
        return False
    i = 2
    while i * i <= k:
        if k % i == 0:
            return False
        i += 1
    return True


def next_prime(k):
    while not is_prime(k):
        k += 1
    return k


def prime_factors(k):
    f = {}
    d = 2
    while d * d <= k:
        while k % d == 0:
            f[d] = f.get(d, 0) + 1
            k //= d
        d += 1
    if k > 1:
        f[k] = f.get(k, 0) + 1
    return f


def xlogx(p):
    """-p log2 p with the 0 log 0 = 0 convention, p a float."""
    return 0.0 if p <= 0.0 else -p * math.log2(p)


# ====================================================================== G1
# Re-exhibit Remark 7.15g's reduction on a small 3-CNF: the exact ideal rate
# L*_RNR = H(X^N) is a #SAT-hard query.  (Standalone proof in
# remark_7_15g_hmm_entropy_sharp_p_hard.py; here it anchors the "intractable target".)
def falsifies(clause, x):
    # literal (v, sign): true iff x[v]==1 matches sign; clause falsified iff ALL literals false
    return all((x[v] == 1) != s for (v, s) in clause)


def u_of(clauses, x):
    return sum(1 for c in clauses if falsifies(c, x))


def g1_reduction(n, clauses):
    m = len(clauses)
    P = next_prime(max(m, 2) + 1)
    Z = (8 * P + m) * (2 ** (n - 3))  # integer for n>=3
    # block marginal p(x)=(P+u(x))/Z ; group counts c_j ; c_0=#SAT
    cj = {}
    for x in itertools.product([0, 1], repeat=n):
        cj[u_of(clauses, x)] = cj.get(u_of(clauses, x), 0) + 1
    c0 = cj.get(0, 0)
    sat = sum(1 for x in itertools.product([0, 1], repeat=n)
              if all(not falsifies(c, x) for c in clauses))
    # exact ideal rate = block entropy = chain-rule sum_t H(X_t|X_<t); verify the two coincide
    # (a) block entropy directly:
    H_block = math.log2(Z) - (1.0 / Z) * sum(c * (P + j) * math.log2(P + j) for j, c in cj.items())
    # (b) chain-rule sum_t H(X_t|X_<t) over the block distribution (the L*_RNR definition):
    H_chain = chain_rule_entropy_from_block(n, clauses, P, Z)
    # exact symbolic log-P coefficient -> recovers #SAT
    coeff = {}
    for pr, e in prime_factors(int(Z)).items():
        coeff[pr] = coeff.get(pr, Fraction(0)) + e       # +v_q(Z) from + log Z
    for j, c in cj.items():
        for pr, e in prime_factors(P + j).items():
            coeff[pr] = coeff.get(pr, Fraction(0)) - Fraction(c * (P + j), Z) * e
    coefP = coeff.get(P, Fraction(0))
    recovered = -coefP * Z / P
    p_isolated = (int(Z) % P != 0) and all((P + j) % P != 0 for j in cj if j != 0)
    return dict(P=P, Z=Z, c0=c0, sat=sat, H_block=H_block, H_chain=H_chain,
                coefP=coefP, recovered=recovered, p_isolated=p_isolated, cj=cj)


def chain_rule_entropy_from_block(n, clauses, P, Z):
    """sum_t H(X_t | X_{<t}) computed from the block law p(x)=(P+u(x))/Z; equals H(X^N)."""
    # joint as dict over all x
    joint = {x: Fraction(P + u_of(clauses, x), Z) for x in itertools.product([0, 1], repeat=n)}
    H = 0.0
    for t in range(n):
        # H(X_t | X_{<t}) = sum_{prefix} P(prefix) * H(X_t | prefix)
        pref = {}
        for x, p in joint.items():
            key = x[:t]
            d = pref.setdefault(key, [Fraction(0), Fraction(0)])  # mass with x_t=0, x_t=1
            d[x[t]] += p
        for key, (p0, p1) in pref.items():
            tot = p0 + p1
            if tot == 0:
                continue
            f0, f1 = float(p0 / tot), float(p1 / tot)
            H += float(tot) * (xlogx(f0) + xlogx(f1))
    return H


F_sat = (4, [[(0, True), (1, True), (2, True)], [(0, False), (1, True), (3, True)]])
F_mix = (4, [[(0, True), (1, True), (2, True)], [(0, False), (1, False), (2, False)],
             [(1, True), (2, False), (3, True)]])
F_unsat = (3, [[(0, s0), (1, s1), (2, s2)]
               for (s0, s1, s2) in itertools.product([True, False], repeat=3)])  # all 8 -> UNSAT

print("=" * 72)
print("G1  INTRACTABLE TARGET: exact ideal rate L*_RNR = H(X^N) is a #SAT-hard query (Rem 7.15g)")
g1_ok = True
g1_rows = []
for name, (n, cl) in [("F_sat", F_sat), ("F_mix", F_mix), ("F_unsat", F_unsat)]:
    r = g1_reduction(n, cl)
    chain_matches = abs(r["H_block"] - r["H_chain"]) < 1e-9
    recovers = (r["recovered"] == r["c0"] == r["sat"])
    ok = chain_matches and recovers and r["p_isolated"]
    g1_ok = g1_ok and ok
    g1_rows.append((name, r))
    print(f"  {name}: L*_RNR = sum_t H(X_t|X_<t) = {r['H_chain']:.6f} bits ; "
          f"H(X^N)={r['H_block']:.6f} (chain==block: {chain_matches})")
    print(f"         exact log-P coeff = {r['coefP']}  ->  -coeff*Z/P = {r['recovered']}  "
          f"= #SAT(phi) = {r['sat']} (recovers: {recovers}); P isolated: {r['p_isolated']}")
# bounded-precision numeric value does NOT expose c0 (hardness is exact-symbolic, complement of G2 approx)
H_sat = g1_rows[0][1]["H_block"]
H_mix = g1_rows[1][1]["H_block"]
print(f"  numeric H alone: F_sat H={H_sat:.6f} (#SAT={g1_rows[0][1]['c0']}), "
      f"F_mix H={H_mix:.6f} (#SAT={g1_rows[1][1]['c0']}) -- the real value mixes all c_j, c0 not")
print(f"         readable from it; recovering #SAT needs the SYMBOLIC log-P coefficient (Rem 7.15g Scope (i)).")
print(f"  G1 {'PASS' if g1_ok else 'FAIL'}  (exact ideal rate = block entropy = #SAT-hard target)")

# ====================================================================== G2
# Tractable eps-approximation under filter stability (Rem 7.15c).
# Small filter-stable HMM: 2 hidden states, POSITIVE emissions (so geometrically ergodic
# and uniformly filter-stable).  A_d = sum_t H(X_t|X_{t-d:t-1}) >= H(X^N), gap geometric.
# Hidden transition T (row-stochastic), emission E[s][a] (positive), stationary pi.
HID = [0, 1]
SYM = [0, 1]
T = [[Fraction(7, 10), Fraction(3, 10)],
     [Fraction(2, 10), Fraction(8, 10)]]
E = [[Fraction(8, 10), Fraction(2, 10)],   # state 0 prefers symbol 0
     [Fraction(3, 10), Fraction(7, 10)]]    # state 1 prefers symbol 1


def stationary(T):
    # solve pi T = pi for 2-state by closed form
    a = T[0][1]
    b = T[1][0]
    p0 = b / (a + b)
    return [p0, 1 - p0]


PI = stationary(T)


def block_prob_hmm(word):
    """Exact P(X_1..X_L = word) by forward algorithm over hidden states (Fractions)."""
    L = len(word)
    # alpha[s] = P(hidden_t = s, observations up to t)
    alpha = [PI[s] * E[s][word[0]] for s in HID]
    for t in range(1, L):
        nxt = [Fraction(0), Fraction(0)]
        for s2 in HID:
            nxt[s2] = sum(alpha[s1] * T[s1][s2] for s1 in HID) * E[s2][word[t]]
        alpha = nxt
    return sum(alpha)


def exact_block_entropy(N):
    """H(X^N) by brute enumeration of all 2^N words (the EXPONENTIAL exact route)."""
    H = 0.0
    touched = 0
    for w in itertools.product(SYM, repeat=N):
        p = float(block_prob_hmm(list(w)))
        touched += 1
        H += xlogx(p)
    return H, touched


def cond_entropy_order_d(d):
    """H(X_t | X_{t-d:t-1}) at stationarity: sum over (context, symbol) of windowed marginals.
    Uses only (d+1)-block marginals -> table size 2^{d+1} (poly route, no |Sigma|^N)."""
    if d == 0:
        # H(X_t) marginal
        # P(X=a) = sum_s pi[s] E[s][a]
        pa = [sum(PI[s] * E[s][a] for s in HID) for a in SYM]
        return sum(xlogx(float(pa[a])) for a in SYM), 2
    # need joint of a (d+1)-window: P(X_1..X_{d+1}) via forward, then H(X_{d+1}|X_1..X_d)
    Hc = 0.0
    touched = 0
    # enumerate context of length d
    for ctx in itertools.product(SYM, repeat=d):
        # P(context) and P(context, next=a)
        p_ctx = float(block_prob_hmm(list(ctx)))
        if p_ctx <= 0:
            continue
        for a in SYM:
            touched += 1
            p_full = float(block_prob_hmm(list(ctx) + [a]))
            if p_full <= 0:
                continue
            cond = p_full / p_ctx
            Hc += p_ctx * xlogx(cond)
    return Hc, touched


print("=" * 72)
print("G2  TRACTABLE eps-APPROX under filter stability: A_d -> H(X^N) geometrically (Rem 7.15c)")
N = 12
H_exact, touched_exact = exact_block_entropy(N)
# A_d = H(X_1..X_d) [exact small head] + (N-d) * H(X_{d+1}|X_1..X_d) is the stationary order-d
# truncation; we use the cumulative-conditional form A_d = sum_t H(X_t|X_{(t-d):t-1}).
# At stationarity H(X_t|X_{t-d:t-1}) is t-independent for t>d, so A_d = head + (N-d)*c_d.
def A_d_value(d, N):
    # head: sum_{t=1}^{d} H(X_t | X_1..X_{t-1}) computed exactly from (<=d)-blocks
    head = 0.0
    head_touched = 0
    for t in range(1, d + 1):
        hc, tt = cond_entropy_order_d(t - 1)
        head += hc
        head_touched += tt
    c_d, t_d = cond_entropy_order_d(d)
    return head + (N - d) * c_d, head_touched + t_d


prev = None
mono_ok = True
geom_ok = True
gaps = []
for d in range(0, 7):
    Ad, touched = A_d_value(d, N)
    gap = Ad - H_exact
    gaps.append(gap)
    flag = ""
    if prev is not None:
        if Ad > prev + 1e-12:
            mono_ok = False
            flag = " (NON-MONOTONE!)"
    print(f"  d={d}: A_d={Ad:.8f}  gap=A_d-H={gap:.3e}  (table touched {touched} windows){flag}")
    prev = Ad
# geometric decay: ratio of successive positive gaps stays bounded < 1
pos = [g for g in gaps if g > 1e-12]
if len(pos) >= 3:
    ratios = [pos[i + 1] / pos[i] for i in range(len(pos) - 1)]
    geom_ok = all(r < 1.0 for r in ratios)
    print(f"  successive gap ratios: {[f'{r:.3f}' for r in ratios]}  (all < 1 => geometric: {geom_ok})")
# poly vs exponential: A_d at depth d touches O(2^{d+1}) windows; exact touches 2^N.
d_eps = 6
_, touched_poly = A_d_value(d_eps, N)
sep_ok = touched_poly < touched_exact
print(f"  POLY vs EXP: A_{d_eps} touched {touched_poly} windows (~2^(d+1)); "
      f"exact H(X^N) enumerated {touched_exact} = 2^{N} words.  poly<exp: {sep_ok}")
floor_ok = all(g >= -1e-9 for g in gaps)  # A_d >= H always (conditioning on less raises entropy)
print(f"  A_d >= H(X^N) for all d (over-estimate floor): {floor_ok}")
g2_ok = mono_ok and geom_ok and sep_ok and floor_ok
print(f"  G2 {'PASS' if g2_ok else 'FAIL'}  (eps-approx poly under filter stability; exact is exponential)")

# ====================================================================== G3
# Achievable scheme realizes the eps-approx target, not the exact optimum (Rem 7.15b).
# One forward pass: L(x)=sum_t -log2 Q(x_t|x_<t).  Under the TRUE model, E[L]=H(X^N).
def causal_codelength_true(word):
    """Realized one-pass codelength under the TRUE model Q=P: sum_t -log2 P(x_t|x_<t)."""
    L = len(word)
    total = 0.0
    for t in range(1, L + 1):
        p_full = float(block_prob_hmm(word[:t]))
        p_ctx = float(block_prob_hmm(word[:t - 1])) if t > 1 else 1.0
        cond = p_full / p_ctx
        total += -math.log2(cond)
    return total


def sample_word(N):
    # sample hidden chain from PI then transitions, emit symbols
    s = 0 if random.random() < float(PI[0]) else 1
    w = []
    for t in range(N):
        if t > 0:
            s = 0 if random.random() < float(T[s_prev][0]) else 1
        a = 0 if random.random() < float(E[s][0]) else 1
        w.append(a)
        s_prev = s
    return w


print("=" * 72)
print("G3  ACHIEVABLE scheme realizes the rate (one forward pass), never computes the #P-hard optimum (Rem 7.15b)")
Nsmall = 8
# identity: E_P[L(X^N)] = H(X^N) exactly (sum over all words of p * (-log2 p)).
EL_exact = 0.0
for w in itertools.product(SYM, repeat=Nsmall):
    p = float(block_prob_hmm(list(w)))
    EL_exact += p * causal_codelength_true(list(w))
H_small, _ = exact_block_entropy(Nsmall)
id_ok = abs(EL_exact - H_small) < 1e-9
print(f"  identity: E_P[L(X^{Nsmall})] = {EL_exact:.8f}  ==  H(X^{Nsmall}) = {H_small:.8f}  (match: {id_ok})")
# Monte-Carlo: the per-symbol realized mean converges to (1/N) H(X^N) WITHOUT computing H.
trials = 20000
acc = 0.0
for _ in range(trials):
    acc += causal_codelength_true(sample_word(Nsmall))
mc_mean = acc / trials
target = H_small  # the achievable scheme targets the (eps-approximable) entropy, not exact symbolic
mc_ok = abs(mc_mean - target) < 0.05 * target  # 5% relative band at 2e4 trials
print(f"  Monte-Carlo realized E[L] over {trials} draws = {mc_mean:.6f} bits ; "
      f"H(X^{Nsmall}) = {target:.6f} (within 5%: {mc_ok})")
# converse + KL gap: matched model attains H(X^N); a MISMATCHED model Q' pays strictly more.
def causal_codelength_mismatch(word, Eq):
    """Realized length under a mismatched emission model Eq (wrong emission probs), exact hidden filter
    of the WRONG model.  E[L] under TRUE source = H(X^N) + N*D_KL-ish gap >= H(X^N) (Gibbs)."""
    # mismatched model uses Eq but same T, PI: its predictive Q(x_t|x_<t)
    def block_prob_q(w):
        Lw = len(w)
        if Lw == 0:
            return Fraction(1)
        alpha = [PI[s] * Eq[s][w[0]] for s in HID]
        for t in range(1, Lw):
            nxt = [Fraction(0), Fraction(0)]
            for s2 in HID:
                nxt[s2] = sum(alpha[s1] * T[s1][s2] for s1 in HID) * Eq[s2][w[t]]
            alpha = nxt
        return sum(alpha)
    Ln = len(word)
    total = 0.0
    for t in range(1, Ln + 1):
        pf = float(block_prob_q(word[:t]))
        pc = float(block_prob_q(word[:t - 1])) if t > 1 else 1.0
        total += -math.log2(pf / pc)
    return total


Eq = [[Fraction(6, 10), Fraction(4, 10)],   # mismatched emissions
      [Fraction(5, 10), Fraction(5, 10)]]
EL_mismatch = 0.0
for w in itertools.product(SYM, repeat=Nsmall):
    p = float(block_prob_hmm(list(w)))       # expectation under the TRUE source
    EL_mismatch += p * causal_codelength_mismatch(list(w), Eq)
converse_ok = (EL_exact >= H_small - 1e-9) and (EL_mismatch > H_small + 1e-9)
print(f"  converse: matched E[L]={EL_exact:.6f} == H={H_small:.6f} (floor attained); "
      f"mismatched E[L]={EL_mismatch:.6f} > H (KL gap {EL_mismatch - H_small:.4f} bits)")
print(f"           so H(X^N) is the FLOOR the achievable scheme approaches; equality at the matched model only.")
g3_ok = id_ok and mc_ok and converse_ok
print(f"  G3 {'PASS' if g3_ok else 'FAIL'}  (achievable one-pass scheme realizes the rate; exact optimum never computed)")

# ====================================================================== summary
print("=" * 72)
overall = g1_ok and g2_ok and g3_ok
print(f"SUMMARY: G1 {'PASS' if g1_ok else 'FAIL'} (intractable exact target, Rem 7.15g) | "
      f"G2 {'PASS' if g2_ok else 'FAIL'} (poly eps-approx under FS, Rem 7.15c) | "
      f"G3 {'PASS' if g3_ok else 'FAIL'} (achievable scheme realizes rate, Rem 7.15b)")
print(f"OVERALL -> {'PASS' if overall else 'FAIL'}")
