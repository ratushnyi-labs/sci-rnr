# Part III dispersion statistics — simulation studies

Empirical program for the operational rate–distortion dispersion results of
Part III (`tex/papers/dispersion/rnr_dispersion.tex`, §7.34). Four
simulation studies, each with executable PASS gates, saved plots (PNG) and
tables (CSV) in `out/`.

## Studies

| study | law (paper source) | script |
|---|---|---|
| (a) | finite-`n` operational dispersion on the BSMS inside Gray: `Var(j_n)/n → V_lossless = p(1-p) log₂²((1-p)/p)`, deficit exactly `V/n` (Remark 7.34b) | `study_a_finite_n_dispersion.py` |
| (b) | third-order recentring: `j_n = i_n − nc` (`c = h₂(D)`), empirical `P(j_n ≤ nR)`, and the `+½ log₂ n` drift of the ball log-mass `G_n = −log₂ P_{Y*}(B_{nD})` over `j_n` (Remark 7.34n′) | `study_b_third_order_recentring.py` |
| (c) | beyond-Gray finite-`n` threshold law `D_c^{(n)} − D_c = K(p)²/n²`, `K(p) = πp(1−2p)^{1/4}/(2(1−p)^{3/2})` (Remarks 7.34o/7.34o′) | `study_c_beyond_gray_law.py` |
| (d) | replica margins `R_A(s;D) − 3/2 > 0` on dense Gray grids, `A = 2..5`, with the binding-corner law `margin ≈ (A−2)/(A−1)·p` (Theorem 7.34m / Remark 7.34m′) | `study_d_replica_margins.py` |

## Running

```sh
# reduced-scale gate (all four studies, < 1 min wall):
/Users/para/.venvs/rnr/bin/python -u run_checks.py

# full-scale configs (longer ladders, more samples, larger blocks):
/Users/para/.venvs/rnr/bin/python -u run_checks.py --full

# individual studies:
/Users/para/.venvs/rnr/bin/python -u study_c_beyond_gray_law.py [--full]
```

Every script prints PASS/FAIL lines per check and an OVERALL verdict
(style of `scripts/verify/`).

## Machinery (`common.py`)

Exact finite-blocklength machinery, no fitting shortcuts:

* Walsh–Hadamard XOR-convolution (`O(N log N)`) for the block
  deconvolution `P_{Y*} = (K^{-1})^{⊗n} P_X`, the tilted kernel
  `M(x) = E_{Y*}[e^{−λ* d_H}]` (hence exact `j_n`), and Hamming-ball masses
  (hence exact `G_n`) — the solver pattern validated against dense
  Blahut–Arimoto in the committed
  `scripts/verify/probe_7_34_beyond_gray_structure.py` (B3) and re-validated
  here against a dense reference (check A0);
* 40-digit transfer-matrix bisection for the alternating binding-word mass
  (study (c); avoids the float64 roundoff bias of the full-WHT bisection at
  `n ≥ 12`, cf. the committed probe's 80-bit control note);
* the pattern-quotient replica builder (4×4 at `A = 2`, 5×5 at `A ≥ 3`) and
  the A-ary Gray threshold as the smallest positive root of the
  relevant-cubic discriminant of the *multiplicity-weighted* 3-class
  alternating transfer — both cross-validated to ≤ 1e-9 against the
  committed `scripts/verify/thm_7_34_route_b_a{3,4,5}_assembly.py`
  (check D0). The collision flag is non-monotone in `D` for `A ≥ 4`, so the
  threshold is located by first-sign-change scan, not bare bisection.

## Headline reduced-scale results

* (a) exact operational `Var(j_n)/n = V_lossless (n−1)/n` to `< 1e-10`
  (n = 6..14, three (p,D) Gray points incl. the endpoint `D = D_c`); Monte
  Carlo ladder to `n = 2048` covers `V_lossless` within CI; deficit fit
  slope `−1.0000`, constant `= V` to 4 digits.
* (b) mean and quantile drift of `G_n − j_n` fit `0.41–0.43 · log₂ n`
  (rising toward the asymptotic `½` with `n`; `O(1)` remainder bounded).
* (c) fitted/predicted `K(p)²` ratios: `1.002 / 0.990 / 0.974` at
  `p = 0.1 / 0.2 / 0.3` (ladder to `n = 24`); log-log slopes `−2.06 / −2.03
  / −1.99`.
* (d) min margin over 960 Gray grid points per alphabet: positive at every
  point (worst: `1.3e-5` at `A = 2`, small-`p` corner, where sharpness is
  expected); corner-law extrapolations `−0.0001 / 0.4998 / 0.6666 / 0.7500`
  vs predicted `0 / 1/2 / 2/3 / 3/4`.

## Honest caveats

* Study (a): at Monte Carlo blocklengths (`n > 14`) `j_n` is evaluated via
  the all-`n` SLB identity `j_n = i_n − n h₂(D)` (proven on Gray in Remark
  7.34b and machine-verified here at the exact-grade rungs); the exact-grade
  rungs are identity-free (full `2^n` solver).
* Study (b): the ladder slopes sit slightly below `½` at these small `n`
  (0.41–0.43, increasing) — the finite-`n` approach to the asymptote; the
  BSMS surprisal is lattice, so the CDF checks carry a lattice band (the
  lattice caveat is stated in Remark 7.34n′ itself). Gates use the band
  `[0.30, 0.70]` (sign + magnitude, separating `½` from `0` and `1`).
* Study (c): the full-WHT float64 threshold is used only as an `n = 10`
  cross-check (binding word = alternating word); at `n = 12`, small `p`, its
  known float64 bias is visible and reported descriptively.
* Study (d) is a *numeric-grade* margin scan (dense grids + corner
  extrapolation), complementing — not replacing — the exact certificates of
  the committed Route B lemmas; grid density is bounded by run time.
* Reduced-scale gate ≈ 10 s wall; `--full` extends ladders/samples
  (largest block `2^24`–`2^25` in study (b) needs ~1–2 GB RAM transiently).
