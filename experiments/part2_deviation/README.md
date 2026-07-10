# Part II deviation-hierarchy statistics — simulation studies

Empirical program for the deviation-hierarchy theorems of Part II
(`tex/papers/random_access/rnr_random_access.tex`): Theorems 7.31–7.33, 7.35,
7.36.  All studies run on **simulated sources** (i.i.d. and 2–3-state Markov
chains with closed-form entropy rate `h` and varentropy rate `V`), fully
seeded.  The exact laws are those already unit-verified in
`scripts/verify/thm_7_3{1,2,3,5,6}_*.py`; this directory turns each into a
statistical study with figures and data files.

## Studies

| Study | Theorem | Question | Script |
|---|---|---|---|
| (a) | 7.31 | Does the empirical overflow probability `P(L_N ≥ N(h+δ))` follow the Legendre rate function `E(R)` (tilted-transfer-matrix CGF for Markov), and is the Markov LDP endpoint the **max cycle mean** `r₊` (not the single-step max)? | `study_a_overflow_ldp.py` |
| (b) | 7.32 | Is a fixed model (wrong **or true**) overflow-exponent-deficient vs `E* = inf{D(Q‖P₀) : H(Q) ≥ R}`, and does the universal mixture cost fit `(d/2)·log₂N`? | `study_b_universality.py` |
| (c) | 7.33 | Does the RA penalty — a constant standardized shift `β = δ/(c√V)` — fade from the exponent (`ln P_RA / ln P → 1`) while staying in the probability (`(1/a_N)·ln(P_RA/P) → β`) along `a_N ~ N^{1/4}`? | `study_c_moderate_fade.py` |
| (d) | 7.35 | Does `max_i L_i` over `m` blocks track `Kh + √(2KV·ln m)` (refined Gumbel centering), and is the validity boundary `ln m = o(K^{1/3})` simulation-visible? | `study_d_query_evt.py` |
| (e) | 7.36 | Does the repair process scale to Brownian motion (KS on increments), with a half-normal net-surplus write buffer, a `sup|W|`-law Lindley high-water-mark, and **one** `V` governing all three buffers (7.29 / 7.35 / 7.36)? | `study_e_buffer_fclt.py` |

## Methods

* **Tilted importance sampling** for deep tails (unbiased): exponential tilt of
  the multinomial law (i.i.d.) and the Perron right-eigenvector tilted kernel
  of the transfer matrix (Markov).  Cross-checked against direct Monte Carlo
  wherever direct MC is feasible, and against exact multinomial type sums in
  study (b).
* **Exact machinery**: Legendre rate functions on a `t`-grid, fundamental-matrix
  Markov CLT variance for `V`, exhaustive simple-cycle enumeration for the max
  cycle mean, DP over paths for the finite-`N` endpoint certificate.
* **Slope fits** of `−log₂P` vs `N` are Bahadur–Rao-corrected (the `½log₂N`
  prefactor is imposed, the exponent is the fitted slope).

## Running

```sh
# reduced scale (< 10 min, seeded, all checks):
/Users/para/.venvs/rnr/bin/python run_checks.py

# full scale (larger ladders / sample sizes):
/Users/para/.venvs/rnr/bin/python run_checks.py --full

# single study:
/Users/para/.venvs/rnr/bin/python study_c_moderate_fade.py [--full]
```

Every run prints `[PASS]/[FAIL]` per check and `OVERALL: PASS/FAIL`
(exit code 0/1), and rewrites `out/*.png` + `out/*.csv`.

## Outputs (`out/`)

* `study_a.png`, `study_a_rates.csv`, `study_a_endpoint.csv`
* `study_b.png`, `study_b_exponents.csv`, `study_b_redundancy.csv`
* `study_c.png`, `study_c_fade.csv`
* `study_d.png`, `study_d_evt.csv`, `study_d_boundary.csv`
* `study_e.png`, `study_e_buffers.csv`

## Dependencies

`numpy`, `scipy`, `matplotlib` (all in `/Users/para/.venvs/rnr`).
