# BUG-009-D - residual-ak-bound-beyond-gray-research

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** BLOCKED
**Complexity:** HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§7.34: Remark 7.34m$''$ rigorous $a_k$ bound / convexity capstone; Remark 7.34o beyond-Gray exact $V_{\mathrm{op}}$)
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** MEDIUM
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-15  **Updated:** 2026-06-15
**Recovered:** N/A

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 009 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-009 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: §7.34 (Remark 7.34m$''$ / Remark 7.34o).

The irreducible RESEARCH residual of BUG-009, of which exactly ONE direction is to
be attacked at a time (waterfall). Two genuinely-open, distinct sub-problems share
the same root obstruction (the optimal output / replica Perron eigenvalue needs a
growing Markov order) and are kept in a SINGLE leaf because each is a deep open
theorem that cannot be honestly split into citable pieces; the in-repo-closeable
parts are already carved into siblings BUG-009-A/B/C:

  (D1) **Rigorous cancellation-aware $a_k$ bound (the convexity capstone).** Prove
  $|a_k(D)|=O(D^{k+1})$ for the Chebyshev coefficients of $g_A(s)=\sum_k a_k(D)
  T_k(\cos s)$ (equivalently the per-harmonic decay of the replica Perron
  eigenvalue $\rho(W_A)$, since the prefactor $|C_s/C_0|^2$ is affine and carries no
  curvature), rigorously enough to make $\mathrm{LB}(D)=4a_2-\sum_{k\ge3}|a_k|
  \tfrac13 k^2(k^2-1)>0$ a THEOREM on the memory bulk — and, ideally, to extend it
  to the near-uniform / memoryless edge $p\to(A-1)/A$ where the crude $T_k''(1)$
  endpoint bound over-counts. This upgrades the all-$D\le D_c$ in-Gray achievability
  from "closed modulo $a_k$" to a full theorem, hence $s=\pi$ binding, hence
  unconditional achievability on the memory regime.

  (D2) **Beyond-Gray exact $V_{\mathrm{op}}(D)$ for $D>D_c$.** Characterize the
  operational dispersion with the optimal non-SLB-tight, memory-aware test channel
  (the deterministic-shift converse and the SLB dual identity both die for $D>D_c$).

FIX TARGET: close ONE of (D1)/(D2), OR pin a partial/weakened-but-honest result on
one of them.

## Clarification / Assumptions
- This leaf is closeable_in_repo = NO in the sense that genuine closure needs new
  mathematics (the memory journals document every standard route as pinned to a
  named blocker). An honest NON-closure ("residual remains, more sharply pinned") is
  an acceptable terminal state and does NOT block shipping siblings A/B/C.
- GOVERNANCE/DOCUMENTATION-class defect in `tex/rnr_coding.tex`. A genuine closure
  would add a new lemma + verify script; otherwise no paper change beyond the
  A/B/C scoping.
- Kept SINGLE (not split) per the project waterfall rule: one open theorem at a
  time; (D1) and (D2) share the growing-Markov-order obstruction but are different
  proofs — whichever is attacked, only ONE is taken at a time.

## Refined Description
**Scope:** Attempt to close (or partially/weakly advance) exactly ONE of:
  - (D1) the rigorous $a_k=O(D^{k+1})$ cancellation-aware bound. The route the
    journal identifies as the only live one: a cancellation-aware coefficient
    recursion from the degree-2 trigonometric-polynomial characteristic equation of
    $W_A$ (the harmonics must cancel the all-forward $O(D^k)$ term of $\rho$), since
    the crude analyticity-strip (Bauer–Fike, $\tau^*\approx\ln(1/D)$) Bernstein bound
    is provably one $D$-power too loose (it gives $a_2\sim D^2$, not the true
    $D^3$) — DEMONSTRATED to FAIL by ~10× at $p=0.2$, do NOT re-attempt the crude
    split.
  - (D2) the beyond-Gray exact $V_{\mathrm{op}}(D)$ via the optimal memory-aware
    test channel. NOTE the established negative (journal): the off-Gray
    converse-ratio $\rho(D)$ has NO elementary closed form (it is the growing-order
    Green–Kubo $\chi''(0)$ rate); do NOT re-hunt one. A genuine advance is a
    characterization of the optimal test channel's dispersion, not a closed form for
    $\rho$.
Any attempt MUST follow the project 3-step post-fix protocol (cold self re-review
+ independent adversarial re-review with status RESOLVED/PARTIAL/UNRESOLVED/NEW +
adversarial disprove "PROVE wrong / Attack failed; result holds") before any commit,
and the "small honest steps" rule (a weakened but proven sub-claim beats a strong
hand-wave). STEALTH: cold-review + amend any subagent auto-commit (strip any
authorship trailers; author Pavlo Ratushnyi).
**Non-goals:** the in-repo scoping ledger (sibling BUG-009-A), the $a_k$-decay /
verification probe (sibling BUG-009-B), and the beyond-Gray map sharpening (sibling
BUG-009-C) — those carry no open math and must be completed first/independently.
Do NOT re-attempt any journalled dead end:
  - the a.c.-Blackwell lead (filter measure singular on all of Gray; category error);
  - the inhomogeneous lattice LLT via DH-2020 (blocked at ellipticity/density-floor),
    Hafouta-2025 re-encoding (LLT is a consumer not a producer; growing-dim cocycle),
    Kloeckner Thm C (homogeneous single test fn only), chaining (provably
    $O(\sqrt n)$), filter-stability-alone (cancellation not forgetting), contour-wide
    KP cluster expansion (Lee–Yang/Fisher zeros on the contour);
  - the crude Bauer–Fike / Bernstein $a_k$ bound (one $D$-power loose, fails by 10×);
  - an elementary closed form for the off-Gray $\rho(D)$ (provably none).
Each is pinned to a named blocker in the BUG-009 memory journals; defeat the blocker
before reusing the route.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-009; sequentially AFTER siblings BUG-009-A and
BUG-009-B (the precise status ledger + the reproducible $a_k$-decay target/falsifier
must exist before attacking (D1); BUG-009-C's map before (D2)) so the attempt starts
from a clean obstacle map.
**Risks / Open Questions:** HIGH research risk — the journals document both
directions as genuinely research-level open (the $a_k$ capstone is a "several-step
analyticity task, not a one-shot"; beyond-Gray needs the growing-order optimal
channel). The most likely terminal state is an honest non-closure (residual remains,
more sharply pinned), acceptable under the project rules and NOT blocking siblings.
**Expected Code Change:** UNCERTAIN (NO if non-closure; a new lemma + verify script
if a sub-case genuinely closes)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | pick ONE direction (D1 or D2) — waterfall, one open theorem at a time |
| ESTIMATION | PENDING | HARD/research; outcome may be honest non-closure |
| DECOMPOSITION | N/A | atomic residual; cannot be honestly split further (shared growing-order obstacle, distinct proofs) |
| DESIGN | PENDING | route must defeat the named blocker (D1: cancellation-aware recursion, not crude Bauer–Fike; D2: optimal memory-aware channel, not a closed-form ρ) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | new lemma + proof IF closed; else no paper change beyond A/B/C scoping |
| TESTING | PENDING | constructive / coefficient-bound or test-channel verify script IF a sub-case closes |
| POST_AUDIT | PENDING | mandatory 3-step protocol (cold + independent adversarial re-review + adversarial disprove) before any commit; amend any auto-commit trailers (stealth) |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 DECOMPOSITION child of BUG-009. This is the irreducible genuinely-open
  research residual; the in-repo-closeable work (scoping, verification, map) is
  carved into siblings A/B/C so the parent's complexity drops to a shippable core
  plus this clearly-flagged residual.
- §12.1.1 (Backlog-File-First): created fresh this turn; no paper change precedes
  this file. `Recovered: N/A`.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- closeable_in_repo = NO: per the BUG-009 memory journals every standard route is
  pinned to a named blocker; genuine closure needs new machinery. Honest
  non-closure is an acceptable terminal state and does not block siblings A/B/C.
- Per project rules (Truth > peremoga, waterfall): one direction at a time; a
  weakened-but-proven result beats a strong hand-wave; "Attack failed; holds"
  required from the adversarial step before any closure commit.
