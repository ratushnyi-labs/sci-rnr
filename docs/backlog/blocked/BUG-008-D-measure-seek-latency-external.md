# BUG-008-D - measure-seek-latency-external

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** BLOCKED (external)
**Complexity:** HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** EXTERNAL (reference implementation + benchmark runs); protocol lives in tex/rnr_experimental_design.tex §2.8 (H10, H11) and §2.9 (H12, H13); claims to confirm/condition: Rem 10.3a 650 ms estimate (tex/rnr_coding.tex ~ll.16112--16121) and the 20--200 ms figures
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
> (Filed gap 008 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-008 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: BUG-008's "ideally a measured seek-latency" branch; the Remark 10.3a wall-clock
estimate ("$K+k=65$ KiB worth of predictor evaluations --- typically on the order of
$650\,\text{ms}$ on a single CPU core ... one forward pass at $\sim 10\,\mu\text{s}$ per
byte"); the abstract/§10.7 "approximately one percent" and the §11/companion 20--200 ms
figures, all flagged as heuristic conjectures, not theorems or measurements.

This is the IRREDUCIBLE external half of BUG-008: actually MEASURING per-query seek latency
(and its independence of archive size $|X|$). The paper's latency numbers — the 650 ms
single-core estimate, the 10 $\mu$s/byte forward-pass rate, the 20--200 ms commodity-hardware
range, the GPU/NPU sub-millisecond claims — are derived from assumed predictor throughputs,
not from a running coder. Confirming them (or refuting/conditioning them) requires building or
obtaining a reference RNR coder with a $W$-bounded predictor, building archives across sizes
(100 MiB to 100 GiB), and benchmarking $T(N,k)$ at uniformly random offsets. The companion
already pre-registers this exact campaign: H10 ($T(N,k)$ independent of $N$; 20--200 ms;
sync-overhead < 2--5%), H11 (LSTM state-snapshot variant 3--10x slower), and the mounted-FS
percentile-latency hypotheses H12/H13. This CANNOT be closed inside this docs/theory repo: it
needs a working implementation, archives, and CPU/GPU/NPU hardware.

FIX TARGET (this leaf): execute the pre-registered seek-latency protocol (H10, H11; and the
mounted-FS H12/H13 if scope extends), produce measured per-query latency with the
size-independence check, then feed the outcome back to the paper — confirm the Rem 10.3a /
§11 latency figures, or refute and downgrade them to measured ranges. Until then the in-repo
posture is: the surfaced `cost(M)` dependence (BUG-008-A), the per-mode granularity table
(BUG-008-B), and the cost-accounting probe (BUG-008-C); the wall-clock numbers stay flagged
as heuristic estimates with a forward pointer to companion §2.8/§2.9.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class file, but the WORK is external: flagged
  closeable_in_repo = NO because a measured seek-latency needs a reference implementation,
  archives, and hardware that do not exist in this repo.
- This is the genuinely-hard residual the §12.7 split isolates so the rest of BUG-008
  (A/B/C) can be honestly closed now. HARD rather than VERY_HARD because the protocol is a
  bounded, already pre-registered latency benchmark (H10/H11, optionally H12/H13) against an
  existing engineering spec — a narrower program than a full multi-corpus ratio campaign —
  yet still irreducibly external. It is intentionally NOT split further: a reference coder
  plus a cross-size latency benchmark is one external program, not a stack of EASY in-repo
  steps.
- Stays open (or moves to a tracking/external-work state) after A/B/C land; it is the
  documented residual, clearly flagged as needing out-of-repo work.
- Non-goal: any in-repo wording / table / accounting probe (those are A/B/C).

## Refined Description
**Scope:** run the pre-registered seek-latency protocol (companion §2.8 H10/H11; optionally
§2.9 H12/H13) and report measured $T(N,k)$ with size-independence; reconcile the Rem 10.3a /
§11 latency estimates with measurement.  **Non-goals:** in-repo scoping / table / accounting
probe (BUG-008-A/B/C); the cost FORMULA verification (already BUG-008-C, which is symbolic,
not wall-clock).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-008; uses the companion's pre-registered H10/H11 (so the run
targets the exact pre-registered hypotheses) and the engineering-spec reference predictor;
does not block A/B/C.
**Risks / Open Questions:** requires a reference implementation (engineering-spec
realization), archive generation at 100 MiB--100 GiB, and CPU/GPU/NPU hardware; results may
refute the 650 ms / 20--200 ms / sub-ms figures and force a downgrade to measured ranges —
that is an acceptable, honest outcome, not a failure of this leaf. The size-independence
claim (the kill-feature) stands or falls on H10's $N$-sweep.
**Expected Code Change:** NO (in THIS repo); external implementation + benchmark are required
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | freeze the latency-protocol scope (H10/H11, optional H12/H13) and hardware set |
| ESTIMATION | DONE | HARD; closeable_in_repo = NO (external work) |
| DECOMPOSITION | N/A | atomic residual leaf (one external benchmark program; not further splittable into EASY in-repo steps) |
| DESIGN | PENDING | reference-coder + latency-benchmark plan (engineering spec + companion §2.8/§2.9) |
| FRONTEND | N/A | no UI |
| BACKEND | N/A (this repo) | measurement is external; no tex/proof change here |
| TESTING | PENDING (external) | measured T(N,k) across sizes + size-independence; percentile latency per hardware path |
| POST_AUDIT | PENDING | reconcile measured latency with Rem 10.3a / §11 estimates; adversarial re-review |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD (external; large)
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at §12.7 decomposition of BUG-008;
  no code/proof change made before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: BLOCKING child leaf of BUG-008 -- the isolated external/measurement residual.
  Flagged closeable_in_repo = NO so the parent's complexity collapses to the shippable
  A/B/C core plus this clearly-marked deferred seek-latency measurement.


## Status (2026-06-15)
IRREDUCIBLE EXTERNAL residual: measured wall-clock seek latency (companion §2.8 H10/H11) with a reference coder across archive sizes/hardware, to confirm/condition the Rem 10.3a 650 ms and §11 20-200 ms / sub-ms predictions. Not closeable in this theory repo; the in-repo cost accounting + per-mode granularity (BUG-008-A/B/C) is done. Parked in blocked/.
