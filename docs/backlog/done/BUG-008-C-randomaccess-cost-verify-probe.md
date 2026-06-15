# BUG-008-C - randomaccess-cost-verify-probe

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** scripts/verify/ (new probe) + CI build job; asserts invariants of tex/rnr_coding.tex Theorem 10.3 cost formula and Remark 10.3a sub-block cost
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
Source: Thm 10.3 cost formula and its polylog-regime conditional clause; Rem 10.3a
sub-block-granular cost; the project rule "every theorem/lemma/proposition has an
executable PASS/FAIL check in scripts/verify/."

BUG-008 scopes the random-access claim (A) and tabulates per-mode granularity (B). Per
project convention each such claim needs an executable probe. This leaf adds the missing
verification probe for the Theorem 10.3 / Remark 10.3a cost accounting. The probe is a
counting / symbolic-numeric model (NOT a neural implementation): it instantiates the
cost-accounting of Theorem 10.3's proof as a deterministic function
$\mathrm{work}(p,k,K,|X|,\mathrm{cost}(M),L(R))$ that counts (i) one seek-index probe of
$O(\log L(R))$ bits, (ii) $\le K+k$ decoded positions of which $\le K$ are warmup, (iii)
$\le (K+k)\cdot\mathrm{cost}(M)$ predictor-evaluation units, and verifies the claimed
relations exactly on small instances.

FIX TARGET (this leaf): write `scripts/verify/thm_10_3_random_access_cost.py` (PASS/FAIL,
included as a CI job) asserting, over a sweep of $(|X|,K,k,p)$ small instances:
1. decoded-position count $= p+k'-P \le K+k$ with $P=K\lfloor p/K\rfloor$ and
   $k'=\min(k,|X|-p)$ (the bound the proof states), including the short-read edge case
   $p+k>|X|$ and the $|X|<K$ single-sync-point case;
2. the per-query cost matches $O(\log L(R)+(K+k)\cdot\mathrm{cost}(M))$ and reduces to
   $\Theta(|X|\cdot\mathrm{cost}(M))$ at $K=\Theta(|X|)$ or $k=\Theta(|X|)$ (the stated
   degeneracy — no better than full sequential decode);
3. the polylog-in-$|X|$ property holds iff $K,k,\mathrm{cost}(M),\log L(R)$ are
   polylog and $L(R)=\mathrm{poly}(|X|)$, and the simultaneous-polylog impossibility
   (polylog per-query AND polylog total sync overhead cannot co-exist) holds against the
   $\Delta_{\mathrm{sync,total}}=\Theta(|X|/K)$ accounting of Lemma 5.1b / §10.7;
4. Remark 10.3a: a read $[p,p+k)$ spans $\lceil k/K\rceil+O(1)$ sub-blocks and costs
   $O(\log L(R)+K\cdot\mathrm{unrank}(K))$ for $k\ll K$, decoding/discarding up to one full
   sub-block per request (sub-block granularity strictly weaker than byte granularity);
5. (if BUG-008-B lands first) a lightweight document-invariant check that every mode row in
   the granularity table is internally consistent (granularity matches the cost formula it
   quotes). This sub-check is optional/soft if B is not yet merged.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class repo, but this leaf produces a runnable artifact in
  `scripts/verify/` — fully closeable in-repo with the project's standard small Python venv
  (`/Users/para/.venvs/rnr/bin/python`); no neural model, no corpus, no external compute.
- The probe is a cost-accounting / combinatorial check (it models the seek + decode-window
  arithmetic), NOT a benchmark; wall-clock seek-latency is explicitly the separate external
  leaf BUG-008-D.
- ABOVE_EASY (not EASY): the probe must encode several distinct regime relations and edge
  cases (short read, $|X|<K$, $K=\Theta(|X|)$ degeneracy, the simultaneous-polylog
  impossibility) and be wired as a CI job following the existing
  `scripts/verify/lemma_5_1*` precedents.
- Non-goal: paper wording (BUG-008-A), the granularity table itself (BUG-008-B), and the
  measured seek-latency (BUG-008-D).

## Refined Description
**Scope:** one `scripts/verify/` PASS/FAIL probe + CI job asserting the Thm 10.3 / Rem 10.3a
cost-accounting invariants and regime boundaries.  **Non-goals:** wording (A); table (B);
wall-clock measurement (D).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-008; soft-depends on BUG-008-A/B for the exact notation and the
optional table-consistency sub-check (probe core does not block on them).
**Risks / Open Questions:** must verify the cost ACCOUNTING (deterministic relations), not a
neural runtime; keep it a counting model so the green/red signal is about the formula, not
about an implementation. Mirror the CI-job wiring of the existing §5.1 verify scripts.
**Expected Code Change:** YES (new scripts/verify file + CI job entry)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | finalize the exact invariant list (1--5) and edge cases |
| ESTIMATION | DONE | ABOVE_EASY (multi-regime probe + CI wiring, single file) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | counting model of seek + (K+k) decode window; regime sweep design |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | author scripts/verify/thm_10_3_random_access_cost.py |
| TESTING | PENDING | the probe itself is the test; PASS on full sweep; add CI job |
| POST_AUDIT | PENDING | adversarial re-review (probe asserts the real bound, not a tautology) |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at §12.7 decomposition of BUG-008;
  no code/proof change made before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
  The scripts/verify probe is the project-standard executable check, not a product runtime.
- §12.7: BLOCKING child leaf of BUG-008 (the closeable-now executable-verification layer).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-008 A/B/C (commit pending). A: cost(M) surfaced as a non-universal, dominant per-position neural-eval term at Abstract/§1.1/§10.7/Thm 10.3 (polylog = in archive params with cost(M) fixed, not wall-clock; aligned w/ Lemma 5.1b). B: 6-row per-mode access-granularity table near Rem 10.3a (byte vs sub-block; per-query cost; wasted decode; sources Thm 10.3/Rem 10.3a/§6.3/§10.7; every cell cited, no new claim). C: scripts/verify/bug_008_random_access_cost.py (decode-window count p+k'-P<=K+k; K=Theta(|X|) degeneracy; simultaneous-polylog impossibility via interior K*; sub-block-granular cost; OVERALL -> PASS) + CI job. Hostile-referee POST_AUDIT clean (table byte-vs-sub-block split verified both directions). Supervisor: shortened the Type-III-C cost cell to a §10.7 pointer; table's 39.81pt overfull is house-style (188 such in the paper, many larger).
