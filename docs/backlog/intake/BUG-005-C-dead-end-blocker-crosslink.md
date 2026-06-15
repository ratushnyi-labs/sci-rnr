# BUG-005-C - dead-end-blocker-crosslink

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§13.4.1/13.4.2/13.4.4 failed-path lists; Remark 4.3d, Lemma 13.4.2a)
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
> (Filed gap 005 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-005 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: §13.4 "What has been tried (failed paths)" lists; OP2/OP4 attempt journals.

The parent gap text carries an explicit warning: "op2/op4 attempt-journals
document dead-end routes -- do not re-attempt without addressing the recorded
blockers." Those journals (op2_attempts, op4_attempts, op4_bitcut_spectral,
op_progress_2026_05/06) record specific blockers that are NOT all reflected in
the paper's §13.4 failed-path prose at the granularity needed to prevent
re-attempts — e.g. (i) the ceiling-discontinuity concrete witness (V=615,
|Σ|+K crossing 2^11) for OP2(a); (ii) the OP2(b) "both-have-poly-witness ⇏
ratio" correction (Lemma 13.4.2a is an UPPER bound, not a closure) and the
PRG-route giving DECISION-only not APX (Theorem 13.4.2b); (iii) the OP4
re-localisation from "row-stochastic F is the obstacle" (now CORRECTED) to the
partial-cube no-robust-NO-promise barrier (Remark 4.3d), plus the three
exhausted reduction families (bit-cut/spectral, cube-host-QAP layout transfer,
max-dilation-1 gap theorem absence). FIX TARGET (this leaf): ensure each
recorded dead-end appears in §13.4 with its precise blocker so the roadmap is a
faithful re-attempt-prevention ledger; correct any §13.4 prose that still states
a superseded obstacle framing.

## Clarification / Assumptions
- This is the IN-REPO DOC-COHERENCE leaf of BUG-005: it does NOT attempt any
  open reduction. It transcribes the memory-journal blockers into the paper's
  attack-vector roadmap and fixes any stale obstacle framing. No new mathematics.
- The corrective is purely prose: align §13.4 / Remark 4.3d / Lemma 13.4.2a
  wording with the validated current understanding so future work does not redo
  the documented dead ends.
- closeable_in_repo = YES: this is editing prose to match already-established
  (and journalled) negative results.
- STEALTH: journal content describes the math/blockers only; no authorship traces
  are introduced into the paper.

## Refined Description
**Scope:** Walk each §13.4 failed-path bullet and confirm it names: (a) the
specific reduction attempted, (b) the precise blocker, (c) (where one exists) the
witnessing negative-control verify script. Add any missing blocker that the
journals record but the paper omits at re-attempt-prevention granularity, in
particular: the OP2(a) ceiling-discontinuity concrete witness; the OP2(b)
"upper-bound ≠ closure" correction and PRG decision-only caveat; the OP4
obstacle re-localisation (correct any residual "row-stochastic F is THE obstacle"
phrasing to the partial-cube no-robust-NO-promise barrier per Remark 4.3d), and
the note that the bit-cut/spectral and layout-transfer routes were tried and
yield no advance. Verify internal consistency (e.g. §6/§4 cross-references to
§13.4 resolve to the corrected framing).
**Non-goals:** the per-sub-question status table (sibling BUG-005-A); verify
script consolidation (sibling BUG-005-B); any new reduction (sibling BUG-005-D).
No change to scripts/verify/.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-005; complements sibling BUG-005-A (A states the
current status; C ensures the FAILED-PATH side of the roadmap is faithful). Run
A and C against a single cold read of §13 to keep them consistent.
**Risks / Open Questions:** must not introduce a NEW overclaim while fixing an
old one — every transcribed blocker must be traceable to a journal entry that is
itself validated (some journal lines are point-in-time; verify against live HEAD).
The op4_bitcut_spectral journal explicitly says its reframing must NOT be added
to the paper as content — only its NEGATIVE conclusion (route exhausted) belongs
in §13.4.
**Expected Code Change:** YES (tex prose edits only; no math, no table)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm which journal blockers belong in §13.4 vs memory-only |
| ESTIMATION | PENDING | ABOVE_EASY: prose alignment of failed-path lists with journals |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | edit §13.4 failed-path bullets + Remark 4.3d / Lemma 13.4.2a framing |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | §13.4 / Remark 4.3d prose edits |
| TESTING | PENDING | none (doc-only); negative-control scripts owned by sibling B |
| POST_AUDIT | PENDING | cold re-read + independent adversarial re-review: each transcribed blocker must survive "Attack failed; holds" framing |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 DECOMPOSITION child of BUG-005. This is the re-attempt-prevention doc leaf:
  it lowers parent complexity by making the paper's failed-path roadmap a faithful
  ledger, so the residual research (sibling BUG-005-D) starts from a clean map of
  dead ends rather than rediscovering them.
- §12.1.1 (Backlog-File-First): created fresh this turn; no paper change precedes
  this file. `Recovered: N/A`.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Per project rules: do not re-attempt the journalled dead ends; transcribe only
  their validated blocker conclusions.
