# SPLIT-PLAN: `tex/rnr_coding.tex` -> three papers (P1 / P2 / P3)

Analysis snapshot: repo `/Users/para/work/rnr`, commit `66c7ebf` (2026-07-04),
`tex/rnr_coding.tex` = 22480 lines (~383 pp). This manifest is designed so that a
follow-up session can execute the split mechanically. All line numbers are
1-indexed into the file at this commit; re-derive them with the scan commands in
Appendix A if the file has moved.

Target papers:

| id | working title | content |
|----|---------------|---------|
| P1 | RNR Coding (core) | taxonomy Type-I/II/III(+A/B/C), 5-type achievability x converse matrix, constructions, hardness (SS1-6), composition/mode selection (7.1x), SS8-9, SS10.1-10.5, SS11.1-11.6, SS12.1-12.7, most of SS13 |
| P2 | Random access and deviation theory for RNR archives | SS7 arc EXCEPT 7.34: Theorems 7.2-7.33 + 7.35-7.36 (random access, 7.15 exact-entropy complexity arc, neural realisations 7.27'/7.30/7.30a, deviation hierarchy), SS10.7-10.8, SS11.7, SS12.8, its share of SS13 |
| P3 | Operational rate-distortion dispersion | the ENTIRE SS7.34 family: Remark 7.34 through Remark 7.34q (Gray region, dispersion, third-order program, beyond-Gray, Route B replica spectral program) |

---

## 1. Heading markup conventions (discovered)

The paper sets `\setcounter{secnumdepth}{-\maxdimen}` and **never uses**
`\section`/`\subsection`/`\paragraph`. There is exactly **one** `\label`
(`eq:G1chain`, line 14372) and its 7 `\eqref` uses are all inside the 7.34
family — no other LaTeX cross-referencing machinery exists. All structure is:

1. **Section / subsection headings**: standalone lines
   `\textbf{N. Title}` / `\textbf{N.M Title}` / `\textbf{N.M.K Title}`
   (plus `\textbf{Abstract}`, `\textbf{Keywords:}...`, `\textbf{References}`).
   CAUTION: two headings contain math with braces and are missed by a naive
   `\textbf\{[^}]*\}$` regex: line 1437 `4.6 Variants of $E_{12}$` and line
   2480 `5.5a The Type-II converse`.
2. **Major result blocks**: `\textbf{Theorem|Lemma|Corollary|Definition|Remark|
   Observation N.Mx (...long title...).}` — the `\textbf{` often spans 2-3
   source lines; primes appear as `$'$`, `$''$`, ... inside the number.
3. **Sub-remark blocks**: `\emph{Remark N.Mx (...).}` at line start (the whole
   7.15a-i arc, 7.25a/b, 7.26a/b, 7.27a/b, 7.28a/b, 7.29a, 7.30, 7.30a, 7.31a,
   7.32a-inline, and — critically — **Remark 7.34 and Remark 7.34a**, which
   open the P3 family in italic, not bold, markup).
4. Mid-statement bold that is NOT a heading: `\textbf{converse}` (14816),
   `\textbf{achievability}` (15633), `\textbf{Step k (...)}` proof steps,
   emphasis like `\textbf{resolved negatively}` (10894). Any mechanical pass
   must anchor on the `(Theorem|Lemma|...) <number>` pattern, not on bare bold.
5. In-text citations are **plain author-year text** ("Kieffer--Yang 2000",
   "(Gray 1970; Zhu--Alajaji)"). The reference list (`\textbf{{[}1{]}}` ...
   `\textbf{{[}135{]}}`, lines 21931-22480) is numbered but the numbers are
   **never cited in the text** — splitting the bibliography is an
   author-surname matching job, not a `\cite`-key job.
6. Document skeleton: preamble lines 1-57 (pandoc-style, only two custom
   macros: `\tightlist` and a `\setminus` remap), `\begin{document}` line 58,
   title block 60-67, Abstract 69-268, Keywords 270-273, `\end{document}` 22480.

---

## 2. Master assignment table (top level)

`end` = line before the next heading. share = of 22480 lines.

| lines | heading | share | paper |
|-------|---------|-------|-------|
| 1-57 | preamble (pandoc skeleton) | 0.3% | SHARED (copy verbatim into all three) |
| 58-68 | `\begin{document}` + title block | 0.0% | SHARED-ADAPT (new subtitle per paper) |
| 69-273 | Abstract + Keywords | 0.9% | REWRITE per paper (see SS6.1) |
| 275-369 | 1. Introduction / 1.1 Motivation | 0.4% | P1 (trim; P2/P3 get fresh intros) |
| 370-427 | 1.2 Position relative to existing work | 0.3% | P1 |
| 428-839 | 1.3 Contributions | 1.8% | SPLIT: bulk P1; bullet 599-604 (7.15 arc) -> P2; bullet 745-800+ ("27-theorem" arc) -> P2 |
| 840-918 | 1.4 Organization | 0.4% | REWRITE per paper (lines 848-896 describe SS7 incl. 7.34) |
| 919-1119 | 2. Preliminaries and notation (2.1-2.7) | 0.9% | SHARED-DUP (all three; trim per paper, SS7) |
| 1120-1249 | 3. The RNR framework (3.1-3.4, Defs 3.1-3.5) | 0.6% | P1 full; P2/P3 condensed recap (SS7) |
| 1250-1803 | 4. Type-I (4.1-4.7, Thm 4.1/4.2/4.3, Rem 4.3a-d) | 2.5% | P1 |
| 1804-4714 | 5. Type-II (5.1-5.8; 5.1a-f cold-context; 5.6a-j TC; 5.7a-d) | 12.9% | P1 (Lemmas 5.1a-e restated in P2 prelims) |
| 4715-6757 | 6. Type-III (6.1-6.5; hardness 6.7/6.8x; 6.9 cluster is in SS13.4) | 9.1% | P1 |
| 6758-6971 | 7. Composition and mode selection: intro + Thm 7.1/7.1b/7.1c + Rem 7.1d + SS7.2 runtime-parameter text | 1.0% | P1 (recommended; see Open Decision D1) |
| 6972-12524 | Theorems 7.2-7.33 (block-level table in SS3) | 24.7% | P2 |
| 12525-15940 | Remark 7.34 ... Remark 7.34q (block-level table in SS3) | 15.2% | **P3 (contiguous!)** |
| 15941-16177 | Theorems 7.35, 7.36 | 1.1% | P2 (note: file order puts these AFTER the P3 range; P2 extraction is 2 disjoint ranges here) |
| 16178-16197 | 8. No universal dominance (Thm 8.1) | 0.1% | P1 (P2 cites companion for Cor 7.31a's use of 8.1) |
| 16198-16271 | 9. Verification and error containment (Thm 9.1) | 0.3% | P1 (P2 restates/cites Thm 9.1; used at 9479, 10150f, 17357) |
| 16272-16745 | 10.1-10.5 practical: bit-exact, amortization (Thm 10.1/10.2, Def 10.1-10.3, Cor 10.2a, Thm 10.6 at 17570 is in 10.8 range) | 2.1% | P1 (DIFFICULT: 10.2's tables quote T7.15/T7.21 figures, lines 16322-16392 — keep, cite companion P2) |
| 16746-17469 | 10.7 Random access (Thm 10.3, Rem 10.3a, Thm 10.4) | 3.2% | P2 |
| 17470-17835 | 10.8 Block-device and filesystem deployment (Thm 10.5, Thm 10.6) | 1.6% | P2 |
| 17836-17944 | 11. Experimental predictions, 11.1-11.6 | 0.5% | P1 (line 17863 cites Thm 7.15 -> companion cite) |
| 17945-17975 | 11.7 Large-scale neural predictors (T7.15 instantiation) | 0.1% | P2 |
| 17976-18049 | 12. Related work, 12.1-12.7 | 0.3% | P1 |
| 18050-18137 | 12.8 Foundations used by the SS7 extensions | 0.4% | P2 (all its content is 7.14/7.15/7.17/7.20/7.21/7.23/7.24; zero 7.34 mentions) |
| 18138-18484 | 13. Status of open problems + 13.1 Resolved problems | 1.5% | SPLIT by paragraph: OP1/5.1a-e narrative (~18356) -> P1; SS7.2-7.28 narrative (18405-18484) -> P2 |
| 18485-19141 | 13.2 Remaining open problems: OP2/OP4 status ledger + items through "compatible Adam" | 2.9% | P1 |
| 19142-19244 | 13.2 (cont.): gamma_LRU item (19143) + Thm 7.27 sub-block-restriction item (19164-19244) | 0.5% | P2 (these are `\item`s of one itemize — split at `\item` boundaries 19142/19164) |
| 19245-19338 | 13.3 Engineering deferrals | 0.4% | SPLIT by item: 7.16-7.18/7.20/7.23/7.24/7.10-7.12 items -> P2; Thm-10.6/hardware items -> P1 |
| 19339-21928 | 13.4 Attack vectors: 13.4.1 OP2(a), 13.4.2 OP2(b), 13.4.3 OP3'', 13.4.4 OP4 (+ Lemma 6.9/6.9a-d cluster 21261-21830), 13.4.5 TC scaling, 13.4.6 method | 11.5% | P1 entirely (zero 7.x-dependency above 7.15; 7.34 never mentioned) |
| 21929-22480 | References ([1]-[135]) + `\end{document}` | 2.5% | SHARED-SUBSET per paper (SS8) |

Resulting sizes: P1 ~11.0k lines (49%), P2 ~7.2k (32%), P3 ~3.4k (15%),
shared/front ~0.8k (4%).

Numbering gaps that are REAL (do not hunt for missing headings): there is no
subsection 4.6->4.7 anomaly (4.6 exists, line 1437); there is no `\textbf{10.6
...}` subsection heading (numbering jumps 10.5 -> 10.7; Theorem 10.6 lives
inside 10.8); SS7 has only two subsection headings (`7.` at 6758, `7.2` at
6930) — everything else in SS7 is theorem blocks.

---

## 3. SS7 block-level table (the P2/P3 frontier)

`end` = next block's start - 1. Assignment: P2 unless marked P3.

| start | block | notes |
|-------|-------|-------|
| 6758 | `7.` section heading + intro | P1 (D1) |
| 6767 | Theorem 7.1 | P1 (composition correctness) |
| 6836 | Theorem 7.1b | P1 |
| 6875 | Theorem 7.1c | P1 |
| 6917 | Remark 7.1d | P1 |
| 6930 | `7.2 Runtime parameter selection` subsection text (ends 6971) | P1 |
| 6972 | Theorem 7.2 | P2 from here on |
| 7077 | Corollary 7.2a |  |
| 7147 | Theorem 7.3 |  |
| 7242 | Theorem 7.4 |  |
| 7331 | Theorem 7.5 |  |
| 7432 | Theorem 7.6 |  |
| 7520 | Theorem 7.6$'$ |  |
| 7660 | Theorem 7.7 (lossy RD mode — stays P2; P3 cites it) |  |
| 7743 | Theorem 7.8 |  |
| 7837 | Theorem 7.9 |  |
| 7941 | Theorem 7.10 |  |
| 8037 | Corollary 7.11 |  |
| 8115 | Theorem 7.12 |  |
| 8224 | Theorem 7.13 |  |
| 8331 | Theorem 7.14 |  |
| 8498 | Theorem 7.15 (+ `\emph` Remarks 7.15a 8681, b 8807, c 8868, d 8931, e 8975, f 9028, g 9087, h 9187, i 9213) |  |
| 9268 | Theorem 7.16 |  |
| 9360 | Theorem 7.17 |  |
| 9485 | Theorem 7.18 |  |
| 9583 | Theorem 7.19 |  |
| 9686 | Theorem 7.20 |  |
| 9823 | Theorem 7.21 |  |
| 9927 | Theorem 7.22 |  |
| 10034 | Theorem 7.23 |  |
| 10158 | Theorem 7.24 |  |
| 10271 | Theorem 7.25 (+ 7.25a 10371, 7.25b 10384) |  |
| 10412 | Theorem 7.26 (+ 7.26a 10534, 7.26b 10564) |  |
| 10576 | Theorem 7.27 (+ Lemma 7.27a 10641, Lemma 7.27b 10675) |  |
| 10870 | Theorem 7.27$'$ |  |
| 10972 | Remark 7.27d |  |
| 11016 | Definition 7.27b |  |
| 11037 | Corollary 7.27c |  |
| 11134 | Corollary 7.27e |  |
| 11288 | Corollary 7.27f |  |
| 11443 | Theorem 7.28 (+ 7.28a 11655, 7.28b 11684) |  |
| 11698 | Theorem 7.29 |  |
| 11781 | Remark 7.29a (`\emph`) |  |
| 11822 | Remark 7.30 (`\emph`; bisection / FIM coding) |  |
| 11900 | Remark 7.30a (`\emph`; neural FV escape) |  |
| 11979 | Theorem 7.31 (+ Corollary 7.31a 12128, uses Thm 8.1) |  |
| 12189 | Theorem 7.32 |  |
| 12331 | Remark 7.32a |  |
| 12413 | Theorem 7.33 (ends 12524) |  |
| **12525** | **Remark 7.34** (`\emph`!) | **P3 starts here** |
| 12572 | Remark 7.34a (`\emph`; Gauss-Markov) | P3 |
| 12650 | Remark 7.34b (contains unnumbered `\textbf{Theorem (the replica spectral inequality...)}` 12930, proof Steps 0-7 13324-13646, unnumbered `\textbf{Theorem (BSMS Gray-region rate-distortion dispersion)}` 13647) | P3 |
| 13806 | Remark 7.34d | P3 |
| 13941 | Lemma 7.34e | P3 |
| 14136 | Remark 7.34e$'$ | P3 |
| 14243 | Lemma 7.34f | P3 |
| 14317 | Lemma 7.34g (+ bold part-head `(IV) Existence below $D_c$` 14357) | P3 |
| 14417 | Lemma 7.34h (the file's only `\label{eq:G1chain}` 14372 sits just above, in 7.34g's display; all 7 `\eqref` uses are within 14302-14525) | P3 |
| 14529 | Theorem 7.34i | P3 |
| 14615 | Lemma 7.34j | P3 |
| 14669 | Lemma 7.34k | P3 |
| 14715 | Lemma 7.34l | P3 |
| 14810 | Theorem 7.34m (statement split by bold `converse` 14816) | P3 |
| 14887 | Remark 7.34m$'$ | P3 |
| 14983 | Remark 7.34m$''$ | P3 |
| 15240 | Remark 7.34m$'''$ (residual ledger TABLE — P3's own open-problem section; includes "witnessing script" column naming verify scripts) | P3 |
| 15342 | Remark 7.34m$''''$ (sharp-floor campaign I) | P3 |
| 15491 | Remark 7.34m$'''''$ (sharp-floor campaign II) | P3 |
| 15617 | Theorem 7.34n (statement split by bold `achievability` 15633) | P3 |
| 15686 | Remark 7.34n$'$ | P3 |
| 15727 | Remark 7.34o | P3 |
| 15794 | Remark 7.34o$'$ | P3 |
| 15855 | Remark 7.34p (deviation-hierarchy transfer; cites 7.31/7.33/7.35/7.36 -> companion) | P3 |
| 15894 | Remark 7.34q (ends **15940**) | P3 |
| 15941 | Theorem 7.35 (ends 16066) | P2 |
| 16067 | Theorem 7.36 (ends 16177) | P2 |

Notes: there is no block numbered 7.34c (gap is real). The whole P3 range is
**one contiguous slice 12525-15940** — the cleanest cut in the paper.
Confirmed: `7.34` appears NOWHERE outside 12525-15940 except SS1 (lines
875-896, Organization) — SS13 never mentions it (P3's open problems are
self-contained in the 7.34m$'''$ ledger).

---

## 4. Difficult ranges (interleaved / rewrite-required)

| lines | issue | resolution |
|-------|-------|-----------|
| 69-273 | Abstract is a P1+P2 blend (cites Lemma 5.1b, SS10.7 sync overhead, Cor 10.2a, Thm 10.1, Thm 8.1, H10(c) of the experimental-design companion). 7.34 is NOT in the abstract. | Write three fresh abstracts. P1: drop the sync-point/SS10.7 sentences (230-252). P2: build from those sentences + 7.15/7.29 material. P3: fresh (source: 7.34b/i/m statements + 875-896). |
| 428-839 | 1.3 Contributions mixes papers at bullet granularity. | Reassign bullets: 599-604 (7.15a-i arc) and 745-800+ (27-theorem bullet) -> P2 intro; everything 7.34 is in 1.4 not 1.3; rest -> P1. |
| 840-918 | 1.4 Organization describes all sections; 848-896 is one giant sentence covering SS7+7.34. | Rewrite per paper. |
| 6930-6971 | Subsection heading `7.2 Runtime parameter selection` collides with `Theorem 7.2` (different objects, same number). If 7.2-text goes to P1 and Theorem 7.2 to P2, both papers have a "7.2". | Acceptable if numbering is kept as-is (D2); flag in each paper's notation note. |
| 16311-16393 | 10.2 amortization regimes: the break-even table hard-codes T7.15 ($\approx$79 MB / 12.0x) and T7.21 ($L^*=\Theta(N^{1/(1+\alpha)})$) figures. | Keep in P1, convert "Theorem 7.15/7.21" to companion citations "[P2, Thm 7.15]". Do NOT strip the numbers. |
| 18150-18484 | 13.1 resolved-problems prose interleaves P1 items (OP1 via 5.1a-e) and P2 items (SS7.2-7.28 narrative, 7.27 saga, "27 theorems verified" 18478-18484). | Split at paragraph level; the SS7 narrative block is 18405-18484. |
| 19132-19244 | 13.2's closing `\item`s: Adam item (19133-19141, P1), gamma_LRU item (19143-19163, P2), 7.27-restriction item (19164-19244, P2). One itemize spans the P1/P2 cut. | Cut between `\item` markers at 19142 and keep list syntax valid in both halves. |
| 19245-19338 | 13.3 engineering deferrals: items reference Theorems 7.16-7.18, 7.20/7.24, 7.23, 7.10-7.12 (P2) but also 10.6/hardware (P1). | Item-level triage (7 items; ~5 to P2, ~2 to P1). |
| 8498-9267 | Theorem 7.15 block embeds the whole 7.15a-i complexity arc as `\emph` remarks; 7.15a leans on SS6.9/Lemma 6.9 (P1, itself located in SS13.4.4!). | Keep whole block in P2. P2 must cite P1 for Lemma 6.9/Remark 6.9a (see SS6.3) — note Lemma 6.9 lives in P1's SS13.4.4, so P1 must keep that cluster (it does). |
| 12930-13646 | Unnumbered bold Theorems inside Remark 7.34b (replica spectral inequality; BSMS dispersion) + proof Steps 0-7. | Travel with 7.34b into P3 untouched; if P3 is renumbered these become numbered theorems (D2). |
| 15941-16177 | 7.35/7.36 sit AFTER the P3 slice; they cite 7.28/7.29/7.31/7.33 only (verified: zero 7.34 references). | P2 extraction = ranges (6972-12524) + (15941-16177); pure concatenation is safe. |

---

## 5. Cross-reference inventory (verified by scan)

LaTeX-level: only `eq:G1chain` (label 14372, 7 eqrefs, all inside P3). Nothing
else to repair at the LaTeX level. Everything below is plain-text mentions.

### 5.1 P3 -> P2 (must become companion citations)

Mentions inside 12525-15940 of non-7.34 results — complete list:
Theorem 7.28 (x4), Theorem 7.29 (x8), Theorem 7.31 (x1), Theorem 7.33 (x1),
Theorem 7.35 (x1), Theorem 7.36 (x1), SS7.15c (x1), SS7.31 (x2), SS7 (x1).
Also first-order inputs referenced by name: Theorem 7.7 (lossy mode),
delta_infty / cold-context notation from Theorem 7.29.

Repair: P3 preliminaries restate (statement-only, with "[P2, Thm 7.28/7.29]"
provenance): the Berry-Esseen CLT (7.28), the second-order dispersion +
random-access interplay (7.29), and the deviation-hierarchy pointers
(7.31/7.33/7.35/7.36 — cite-only, used in Remarks 7.34/7.34p).

### 5.2 P3 -> P1

None by theorem number (verified: no Theorem/Lemma 2.x-6.x mentions inside
12525-15940). P3 needs only the generic framework vocabulary (repair length
L^rep, sub-block K, conditionally-advantageous predicate) — covered by the
condensed recap (SS7.3).

### 5.3 P2 -> P1 (must become companion citations or restated prelims)

Complete mention counts inside P2 ranges:
- Cold-context lemmas: Lemma 5.1b (x9+3), 5.1c (x2+2), 5.1d (x7+1), 5.1e (x8),
  5.1a (x2), Corollary 5.1f (x1) — load-bearing for Theorems 7.2/7.3/7.29 and
  SS10.7. RESTATE statements in P2 preliminaries.
- TC lemmas: Lemma 5.6b/c/d/e/f (x1/3/1/3/1) — cite companion.
- SS6.9 cluster: Lemma 6.9 (x6), SS6.9 (x6), Remark 6.9a (x2), Corollary 6.9b
  (x2) — used by Remark 7.15a arc. Cite companion (they live in P1 SS13.4.4).
- Theorems 4.2, 5.3, 5.5, 6.1, 7.1 (x1 each in 6972+ range), Theorem 5.4 (x1 in
  10.7), SS6.3 strategies (i)/(ii)/(iii) (x9 in 10.7 — P2 must include a short
  recap of the three Type-III-B coding strategies or the SS10.7 tables dangle).
- Theorem 8.1 (used by Corollary 7.31a, 12129/12148), Theorem 9.1 (used at
  9479, 10150-10151, 17357) — restate statements in P2 prelims.
- SS2.5/SS2.6 (MDL, bits-back) mentioned once each — covered by shared prelims.

### 5.4 P1 -> P2/P3 (forward references out of the core)

- Abstract lines 230-252 (SS10.7 sync-overhead discussion) -> moves to P2.
- 1.3 bullets 599-604, 745-800+ -> move to P2 intro.
- 1.4 lines 848-896 -> rewritten (P2+P3 trailers become one "companion papers"
  paragraph).
- SS10.2 tables 16322-16392: T7.15/T7.21 figures -> "[P2]" citations, keep.
- SS11 line 17863 (predictions table "derives from Theorem 7.15") -> "[P2]".
- SS13.1 narrative 18405-18484 -> moves to P2's own status section.
- SS12.8 (18050-18137) -> moves to P2 related work.
- Zero 7.34 references anywhere in P1-bound text (verified) — P1 needs only a
  single "companion paper III" sentence in the rewritten 1.4.

### 5.5 Companion-citation mechanics (recommended)

Because references are plain text, the repair is a set of literal string
substitutions applied ONLY to the receiving paper's extracted body, e.g. in P3:
`Theorem 7.29` -> `Theorem 7.29 of [RNR-RA]`  (first occurrence per block;
subsequent occurrences can stay bare since the numbering is preserved by D2).
Add to each paper's reference list: entries for the two companions, e.g.
`[C1] Ratushnyi, P. (2026). RNR Coding I: ... (companion paper).` Cross-file
integrity is then still checkable by `scripts/verify/crossref_integrity.py`
after adding the new tex files to its file set (it already resolves against
the UNION of all documents' defined blocks).

---

## 6-7. Preliminaries duplication list (what each paper must carry)

SHARED VERBATIM (all three): preamble lines 1-57; title-block skeleton 58-68
(subtitle per paper); SS2.1 lossless problem + Definition 2.1
(conditionally-advantageous predicate, lines 921-955) — the evaluation frame
every paper's claims are scoped by; SS2.7 notation summary (1094-1119).

P1 keeps all of SS2 (919-1119) and SS3 (1120-1249) as-is.

P2 additionally needs:
- Condensed SS3 recap: Definitions 3.1-3.5 (five types, statements only) +
  SS3.3 repair-length decomposition (L^rep, L(Theta), cost(M)); ~1 page.
- SS2.2 side-information, SS2.3 entropy coding, SS2.5 MDL, SS2.6 bits-back
  (7.6'/7.15 use them) — take from 956-1093.
- Restated statements: Lemmas 5.1a-e (+ Cor 5.1f pointer), Theorem 4.2
  (Type-I converse baseline), Theorem 8.1, Theorem 9.1, Type-III-B strategy
  list (i)/(ii)/(iii) from SS6.3.
- Definition of sub-block architecture / sync points: Definitions 10.2, 10.3
  (16783, 16797) MOVE to P2 (they sit in P1's 10.5 range but are the
  random-access vocabulary; P1 keeps a one-line pointer) — or duplicate.
  Flag as Open Decision D5.

P3 additionally needs:
- Condensed RNR recap (half page: repair principle, sub-block K, L^rep).
- Restated: Theorem 7.28 statement, Theorem 7.29 statement (with delta_infty
  cold-context definition), pointers to 7.31/7.33/7.35/7.36.
- The lossy operational model: Theorem 7.7's byte-granular RD setup
  (statement only).
- d-tilted information / dispersion definitions: already self-contained
  inside 7.34b/7.34d text (verified — the blocks define what they use);
  no extra import needed beyond notation for V_op, V_lossless.

---

## 8. Bibliography plan

Master list: 135 numbered entries (21931-22478), alphabetical, never cited by
number. Split method: surname+year text match per paper body (script in
Appendix A generated the current pass; heuristic — accents and multi-author
cites need a manual pass).

Heuristic subsets (entry numbers, current pass):
- P1 (86): 1-14, 16-25, 27-29, 31-35, 39-45, 47-50, 52, 58, 59, 69-71, 78, 83,
  84, 87, 88, 93-111, 113, 114, 117, 121, 125, 128-135
- P2 (85): 2-4, 7-10, 16, 21, 22, 25, 27, 28, 31, 32, 40, 42-44, 49, 52-59,
  61, 63-96, 99, 101, 108-111, 113-127, 130
- P3 (12): 3, 31, 52, 56, 75, 81, 88, 100, 112, 123, 125, 130
- Unmatched by the heuristic (manual assignment needed): 15 (Gallager 1962),
  26 (Reed 1960), 30 (Sahni 1976), 36-38 (parse failures — inspect), 46
  (Studeny 1998), 51 (Ar... 2009), 60 (Csiszar 2011), 62 (Deletang 2024 — P2:
  cited as "Delétang et al. 2024" in 7.15/11.7), plus all accented surnames
  (Verdu, Csiszar, Renyi forms) which the match missed.

**P3 gap — MISSING ENTRIES.** The 7.34 family cites ~25 works inline that are
NOT in the master reference list at all (verified by scanning entries for the
surnames): Kostina--Verdu 2012 (x4), Kostina 2017 (x3), Kontoyiannis--Verdu
2014, Gray 1971 (only Gray 1970-ish entries 73/100/131 exist — check which),
Zhu--Alajaji, Hafouta 2020/2025, Dolgopyat--Sarig 2023, Dolgopyat--Hafouta
2022, Krishnamachari 2026 (x8), Tasci--Kostina 2026, Eswaran--Gastpar 2022,
Elshafiy--Rose 2022, Paulin 2015 (=123, exists), Kontorovich--Ramanan 2008
(=112, exists), Diaconis--Saloff-Coste 1996, Kingman 1961, Nussbaum 1986,
Hennion 1997, Lalley 1989, Hipp 1983, Lahiri 1993, Avram--Berger 1985,
Anantharam--Borkar 2019 (=107, exists), Aaronson--Denker 2001, Cohen 1981,
Kato 1982. **P3's reference list must be composed largely from scratch** —
budget a dedicated pass; the same is true to a smaller degree for P2 (check
the 7.31-7.36 blocks for inline-only citations the same way).

---

## 9. scripts/verify mapping (234 files)

Filename convention maps scripts to papers nearly 1:1. Rules:

| pattern | paper | count (approx) |
|---------|-------|-------|
| `thm_4_*`, `thm_5_*`, `thm_6_*`, `lemma_5_*`, `lemma_6_*`, `cor_6_9a*`, `remark_4_3*`, `remark_6_8m*`, `remark_6_9*`, `observation_6_8f*`, `prop_5_7*`, `op1_5_*`, `op2a_*`, `op3_*`, `op4_*`, `empirical_tc_*`, `defs_3_round_trip`, `thm_10_1_*`, `thm_7_1bc_*` (7.1b/c stay P1 per D1) | P1 | ~75 |
| `thm_7_*` (7_2..7_33, 7_35, 7_36; NOT 7_34*), `remark_7_15*`, `remark_7_29a*`, `remark_7_30a*`, `remark_7_32a*`, `cor_7_11*`, `cor_7_27f*`, `thm_7_30_bidirectional*`, `bug_008_random_access_cost` (Thm 10.3) | P2 | ~55 |
| ALL `*7_34*`: `lemma_7_34*`, `probe_7_34*`, `remark_7_34*`, `theorem_7_34i*`, `thm_7_34_*`, plus `bug_009_ak_decay_convexity` (the a_k/convexity residual) — includes the 17 currently-untracked `probe_7_34_route_b_*` / `lemma_7_34m_curvature_routeE.py` files in git status | P3 | ~100 |
| `crossref_integrity.py` — resolves refs against the UNION of all four docs; MUST be extended with the three new tex files | ALL | 1 |
| `bug_001_ratio_claim_conditional.py`, `bug_010_precondition_invariant.py` — document-invariant probes that PARSE rnr_coding.tex (abstract/SS1.1/SS10.2.1/SS11 cross-links); they hard-break when those sections move | ALL (retarget) | 2 |
| `bug_002/003/004/005` (accounting, sync-overhead, ideal-vs-achievable, status-index — parse SS10.2/abstract/SS13.2) | P1 (bug_003 -> P2 if abstract sync text moves there) | 4 |

15 spot-checks performed (docstring section references confirm the rules):
bug_003 (abstract sync-overhead -> P2-leaning), bug_008 (Thm 10.3 -> P2),
bug_009 (BUG-009-D a_k residual -> P3), crossref_integrity (all four docs),
defs_3_round_trip (SS3 -> P1), empirical_tc_2d_ising (SS13.4.5 -> P1),
thm_10_1 (SS10.5 -> P1), thm_7_1bc (7.1b/c -> P1 per D1), op1_5 (SS13.4.4/5.7
-> P1), bug_001 (abstract/SS11 invariants -> retarget), bug_005 (SS13.2 ledger
-> P1), bug_010 (abstract/SS10.2.1 -> retarget), thm_7_34_lossy_rnr_dispersion
(Remark 7.34 -> P3), remark_7_32a (7.32a -> P2), probe_7_34_routeB_adversarial
(BUG-009-D -> P3).

CI note: `.github/workflows/build.yml` builds the PDFs and runs the verify
job; after the split it needs three build targets. CI is NOT observable from
this account (`gh` mismatch) — verify locally with the Docker
`pandoc/latex:latest` flow before pushing.

---

## 10. Execution recipe (mechanical)

Order: **P3 first** (cleanest cut, one contiguous range), then P2, then P1
(P1 = what remains + rewrites). Keep `tex/rnr_coding.tex` untouched until all
three compile (it is the reference for `crossref_integrity.py` during the
transition).

Step 0 — freeze. Tag the pre-split commit. Do all work on a branch.

Step 1 — extract bodies (python, ranges are inclusive):

```python
# /Users/para/.venvs/rnr/bin/python
SRC = 'tex/rnr_coding.tex'
lines = open(SRC).read().split('\n')          # lines[i-1] == line i
def take(*ranges):
    return '\n'.join('\n'.join(lines[a-1:b]) for a, b in ranges)

P3_BODY = take((12525, 15940))
P2_BODY = take((6972, 12524), (15941, 16177),          # theorem arc
               (16746, 17835),                          # 10.7 + 10.8
               (17945, 17975),                          # 11.7
               (18050, 18137),                          # 12.8
               (18405, 18484), (19142, 19244))          # 13.1/13.2 P2 items
P1_BODY = take((275, 6971), (16178, 16745),
               (17836, 17944), (17976, 18049),
               (18138, 18404), (18485, 19141), (19245, 21928))
# NOTE: (19245,19338) 13.3 needs the item-level triage of SS4 before final cut;
# (18138,18404)/(18485,19141) carry the SS13 headings for P1.
PREAMBLE = take((1, 57)); REFS = take((21929, 22479))
```

Step 2 — assemble each paper:
`PREAMBLE + \begin{document} + [new title block] + [new abstract+keywords] +
[new SS1 intro per SS4] + [prelims per SS6-7] + BODY + [status/open-problems
section] + [references subset per SS8] + \end{document}`.
Suggested files: `tex/rnr_core.tex`, `tex/rnr_random_access.tex`,
`tex/rnr_rd_dispersion.tex` (final names: Open Decision D7).

Step 3 — companion-citation substitutions (apply to extracted bodies only):
- In P3 body: prepend restated prelims (SS7); first mention per block of
  `Theorem 7.28|7.29|7.31|7.33|7.35|7.36` gets ` of [RNR-RA]`; `SS7.15c` ->
  `[RNR-RA, SS7.15]`.
- In P2 body: `Lemma 5.1[a-e]` -> keep (restated in P2 prelims, note "proved
  in [RNR-I]"); `Lemma 6.9|Remark 6.9a|Corollary 6.9b|SS6.9` -> `[RNR-I]`;
  `Theorem (4.2|5.3|5.4|5.5|6.1|8.1|9.1)` -> restated-or-companion per SS5.3;
  `SS6.3` (in 10.7) -> point at the strategy recap.
- In P1 body: SS10.2 table + SS11 mentions of `Theorem 7.15|7.21` ->
  `[RNR-RA]`; drop/replace the moved 1.3 bullets and 1.4 SS7 sentence.

Step 4 — compile: Docker `pandoc/latex:latest` (the repo's local check flow),
all three PDFs. Expect zero undefined references by construction (only
`eq:G1chain`, internal to P3). Grep each PDF log for `Missing character` /
overfull regressions only.

Step 5 — verify guards:
- Update `crossref_integrity.py` file list (add three new tex files; keep
  union semantics).
- Retarget `bug_001`/`bug_010` (and `bug_003` if the abstract sync text moved)
  to the new file(s).
- Run the full `scripts/verify` suite locally with
  `/Users/para/.venvs/rnr/bin/python`.
- Update `.github/workflows/build.yml` with the new build targets.

Step 6 — only after all three compile and guards pass: decide the fate of
`rnr_coding.tex` (keep as frozen monolith vs delete; Open Decision D8) and
commit (stealth rules: no tool references, author Pavlo Ratushnyi, no
Co-Authored-By trailers).

---

## 11. Open decisions for the author

- **D1 — SS7 lead (6758-6971).** Recommended: Theorem 7.1/7.1b/7.1c/Remark
  7.1d + the 7.2 runtime-parameter text go to **P1** (they are mode-selection/
  composition, i.e. core; 7.1b cites Type-III-C machinery). Alternative: put
  them in P2 so "Theorem 7.x" numbering is uninterrupted there. Affects
  `thm_7_1bc_rnr_hierarchy_limits.py` assignment.
- **D2 — numbering.** Recommended: KEEP original numbers in all three papers
  ("numbering follows the series; Part I contains SS1-6, ..."). All references
  are plain text; renumbering = thousands of hand-edits with high silent-error
  risk. Cost: P2 starts at Theorem 7.2, P3's results are 7.34x — cosmetically
  unusual, mechanically safe. If renumbering is insisted on, do it AFTER the
  split compiles, one paper at a time, with a regex ledger.
- **D3 — SS13.1/13.3 item triage.** The paragraph/item boundaries are listed
  in SS4; someone must confirm each item's destination (est. 30 min).
- **D4 — SS10.2 break-even tables** stay in P1 with companion citations
  (recommended) or move to P2. Keeping them preserves P1's "when does the
  model pay for itself" story.
- **D5 — Definitions 10.2/10.3 (W-bounded predictor, sync point)** currently
  in P1's 10.5 range but they are P2 vocabulary: move or duplicate?
  Recommended: duplicate (they are 15 lines total; both papers read cleanly).
- **D6 — bibliography.** Approve the per-paper subsets (SS8), assign the 10
  unmatched entries, and budget the from-scratch P3 reference-list pass
  (~25 new entries for inline-only citations).
- **D7 — file names and titles** for the three papers; and whether the series
  gets an umbrella name ("RNR Coding I/II/III").
- **D8 — the monolith.** Keep `rnr_coding.tex` in-tree (frozen, CI-excluded)
  as the historical record, or remove after the split lands? The doc-invariant
  probes and the three companion docs (`rnr_summary`, `rnr_engineering_spec`,
  `rnr_experimental_design`) reference "the main paper" — their pointers need
  a sweep either way (they cite results by number, which D2 preserves).
- **D9 — the three companion docs**: update their cross-links to point at the
  correct member of the new trio (mechanical once D7 is fixed;
  `crossref_integrity.py` will catch stragglers).

---

## Appendix A — regeneration commands

All scans used `/Users/para/.venvs/rnr/bin/python`. Re-derive the heading
table after any upstream edit with:

```python
import re
lines = open('tex/rnr_coding.tex').read().split('\n')
sec   = re.compile(r'^\\textbf\{(\d+(?:\.\d+)*[a-z]?\.?\s+.*)\}\s*$')       # sections (also catches 4.6/5.5a w/ math)
blk   = re.compile(r'^\s*\\textbf\{(Theorem|Lemma|Remark|Corollary|Definition|Observation)\s+([0-9]+(?:\.[0-9]+)?[a-z]?(?:\$[^$]*\$)?)')
sub   = re.compile(r'^\s*\\emph\{(Remark|Lemma)\s+([0-9]+(?:\.[0-9]+)*[a-z]?)')
for i, l in enumerate(lines, 1):
    for p in (sec, blk, sub):
        m = p.match(l)
        if m: print(i, m.group(0)[:90]); break
```

Boundary invariants to re-check before cutting (cheap, high-value):
1. `grep -c '7\.34'` on lines outside 12525-15940 and 275-918 must be 0.
2. No `\label`/`\ref` besides `eq:G1chain` (14372; uses 14302-14525).
3. `\textbf{Theorem 7.35` immediately follows line 15940's block.
4. The `\emph{Remark 7.34 (` opener at 12525 (NOT bold — easy to miss).
