# Three-paper split of `rnr_coding.tex`

`../rnr_coding.tex` (the monolith, 22480 lines) remains the source of truth.
These are **mechanical first-cut copies** produced per `docs/backlog/SPLIT-PLAN.md`,
each keeping the shared preamble + `\begin/\end{document}` + reference list and
deleting the line ranges owned by the other two papers. No content was
rewritten; nothing was renumbered.

| Folder | File | Scope | Body |
|---|---|---|---|
| `core/` | `rnr_core.tex` | Paper I — the RNR framework | §1–6 (Type-I/II/III taxonomy, the five-type achievability×converse matrix, constructions, hardness), §7 lead, §8–9, and §10–13 |
| `random_access/` | `rnr_random_access.tex` | Paper II — random access & deviation theory | §7.2–7.33 + §7.35–7.36, §10.7–10.8 |
| `dispersion/` | `rnr_dispersion.tex` | Paper III — operational RD dispersion | the §7.34 family (Remarks 7.34–7.34q) incl. the Route B replica spectral theorem |

Line ranges kept (1-indexed, against the monolith):

- **core**: 1–6971, 16178–16745, 17836–22480
- **random_access**: 1–273, 919–1119, 6758–12524, 15941–16177, 16746–17835, 21929–22480
- **dispersion**: 1–273, 919–1119, 6758–6971, 12525–15940, 21929–22480

Structural checks passed on all three: exactly one `document` environment,
balanced braces, balanced `description/itemize/enumerate/align/…` environments
(no cut landed mid-environment).

## Still to polish (not blocking the cut)

- **Front matter**: each copy currently carries the monolith's title/abstract/keywords
  (they mention all three types). Rewrite per paper.
- **Section numbering**: literal (`\textbf{8. …}`); after the cuts the numbers have
  gaps (e.g. random_access jumps §2 → §7). Renumber, or keep the original series
  numbering with a note (the manifest recommends keeping the series).
- **Cross-references**: plain-text theorem numbers that now point across papers
  (core still names 7.34/7.35; dispersion names 7.29/7.31/7.33/7.35/7.36) should
  become "Paper II, Theorem 7.29"-style companion citations.
- **References**: the full `[1]`–`[135]` list is duplicated into each copy; prune to
  each paper's cited subset. ~25 works cited by name in the §7.34 body
  (Kostina–Verdú 2012, Tian–Kostina, Hafouta, Dolgopyat–Sarig, …) are absent from
  the master list and need proper entries in the dispersion paper.
- **Coarse tail assignment**: to avoid item-level surgery, `core/` still retains a
  few P2-owned tail items (§11.7, §12.8, some §13 items); move them to
  `random_access/` in polish.
- **verify scripts**: `scripts/verify/` partitions ~1:1 by filename; retarget the
  document-invariant probes (`bug_001`, `bug_010`, `bug_003`, `crossref_integrity`)
  to the three files.
