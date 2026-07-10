#!/usr/bin/env python3
r"""
Cross-reference integrity guard for the RNR manuscripts, in two groups.

GROUP A -- the monolith rnr_coding.tex and the three companion docs (rnr_summary /
rnr_engineering_spec / rnr_experimental_design). A reference in any of these four
resolves against the UNION of every Group-A document's defined blocks, because the
companions legitimately cite the main paper's results.

GROUP B -- the three-paper split (tex/papers/core/rnr_core.tex,
tex/papers/random_access/rnr_random_access.tex,
tex/papers/dispersion/rnr_dispersion.tex). Each split paper is checked as an
INDEPENDENT document: a reference must resolve to a block defined IN THE SAME FILE
(prerequisite restatements count -- they are real \textbf{...} headers), OR be a
companion citation, i.e. carry an [RNR-I]/[RNR-II]/[RNR-III] tag within +/-2 lines
of the reference (the split's convention for results proved in a companion Part).
A bare cross-Part reference with neither is a phantom for THAT paper -- exactly the
defect class the 2026-07-07 split audit found (AUDIT-SPLIT-2026-07-07.md items
1, 10-13, 18-22).

The paper numbers its structure MANUALLY -- named blocks are
\textbf{Definition N.M (...)} / \textbf{Theorem N.M (...)} / \textbf{Lemma ...} /
\textbf{Proposition ...} / \textbf{Corollary ...} / \textbf{Remark ...} (a few via
\emph{...}). Cross-references are therefore PROSE ("see Theorem 8.1", "(Remark
7.34o)") and NOT LaTeX \ref/\label, so a pointer to a non-existent object compiles
with NO warning. Two such dangling pointers were once present (the phantom "Remark
7.34c"; "Definition 2.1" cited six times before any such header existed). This
probe makes that failure mode mechanical.

CHECK -- phantom named-block references. For every prose reference "<Kind>
N.M[suffix]" it checks that SOME named block numbered N.M[suffix] exists. The
match is KIND-AGNOSTIC: the paper sometimes calls a "theorem"-titled Remark a
Theorem in a list ("Theorems 7.34b/i/m/n", where 7.34b/i are Remarks) -- the
NUMBER is what must resolve to a real header, not the kind word. A reference whose
number matches no block, and which is not an external citation, is a phantom and
fails the probe.

False-positive controls (so the gate is trustworthy):
  - Primes are matched ONLY as $'$/$''$/$'''$ (math mode); a bare apostrophe is an
    English possessive ("Theorem 5.4's"), so the number is taken as "5.4".
  - A block number is two-part (N.M[letters][primes]) and must NOT be followed by
    a further ".K" -- a three-part number is the signature of an EXTERNAL citation
    (Cover & Thomas 2006 Theorem 10.3.1) and is skipped entirely.
  - A reference preceded (within ~48 chars) by an author/year citation cue (a
    19xx/20xx year, "et al", "&", a known external author, "ed.") is external and
    is skipped.
  - A small, documented ALLOWLIST covers the residue: external theorems cited
    without a nearby cue, and objects DELIBERATELY mentioned as retracted.

Section-reference (§N.M) integrity is deliberately OUT OF SCOPE: §-subsection
header conventions in this paper vary (some groupings carry no \textbf{N.M Title}
header) and many §-refs are external (Wyner--Ziv 1976 §3.5), which would make a
§-gate untrustworthy. The two real defects this guard targets were both named-
block references.

Stdlib only. Prints one PASS/FAIL line per document and "OVERALL -> PASS" on
success; exits nonzero if ANY document fails.
Run:  /Users/para/.venvs/rnr/bin/python scripts/verify/crossref_integrity.py
"""

import bisect
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
TEX = os.path.join(REPO, "tex")
# Group A: monolith + companions. A reference in any one resolves against the
# UNION of every Group-A document's defined blocks, because the companion docs
# legitimately cross-reference the main paper (e.g. rnr_experimental_design.tex
# cites "Theorem 7.15" of rnr_coding.tex). A number that resolves in no Group-A
# document is a phantom.
DOCS = [
    os.path.join(TEX, "rnr_coding.tex"),
    os.path.join(TEX, "rnr_summary.tex"),
    os.path.join(TEX, "rnr_engineering_spec.tex"),
    os.path.join(TEX, "rnr_experimental_design.tex"),
]

# Group B: the three-paper split. Each is INDEPENDENT -- own defined blocks only,
# plus the [RNR-x] companion-citation convention (tag within +/-2 lines).
PAPERS = os.path.join(TEX, "papers")
SPLIT_DOCS = [
    os.path.join(PAPERS, "core", "rnr_core.tex"),
    os.path.join(PAPERS, "random_access", "rnr_random_access.tex"),
    os.path.join(PAPERS, "dispersion", "rnr_dispersion.tex"),
]

# Companion-citation tag: [RNR-I] / [RNR-II] / [RNR-III], inline or in the
# LaTeX-escaped reference-list form {[}RNR-II{]}. Case-sensitive so that prose
# like "RNR-independent" / "RNR-specific" does not match.
COMPANION_TAG = re.compile(r"RNR-I{1,3}(?![A-Za-z])")
COMPANION_WINDOW = 2  # lines on each side

KINDS = ["Definition", "Theorem", "Lemma", "Proposition", "Corollary", "Remark",
         "Conjecture", "Observation"]

# Internal block number: N.M, optional letter suffix, optional math-mode primes,
# and NOT the prefix of a three-part (external) number.
BNUM = r"[0-9]+\.[0-9]+[a-z]*(?:\$'+\$)?(?!\.[0-9])"

DEF_BLOCK = re.compile(r"\\(?:textbf|emph)\{(?:" + "|".join(KINDS) + r")~?\s+(" + BNUM + r")")
REF_BLOCK = re.compile(r"(" + "|".join(KINDS) + r")s?~?\s+(" + BNUM + r")")

# Citation cue in the ~48 chars preceding a reference => external, skip.
EXTERNAL_CUE = re.compile(
    r"(?:19|20)\d\d"                              # a year
    r"|et al|&|\bed\.|2nd ed|, ?p\.|pp\.|Press"   # bibliographic cues
    r"|Cover|Thomas|Gray 1|Muirhead|Ibragimov|Lezaud|Dembo|Zeitouni"
    r"|Aldous|Wyner|Ziv|Slepian|Charikar|Watanabe|Paulin|Bryc|Petrov"
    r"|Tian|Kostina|Polyanskiy|Barron|LPW|Levin|Peres|Bobkov|Marton"
)
PRECEDING = 48

# Documented residue. Each entry is a normalized block number with a reason it
# legitimately has no internal header. Kept tiny and justified so the gate stays
# meaningful.
ALLOWLIST_BLOCK_NUM = {
    "1.1": "external: Lezaud 1998 Thm 1.1 (concentration), back-referenced without a nearby cue; no internal block 1.1",
    "3.6": "external: a beta-mixing corollary cited with a DOI; no internal block 3.6",
    "6.8c": "retracted: Lemma 6.8c was removed in an earlier draft and is mentioned only historically",
    "7.3a": "retracted: Theorem 7.3a was an attempted closure, retracted; mentioned only historically",
    "4.4": "deleted: Theorem 4.4 (an earlier APX-hardness claim with a broken QAP-Hamming reduction) was removed; rnr_summary.tex documents the deletion and the resulting open E_12 inapproximability",
}


def read(path):
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def norm(s):
    return re.sub(r"[$~\\\s]", "", s).lower()


def strip_comments(text):
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("%"))


def defined_block_nums(text):
    return {norm(m.group(1)) for m in DEF_BLOCK.finditer(text)}


def phantom_block_refs(text, block_nums, companion_aware=False):
    """Distinct unresolved (kind, number) references in `text`.

    companion_aware=True enables the split-paper convention: a reference is
    satisfied when an [RNR-x] companion tag sits within COMPANION_WINDOW lines.
    """
    newlines = [i for i, ch in enumerate(text) if ch == "\n"]

    def line_of(pos):
        return bisect.bisect_right(newlines, pos)

    tag_lines = set()
    if companion_aware:
        for m in COMPANION_TAG.finditer(text):
            tag_lines.add(line_of(m.start()))

    out = {}
    for m in REF_BLOCK.finditer(text):
        num = norm(m.group(2))
        if num in block_nums or num in ALLOWLIST_BLOCK_NUM:
            continue
        if EXTERNAL_CUE.search(text[max(0, m.start() - PRECEDING):m.start()]):
            continue
        if companion_aware:
            ln = line_of(m.start())
            if any(t in tag_lines
                   for t in range(ln - COMPANION_WINDOW, ln + COMPANION_WINDOW + 1)):
                continue
        ctx = text[max(0, m.start() - 22):m.start() + 46].replace("\n", " ")
        out.setdefault((m.group(1), num), ctx)
    return sorted(out.items())


def report_doc(name, dangling):
    print(f"\n-- {name}: phantom named-block references --")
    if dangling:
        for (kind, num), ctx in dangling:
            print(f"  [PHANTOM] {kind} {num}   ...{ctx.strip()}...")
    else:
        print("  (none -- every block reference resolves)")


def main():
    print("Cross-reference integrity: prose Definition/Theorem/Lemma/Proposition/")
    print("Corollary/Remark pointers resolve to a numbered block (kind-agnostic).")
    print("Group A (monolith + companions): refs resolve against the union.")
    print("Group B (three-paper split): each paper INDEPENDENT -- own blocks or")
    print("an [RNR-x] companion tag within +/-%d lines." % COMPANION_WINDOW)
    print("=" * 74)

    verdicts = []  # (doc name, phantom count)

    # -- Group A: union resolution ------------------------------------------
    texts_a = {os.path.basename(p): strip_comments(read(p)) for p in DOCS}
    union_nums = set()
    for t in texts_a.values():
        union_nums |= defined_block_nums(t)
    print(f"  Group A defined: {len(union_nums)} distinct block numbers "
          f"(union of {len(texts_a)} docs)")
    print(f"  allowlisted (external/retracted, documented): "
          f"{sorted(ALLOWLIST_BLOCK_NUM)}")

    for name, text in texts_a.items():
        dangling = phantom_block_refs(text, union_nums)
        report_doc(name, dangling)
        verdicts.append((name, len(dangling)))

    # -- Group B: each split paper independent ------------------------------
    texts_b = {os.path.basename(p): strip_comments(read(p)) for p in SPLIT_DOCS}
    for name, text in texts_b.items():
        own_nums = defined_block_nums(text)
        dangling = phantom_block_refs(text, own_nums, companion_aware=True)
        print(f"\n  {name}: {len(own_nums)} own block numbers")
        report_doc(name, dangling)
        verdicts.append((name, len(dangling)))

    # -- Self-tests (non-vacuity) -------------------------------------------
    # (a) Union mode: a fabricated reference to an absent block must be caught.
    probe_a = next(iter(texts_a.values())) + "\n\nCites nonexistent Theorem 99.99 here.\n"
    caught_a = any(num == "99.99"
                   for (_, num), _ in phantom_block_refs(probe_a, union_nums))
    # (b) Split mode: a fabricated BARE reference must be caught ...
    split_text = texts_b[os.path.basename(SPLIT_DOCS[-1])]
    split_nums = defined_block_nums(split_text)
    probe_b = split_text + "\n\nCites nonexistent Theorem 99.99 here.\n"
    caught_b = any(num == "99.99"
                   for (_, num), _ in phantom_block_refs(probe_b, split_nums,
                                                         companion_aware=True))
    # (c) ... while the same reference WITH an adjacent companion tag is satisfied.
    probe_c = split_text + "\n\nUses companion Theorem 99.98 of [RNR-II] here.\n"
    tag_ok = not any(num == "99.98"
                     for (_, num), _ in phantom_block_refs(probe_c, split_nums,
                                                           companion_aware=True))
    print("\n-- self-tests (non-vacuity) --")
    print(f"  [{'PASS' if caught_a else 'FAIL'}] union mode: fabricated bare "
          f"'Theorem 99.99' is {'flagged' if caught_a else 'NOT flagged -- vacuous!'}")
    print(f"  [{'PASS' if caught_b else 'FAIL'}] split mode: fabricated bare "
          f"'Theorem 99.99' is {'flagged' if caught_b else 'NOT flagged -- vacuous!'}")
    print(f"  [{'PASS' if tag_ok else 'FAIL'}] split mode: fabricated "
          f"'Theorem 99.98 of [RNR-II]' is "
          f"{'accepted as companion-cited' if tag_ok else 'wrongly flagged!'}")

    # -- Per-document verdicts ----------------------------------------------
    self_ok = caught_a and caught_b and tag_ok
    print("\n" + "=" * 74)
    print("Per-document verdicts:")
    for name, n in verdicts:
        status = "PASS" if n == 0 else "FAIL"
        detail = "" if n == 0 else f" -- {n} phantom reference{'s' if n != 1 else ''}"
        print(f"  [{status}] {name}{detail}")
    print(f"  [{'PASS' if self_ok else 'FAIL'}] self-tests")

    ok = self_ok and all(n == 0 for _, n in verdicts)
    print("=" * 74)
    if ok:
        print("OVERALL -> PASS")
        return 0
    print("OVERALL -> FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
