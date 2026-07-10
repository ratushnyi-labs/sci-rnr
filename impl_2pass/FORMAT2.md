# RNR1 two-pass (train-then-freeze) archive format — lab prototype

Extension of `impl/FORMAT.md` (reference container, version 1.0.0) by a
frozen-model section.  Version **1.1.0** (MINOR bump per R-13.1: a 1.0.0
reader refuses 1.1.0 archives per R-13.2, as required — it cannot decode
them; a 1.0.0 archive is still what the fallback path emits, bit-identical
to `rnr1.pack`).

## Mode flag and fallback rule

Header `flags` bit 1 (`0x0002`) is the **two-pass mode flag**.  The
encoder (`two_pass.pack2`) always produces the single-pass reference
archive as well and emits whichever is smaller in total bytes (ties →
single-pass).  Consequently:

* a **two-pass** archive has version 1.1.0, flag bit 1 set,
  `model_hash` = the two-pass canonical description hash (below);
* a **fallback** archive is byte-identical to the reference coder's
  output (version 1.0.0, flag bit 1 clear) — two-pass mode never loses
  on archive size, by construction.

## Layout (two-pass mode)

```
+--------------------------------------+
| Header (48 bytes, as FORMAT.md)      |  version 1.1.0, flags bit 1 set
+--------------------------------------+
| Model section: uint32 LE length,     |
|                then the model blob   |
+--------------------------------------+
| Seek index (m x 32, as FORMAT.md)    |
+--------------------------------------+
| Repair stream (as FORMAT.md)         |
+--------------------------------------+
```

Everything except the inserted model section is unchanged from
FORMAT.md: same header struct, same seek-index entries (byte offsets,
bit offsets, store-mode flag, per-sub-block hashes), same per-sub-block
arithmetic-coded blobs with byte-aligned flush, same store-mode escape,
same `v = H_v(X)` verification.

`model_hash` covers the canonical two-pass description string
`rnr-type1-2pass;base=<reference canonical hash for W>;mode=train-freeze;
prune-v1;serial=varint-v1`, pinning the base Type-I family (including
the E12 table via the reference hash), the mode, the prune rule and the
serialization version.  The concrete frozen counts are carried by the
archive itself and are integrity-covered by `v`/sub-block hashes.

## Model blob

All integers are unsigned LEB128 varints.

```
uint8  W
for order o = 0 .. W:
  varint n_contexts[o]
  n_contexts entries, sorted by ascending context key:
    varint key_delta    # first entry: the key itself; subsequent
                        # entries: key - prev_key (strictly positive)
    varint n_nonzero    # 1..256
    n_nonzero pairs, ascending byte value:
      varint byte_delta # first: the byte value; subsequent: strictly
                        # positive difference
      varint count      # >= 1; per-context total (= sum of counts,
                        # not stored) < 2^16
```

A context key packs the last `o` bytes with the most recent byte in the
low 8 bits (the reference `_ctx_key` convention).  Readers MUST reject
non-ascending keys/bytes, out-of-range values, totals >= 2^16, and
trailing bytes (strict parse).

## Coding with the frozen model (PASS 2)

Frequencies are the reference mapping: `freq[b] = 1 + 32*count[b]`,
`total = 256 + 32*sum(count)`, so the probability floor (eta >= 2^-24)
and the arithmetic-coder total bound (< 2^24) of FORMAT.md hold.

At position `i` of sub-block `j`, the distribution comes from the
longest **kept** context among the last `min(W, i - jK)` bytes, falling
through pruned contexts down to order 0 (always kept).  No adaptation,
no reset: the model carries no cross-sub-block state, so Definition
10.3(a) holds trivially and random access is identical to FORMAT.md
(`read(p, k)` decodes at most `K + k` positions).  The class/payload
partition and the arithmetic coder are exactly the reference ones.

Note: coding contexts are truncated at sync points (random access needs
it) while PASS 1 training ran over the whole input without resets; the
frozen tables are a model, not a replay, so the mismatch affects rate
only for the first <= W-1 positions of each sub-block, never
correctness.

## PASS 1 training and the prune rule (prune-v1)

Training: the reference adaptive integer update rules (+1 per observed
byte at every order `0..min(W, i)`; halve counts by integer shift when
a context total reaches 2^16), applied once over the whole input with
no sync resets.  Deterministic.

Prune (encoder-side model selection): order 0 is always kept; for
orders `o = 1..W` in increasing order, a context with counts `c[b]`,
total `T` is kept iff

```
sum_b c[b]*( log2((256+32*T_a)/(1+32*a[b])) - log2((256+32*T)/(1+32*c[b])) )
      >  8 * serialized_entry_bytes
```

where `(a, T_a)` is the nearest kept ancestor (longest kept suffix at a
lower order; exists because order 0 is kept) and
`serialized_entry_bytes` is the entry's actual varint cost with the key
costed as a full varint (an upper bound on its delta encoding).  This
is the estimated inline-coding saving versus backing off, against the
model-header cost.

Documented approximations (encoder-side only; decoding is exact):

* the saving charges all `T` training occurrences to the context, but
  occurrences also matched by a kept longer context are coded by that
  longer context in PASS 2 — usage of shorter contexts is
  overestimated; the whole-archive fallback comparison bounds any loss;
* the decision uses float `log2`; same-platform encoder determinism
  holds (fixed iteration order), but cross-platform bit-reproducibility
  of the *encoder* would need a fixed-point log (not implemented in
  this prototype).  Decoders are unaffected: the model ships in the
  archive.
