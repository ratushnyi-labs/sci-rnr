# RNR1 archive format — reference Type-I coder

Container format for the reference RNR Type-I coder (`rnr1.py`).
Version 1.0.0 (semantic versioning per spec R-13.1: MAJOR breaks
compatibility, MINOR is backward-compatible, PATCH is bidirectional).

All multi-byte integers are little-endian. All offsets are absolute
within the section stated.

## Layout

```
+--------------------+
| Header  (48 bytes) |
+--------------------+
| Seek index (m x 28)|
+--------------------+
| Repair stream      |
+--------------------+
```

## Header (48 bytes)

| offset | size | field        | meaning                                             |
|-------:|-----:|--------------|-----------------------------------------------------|
| 0      | 4    | magic        | ASCII `RNR1`                                        |
| 4      | 3    | version      | major, minor, patch (1, 0, 0)                       |
| 7      | 1    | reserved     | 0                                                   |
| 8      | 2    | flags        | bit 0: per-sync-point hashes present (R-11.1.1)     |
| 10     | 2    | W            | predictor context order (Definition 10.2)           |
| 12     | 4    | K            | sync-point spacing in source bytes (Definition 10.3)|
| 16     | 8    | N            | uncompressed source length in bytes                 |
| 24     | 8    | model_hash   | truncated SHA-256 of the canonical model description (R-13.3) |
| 32     | 8    | v            | verification value H_v(X): truncated SHA-256 of the source (Definition 3.1) |
| 40     | 8    | m            | sub-block count, `ceil(N/K)` (uint64 length prefix, R-8.2.1) |

A reader MUST refuse an archive whose `model_hash` does not match its
own model description hash (R-13.3 hard failure), and MUST refuse a
different MAJOR or a newer MINOR version (R-13.2).

The model description covered by `model_hash` is the canonical string
`rnr-type1-code12a;predictor=ngram;W=<W>;scale=32;cap=65536;cand=15;
floor=1;classes=c0,c1,c2,cand,esc;e12=<hash of E12 table>`; it pins the
predictor family, the integer count parameters, the class alphabet, and
the exact E12 table.

## Seek index (m entries of 32 bytes, R-8.1.1 / R-8.2.1)

| offset | size | field             | meaning                                        |
|-------:|-----:|-------------------|------------------------------------------------|
| 0      | 8    | byte_offset       | source position of sync point j (uniform: jK)  |
| 8      | 8    | repair_bit_offset | bit offset of sub-block j in the repair stream |
| 16     | 4    | snapshot_length   | always 0 (W-bounded predictor, R-8.1.1)        |
| 20     | 4    | mode              | 0 = arithmetic-coded, 1 = raw store mode       |
| 24     | 8    | subblock_hash     | truncated SHA-256 of source bytes [jK, jK+len) (R-11.1.1) |

Entries are sorted by `byte_offset` ascending; lookup uses binary
search for the largest `byte_offset <= p` (R-8.2.2). Spacing is the
uniform variant of Definition 10.3 (`byte_offset = jK` is redundant
but stored explicitly to keep the record layout of R-8.1.1). Because
every sub-block is flushed to a byte boundary, `repair_bit_offset` is
always a multiple of 8. The index is stored uncompressed (compression
is a MAY, R-8.2.3).

## Repair stream

Concatenation of the m per-sub-block arithmetic-coder blobs. Sub-block
j covers source bytes `[jK, min((j+1)K, N))`. Each sub-block is a
self-contained restart point (Definition 10.3):

* (a) predictor reset — context history and all adaptive count tables
  are reset to the fixed initial state Q_init (empty tables, uniform
  distribution), independent of bytes `[0, jK)`;
* (b) entropy-coder flush — the arithmetic coder terminates (2 + carry
  bits) and pads to a byte boundary; the decoder starts from a fresh
  coder state at the recorded bit offset.

The sub-block length is implied by `N` and `K`; no in-stream length
field or end-of-block symbol is needed (the decoder decodes exactly
`min(K, N - jK)` positions).

Fail-safe store mode: if the arithmetic-coded blob of a sub-block is
not smaller than the raw sub-block, the encoder stores the sub-block
raw and sets `mode = 1` in its seek-index entry. This bounds the
archive at 8 bits/byte plus header/index overhead on incompressible or
adversarial input (the store-mode upper bound of Theorem 7.22 /
hypothesis H16), while keeping every sub-block independently decodable
(a raw sub-block is trivially so).

## Per-position coding (Type-I, Definition 3.1, Sections 4.2–4.3)

For each source byte x with predictor distribution `freq[0..255]`
(integers, `freq[b] >= 1`, `total < 2^24`):

1. `b* = argmax(freq)` (ties broken toward the smallest byte value);
   the predicted codeword is `c_hat = E12(b*)` where E12 is the
   systematic shortened-Hamming-like (12,8) map of Section 4.2 (data
   bits at codeword positions 3,5,6,7,9,10,11,12; even parity at
   1,2,4,8). The error mask is `e = E12(x) XOR c_hat`.
2. The class of x is determined with precedence C0 > C1 > C2 > CAND > ESC:
   * C0 — `e = 0` (x = b*);
   * C1 — Hamming weight(e) = 1;
   * C2 — Hamming weight(e) = 2;
   * CAND — x is one of the 15 next-most-probable bytes after b*
     (rank order: freq descending, byte value ascending on ties);
   * ESC — everything else.
3. The class symbol is arithmetic-coded with the exact partition sums
   of `freq` as class frequencies; the within-class payload (rank for
   CAND, byte in ascending order for the other classes) is coded with
   the conditional frequencies `freq` restricted to the class. By the
   chain rule the pair costs exactly `-log2(freq[x]/total)` per byte
   up to arithmetic-coder quantization, realizing the ideal
   class-and-payload chain of Theorem 4.2(A) with the model law in
   place of the true law.

Note on C1/C2: E12-A is linear, so the XOR of two valid codewords is a
codeword and has Hamming weight 0 or >= 3 (minimum distance 3). With a
byte-aligned argmax predictor both c and c_hat are valid codewords, so
classes C1 and C2 are structurally empty and carry zero probability
mass (they cost nothing). They are retained in the class alphabet for
format compatibility with codeword-level predictors, whose raw 12-bit
predictions can produce weight-1/2 masks (Section 4.3 taxonomy).

## Predictor (order-W byte context model)

Integer-deterministic order-W context model:

* `dist()` — from the longest available context `min(W, i - jK)` down
  to order 0, the first context with at least one observation supplies
  `freq[b] = 1 + 32 * count[ctx][b]`, `total = 256 + 32 * count_total`;
  if none, the distribution is uniform (`freq = 1`). The `+1` floor
  guarantees a per-byte probability floor `eta >= 2^-24` (Section 10.7
  bounded-cost assumption; worst-case 24 bits/byte).
* `update(x)` — after coding each byte, `count[ctx][x] += 1` for every
  order 0..min(W, available); when a context total reaches 2^16 all its
  counts are halved (integer shift). This is a deterministic integer
  state update in the sense of Definition 10.1 (compatible online
  optimizer), covered by Theorem 10.2.

W-boundedness: the context *argument* of the prediction is the last W
bytes only (Definition 10.2 semantics). The adaptive count state is a
deterministic function of the bytes coded since the last sync point and
is reset at every sync point, so the property that Definition 10.2 is
used for in Theorem 10.3 — decodability of any sub-block from Q_init
with no bytes before the sync point — holds verbatim. A strictly
static-table reading of Definition 10.2 (frozen counts shipped as the
model description) is a trivial specialization: freeze the tables and
skip `update`.

## Arithmetic coder

Binary (bit-output) integer arithmetic coder, 32-bit registers,
Witten–Neal–Cleary renormalization with carry handling via pending-bit
counting. Cumulative frequency totals are bounded by `2^24 < 2^30`, so
`range // total >= 1` for every coded symbol and all intermediates fit
in 64 bits (R-3.2/R-5.2 style overflow safety; the reference uses
Python arbitrary-precision integers, with assertions enforcing the
bounds). Termination emits 2 bits plus pending carries, then pads to a
byte boundary — the O(1)-bits-per-sync-point flush of Definition 10.3.

## Random access (Theorem 10.3)

`read(p, k)`: binary-search the seek index for the latest sync point
`P <= p`, decode forward from the sub-block start, discard the warm-up
prefix `[P, p)` (at most K-1 bytes), output `[p, p + k')` with
`k' = min(k, N - p)` (short read at end of source). Reads spanning
sub-block boundaries decode each covering sub-block independently.
Total decoded positions are at most `K + k`. With `verify=True` the
covering sub-blocks are decoded in full so their R-11.1.1 hashes can be
checked (decode window at most `2K + k`); the default leaves partial
sub-blocks unverified, relying on archive-level trust (R-11.3
documented mode).

## Verification

* Archive level: `v = H_v(X)` (truncated SHA-256, Definition 3.1) is
  checked after full sequential decode.
* Sub-block level: 8-byte truncated SHA-256 per sync point (R-11.1.1),
  checked whenever a sub-block is decoded in full (sequential decode
  always; random access with `verify=True`).

## CLI

```
rnr1.py pack   <input> <output.rnr> [--W 3] [--K 65536]
rnr1.py unpack <archive.rnr> <output> [--no-verify]
rnr1.py read   <archive.rnr> --pos P --len k [--out FILE] [--verify]
rnr1.py info   <archive.rnr>
```
