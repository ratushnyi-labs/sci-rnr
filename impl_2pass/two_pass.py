#!/usr/bin/env python3
"""two_pass.py -- Two-pass (train-then-freeze) coding mode for the RNR
Type-I coder.  Lab prototype.

Mode summary
------------
PASS 1 trains the order-0..W count tables over the WHOLE input with the
same integer update rules as the reference adaptive path (impl/rnr1.py
NGramPredictor.update: +1 per observation at every order 0..min(W,i);
halve counts by integer shift when a context total reaches 2^16), but
with NO sync-point resets -- the model is global.  The trained tables
are PRUNED to the profitable contexts (rule below), SERIALIZED into a
model section of the archive (format below and in FORMAT2.md), and
PASS 2 codes every sub-block with the FROZEN tables: no adaptation, no
cold start.  Sync semantics and random access are preserved -- a frozen
model carries no cross-block state, so a sub-block decodes from the
shared model plus its own bytes only.  The class/payload partition, the
arithmetic coder, the seek index, store mode and the verification
values are exactly those of the reference coder.

Fallback guarantee: pack2 also produces the single-pass reference
archive and emits whichever is smaller (ties -> single-pass), so the
two-pass mode never loses on total archive size.  A single-pass
fallback archive is byte-identical to rnr1.pack output (version 1.0.0,
two-pass flag clear); a two-pass archive carries format version 1.1.0
and header flag bit 1 (the mode flag).

Prune rule (prune-v1)
---------------------
Order-0 is always kept.  For orders o = 1..W in increasing order, and
each trained context ctx (counts c[b], total T = sum c[b]):

  bits_self = sum_b c[b] * (log2(256 + 32*T)   - log2(1 + 32*c[b]))
  bits_anc  = sum_b c[b] * (log2(256 + 32*T_a) - log2(1 + 32*a[b]))

where (a, T_a) is the NEAREST KEPT ANCESTOR of ctx: the longest kept
context among the suffixes of ctx at orders o-1, o-2, ..., 0 (kept-ness
of shorter orders is already decided because orders are processed in
increasing order; order 0 is always kept, so an ancestor always
exists).  ctx is kept iff

  bits_anc - bits_self  >  8 * serialized_entry_bytes(ctx)

with serialized_entry_bytes computed from the actual varint encoding
(the context key is costed as a full varint, an upper bound on its
delta encoding, so pruning errs slightly toward dropping).  This is the
estimated inline-coding saving of the context against what the coder
would fall back to without it, versus its cost in the model header.

Known approximation (documented, encoder-side only): bits_anc/bits_self
charge ALL T training occurrences of ctx to ctx, but occurrences also
matched by a kept LONGER context are actually coded by that longer
context in pass 2.  The rule therefore overestimates usage of shorter
contexts; the fallback comparison bounds any resulting loss at the
whole-archive level.

Backoff in pass 2 mirrors the reference predictor with "has a table"
replaced by "kept in the model": at each position the distribution
comes from the longest kept context among the last min(W, i - jK)
bytes, falling through pruned contexts down to order 0.  Encoder and
decoder share the serialized model, so the backoff is reproducible.

Archive layout (two-pass mode; single-pass fallback == FORMAT.md)
-----------------------------------------------------------------
  Header (48 bytes, HEADER_FMT of FORMAT.md; version = 1.1.0,
          flags bit 1 = two-pass mode, model_hash = two-pass canonical
          description hash)
  Model section:  uint32 LE byte length, then the model blob
  Seek index:     m x 32 bytes, unchanged
  Repair stream:  unchanged (frozen-model arithmetic-coded sub-blocks,
                  per-sub-block store-mode escape preserved)

Model blob (all integers LEB128 varints):
  uint8 W
  for order o = 0..W:
    varint n_contexts[o]
    n_contexts entries, sorted by ascending context key:
      varint key_delta   (first entry: the key; later: key - prev_key,
                          strictly positive)
      varint n_nonzero   (1..256)
      n_nonzero (byte, count) pairs, ascending byte value:
        varint byte_delta (first: byte value; later: strictly positive)
        varint count      (1 <= count, context total < 2^16)

Totals are not stored; total = sum of counts (invariant of the update
rule).  Frequencies in pass 2 are freq[b] = 1 + 32*count[b] out of
total_f = 256 + 32*total, exactly the reference mapping, so the coder
floor (eta >= 2^-24) and the total < 2^24 arithmetic-coder bound hold.

Honest caveats
--------------
* The prune decision uses float log2 (encoder-side model selection
  only; the archive carries the resulting model explicitly, so decoding
  and the R-13.3-style model hash are unaffected).  Same-platform
  encoder determinism holds (fixed iteration order, pure functions);
  cross-platform bit-reproducibility of the ENCODER would additionally
  need a fixed-point log, which this prototype does not implement.
* Training context crosses sub-block boundaries (global pass), while
  pass-2 coding contexts are truncated to the current sub-block for the
  first <= W-1 positions after each sync point (random access needs
  this).  The frozen tables are a model, not a replay, so this small
  train/code mismatch affects rate only, never correctness.
* Pure-Python throughput is reference-grade (slow); use for gates and
  small measurements, not campaigns (see measure_smoke.sh).

CLI
---
  two_pass.py pack   <input> <output.rnr> [--W 3] [--K 65536]
                     [--force {single,two}]
  two_pass.py unpack <archive.rnr> <output> [--no-verify]
  two_pass.py read   <archive.rnr> --pos P --len k [--out F] [--verify]
  two_pass.py info   <archive.rnr>
"""

import argparse
import hashlib
import math
import os
import struct
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "impl"))
import rnr1  # noqa: E402  (normative reference coder)

# ---------------------------------------------------------------------------
# Format constants (see FORMAT2.md; base constants come from the reference)
# ---------------------------------------------------------------------------

VERSION2 = (1, 1, 0)          # MINOR bump: adds the model section (R-13.1)
FLAG_TWO_PASS = 0x0002        # header flags bit 1: two-pass mode flag
PRUNE_RULE = "prune-v1"

COUNT_SCALE = rnr1.COUNT_SCALE
COUNT_CAP = rnr1.COUNT_CAP
DEFAULT_W = rnr1.DEFAULT_W
DEFAULT_K = rnr1.DEFAULT_K


def model_description_hash2(W):
    """Truncated SHA-256 of the two-pass canonical model description.

    Pins the base Type-I model family (via the reference canonical
    string hash, which includes the E12 table), the train-then-freeze
    mode, the prune rule and the model serialization version.
    """
    desc = (
        "rnr-type1-2pass;base=%s;mode=train-freeze;%s;serial=varint-v1"
        % (rnr1.model_description_hash(W).hex(), PRUNE_RULE)
    ).encode("ascii")
    return hashlib.sha256(desc).digest()[:8]


# ---------------------------------------------------------------------------
# Varints (LEB128, unsigned)
# ---------------------------------------------------------------------------

def _varint(v):
    if v < 0:
        raise ValueError("varint of negative value")
    out = bytearray()
    while True:
        b = v & 0x7F
        v >>= 7
        if v:
            out.append(b | 0x80)
        else:
            out.append(b)
            return bytes(out)


def _varint_len(v):
    n = 1
    while v >= 0x80:
        v >>= 7
        n += 1
    return n


class _Reader:
    __slots__ = ("buf", "pos")

    def __init__(self, buf):
        self.buf = buf
        self.pos = 0

    def varint(self):
        shift = 0
        val = 0
        while True:
            if self.pos >= len(self.buf):
                raise ValueError("truncated varint in model section")
            b = self.buf[self.pos]
            self.pos += 1
            val |= (b & 0x7F) << shift
            if not (b & 0x80):
                return val
            shift += 7
            if shift > 63:
                raise ValueError("overlong varint in model section")

    def eof(self):
        return self.pos == len(self.buf)


# ---------------------------------------------------------------------------
# PASS 1: global training (reference integer update rules, no resets)
# ---------------------------------------------------------------------------

def train_counts(data, W):
    """Train order-0..W count tables over the whole input.

    Same integer update rules as rnr1.NGramPredictor.update (+1 at every
    order 0..min(W, i) per position; halve-and-drop-zeros when a context
    total reaches COUNT_CAP), with no sync resets.  Sparse tables:
    tables[o] maps context key -> [dict byte -> count, total], with the
    invariant total == sum(counts) (kept by the update rule).
    Deterministic: pure integer state, fixed left-to-right order.
    """
    tables = [dict() for _ in range(W + 1)]
    masks = [(1 << (8 * o)) - 1 for o in range(W + 1)]
    wmask = masks[W]
    hk = 0        # rolling key of the last W bytes (newest in low bits)
    nseen = 0
    cap = COUNT_CAP
    for x in data:
        avail = nseen if nseen < W else W
        for o in range(avail + 1):
            tab = tables[o]
            key = hk & masks[o]
            ent = tab.get(key)
            if ent is None:
                ent = [{}, 0]
                tab[key] = ent
            c = ent[0]
            c[x] = c.get(x, 0) + 1
            tot = ent[1] + 1
            if tot >= cap:
                nc = {}
                for b, v in c.items():
                    v >>= 1
                    if v:
                        nc[b] = v
                ent[0] = nc
                ent[1] = sum(nc.values())
            else:
                ent[1] = tot
        hk = ((hk << 8) | x) & wmask
        nseen += 1
    return tables


# ---------------------------------------------------------------------------
# Pruning (prune-v1; rule documented in the module docstring / FORMAT2.md)
# ---------------------------------------------------------------------------

def _entry_cost_bits(key, cnts):
    """Serialized size of one model entry in bits (full-key upper bound)."""
    nbytes = _varint_len(key) + _varint_len(len(cnts))
    prev = 0
    first = True
    for b in sorted(cnts):
        nbytes += _varint_len(b if first else b - prev) + _varint_len(cnts[b])
        first = False
        prev = b
    return 8 * nbytes


def prune_tables(tables, W):
    """Keep a context iff its estimated coding saving over its nearest
    kept ancestor exceeds its serialized cost.  Returns (kept, stats).

    Orders are processed in increasing order so ancestor kept-ness is
    final when a context is judged; order 0 is always kept.  Float log2
    is used for the estimate (encoder-side model selection only; see
    module docstring caveat).  Iteration over sorted keys keeps the
    procedure deterministic.
    """
    kept = [dict() for _ in range(W + 1)]
    kept[0] = dict(tables[0])
    lg_cache = {}
    log2 = math.log2

    def lg(v):
        r = lg_cache.get(v)
        if r is None:
            r = log2(v)
            lg_cache[v] = r
        return r

    for o in range(1, W + 1):
        ko = kept[o]
        amasks = [(1 << (8 * ao)) - 1 for ao in range(o)]
        for key in sorted(tables[o]):
            cnts, T = tables[o][key]
            if T <= 0:
                continue
            anc = None
            for ao in range(o - 1, -1, -1):
                aent = kept[ao].get(key & amasks[ao])
                if aent is not None and aent[1] > 0:
                    anc = aent
                    break
            acnts, aT = anc  # order 0 is kept and non-empty, so anc exists
            lt_self = lg(256 + COUNT_SCALE * T)
            lt_anc = lg(256 + COUNT_SCALE * aT)
            saving = 0.0
            for b, cb in cnts.items():
                fs = 1 + COUNT_SCALE * cb
                fa = 1 + COUNT_SCALE * acnts.get(b, 0)
                saving += cb * ((lt_anc - lg(fa)) - (lt_self - lg(fs)))
            if saving > _entry_cost_bits(key, cnts):
                ko[key] = [cnts, T]
    stats = {
        "contexts_trained": [len(t) for t in tables],
        "contexts_kept": [len(k) for k in kept],
    }
    return kept, stats


# ---------------------------------------------------------------------------
# Model serialization (varint, context-sorted; layout in FORMAT2.md)
# ---------------------------------------------------------------------------

def serialize_model(kept, W):
    out = bytearray()
    out.append(W)
    for o in range(W + 1):
        keys = sorted(kept[o])
        out += _varint(len(keys))
        prev = -1
        for key in keys:
            out += _varint(key if prev < 0 else key - prev)
            prev = key
            cnts = kept[o][key][0]
            bs = sorted(cnts)
            out += _varint(len(bs))
            pb = -1
            for b in bs:
                out += _varint(b if pb < 0 else b - pb)
                pb = b
                out += _varint(cnts[b])
    return bytes(out)


def deserialize_model(blob, W_expected):
    """Strict parse of the model blob; validates ordering and bounds."""
    if len(blob) < 1:
        raise ValueError("empty model section")
    W = blob[0]
    if W != W_expected:
        raise ValueError("model W mismatch (header %d, model %d)"
                         % (W_expected, W))
    r = _Reader(blob)
    r.pos = 1
    kept = []
    for o in range(W + 1):
        maxkey = (1 << (8 * o)) - 1
        nctx = r.varint()
        tab = {}
        prev = -1
        for _ in range(nctx):
            d = r.varint()
            if prev >= 0 and d == 0:
                raise ValueError("model context keys not strictly ascending")
            key = d if prev < 0 else prev + d
            if key > maxkey:
                raise ValueError("context key out of range at order %d" % o)
            prev = key
            nnz = r.varint()
            if not 1 <= nnz <= 256:
                raise ValueError("bad nonzero-count in model entry")
            cnts = {}
            pb = -1
            total = 0
            for _ in range(nnz):
                db = r.varint()
                if pb >= 0 and db == 0:
                    raise ValueError("model byte values not strictly ascending")
                b = db if pb < 0 else pb + db
                if b > 255:
                    raise ValueError("byte value out of range in model entry")
                pb = b
                c = r.varint()
                if c < 1 or c >= COUNT_CAP:
                    raise ValueError("bad count in model entry")
                cnts[b] = c
                total += c
            if total >= COUNT_CAP:
                raise ValueError("context total exceeds cap")
            tab[key] = [cnts, total]
        kept.append(tab)
    if not r.eof():
        raise ValueError("trailing bytes in model section")
    return kept


# ---------------------------------------------------------------------------
# Frozen model + PASS 2 coding (no adaptation, no cold start)
# ---------------------------------------------------------------------------

class _Ctx:
    """Cached coding tables of one frozen context (pure functions of it).

    Storage is kept small (int32/uint8, candidate-class cache only) so
    the per-context cache stays modest even with tens of thousands of
    kept contexts; the rarely-hit non-candidate payload classes are
    recomputed on the fly.
    """

    __slots__ = ("freq", "total", "bstar", "cand", "classid", "clsb",
                 "ccum", "_cand_tab")

    def __init__(self, cnts, total):
        freq = np.ones(256, dtype=np.int64)
        for b, c in cnts.items():
            freq[b] += COUNT_SCALE * c
        total_f = 256 + COUNT_SCALE * total
        assert total_f < rnr1.PROB_TOTAL_BOUND
        bstar, cand, classid, class_freq = rnr1._partition(freq)
        self.freq = freq.astype(np.int32)      # values < 2^22
        self.total = int(total_f)
        self.bstar = int(bstar)
        self.cand = np.ascontiguousarray(cand)
        self.classid = classid.astype(np.uint8)
        self.clsb = self.classid.tobytes()
        cc = np.cumsum(class_freq)
        self.ccum = [0] + [int(v) for v in cc]
        self._cand_tab = None

    def cand_tab(self):
        """Candidate-class payload tables (rank order), cached."""
        if self._cand_tab is None:
            m = self.cand[self.classid[self.cand] == rnr1.CLS_CAND]
            mc = np.cumsum(self.freq[m])
            mcum = [0] + [int(v) for v in mc]
            pos = {int(b): i for i, b in enumerate(m)}
            self._cand_tab = ([int(b) for b in m], mcum, pos)
        return self._cand_tab

    def other_tab(self, cls):
        """Payload tables of a non-candidate class (ascending byte
        order, so members are sorted); computed fresh, not cached."""
        m = np.nonzero(self.classid == cls)[0]
        mc = np.cumsum(self.freq[m])
        return m, mc


class FrozenModel:
    """Frozen order-0..W tables with longest-kept-context backoff.

    ctx(hist_key, avail) returns the coding context for a position whose
    last min(avail, W) bytes are the low bytes of hist_key: the longest
    order o <= min(avail, W) whose context key is kept, falling through
    pruned contexts down to order 0 (always kept for non-empty input).
    No state is mutated by coding, so sub-blocks share the model freely.
    """

    __slots__ = ("W", "kept", "mask", "_masks", "_cache", "_uniform")

    def __init__(self, kept, W):
        self.W = W
        self.kept = kept
        self._masks = [(1 << (8 * o)) - 1 for o in range(W + 1)]
        self.mask = self._masks[W]
        self._cache = {}
        self._uniform = None

    def ctx(self, hist_key, avail):
        top = avail if avail < self.W else self.W
        for o in range(top, -1, -1):
            key = hist_key & self._masks[o]
            ck = (o, key)
            got = self._cache.get(ck)
            if got is not None:
                return got
            ent = self.kept[o].get(key)
            if ent is not None and ent[1] > 0:
                c = _Ctx(ent[0], ent[1])
                self._cache[ck] = c
                return c
        # Only reachable for an empty model (never emitted by pack2);
        # uniform distribution keeps the reader total.
        if self._uniform is None:
            self._uniform = _Ctx({}, 0)
        return self._uniform


def encode_subblock_frozen(block, model):
    """PASS 2 encoding of one sub-block with the frozen model."""
    enc = rnr1.ArithmeticEncoder()
    W = model.W
    mask = model.mask
    ctxf = model.ctx
    hk = 0
    avail = 0
    for x in block:
        c = ctxf(hk, avail)
        cls = c.clsb[x]
        ccum = c.ccum
        enc.encode(ccum[cls], ccum[cls + 1], c.total)
        if cls == rnr1.CLS_CAND:
            mlist, mcum, pos = c.cand_tab()
            i = pos[x]
            enc.encode(mcum[i], mcum[i + 1], mcum[-1])
        elif cls != rnr1.CLS_0:
            m, mc = c.other_tab(cls)
            i = int(np.searchsorted(m, x))
            lo = 0 if i == 0 else int(mc[i - 1])
            enc.encode(lo, int(mc[i]), int(mc[-1]))
        hk = ((hk << 8) | x) & mask
        if avail < W:
            avail += 1
    return enc.finish()


def decode_subblock_frozen(blob, nbytes, model):
    """Decode the first nbytes positions of a frozen-coded sub-block."""
    from bisect import bisect_right
    dec = rnr1.ArithmeticDecoder(blob)
    out = bytearray(nbytes)
    W = model.W
    mask = model.mask
    ctxf = model.ctx
    hk = 0
    avail = 0
    for i in range(nbytes):
        c = ctxf(hk, avail)
        ccum = c.ccum
        t = dec.target(c.total)
        cls = bisect_right(ccum, t) - 1
        dec.consume(ccum[cls], ccum[cls + 1], c.total)
        if cls == rnr1.CLS_0:
            x = c.bstar
        elif cls == rnr1.CLS_CAND:
            mlist, mcum, pos = c.cand_tab()
            t2 = dec.target(mcum[-1])
            k = bisect_right(mcum, t2) - 1
            dec.consume(mcum[k], mcum[k + 1], mcum[-1])
            x = mlist[k]
        else:
            m, mc = c.other_tab(cls)
            t2 = dec.target(int(mc[-1]))
            k = int(np.searchsorted(mc, t2, side="right"))
            lo = 0 if k == 0 else int(mc[k - 1])
            dec.consume(lo, int(mc[k]), int(mc[-1]))
            x = int(m[k])
        out[i] = x
        hk = ((hk << 8) | x) & mask
        if avail < W:
            avail += 1
    return bytes(out)


# ---------------------------------------------------------------------------
# Container: pack2 / unpack2 / random access
# ---------------------------------------------------------------------------

def build_frozen_model(data, W):
    """PASS 1 + prune + serialize.  Returns (model, blob, prune_stats).

    Depends on (data, W) only -- callers may reuse the result across
    different sync spacings K.
    """
    tables = train_counts(data, W)
    kept, stats = prune_tables(tables, W)
    blob = serialize_model(kept, W)
    return FrozenModel(kept, W), blob, stats


def _pack_two_pass(data, W, K, model, model_blob):
    """Assemble the two-pass archive (header + model + index + stream)."""
    n = len(data)
    m = (n + K - 1) // K
    blobs = []
    entries = []
    repair_bit_offset = 0
    for j in range(m):
        block = data[j * K: (j + 1) * K]
        blob = encode_subblock_frozen(block, model)
        mode = rnr1.MODE_CODED
        if len(blob) >= len(block):
            blob = block          # store-mode escape, unchanged semantics
            mode = rnr1.MODE_RAW
        entries.append(struct.pack(
            rnr1.INDEX_ENTRY_FMT,
            j * K, repair_bit_offset, 0, mode, rnr1.h_v(block)))
        blobs.append(blob)
        repair_bit_offset += 8 * len(blob)
    header = struct.pack(
        rnr1.HEADER_FMT,
        rnr1.MAGIC,
        bytes(VERSION2),
        0,
        rnr1.FLAG_SUBBLOCK_HASHES | FLAG_TWO_PASS,
        W,
        K,
        n,
        model_description_hash2(W),
        rnr1.h_v(data),
        m,
    )
    return (header + struct.pack("<I", len(model_blob)) + model_blob
            + b"".join(entries) + b"".join(blobs))


def pack2(data, W=DEFAULT_W, K=DEFAULT_K, single_raw=None, force=None,
          frozen=None, return_info=False):
    """Encode data; emit two-pass archive iff strictly smaller (fallback
    rule: two-pass must never lose on total archive size).

    single_raw : optional precomputed rnr1.pack(data, W, K) bytes (e.g.
                 from the byte-identity-gated fast engine) used for the
                 size comparison and as the fallback archive.  Caller is
                 responsible for it matching the reference output.
    force      : None (size comparison, default) | 'single' | 'two'
                 (testing/measurement knob; 'two' skips the single-pass
                 encode entirely).
    frozen     : optional (model, blob) from build_frozen_model(data, W)
                 to amortize PASS 1 across several K.
    """
    if W < 1 or W > 8:
        raise ValueError("W must be in [1, 8]")
    if K < 256:
        raise ValueError("K must be >= 256")
    n = len(data)
    info = {"n": n, "W": W, "K": K, "single_bytes": None,
            "two_pass_bytes": None, "model_bytes": None, "mode": None}

    two_raw = None
    if n > 0 and force != "single":
        if frozen is None:
            model, blob, stats = build_frozen_model(data, W)
        else:
            model, blob = frozen
            stats = {}
        two_raw = _pack_two_pass(data, W, K, model, blob)
        info.update(stats)
        info["model_bytes"] = len(blob)
        info["model_bpb_share"] = 8.0 * len(blob) / n
        info["two_pass_bytes"] = len(two_raw)

    if force == "two":
        if two_raw is None:
            raise ValueError("cannot force two-pass mode on empty input")
        chosen, mode = two_raw, "two-pass"
    else:
        if single_raw is None:
            single_raw = rnr1.pack(data, W=W, K=K)
        info["single_bytes"] = len(single_raw)
        if force == "single" or two_raw is None \
                or len(two_raw) >= len(single_raw):
            chosen, mode = single_raw, "single-pass"
        else:
            chosen, mode = two_raw, "two-pass"
    info["mode"] = mode
    info["archive_bytes"] = len(chosen)
    if return_info:
        return chosen, info
    return chosen


class TwoPassArchive(rnr1.Archive):
    """Parsed two-pass archive.  Inherits the seek/read machinery of the
    reference Archive (identical index and sub-block semantics); only
    parsing and per-sub-block decoding differ (frozen model, no reset
    needed -- the model is stateless across sub-blocks)."""

    def __init__(self, raw):  # noqa: super().__init__ intentionally not called
        if len(raw) < rnr1.HEADER_SIZE + 4:
            raise ValueError("truncated archive: no header/model section")
        (magic, version, _reserved, self.flags, self.W, self.K, self.n,
         self.model_hash, self.v, self.m) = struct.unpack_from(
            rnr1.HEADER_FMT, raw, 0)
        if magic != rnr1.MAGIC:
            raise ValueError("bad magic: not an RNR1 archive")
        if version[0] != VERSION2[0] or version[1] != VERSION2[1]:
            raise ValueError("unsupported two-pass archive version %s"
                             % (tuple(version),))
        if not (self.flags & FLAG_TWO_PASS):
            raise ValueError("two-pass flag not set")
        if self.model_hash != model_description_hash2(self.W):
            raise ValueError("model description hash mismatch")
        if self.m != (self.n + self.K - 1) // self.K:
            raise ValueError("inconsistent sub-block count")
        off = rnr1.HEADER_SIZE
        (mlen,) = struct.unpack_from("<I", raw, off)
        off += 4
        if off + mlen > len(raw):
            raise ValueError("truncated model section")
        self.model_bytes = mlen
        self.model = FrozenModel(deserialize_model(raw[off:off + mlen],
                                                   self.W), self.W)
        off += mlen
        self.entries = []
        for _ in range(self.m):
            self.entries.append(
                struct.unpack_from(rnr1.INDEX_ENTRY_FMT, raw, off))
            off += rnr1.INDEX_ENTRY_SIZE
        self.repair = raw[off:]
        for j in range(1, self.m):
            if self.entries[j][0] <= self.entries[j - 1][0]:
                raise ValueError("seek index not sorted")

    def decode_block(self, j, upto=None, verify=True):
        blk_len = self._block_len(j)
        take = blk_len if upto is None else min(upto, blk_len)
        if self.entries[j][3] == rnr1.MODE_RAW:
            out = self._block_blob(j)[:take]
        else:
            out = decode_subblock_frozen(self._block_blob(j), take, self.model)
        if verify and take == blk_len and (self.flags & rnr1.FLAG_SUBBLOCK_HASHES):
            if rnr1.h_v(out) != self.entries[j][4]:
                raise ValueError("sub-block %d hash mismatch" % j)
        return out


def open_archive(raw):
    """Parse either archive flavor (mode flag in the header decides)."""
    if len(raw) < rnr1.HEADER_SIZE:
        raise ValueError("truncated archive: no header")
    magic, _version, _res, flags = struct.unpack_from("<4s3sBH", raw, 0)
    if magic != rnr1.MAGIC:
        raise ValueError("bad magic: not an RNR1 archive")
    if flags & FLAG_TWO_PASS:
        return TwoPassArchive(raw)
    return rnr1.Archive(raw)


def unpack2(raw, verify=True):
    """Full sequential decode of either archive flavor."""
    return open_archive(raw).unpack(verify=verify)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cmd_pack(args):
    with open(args.input, "rb") as f:
        data = f.read()
    raw, info = pack2(data, W=args.W, K=args.K, force=args.force,
                      return_info=True)
    with open(args.output, "wb") as f:
        f.write(raw)
    extra = ""
    if info["model_bytes"] is not None:
        extra = ", model %d bytes (%.4f bpb)" % (
            info["model_bytes"], 8.0 * info["model_bytes"] / max(1, len(data)))
    print("packed %d bytes -> %d bytes (%.4f bpb), mode=%s%s"
          % (len(data), len(raw), 8.0 * len(raw) / max(1, len(data)),
             info["mode"], extra))


def _cmd_unpack(args):
    with open(args.input, "rb") as f:
        raw = f.read()
    data = unpack2(raw, verify=not args.no_verify)
    with open(args.output, "wb") as f:
        f.write(data)
    print("unpacked %d bytes -> %d bytes (verified: %s)"
          % (len(raw), len(data), "no" if args.no_verify else "yes"))


def _cmd_read(args):
    with open(args.archive, "rb") as f:
        arc = open_archive(f.read())
    data, stats = arc.read(args.pos, args.len, verify=args.verify)
    if args.out:
        with open(args.out, "wb") as f:
            f.write(data)
    else:
        sys.stdout.buffer.write(data)
        sys.stdout.buffer.flush()
    print("\nread %d byte(s) at %d: warmup=%d decoded=%d"
          % (len(data), args.pos, stats["warmup_bytes"],
             stats["decoded_bytes"]), file=sys.stderr)


def _cmd_info(args):
    with open(args.archive, "rb") as f:
        raw = f.read()
    arc = open_archive(raw)
    two = isinstance(arc, TwoPassArchive)
    print("RNR1 archive: mode=%s flags=0x%04x"
          % ("two-pass (v1.1.0)" if two else "single-pass (v1.0.0)",
             arc.flags))
    print("  source length : %d bytes" % arc.n)
    print("  W (order)     : %d" % arc.W)
    print("  K (sync)      : %d bytes" % arc.K)
    print("  sub-blocks    : %d" % arc.m)
    print("  model hash    : %s" % arc.model_hash.hex())
    print("  verification v: %s" % arc.v.hex())
    if two:
        kept = [len(t) for t in arc.model.kept]
        print("  model section : %d bytes (%.4f bpb), kept contexts %s"
              % (arc.model_bytes,
                 8.0 * arc.model_bytes / max(1, arc.n), kept))
    print("  archive size  : %d bytes (%.4f bpb)"
          % (len(raw), 8.0 * len(raw) / max(1, arc.n)))
    print("  repair stream : %d bytes (%.4f bpb)"
          % (len(arc.repair), 8.0 * len(arc.repair) / max(1, arc.n)))


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="two_pass",
        description="Two-pass (train-then-freeze) RNR Type-I coder with "
                    "single-pass fallback (never loses on size).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("pack", help="encode (two-pass vs single-pass, "
                                    "smaller wins)")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--W", type=int, default=DEFAULT_W)
    p.add_argument("--K", type=int, default=DEFAULT_K)
    p.add_argument("--force", choices=["single", "two"], default=None,
                   help="testing knob: skip the size comparison")
    p.set_defaults(func=_cmd_pack)

    p = sub.add_parser("unpack", help="full sequential decode (either mode)")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--no-verify", action="store_true")
    p.set_defaults(func=_cmd_unpack)

    p = sub.add_parser("read", help="random-access read (either mode)")
    p.add_argument("archive")
    p.add_argument("--pos", type=int, required=True)
    p.add_argument("--len", type=int, required=True)
    p.add_argument("--out", default=None)
    p.add_argument("--verify", action="store_true")
    p.set_defaults(func=_cmd_read)

    p = sub.add_parser("info", help="print archive header summary")
    p.add_argument("archive")
    p.set_defaults(func=_cmd_info)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
