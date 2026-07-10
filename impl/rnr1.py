#!/usr/bin/env python3
"""
rnr1.py -- Reference RNR Type-I coder (Code12 E12-A, order-W byte context
predictor, binary arithmetic coder, sub-block/sync architecture).

Scope
-----
Correctness-first reference implementation of the RNR Type-I stream-repair
coder (main paper Definition 3.1, Sections 4.2-4.4) with the random-access
sub-block architecture (Definitions 10.2/10.3, Theorem 10.3) and the
integer-determinism conformance requirements of the engineering spec
(R-3.1..R-3.8, R-5.x, R-8.x, R-11.1, R-13.x).  Archive layout is documented
in FORMAT.md next to this file.

Integer determinism (R-3.1..R-3.8)
----------------------------------
The entire coding path (predictor inference, class/payload partitioning,
arithmetic coding, adaptive count updates) uses integer arithmetic only:
Python arbitrary-precision ints and numpy int64/uint8 arrays.  There are no
floating-point operations anywhere in the coding path (R-3.1); floats are
permitted only in offline analysis tooling (R-3.6), which lives in
run_checks.py, not here.  All reductions are exact integer sums, hence
order-independent with a well-defined value; the implementation nevertheless
uses fixed sequential/cumulative orders (np.cumsum, np.add.at, explicit
left-to-right loops) so the evaluation order is also syntactically fixed
(R-3.4, R-5.4).  Accumulator magnitudes are bounded by construction
(total frequency < 2^24, see PROB_TOTAL_BOUND) and guarded by assertions,
far below int64 range (R-3.2, R-5.2).  No randomness (R-3.5), no runtime
approximations (R-3.3), no FMA (R-3.8).

Coder architecture (Type-I, Definition 3.1 / Section 4.2-4.3)
-------------------------------------------------------------
Each source byte x is mapped through the injective systematic Code12 map
E12-A to a 12-bit codeword c = E12(x).  The predictor deterministically
predicts a codeword c_hat = E12(b*), where b* is the argmax byte of the
predictor's integer distribution (deterministic tie-break: smallest byte).
The error mask is e = c XOR c_hat, and the repair symbol is coded as a
class + payload pair per Section 4.3:

  class 0     : e == 0  (prediction exact)
  class 1     : Hamming weight(e) == 1
  class 2     : Hamming weight(e) == 2
  class CAND  : x is the k-th entry (k = 1..15) of the predicted candidate
                list (next-most-probable bytes after b*)
  class ESC   : everything else (payload = the byte, coded under the
                predictor's conditional distribution restricted to ESC)

Class probabilities and within-class payload probabilities are derived
from the same integer predictor distribution by exact partition sums, so
by the chain rule the two-symbol code costs exactly
-log2(freq[x]/total) per byte, up to arithmetic-coder quantization.
Note: because E12-A is linear and both c and c_hat are valid codewords,
e is itself a codeword and its weight is 0 or >= 3 (minimum distance 3);
classes 1 and 2 are therefore structurally empty in this byte-aligned
reference coder and carry zero probability mass.  They are retained in
the class alphabet for format compatibility with codeword-level
predictors (see FORMAT.md).

Random access (Definitions 10.2/10.3, Theorem 10.3)
---------------------------------------------------
Sync points at uniform spacing K: at each sub-block boundary the predictor
context/adaptive state is reset to Q_init (empty counts) and the arithmetic
coder is flushed and byte-aligned, so each sub-block is independently
decodable.  A seek index in the header maps sync point j to its repair-
stream bit offset (R-8.1.1/R-8.2.1) plus an 8-byte truncated SHA-256 of
the sub-block's source bytes (R-11.1.1).  read(p, k) decodes only the
covering sub-block(s), discarding at most K-1 warm-up bytes, i.e. at most
K + k decoded positions total (Theorem 10.3 access pattern).

Predictor (Definition 10.2 semantics)
-------------------------------------
The concrete predictor is an order-W byte context model: the distribution
at position i is read from the count table addressed by the last
min(W, i - P) bytes only (P = latest sync point), backing off to shorter
contexts when the addressed context has no observations.  The context
argument is W-bounded by construction.  The count tables adapt (increment
by 1 per observed byte) deterministically as coding proceeds -- this is
the compatible-online-optimizer path (Definition 10.1, Theorem 10.2), an
integer state update -- and the whole adaptive state is reset at every
sync point, so the Definition 10.3(a) requirement (predictor state at a
sync point is a fixed Q_init independent of prior bytes) and the
Theorem 10.3 decode-window bound hold verbatim.

Probability floor: every byte has frequency >= 1 out of a total < 2^24,
so the coder-level floor of Section 10.7 holds with eta >= 2^-24 (bounded
worst-case cost of 24 bits/byte per position).  Independently, a per-sub-
block fail-safe store mode (raw copy when the coded blob is not smaller,
Theorem 7.22 semantics) bounds the whole archive at 8 bits/byte plus
header/index overhead on adversarial input.
"""

import argparse
import hashlib
import struct
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Format constants (see FORMAT.md)
# ---------------------------------------------------------------------------

MAGIC = b"RNR1"
VERSION = (1, 0, 0)  # semantic versioning, R-13.1
FLAG_SUBBLOCK_HASHES = 0x0001  # R-11.1.1 per-sync-point hashes present

DEFAULT_W = 3
DEFAULT_K = 65536

# Predictor integer parameters (part of the model description, R-13.3)
COUNT_SCALE = 32          # counts contribute SCALE per observation
COUNT_CAP = 1 << 16       # per-context total cap; halve on reaching it
NUM_CANDIDATES = 15       # candidate-list length (ranks 1..15 after argmax)
PROB_TOTAL_BOUND = 1 << 24  # invariant: arithmetic-coder totals stay below

# Class alphabet (Section 4.3 of the main paper)
CLS_0, CLS_1, CLS_2, CLS_CAND, CLS_ESC = 0, 1, 2, 3, 4
NUM_CLASSES = 5

HEADER_FMT = "<4s3sB H H I Q 8s 8s Q"  # see FORMAT.md
HEADER_SIZE = struct.calcsize(HEADER_FMT)
# byte_offset, repair_bit_offset, snapshot_len, mode, hash8
INDEX_ENTRY_FMT = "<QQII8s"
INDEX_ENTRY_SIZE = struct.calcsize(INDEX_ENTRY_FMT)

# Sub-block modes: fail-safe store-mode escape bounds the worst case near
# 8 bits/byte on incompressible/adversarial input (Theorem 7.22 semantics).
MODE_CODED = 0
MODE_RAW = 1


# ---------------------------------------------------------------------------
# Code12 E12-A: shortened-Hamming-like systematic (12,8) map, Section 4.2
# ---------------------------------------------------------------------------

def _build_e12a():
    """Build the E12-A lookup table exactly per Section 4.2 of the main paper.

    1-based codeword positions c1..c12; data bits b0..b7 at positions
    3,5,6,7,9,10,11,12; even parity at positions 1,2,4,8.
    Bit i of the returned integer holds c_{i+1} (i.e. c1 is bit 0).
    """
    table = np.zeros(256, dtype=np.int64)
    for b in range(256):
        bits = [(b >> j) & 1 for j in range(8)]  # b0..b7
        c = [0] * 13  # 1-based
        c[3], c[5], c[6], c[7] = bits[0], bits[1], bits[2], bits[3]
        c[9], c[10], c[11], c[12] = bits[4], bits[5], bits[6], bits[7]
        c[1] = c[3] ^ c[5] ^ c[7] ^ c[9] ^ c[11]
        c[2] = c[3] ^ c[6] ^ c[7] ^ c[10] ^ c[11]
        c[4] = c[5] ^ c[6] ^ c[7] ^ c[12]
        c[8] = c[9] ^ c[10] ^ c[11] ^ c[12]
        word = 0
        for i in range(1, 13):
            word |= c[i] << (i - 1)
        table[b] = word
    return table


E12 = _build_e12a()

# WT[a, b] = Hamming weight of E12[a] XOR E12[b]; precomputed offline
# (integer-only construction, deterministic).
_x = np.bitwise_xor.outer(E12, E12)
WT = np.zeros((256, 256), dtype=np.uint8)
for _bit in range(12):
    WT += ((_x >> _bit) & 1).astype(np.uint8)
del _x


def model_description_hash(W):
    """Truncated SHA-256 of the canonical model description (R-13.3)."""
    desc = (
        "rnr-type1-code12a;predictor=ngram;W=%d;scale=%d;cap=%d;cand=%d;"
        "floor=1;classes=c0,c1,c2,cand,esc;e12=%s"
        % (
            W,
            COUNT_SCALE,
            COUNT_CAP,
            NUM_CANDIDATES,
            hashlib.sha256(E12.astype("<i8").tobytes()).hexdigest()[:16],
        )
    ).encode("ascii")
    return hashlib.sha256(desc).digest()[:8]


def h_v(data):
    """Verification value H_v (Definition 3.1): truncated SHA-256."""
    return hashlib.sha256(data).digest()[:8]


# ---------------------------------------------------------------------------
# Binary arithmetic coder: 32-bit registers, bit renormalization,
# carry handling via pending-bit counting (integer-only path).
# ---------------------------------------------------------------------------

_AC_BITS = 32
_AC_TOP = 1 << _AC_BITS
_AC_MASK = _AC_TOP - 1
_AC_HALF = _AC_TOP >> 1
_AC_Q1 = _AC_TOP >> 2
_AC_Q3 = _AC_HALF + _AC_Q1


class ArithmeticEncoder:
    """Witten-Neal-Cleary style integer arithmetic encoder (bit output)."""

    __slots__ = ("low", "high", "pending", "_buf", "_cur", "_ncur")

    def __init__(self):
        self.low = 0
        self.high = _AC_MASK
        self.pending = 0
        self._buf = bytearray()
        self._cur = 0
        self._ncur = 0

    def _put(self, bit):
        self._cur = (self._cur << 1) | bit
        self._ncur += 1
        if self._ncur == 8:
            self._buf.append(self._cur)
            self._cur = 0
            self._ncur = 0

    def _emit(self, bit):
        self._put(bit)
        inv = 1 - bit
        while self.pending:
            self._put(inv)
            self.pending -= 1

    def encode(self, cum_lo, cum_hi, total):
        """Encode a symbol occupying [cum_lo, cum_hi) of [0, total)."""
        assert 0 <= cum_lo < cum_hi <= total < PROB_TOTAL_BOUND
        rng = self.high - self.low + 1
        self.high = self.low + (rng * cum_hi) // total - 1
        self.low = self.low + (rng * cum_lo) // total
        while True:
            if self.high < _AC_HALF:
                self._emit(0)
            elif self.low >= _AC_HALF:
                self._emit(1)
                self.low -= _AC_HALF
                self.high -= _AC_HALF
            elif self.low >= _AC_Q1 and self.high < _AC_Q3:
                self.pending += 1
                self.low -= _AC_Q1
                self.high -= _AC_Q1
            else:
                break
            self.low <<= 1
            self.high = (self.high << 1) | 1
            assert self.high <= _AC_MASK

    def finish(self):
        """Flush (Definition 10.3(b)): terminate and byte-align the block."""
        self.pending += 1
        self._emit(0 if self.low < _AC_Q1 else 1)
        while self._ncur:
            self._put(0)
        return bytes(self._buf)


class ArithmeticDecoder:
    """Mirror-image integer arithmetic decoder; reads 0 past end of block."""

    __slots__ = ("low", "high", "value", "_data", "_bitpos")

    def __init__(self, data):
        self.low = 0
        self.high = _AC_MASK
        self.value = 0
        self._data = data
        self._bitpos = 0
        for _ in range(_AC_BITS):
            self.value = (self.value << 1) | self._get()

    def _get(self):
        i = self._bitpos >> 3
        bit = 0
        if i < len(self._data):
            bit = (self._data[i] >> (7 - (self._bitpos & 7))) & 1
        self._bitpos += 1
        return bit

    def target(self, total):
        """Scaled position of the coded point inside the current interval."""
        rng = self.high - self.low + 1
        t = ((self.value - self.low + 1) * total - 1) // rng
        assert 0 <= t < total
        return t

    def consume(self, cum_lo, cum_hi, total):
        """Advance past a symbol occupying [cum_lo, cum_hi) of [0, total)."""
        rng = self.high - self.low + 1
        self.high = self.low + (rng * cum_hi) // total - 1
        self.low = self.low + (rng * cum_lo) // total
        while True:
            if self.high < _AC_HALF:
                pass
            elif self.low >= _AC_HALF:
                self.low -= _AC_HALF
                self.high -= _AC_HALF
                self.value -= _AC_HALF
            elif self.low >= _AC_Q1 and self.high < _AC_Q3:
                self.low -= _AC_Q1
                self.high -= _AC_Q1
                self.value -= _AC_Q1
            else:
                break
            self.low <<= 1
            self.high = (self.high << 1) | 1
            self.value = (self.value << 1) | self._get()


# ---------------------------------------------------------------------------
# Predictor: order-W adaptive byte context model (integer path)
# ---------------------------------------------------------------------------

class NGramPredictor:
    """Order-W byte context predictor with deterministic adaptive counts.

    Prediction context is the last min(W, bytes-since-sync) bytes only
    (Definition 10.2 semantics for the context argument); adaptive count
    state is a deterministic integer function of the bytes coded since the
    last sync point and is reset to Q_init (empty) at every sync point
    (Definition 10.3(a), Definition 10.1 / Theorem 10.2 adaptive path).

    dist() returns (freq, total): int64 frequency vector over the 256
    bytes, freq[b] >= 1 (probability floor), total = freq.sum() < 2^24.
    """

    __slots__ = ("W", "tables", "hist")

    def __init__(self, W):
        self.W = W
        self.reset()

    def reset(self):
        """Reset to Q_init: empty count tables, empty context history."""
        self.tables = [dict() for _ in range(self.W + 1)]
        self.hist = []

    def _ctx_key(self, order):
        h = self.hist
        key = 0
        for b in h[len(h) - order:]:
            key = (key << 8) | b
        return key

    def dist(self):
        """Integer distribution from the longest non-empty context."""
        avail = min(len(self.hist), self.W)
        for order in range(avail, -1, -1):
            ent = self.tables[order].get(self._ctx_key(order))
            if ent is not None and ent[1] > 0:
                freq = 1 + COUNT_SCALE * ent[0]
                total = 256 + COUNT_SCALE * ent[1]
                assert total < PROB_TOTAL_BOUND
                return freq, total
        return np.ones(256, dtype=np.int64), 256

    def update(self, byte):
        """Deterministic integer count update (Definition 10.1 semantics)."""
        avail = min(len(self.hist), self.W)
        for order in range(avail + 1):
            key = self._ctx_key(order)
            ent = self.tables[order].get(key)
            if ent is None:
                ent = [np.zeros(256, dtype=np.int64), 0]
                self.tables[order][key] = ent
            ent[0][byte] += 1
            ent[1] += 1
            if ent[1] >= COUNT_CAP:
                ent[0] >>= 1
                ent[1] = int(ent[0].sum())
        self.hist.append(byte)
        if len(self.hist) > self.W:
            del self.hist[: len(self.hist) - self.W]


# ---------------------------------------------------------------------------
# Type-I class/payload partition (Section 4.3)
# ---------------------------------------------------------------------------

def _partition(freq):
    """Deterministic Type-I partition derived from the predictor output.

    Returns (bstar, cand, classid, class_freq):
      bstar      : argmax byte (ties -> smallest byte value)
      cand       : candidate list, ranks 1..NUM_CANDIDATES (by freq desc,
                   byte value asc on ties; deterministic via stable sort)
      classid    : int64[256], class of each byte with precedence
                   C0 > C1 > C2 > CAND > ESC (Section 4.3 listing order)
      class_freq : int64[NUM_CLASSES], exact partition sums of freq
    """
    order = np.argsort(-freq, kind="stable")
    bstar = int(order[0])
    cand = order[1: 1 + NUM_CANDIDATES]
    wt_row = WT[bstar]
    classid = np.full(256, CLS_ESC, dtype=np.int64)
    classid[cand] = CLS_CAND
    classid[wt_row == 2] = CLS_2
    classid[wt_row == 1] = CLS_1
    classid[bstar] = CLS_0
    class_freq = np.zeros(NUM_CLASSES, dtype=np.int64)
    np.add.at(class_freq, classid, freq)
    return bstar, cand, classid, class_freq


def _class_members(cls, cand, classid):
    """Canonical within-class payload ordering (decoder-reproducible).

    CAND: candidate-list rank order (Section 4.3 'k-th in candidate list');
    other classes: ascending byte value.
    """
    if cls == CLS_CAND:
        return cand[classid[cand] == CLS_CAND]
    return np.nonzero(classid == cls)[0]


# ---------------------------------------------------------------------------
# Sub-block encode / decode
# ---------------------------------------------------------------------------

def encode_subblock(block, W):
    """Encode one sub-block (fresh predictor state, fresh coder state)."""
    pred = NGramPredictor(W)
    enc = ArithmeticEncoder()
    for x in block:
        freq, total = pred.dist()
        bstar, cand, classid, class_freq = _partition(freq)
        cls = int(classid[x])
        ccum = np.concatenate(([0], np.cumsum(class_freq)))
        enc.encode(int(ccum[cls]), int(ccum[cls + 1]), total)
        if cls != CLS_0:
            members = _class_members(cls, cand, classid)
            mcum = np.concatenate(([0], np.cumsum(freq[members])))
            idx = int(np.nonzero(members == x)[0][0])
            enc.encode(int(mcum[idx]), int(mcum[idx + 1]), int(mcum[-1]))
        pred.update(x)
    return enc.finish()


def decode_subblock(blob, nbytes, W):
    """Decode the first nbytes positions of a sub-block blob."""
    pred = NGramPredictor(W)
    dec = ArithmeticDecoder(blob)
    out = bytearray(nbytes)
    for i in range(nbytes):
        freq, total = pred.dist()
        bstar, cand, classid, class_freq = _partition(freq)
        ccum = np.concatenate(([0], np.cumsum(class_freq)))
        t = dec.target(total)
        cls = int(np.searchsorted(ccum, t, side="right")) - 1
        dec.consume(int(ccum[cls]), int(ccum[cls + 1]), total)
        if cls == CLS_0:
            x = bstar
        else:
            members = _class_members(cls, cand, classid)
            mcum = np.concatenate(([0], np.cumsum(freq[members])))
            t = dec.target(int(mcum[-1]))
            idx = int(np.searchsorted(mcum, t, side="right")) - 1
            dec.consume(int(mcum[idx]), int(mcum[idx + 1]), int(mcum[-1]))
            x = int(members[idx])
        out[i] = x
        pred.update(x)
    return bytes(out)


# ---------------------------------------------------------------------------
# Container: pack / unpack / random-access read (FORMAT.md)
# ---------------------------------------------------------------------------

def pack(data, W=DEFAULT_W, K=DEFAULT_K):
    """Encode data into an RNR1 archive (bytes)."""
    if W < 1 or W > 8:
        raise ValueError("W must be in [1, 8]")
    if K < 256:
        raise ValueError("K must be >= 256")
    n = len(data)
    m = (n + K - 1) // K
    blobs = []
    entries = []
    repair_bit_offset = 0
    for j in range(m):
        block = data[j * K: (j + 1) * K]
        blob = encode_subblock(block, W)
        mode = MODE_CODED
        if len(blob) >= len(block):
            # Fail-safe store mode (Theorem 7.22 semantics): never pay more
            # than raw for a sub-block of incompressible/adversarial data.
            blob = block
            mode = MODE_RAW
        entries.append(
            struct.pack(
                INDEX_ENTRY_FMT,
                j * K,               # byte_offset in source (uniform: jK)
                repair_bit_offset,   # bit offset in repair stream
                0,                   # snapshot_length: 0, W-bounded (R-8.1.1)
                mode,                # MODE_CODED or MODE_RAW
                h_v(block),          # per-sync-point hash (R-11.1.1)
            )
        )
        blobs.append(blob)
        repair_bit_offset += 8 * len(blob)
    header = struct.pack(
        HEADER_FMT,
        MAGIC,
        bytes(VERSION),
        0,                       # reserved
        FLAG_SUBBLOCK_HASHES,    # flags
        W,
        K,
        n,
        model_description_hash(W),
        h_v(data),               # archive verification value v (Def 3.1)
        m,
    )
    return header + b"".join(entries) + b"".join(blobs)


class Archive:
    """Parsed RNR1 archive with sequential and random-access decoding."""

    def __init__(self, raw):
        if len(raw) < HEADER_SIZE:
            raise ValueError("truncated archive: no header")
        (
            magic,
            version,
            _reserved,
            self.flags,
            self.W,
            self.K,
            self.n,
            self.model_hash,
            self.v,
            self.m,
        ) = struct.unpack_from(HEADER_FMT, raw, 0)
        if magic != MAGIC:
            raise ValueError("bad magic: not an RNR1 archive")
        if version[0] != VERSION[0] or version[1] > VERSION[1]:
            # R-13.2: same MAJOR required; refuse newer MINOR.
            raise ValueError("unsupported archive format version %s" % (tuple(version),))
        if self.model_hash != model_description_hash(self.W):
            # R-13.3: predictor hash mismatch is a hard failure.
            raise ValueError("model description hash mismatch")
        expect_m = (self.n + self.K - 1) // self.K
        if self.m != expect_m:
            raise ValueError("inconsistent sub-block count")
        self.entries = []
        off = HEADER_SIZE
        for _ in range(self.m):
            self.entries.append(struct.unpack_from(INDEX_ENTRY_FMT, raw, off))
            off += INDEX_ENTRY_SIZE
        self.repair = raw[off:]
        # R-8.2.1: index sorted by byte_offset ascending.
        for j in range(1, self.m):
            if self.entries[j][0] <= self.entries[j - 1][0]:
                raise ValueError("seek index not sorted")

    def _find_sync(self, p):
        """Binary search for the latest sync point <= p (R-8.2.2)."""
        lo, hi = 0, self.m - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.entries[mid][0] <= p:
                lo = mid
            else:
                hi = mid - 1
        return lo

    def _block_blob(self, j):
        bit_off = self.entries[j][1]
        assert bit_off % 8 == 0  # blocks are byte-aligned by the flush
        start = bit_off // 8
        end = self.entries[j + 1][1] // 8 if j + 1 < self.m else len(self.repair)
        return self.repair[start:end]

    def _block_len(self, j):
        return min(self.K, self.n - j * self.K)

    def decode_block(self, j, upto=None, verify=True):
        """Decode sub-block j (optionally only its first upto bytes)."""
        blk_len = self._block_len(j)
        take = blk_len if upto is None else min(upto, blk_len)
        if self.entries[j][3] == MODE_RAW:
            out = self._block_blob(j)[:take]
        else:
            out = decode_subblock(self._block_blob(j), take, self.W)
        if verify and take == blk_len and (self.flags & FLAG_SUBBLOCK_HASHES):
            if h_v(out) != self.entries[j][4]:
                raise ValueError("sub-block %d hash mismatch" % j)
        return out

    def unpack(self, verify=True):
        """Full sequential decode; checks H_v(X) = v (Definition 3.1)."""
        out = b"".join(self.decode_block(j, verify=verify) for j in range(self.m))
        if verify and h_v(out) != self.v:
            raise ValueError("archive verification failed: H_v mismatch")
        return out

    def read(self, p, k, verify=False):
        """Random access: return source bytes [p, p+k') with
        k' = min(k, n-p), decoding only the covering sub-block(s).

        Returns (data, stats) where stats reports warmup_bytes (discarded
        decoded positions before p; always <= K-1) and decoded_bytes
        (total decoded positions; always <= K + k, Theorem 10.3).
        With verify=True, covering sub-blocks are decoded in full so their
        R-11.1.1 hashes can be checked (decode window then <= 2K + k).
        """
        if p < 0 or p > self.n:
            raise ValueError("position out of range")
        kprime = min(k, self.n - p)
        if kprime <= 0:
            return b"", {"warmup_bytes": 0, "decoded_bytes": 0}
        j0 = self._find_sync(p)
        j1 = self._find_sync(p + kprime - 1)
        pieces = []
        decoded = 0
        for j in range(j0, j1 + 1):
            base = self.entries[j][0]
            if j < j1 or verify:
                upto = None  # full block
            else:
                upto = (p + kprime) - base
            blk = self.decode_block(j, upto=upto, verify=verify)
            decoded += len(blk)
            lo = max(p, base) - base
            hi = min(p + kprime, base + len(blk)) - base
            pieces.append(blk[lo:hi])
        stats = {"warmup_bytes": p - self.entries[j0][0], "decoded_bytes": decoded}
        return b"".join(pieces), stats


def unpack(raw, verify=True):
    return Archive(raw).unpack(verify=verify)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cmd_pack(args):
    with open(args.input, "rb") as f:
        data = f.read()
    raw = pack(data, W=args.W, K=args.K)
    with open(args.output, "wb") as f:
        f.write(raw)
    print("packed %d bytes -> %d bytes (%.4f bpb), %d sub-block(s)"
          % (len(data), len(raw),
             8.0 * len(raw) / max(1, len(data)), (len(data) + args.K - 1) // args.K))


def _cmd_unpack(args):
    with open(args.input, "rb") as f:
        raw = f.read()
    data = unpack(raw, verify=not args.no_verify)
    with open(args.output, "wb") as f:
        f.write(data)
    print("unpacked %d bytes -> %d bytes (verified: %s)"
          % (len(raw), len(data), "no" if args.no_verify else "yes"))


def _cmd_read(args):
    with open(args.archive, "rb") as f:
        arc = Archive(f.read())
    data, stats = arc.read(args.pos, args.len, verify=args.verify)
    if args.out:
        with open(args.out, "wb") as f:
            f.write(data)
    else:
        sys.stdout.buffer.write(data)
        sys.stdout.buffer.flush()
    print("\nread %d byte(s) at %d: warmup=%d decoded=%d"
          % (len(data), args.pos, stats["warmup_bytes"], stats["decoded_bytes"]),
          file=sys.stderr)


def _cmd_info(args):
    with open(args.archive, "rb") as f:
        raw = f.read()
    arc = Archive(raw)
    print("RNR1 archive: version=%d.%d.%d flags=0x%04x" % (*VERSION, arc.flags))
    print("  source length : %d bytes" % arc.n)
    print("  W (order)     : %d" % arc.W)
    print("  K (sync)      : %d bytes" % arc.K)
    print("  sub-blocks    : %d" % arc.m)
    print("  model hash    : %s" % arc.model_hash.hex())
    print("  verification v: %s" % arc.v.hex())
    print("  archive size  : %d bytes (%.4f bpb)"
          % (len(raw), 8.0 * len(raw) / max(1, arc.n)))
    print("  repair stream : %d bytes (%.4f bpb)"
          % (len(arc.repair), 8.0 * len(arc.repair) / max(1, arc.n)))


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="rnr1",
        description="Reference RNR Type-I coder (Code12 E12-A, order-W "
                    "context predictor, arithmetic coding, sync sub-blocks).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("pack", help="encode a file into an RNR1 archive")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--W", type=int, default=DEFAULT_W, help="context order")
    p.add_argument("--K", type=int, default=DEFAULT_K, help="sync spacing (bytes)")
    p.set_defaults(func=_cmd_pack)

    p = sub.add_parser("unpack", help="full sequential decode")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--no-verify", action="store_true")
    p.set_defaults(func=_cmd_unpack)

    p = sub.add_parser("read", help="random-access read of k bytes at position p")
    p.add_argument("archive")
    p.add_argument("--pos", type=int, required=True)
    p.add_argument("--len", type=int, required=True)
    p.add_argument("--out", default=None)
    p.add_argument("--verify", action="store_true",
                   help="decode covering sub-blocks fully and check hashes")
    p.set_defaults(func=_cmd_read)

    p = sub.add_parser("info", help="print archive header summary")
    p.add_argument("archive")
    p.set_defaults(func=_cmd_info)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
