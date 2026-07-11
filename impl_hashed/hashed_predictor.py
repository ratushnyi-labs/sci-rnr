#!/usr/bin/env python3
"""hashed_predictor.py -- MEMORY-BOUNDED hashed context-mixing predictor
for the RNR Type-I reference coder.

Motivation
----------
The exact-dictionary context models of impl_cm/cm_predictor.py and
impl/rnr1.py store one growable table entry per *distinct* observed
context, so the resident set grows ~linearly with the number of bytes
coded and with the context order W.  A high order W on a large file
therefore OOMs.  This predictor decouples the resident set from both the
file size and (nearly) from W by replacing every growable per-context
table with a FIXED-SIZE, direct-mapped HASHED array that is allocated
once and never grows.

Design (bit-level hashed context mixing, lpaq/zpaq family)
---------------------------------------------------------
A byte is decomposed MSB-first over the 255 internal nodes of the binary
alphabet tree (identical decomposition to cm_predictor).  For every
context order o = 0..W and every tree node j = 1..255 the model keeps a
tiny non-stationary bit counter -- two small integer counts (n0, n1) --
addressed by a deterministic integer hash of (o, last-o-bytes, j) into a
direct-mapped array of 2^HBITS buckets.  Each bucket holds

    n0  uint8   count of 0-bits seen in this (order, context, node)
    n1  uint8   count of 1-bits
    chk uint8   truncated hash "check byte" (collision confirmation)

so a bucket is 3 bytes and one order's table is 3 * 2^HBITS bytes,
INDEPENDENT of the file size.  Peak resident set is thus fixed at
approximately 3 * (W+1) * 2^HBITS bytes plus the input buffer, decoupled
from N.  (A per-bucket 256-entry frequency vector -- the cm_predictor
representation -- would be 512 * 2^HBITS bytes per order and is exactly
what makes HBITS=22 infeasible; the bit-level counter is what lets the
2^22-bucket example fit.)

The node bit-probability is the KT estimator p1 = (2 n1 + 1) /
(2 (n0 + n1 + 1)) in 12-bit fixed point (same as cm_predictor); the W+1
per-order stretches are mixed with 16.16 fixed-point per-node weights,
squashed, and the 256 leaf probabilities are exact fixed-point products
down the tree -- byte-for-byte the same mixing/leaf pipeline as
cm_predictor, only the *source* of (n0, n1) changed from exact subtree
sums to hashed bit counters.

Collision policy (documented, deterministic)
--------------------------------------------
Direct-mapped, one check byte per bucket.
  * READ  : if the stored check byte does not equal the context's check
            byte, the bucket is treated as EMPTY (n0 = n1 = 0, i.e. the
            neutral bit probability 1/2, stretch 0).  A colliding context
            therefore contributes nothing rather than wrong evidence.
  * WRITE : on a check-byte mismatch the bucket is EVICTED -- overwritten
            with the check byte of the current context and the single new
            observation (n0/n1 = 1).  On a match the counter is
            incremented with cap-halving.
Both encoder and decoder replay identical bytes and identical history, so
they perform identical evictions and stay in lock-step; the model is
never shipped in the archive (fully adaptive, rebuilt left-to-right).

Determinism / integer arithmetic (R-3.1..R-3.8)
-----------------------------------------------
The whole inference path is Python ints and numpy integers (uint8/int64/
uint64).  The 64-bit hash uses uint64 wraparound multiply/xor/shift
(splitmix64), which is bit-for-bit reproducible and matches a C port.
No floats, no randomness, no runtime approximation.  SQUASH/STRETCH,
NODE_PATH, BIT_PATH and the mixer step are imported unchanged from
cm_predictor so the two predictors share one audited non-linearity.
"""

import sys
from pathlib import Path

import numpy as np

# Self-contained import: the shared integer machinery lives in impl_cm.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "impl_cm"))

import cm_predictor as _cm  # noqa: E402
from cm_predictor import (  # noqa: E402  shared, audited integer machinery
    SQUASH, STRETCH, NODE_PATH, BIT_PATH,
    PSCALE_BITS, PSCALE, STRETCH_MAX,
    WEIGHT_FIX_BITS, WEIGHT_CLAMP, LR_NUM, LR_ROUND, ACC_INIT,
)

# ---------------------------------------------------------------------------
# Hashed-model constants (part of the model description)
# ---------------------------------------------------------------------------

DEFAULT_HBITS = 22               # 2^22 = 4,194,304 buckets per order
BIT_COUNT_CAP = 250              # cap-halving threshold for n0/n1 (uint8-safe)

# splitmix64 constants (deterministic, C-portable via uint64 wraparound)
_SM_ADD = np.uint64(0x9E3779B97F4A7C15)
_SM_M1 = np.uint64(0xBF58476D1CE4E5B9)
_SM_M2 = np.uint64(0x94D049BB133111EB)
_S30, _S27, _S31 = np.uint64(30), np.uint64(27), np.uint64(31)
_NODE_MULT = np.uint64(0x100000001B3)     # FNV-style node spreader
_U64 = np.uint64

# node indices 1..255 as a uint64 vector (the 255 internal tree nodes)
_NODES_U64 = np.arange(1, 256, dtype=np.uint64)


def _splitmix64_scalar(x):
    """Pure-python splitmix64 on a python int (deterministic, mod 2^64)."""
    m = (1 << 64) - 1
    z = (x + 0x9E3779B97F4A7C15) & m
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & m
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & m
    z = z ^ (z >> 31)
    return z & m


def _splitmix64_vec(x):
    """splitmix64 over a uint64 numpy array (matches the scalar version)."""
    z = x + _SM_ADD
    z = (z ^ (z >> _S30)) * _SM_M1
    z = (z ^ (z >> _S27)) * _SM_M2
    z = z ^ (z >> _S31)
    return z


class HashedCMPredictor:
    """Bit-level hashed context-mixing predictor (rnr1 Predictor interface).

    Interface parity with impl/rnr1.py NGramPredictor / cm_predictor:
        dist()      -> (freq, total)  int64 freq[256] >= 1, total < 2^24
        update(sym) -> deterministic integer state update
        reset()     -> Q_init: zeroed hashed tables, initial weights, empty
                       history (the buckets are allocated ONCE and zeroed
                       in place on reset, so the resident set is fixed).

    Parameters
    ----------
    W     : context order (>= 1); orders 0..W are mixed.
    HBITS : log2 of the per-order bucket count (fixed memory knob).
    """

    __slots__ = ("W", "HBITS", "MASK", "N0", "N1", "CHK", "hist", "_rows",
                 "weights", "freeze", "_st", "_pmix_nodes", "_idx", "_chkv",
                 "_avail")

    def __init__(self, W, HBITS=DEFAULT_HBITS):
        if W < 1:
            raise ValueError("W must be >= 1")
        self.W = int(W)
        self.HBITS = int(HBITS)
        self.MASK = _U64((1 << self.HBITS) - 1)
        self.N0 = None
        self._rows = np.arange(self.W + 1, dtype=np.int64)[:, None]
        self.freeze = False       # multi-pass: freeze mixer weights at code
        self.reset()

    # -- lifecycle ---------------------------------------------------------

    def _alloc_or_zero(self):
        """Fixed 2D tables (M orders x 2^HBITS buckets), uint8 counters.

        Contiguous per-order rows let dist()/update() gather across all
        orders in a single vectorized numpy op (no per-order Python loop),
        while the resident set stays 3 * (W+1) * 2^HBITS bytes -- fixed by
        (W, HBITS) alone, independent of the file size."""
        size = 1 << self.HBITS
        M = self.W + 1
        if self.N0 is None or self.N0.shape != (M, size):
            self.N0 = np.zeros((M, size), dtype=np.uint8)
            self.N1 = np.zeros((M, size), dtype=np.uint8)
            self.CHK = np.zeros((M, size), dtype=np.uint8)
        else:
            self.N0[:] = 0
            self.N1[:] = 0
            self.CHK[:] = 0

    def reset(self):
        """Reset to Q_init: zero the fixed hashed tables and mixer weights.

        The buckets are allocated ONCE and zeroed in place on subsequent
        resets, so the resident set is fixed regardless of file size.
        """
        self._alloc_or_zero()
        self.hist = []
        winit = (1 << WEIGHT_FIX_BITS) // (self.W + 1)
        self.weights = np.full((self.W + 1, 256), winit, dtype=np.int64)
        self._clear_scratch()

    def reset_counts(self):
        """Zero the hashed counts + history but KEEP the mixer weights.

        Used by multi-pass training: after epoch weight adaptation, the
        coded pass restarts the counts from empty (which the decoder can
        reproduce) while retaining the converged weights."""
        self._alloc_or_zero()
        self.hist = []
        self._clear_scratch()

    def _clear_scratch(self):
        self._st = None
        self._pmix_nodes = None
        self._idx = None
        self._chkv = None
        self._avail = 0

    def load_weights(self, w):
        """Install a fixed mixer-weight matrix (semi-static coded pass)."""
        w = np.asarray(w, dtype=np.int64)
        assert w.shape == (self.W + 1, 256)
        self.weights = w.copy()

    # -- context addressing ------------------------------------------------

    def _ctx_key(self, order):
        h = self.hist
        key = 0
        for b in h[len(h) - order:]:
            key = (key << 8) | b
        return key

    def _addr_grid(self, avail):
        """Vector (idx, chk) of shape (avail+1, 255) for the current context.

        Row o addresses tree nodes 1..255 under (order o, last-o bytes):
        a single splitmix64 over the (order, node) grid, so the whole
        address computation is one vectorized op."""
        bases = np.empty(avail + 1, dtype=np.uint64)
        for o in range(avail + 1):
            bases[o] = _U64(self._order_base(o, self._ctx_key(o)))
        h = _splitmix64_vec(bases[:, None] + _NODES_U64[None, :] * _NODE_MULT)
        idx = (h & self.MASK).astype(np.int64)          # (avail+1, 255)
        chkv = (h >> _U64(56)).astype(np.uint8)
        return idx, chkv

    def _order_base(self, order, ctx_key):
        """Scalar splitmix seed for (order, context), shared by all nodes."""
        return _splitmix64_scalar((order * 0x1000193) ^ ctx_key
                                  ^ ((order + 1) << 40))

    # -- inference ---------------------------------------------------------

    def dist(self):
        """Integer distribution by logistic mixing of hashed orders 0..W."""
        M = self.W + 1
        avail = min(len(self.hist), self.W)
        idx, chkv = self._addr_grid(avail)              # (avail+1, 255)
        rows = self._rows[:avail + 1]
        valid = self.CHK[rows, idx] == chkv
        n0 = self.N0[rows, idx].astype(np.int64)
        n1 = self.N1[rows, idx].astype(np.int64)
        S0 = np.zeros((M, 255), dtype=np.int64)
        S1 = np.zeros((M, 255), dtype=np.int64)
        S0[:avail + 1] = np.where(valid, n0, 0)         # unavailable orders
        S1[:avail + 1] = np.where(valid, n1, 0)         # -> neutral (0,0)
        # KT estimator per node, 12-bit fixed point, clamp to [1, 4095].
        p1 = ((2 * S1 + 1) << PSCALE_BITS) // (2 * (S0 + S1 + 1))
        np.clip(p1, 1, PSCALE - 1, out=p1)
        st = STRETCH[p1]                                # (M, 255)
        mixed = (self.weights[:, 1:] * st).sum(axis=0) >> WEIGHT_FIX_BITS
        np.clip(mixed, -STRETCH_MAX, STRETCH_MAX, out=mixed)
        pmix_nodes = np.empty(256, dtype=np.int64)
        pmix_nodes[0] = 0
        pmix_nodes[1:] = SQUASH[mixed + STRETCH_MAX]
        # Leaf probabilities: exact fixed-point products down the tree.
        acc = np.full(256, ACC_INIT, dtype=np.int64)
        for lvl in range(8):
            pn = pmix_nodes[NODE_PATH[:, lvl]]
            pb = (PSCALE - pn) + BIT_PATH[:, lvl] * (2 * pn - PSCALE)
            acc = (acc * pb) >> PSCALE_BITS
        freq = acc + 1
        total = int(freq.sum())
        assert total < (1 << 24)                        # PROB_TOTAL_BOUND
        self._st = st
        self._pmix_nodes = pmix_nodes
        self._idx = idx                                 # (avail+1, 255)
        self._chkv = chkv
        self._avail = avail
        return freq, total

    def update(self, byte):
        """Mixer gradient step + hashed bit-counter update (cap-halving).

        Vectorized over the 8 tree nodes on the byte's path (and, for the
        counters, per order).  Duplicate bucket indices among the 8 path
        nodes (a same-context hash collision, ~8^2/2^HBITS, astronomically
        rare) resolve last-write-wins -- deterministic across encoder and
        decoder since both run this identical vectorized path.
        """
        byte = int(byte)
        path = NODE_PATH[byte]                           # (8,) node indices
        bits = BIT_PATH[byte]                            # (8,) branch bits
        if self._st is not None and not self.freeze:     # mixer gradient step
            pm = self._pmix_nodes[path]
            err7 = ((bits << PSCALE_BITS) - pm) * LR_NUM
            stp = self._st[:, path - 1]                  # (M, 8)
            dw = (stp * err7 + LR_ROUND) >> WEIGHT_FIX_BITS
            self.weights[:, path] += dw
            np.clip(self.weights, -WEIGHT_CLAMP, WEIGHT_CLAMP,
                    out=self.weights)
        avail = self._avail
        sel = path - 1                                   # node -> row in 0..254
        idx8 = self._idx[:, sel]                         # (avail+1, 8) buckets
        chk8 = self._chkv[:, sel]                        # (avail+1, 8) checks
        rows = self._rows[:avail + 1]
        hit = self.CHK[rows, idx8] == chk8
        n0v = self.N0[rows, idx8].astype(np.int64)
        n1v = self.N1[rows, idx8].astype(np.int64)
        n0v = np.where(hit, n0v, 0) + (1 - bits)         # evict-miss + observe
        n1v = np.where(hit, n1v, 0) + bits
        over = (n0v + n1v) >= BIT_COUNT_CAP
        n0v = np.where(over, n0v >> 1, n0v)
        n1v = np.where(over, n1v >> 1, n1v)
        self.CHK[rows, idx8] = chk8
        self.N0[rows, idx8] = n0v.astype(np.uint8)
        self.N1[rows, idx8] = n1v.astype(np.uint8)
        self._clear_scratch()
        self.hist.append(byte)
        if len(self.hist) > self.W:
            del self.hist[: len(self.hist) - self.W]

    # -- introspection -----------------------------------------------------

    def table_bytes(self):
        """Fixed resident bytes of the hashed tables (n0+n1+chk)."""
        return 3 * (self.W + 1) * (1 << self.HBITS)


def _self_test():
    """Determinism + monotone-learning smoke on a tiny repeating pattern."""
    import math
    p = HashedCMPredictor(4, HBITS=16)
    data = (b"the quick brown fox jumps over the lazy dog. " * 60)
    bits = 0.0
    for x in data:
        freq, total = p.dist()
        bits += -math.log2(freq[x] / total)
        p.update(x)
    bpb = bits / len(data)
    # Two fresh predictors must produce identical CE on identical input.
    q = HashedCMPredictor(4, HBITS=16)
    bits2 = 0.0
    for x in data:
        f, t = q.dist()
        bits2 += -math.log2(f[x] / t)
        q.update(x)
    assert abs(bits - bits2) < 1e-9, "non-deterministic"
    assert bpb < 3.0, "failed to learn a trivial pattern (%.3f bpb)" % bpb
    return bpb


if __name__ == "__main__":
    bpb = _self_test()
    print("hashed CM self-test PASS; repeating-text model cost %.3f bpb" % bpb)
