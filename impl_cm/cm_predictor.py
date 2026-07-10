#!/usr/bin/env python3
"""cm_predictor.py -- integer-deterministic context-mixing (CM) predictor
for the RNR Type-I reference coder (H1 instrument).

Architecture
------------
Logistic mixing of order-0..W adaptive byte context models (the classic
context-mixing construction: binary decomposition of the byte, per-node
probabilities mixed in the stretch = logit domain, gradient-descent
weight adaptation), implemented entirely in fixed-point integer
arithmetic behind the exact Predictor interface of impl/rnr1.py:

    dist()      -> (freq, total)  int64 frequency table over 256 bytes,
                   freq[b] >= 1, total < 2^24 (PROB_TOTAL_BOUND)
    update(sym) -> deterministic integer state update
    reset()     -> Q_init (Definition 10.3(a)): empty tables, initial
                   weights, empty history

Per position the byte is decomposed MSB-first over the 255 internal
nodes of the binary alphabet tree.  For each model order o = 0..W the
node bit-probability is the KT estimator over exact subtree count sums
p1 = (2*S1 + 1) / (2*(S0 + S1 + 1)) in 12-bit fixed point; each order's
node probability is mapped to the stretch domain by a lookup table, the
W+1 stretches are mixed with per-node 16.16 fixed-point weights, the mix
is mapped back by a table-based squash, and the 256 leaf probabilities
are exact fixed-point products down the tree, emitted as an integer
frequency table (scale 2^22, +1 floor per symbol).

After coding a byte, the mixer weights of the 8 nodes on its path are
updated by integer gradient descent
    w += (stretch * 7 * ((bit<<12) - p_mix) + 2^15) >> 16
(the standard integer CM mixer step) with symmetric clamping, and the
order-0..W count tables are incremented with cap-halving identical in
semantics to impl/rnr1.py's NGramPredictor.

Engineering-spec conformance (R-3.1..R-3.8)
-------------------------------------------
R-3.1 integer arithmetic : the whole inference path is Python ints and
      numpy int64; no floats anywhere in dist()/update()/reset().
R-3.2 overflow avoidance : declared bounds -- counts < 2^16 per context
      (cap-halving), subtree sums < 2^24, KT numerator < 2^37, stretch
      in [-2047, 2047], weights clamped to +/-2^20, mixer dot product
      |sum| <= (W+1)*2047*2^20 < 2^42, leaf accumulator <= 2^22 * 4095
      < 2^34; all far inside int64, asserted where cheap.
R-3.3 non-linearities    : squash and stretch are deterministic
      integer-to-integer lookup tables (SQUASH, STRETCH) built once at
      import from integer anchor constants; no runtime approximation.
R-3.4 reduction order    : all reductions are exact integer sums (order
      independent with a well-defined value); the implementation uses
      fixed numpy axis-sum / pairwise tree orders, so the evaluation
      order is also syntactically fixed.
R-3.5 determinism        : no randomness anywhere; deterministic init.
R-3.6 initialization     : even table construction is integer-only, so
      nothing float-valued exists at parse time.
R-3.7 layer ordering     : the fixed pipeline counts -> subtree sums ->
      KT -> stretch -> mix -> squash -> leaf products is hard-coded.
R-3.8 no FMA             : plain integer multiply and add only.

The model description (all constants below plus the table anchors) is
hashed into the archive header by impl_cm/cm_coder.py, mirroring
impl/rnr1.py's R-13.3 mechanism.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Model constants (part of the model description)
# ---------------------------------------------------------------------------

PSCALE_BITS = 12                 # probabilities are 12-bit: [1, 4095] of 4096
PSCALE = 1 << PSCALE_BITS
STRETCH_MAX = 2047               # stretch/squash domain [-2047, 2047]
WEIGHT_FIX_BITS = 16             # mixer weights are 16.16 fixed point
WEIGHT_CLAMP = 1 << 20           # symmetric weight bound (overflow guard)
LR_NUM = 7                       # integer mixer learning-rate numerator
LR_ROUND = 1 << 15               # rounding constant of the >>16 update
ACC_INIT = 1 << 22               # leaf-product scale; total < 2^22 + 256
COUNT_CAP = 1 << 16              # per-context total cap (halve on reach),
                                 # same value as impl/rnr1.py

# Integer anchor table of the logistic squash (the classic 33-point
# integer CM anchor set; linear interpolation between anchors).
_SQUASH_ANCHORS = (
    1, 2, 3, 6, 10, 16, 27, 45, 73, 120, 194, 310, 488, 747, 1101, 1546,
    2047, 2549, 2994, 3348, 3607, 3785, 3901, 3965, 4023, 4050, 4068,
    4079, 4085, 4089, 4092, 4093, 4094,
)


def _squash_int(d):
    """Integer logistic: stretch domain -> probability in [1, 4095]."""
    if d >= STRETCH_MAX:
        return PSCALE - 1
    if d <= -STRETCH_MAX:
        return 1
    w = d & 127
    i = (d >> 7) + 16
    v = (_SQUASH_ANCHORS[i] * (128 - w) + _SQUASH_ANCHORS[i + 1] * w + 64) >> 7
    return min(max(v, 1), PSCALE - 1)


def _build_tables():
    """Build SQUASH (indexed by x + 2047) and STRETCH (indexed by p).

    Integer-only construction (R-3.6): STRETCH is the monotone inverse
    of SQUASH obtained by a single left-to-right sweep.
    """
    squash = np.zeros(2 * STRETCH_MAX + 1, dtype=np.int64)
    for x in range(-STRETCH_MAX, STRETCH_MAX + 1):
        squash[x + STRETCH_MAX] = _squash_int(x)
    stretch = np.zeros(PSCALE, dtype=np.int64)
    prev = 0
    for x in range(-STRETCH_MAX, STRETCH_MAX + 1):
        v = _squash_int(x)
        if v >= prev:
            stretch[prev:v + 1] = x
            prev = v + 1
    stretch[prev:] = STRETCH_MAX
    return squash, stretch


SQUASH, STRETCH = _build_tables()

# Byte tree paths, MSB first: heap node indices 1..255 and branch bits.
NODE_PATH = np.zeros((256, 8), dtype=np.int64)
BIT_PATH = np.zeros((256, 8), dtype=np.int64)
for _b in range(256):
    _i = 1
    for _l in range(8):
        _bit = (_b >> (7 - _l)) & 1
        NODE_PATH[_b, _l] = _i
        BIT_PATH[_b, _l] = _bit
        _i = 2 * _i + _bit


class CMPredictor:
    """Context-mixing predictor with the impl/rnr1.py Predictor interface.

    Same context semantics as NGramPredictor: the order-o context is the
    last o bytes coded since the most recent reset (sync point); orders
    with fewer than o bytes of history, or with no observations in the
    addressed context, contribute the neutral probability 1/2 (stretch
    0), so the mixer's fixed reduction is well-defined at every position.
    """

    __slots__ = ("W", "tables", "hist", "weights", "_st", "_pmix_nodes")

    def __init__(self, W):
        self.W = W
        self.reset()

    def reset(self):
        """Reset to Q_init: empty counts, initial weights, empty history."""
        self.tables = [dict() for _ in range(self.W + 1)]
        self.hist = []
        # Deterministic init: weights sum to ~1.0 in 16.16 fixed point.
        winit = (1 << WEIGHT_FIX_BITS) // (self.W + 1)
        self.weights = np.full((self.W + 1, 256), winit, dtype=np.int64)
        self._st = None
        self._pmix_nodes = None

    def _ctx_key(self, order):
        h = self.hist
        key = 0
        for b in h[len(h) - order:]:
            key = (key << 8) | b
        return key

    def dist(self):
        """Integer distribution by logistic mixing of orders 0..W."""
        M = self.W + 1
        counts = np.zeros((M, 256), dtype=np.int64)
        avail = min(len(self.hist), self.W)
        for o in range(avail + 1):
            ent = self.tables[o].get(self._ctx_key(o))
            if ent is not None:
                counts[o] = ent[0]
        # Exact subtree sums over the heap-layout byte tree (fixed
        # pairwise reduction order, R-3.4). Leaf of byte b is node 256+b.
        sums = np.zeros((M, 512), dtype=np.int64)
        sums[:, 256:] = counts
        size = 128
        while size >= 1:
            hi = 2 * size
            sums[:, size:hi] = sums[:, hi:2 * hi:2] + sums[:, hi + 1:2 * hi:2]
            size >>= 1
        S0 = sums[:, 2::2]          # node i in 1..255: left-child sums
        S1 = sums[:, 3::2]          # node i in 1..255: right-child sums
        # KT estimator, 12-bit fixed point, then clamp into [1, 4095].
        p1 = ((2 * S1 + 1) << PSCALE_BITS) // (2 * (S0 + S1 + 1))
        np.clip(p1, 1, PSCALE - 1, out=p1)
        st = STRETCH[p1]            # (M, 255) stretches in [-2047, 2047]
        mixed = (self.weights[:, 1:] * st).sum(axis=0) >> WEIGHT_FIX_BITS
        np.clip(mixed, -STRETCH_MAX, STRETCH_MAX, out=mixed)
        pmix_nodes = np.empty(256, dtype=np.int64)
        pmix_nodes[0] = 0           # heap index 0 unused
        pmix_nodes[1:] = SQUASH[mixed + STRETCH_MAX]
        # Leaf probabilities: exact fixed-point products down the tree.
        acc = np.full(256, ACC_INIT, dtype=np.int64)
        for lvl in range(8):
            pn = pmix_nodes[NODE_PATH[:, lvl]]
            pb = (PSCALE - pn) + BIT_PATH[:, lvl] * (2 * pn - PSCALE)
            acc = (acc * pb) >> PSCALE_BITS
        freq = acc + 1
        total = int(freq.sum())
        assert total < (1 << 24)    # coder invariant (PROB_TOTAL_BOUND)
        self._st = st
        self._pmix_nodes = pmix_nodes
        return freq, total

    def update(self, byte):
        """Mixer gradient step at the 8 path nodes + count updates."""
        if self._st is not None:
            path = NODE_PATH[byte]
            bits = BIT_PATH[byte]
            pm = self._pmix_nodes[path]
            err7 = ((bits << PSCALE_BITS) - pm) * LR_NUM
            stp = self._st[:, path - 1]
            dw = (stp * err7 + LR_ROUND) >> WEIGHT_FIX_BITS
            self.weights[:, path] += dw
            np.clip(self.weights, -WEIGHT_CLAMP, WEIGHT_CLAMP,
                    out=self.weights)
            self._st = None
            self._pmix_nodes = None
        avail = min(len(self.hist), self.W)
        for o in range(avail + 1):
            key = self._ctx_key(o)
            ent = self.tables[o].get(key)
            if ent is None:
                ent = [np.zeros(256, dtype=np.int64), 0]
                self.tables[o][key] = ent
            ent[0][byte] += 1
            ent[1] += 1
            if ent[1] >= COUNT_CAP:
                ent[0] >>= 1
                ent[1] = int(ent[0].sum())
        self.hist.append(byte)
        if len(self.hist) > self.W:
            del self.hist[: len(self.hist) - self.W]


def table_self_test():
    """Integer table sanity: ranges, monotonicity, inverse consistency."""
    assert SQUASH.shape == (2 * STRETCH_MAX + 1,)
    assert STRETCH.shape == (PSCALE,)
    assert int(SQUASH.min()) >= 1 and int(SQUASH.max()) <= PSCALE - 1
    assert np.all(np.diff(SQUASH) >= 0), "squash not monotone"
    assert np.all(np.diff(STRETCH) >= 0), "stretch not monotone"
    assert int(STRETCH.min()) >= -STRETCH_MAX
    assert int(STRETCH.max()) <= STRETCH_MAX
    # Inverse consistency: squash(stretch(p)) stays close to p, and the
    # neutral point maps to stretch ~ 0.
    ps = np.arange(1, PSCALE)
    back = SQUASH[STRETCH[ps] + STRETCH_MAX]
    assert int(np.abs(back - ps).max()) <= 64, "squash/stretch inverse drift"
    assert abs(int(STRETCH[PSCALE // 2])) <= 1
    return True


if __name__ == "__main__":
    table_self_test()
    # Micro sanity: the predictor should learn a repeating pattern.
    p = CMPredictor(3)
    data = (b"the quick brown fox jumps over the lazy dog. " * 40)
    import math
    bits = 0.0
    for x in data:
        freq, total = p.dist()
        bits += -math.log2(freq[x] / total)
        p.update(x)
    bpb = bits / len(data)
    print("table self-test PASS; repeating-text model cost %.3f bpb" % bpb)
    assert bpb < 3.0, "CM predictor failed to learn a trivial pattern"
