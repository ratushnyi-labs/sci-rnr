#!/usr/bin/env python3
"""Shared infrastructure for the H5/H6/H7 campaign (Type-III modes).

Wires the committed instruments together (never rebuilds them):
  impl/rnr1.py         -- Type-I reference coder (ngram predictor)
  impl_cm/cm_coder.py  -- same container with the CM predictor seam
  bench/               -- frozen methodology (metrics, results store)
  data/loader.py       -- manifest-checked corpus access

All statistics here are SECONDARY relative to the frozen H1-H18
primary families (bench/METHODS.md section 1 allows labeled secondary
procedures); every reported interval is labeled with its method.
"""

from __future__ import annotations

import bz2
import math
import sys
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
for sub in ("bench", "data", "impl", "impl_cm"):
    p = str(ROOT / sub)
    if p not in sys.path:
        sys.path.insert(0, p)

import loader          # noqa: E402
import metrics         # noqa: E402
import results         # noqa: E402
import rnr1            # noqa: E402
import cm_coder        # noqa: E402

OUT = HERE / "out"
OUT.mkdir(exist_ok=True)

EXPERIMENT_ID = "h5_h7_modes"


def derive_seed(tag: str, index: int) -> int:
    """Frozen seed policy (METHODS.md section 4)."""
    return metrics.derive_seed(EXPERIMENT_ID, index, tag=tag)


def corpus_sha(name: str) -> str:
    return loader.entry(name)["sha256"]


# ---------------------------------------------------------------------------
# Coding modes (H5).  Each mode maps bytes -> bit cost of a standalone,
# cold-start blob (bare-blob framing, uniform across modes; container
# and seek-index bytes are excluded uniformly, selector flag bits are
# accounted separately by the harness).
# ---------------------------------------------------------------------------

def _lzma_bits(data: bytes) -> int:
    import lzma
    return 8 * len(lzma.compress(data, preset=9 | lzma.PRESET_EXTREME))


def _zstd_bits(data: bytes) -> int:
    import zstandard
    return 8 * len(zstandard.ZstdCompressor(level=22).compress(data))


def _brotli_bits(data: bytes) -> int:
    import brotli
    return 8 * len(brotli.compress(data, quality=11))


def _ngram_bits(data: bytes, W: int = 3) -> int:
    blob = rnr1.encode_subblock(bytes(data), W)
    return 8 * min(len(blob), len(data))   # store-mode escape, Thm 7.22


def _cm_bits(data: bytes, W: int = 3) -> int:
    with cm_coder.cm_patched():
        blob = rnr1.encode_subblock(bytes(data), W)
    return 8 * min(len(blob), len(data))


MODES = {
    "store":      lambda d: 8 * len(d),
    "gzip9":      lambda d: 8 * len(zlib.compress(bytes(d), 9)),
    "bz2-9":      lambda d: 8 * len(bz2.compress(bytes(d), 9)),
    "lzma-9e":    _lzma_bits,
    "zstd-22":    _zstd_bits,
    "brotli-11":  _brotli_bits,
    "rnr1-ngram": _ngram_bits,
    "rnr1-cm":    _cm_bits,
}
MODE_NAMES = list(MODES)
MODE_FLAG_BITS = math.ceil(math.log2(len(MODES)))


def rnr_roundtrip_ok(data: bytes, W: int = 3, cm: bool = False) -> bool:
    """Byte-exact round-trip of the rnr sub-block coders on `data`."""
    if cm:
        with cm_coder.cm_patched():
            blob = rnr1.encode_subblock(bytes(data), W)
            if len(blob) >= len(data):
                return True  # store escape: trivially exact
            return rnr1.decode_subblock(blob, len(data), W) == bytes(data)
    blob = rnr1.encode_subblock(bytes(data), W)
    if len(blob) >= len(data):
        return True
    return rnr1.decode_subblock(blob, len(data), W) == bytes(data)


# ---------------------------------------------------------------------------
# Conditional mutual information I(H;L|C), Miller-Madow corrected (H6).
# Generalizes experiments/preconditions/precond_lib.nibble_cmi to an
# arbitrary integer context array.
# ---------------------------------------------------------------------------

def _entropy_from_counts(counts: np.ndarray, n: int) -> float:
    nz = counts[counts > 0].astype(float)
    return float(-(nz / n * np.log2(nz / n)).sum())


def cmi_mm(c: np.ndarray, s: np.ndarray, c_card: int) -> dict:
    """I(H;L|C) with H/L the nibbles of byte array s, context array c.

    Plug-in over the empirical joint law with Miller-Madow correction
    using observed supports (bias ~ (K_CHL - K_CH - K_CL + K_C)/(2 n ln 2),
    validated on an i.i.d. control in run_checks.py).
    Returns {"raw", "mm", "n", "k_joint"}.
    """
    c = np.asarray(c, dtype=np.int64)
    s = np.asarray(s, dtype=np.int64)
    assert c.shape == s.shape and c.min() >= 0 and c.max() < c_card
    n = s.size
    cnt_c = np.bincount(c, minlength=c_card)
    cnt_j = np.bincount(c * 256 + s, minlength=c_card * 256)
    cnt_ch = np.bincount(c * 16 + (s >> 4), minlength=c_card * 16)
    cnt_cl = np.bincount(c * 16 + (s & 15), minlength=c_card * 16)
    hc = _entropy_from_counts(cnt_c, n)
    hx_c = _entropy_from_counts(cnt_j, n) - hc
    hh_c = _entropy_from_counts(cnt_ch, n) - hc
    hl_c = _entropy_from_counts(cnt_cl, n) - hc
    raw = hh_c + hl_c - hx_c
    k_c = int((cnt_c > 0).sum())
    k_j = int((cnt_j > 0).sum())
    k_ch = int((cnt_ch > 0).sum())
    k_cl = int((cnt_cl > 0).sum())
    bias = (k_j - k_ch - k_cl + k_c) / (2.0 * n * math.log(2.0))
    return {"raw": float(raw), "mm": float(raw - bias), "n": int(n),
            "k_joint": k_j}


def cmi_perm_null(c, s, c_card, seed) -> float:
    """Permutation null for I(H;L|C): permute the LOW nibbles among
    positions sharing the same context value (seeded).  This generates
    data from the null H (independent of) L | C with the exact stratum
    sizes and both (C,H) and (C,L) marginals preserved, so the plug-in
    MM statistic of the permuted data measures precisely the residual
    small-sample bias to subtract (the standard conditional-permutation
    debias).  Note a cyclic-shift null is NOT valid here: it destroys
    the context relevance but keeps the within-byte H-L coupling."""
    c = np.asarray(c, dtype=np.int64)
    s = np.asarray(s, dtype=np.int64)
    rng = np.random.default_rng(seed)
    order = np.argsort(c, kind="stable")
    l_sorted = (s & 15)[order]
    cs = c[order]
    bounds = np.flatnonzero(np.diff(cs)) + 1
    start = 0
    for end in list(bounds) + [order.size]:
        g = end - start
        if g > 1:
            l_sorted[start:end] = l_sorted[start + rng.permutation(g)]
        start = end
    l_null = np.empty_like(l_sorted)
    l_null[order] = l_sorted
    s_null = (s & ~np.int64(15)) | l_null
    return cmi_mm(c, s_null, c_card)["mm"]


def cmi_block_bootstrap(c, s, c_card, seed, n_resamples=1000,
                        chunk=2048, alpha=0.05) -> dict:
    """Percentile CI on the Miller-Madow CMI by CHUNK bootstrap
    (contiguous chunks resampled with replacement; respects short-range
    dependence better than i.i.d. position resampling). Secondary,
    labeled; method recorded in the result."""
    c = np.asarray(c, dtype=np.int64)
    s = np.asarray(s, dtype=np.int64)
    n = s.size
    nch = max(2, n // chunk)
    bounds = [(i * n // nch, (i + 1) * n // nch) for i in range(nch)]
    rng = np.random.default_rng(seed)
    point = cmi_mm(c, s, c_card)["mm"]
    boots = np.empty(n_resamples)
    for b in range(n_resamples):
        take = rng.integers(0, nch, size=nch)
        idx = np.concatenate([np.arange(bounds[t][0], bounds[t][1])
                              for t in take])
        boots[b] = cmi_mm(c[idx], s[idx], c_card)["mm"]
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return {"point": float(point), "lo": float(lo), "hi": float(hi),
            "method": "chunk-percentile", "B": int(n_resamples),
            "chunks": int(nch)}


def cm_argmax_stream(data: bytes, W: int = 3) -> np.ndarray:
    """CM-predictor context instrument: at each position i, the
    predictor's argmax byte BEFORE observing data[i] (deterministic
    integer path; a decoder-reproducible function of the past)."""
    from cm_predictor import CMPredictor
    pred = CMPredictor(W)
    out = np.empty(len(data), dtype=np.int64)
    for i, x in enumerate(data):
        freq, _total = pred.dist()
        out[i] = int(np.argmax(freq))   # ties -> smallest byte (argmax)
        pred.update(x)
    return out


# ---------------------------------------------------------------------------
# Ratio bootstrap for H5 regret (secondary, labeled)
# ---------------------------------------------------------------------------

def ratio_bootstrap(num, den, seed, n_resamples=10000, alpha=0.05) -> dict:
    """CI on sum(num)/sum(den) by resampling paired per-tile values."""
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    n = num.size
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_resamples, n))
    boots = num[idx].sum(axis=1) / den[idx].sum(axis=1)
    point = float(num.sum() / den.sum())
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return {"point": point, "lo": float(lo), "hi": float(hi),
            "method": "paired-tile-percentile", "B": int(n_resamples)}


# ---------------------------------------------------------------------------
# Held-out cross-entropy ladder (H7 reference rate; honest UPPER bound
# on the order-k entropy rate, same construction as the preconditions
# campaign: counts from the first half, add-one evaluation on the
# second half)
# ---------------------------------------------------------------------------

def holdout_ce(x: np.ndarray, order: int) -> float:
    x = np.asarray(x, dtype=np.int64)
    half = x.size // 2
    train, test = x[:half], x[half:]

    def keys(arr, k):
        if k == 0:
            return np.zeros(arr.size, dtype=np.int64), arr
        ctx = np.zeros(arr.size - k, dtype=np.int64)
        for j in range(k):
            ctx = ctx * 256 + arr[j:arr.size - k + j]
        return ctx, arr[k:]

    tctx, tsym = keys(train, order)
    code = tctx * 256 + tsym
    uniq, cnt = np.unique(code, return_counts=True)
    uctx, ccnt = np.unique(tctx, return_counts=True)
    pair = dict(zip(uniq.tolist(), cnt.tolist()))
    ctxc = dict(zip(uctx.tolist(), ccnt.tolist()))
    ectx, esym = keys(test, order)
    bits = 0.0
    for cv, sv in zip(ectx.tolist(), esym.tolist()):
        num = pair.get(cv * 256 + sv, 0) + 1
        den = ctxc.get(cv, 0) + 256
        bits += -math.log2(num / den)
    return bits / esym.size


def entropy_rate_estimate(x: np.ndarray, orders=(0, 1, 2, 3)) -> dict:
    vals = {k: holdout_ce(x, k) for k in orders}
    best = min(vals, key=vals.get)
    return {"per_order": vals, "h_hat": vals[best], "order": best}
