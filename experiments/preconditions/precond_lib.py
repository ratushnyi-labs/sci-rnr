#!/usr/bin/env python3
"""Measurement library for the corpus-precondition statistics campaign.

Instruments the qualifying data-class precondition of the main paper
(rnr_core.tex section 10.2.1, inequality (10.1)) on the local corpus
suite, using the reference Type-I coder (impl/rnr1.py) as the
side-information instrument:

  (a) conditional-advantage ratio  E[L(Enc(X))] / E[L(B(X))]
      (Definition 2.1 as a ratio), with the (10.1) term breakdown
      (a) H(X|Y) proxy = ideal predictor cross-entropy,
      (b) L(M)/mN      = 0 by construction for the counting predictor
                         (only the 8-byte model-description hash is
                         charged; the neural L(M) term is NOT exercised),
      (c) rho          = realized repair rate minus the ideal rate,
      (d) nu           = verification bytes (per-sub-block hashes + v),
      (e) mu           = header/index metadata;
  (b) stationarity / mixing indicators (per-third entropy ladder,
      byte-stream autocorrelation at structural lags);
  (c) per-corpus entropy-rate ladder (plug-in and held-out
      order-0..4 conditional empirical entropy) and the Type-III-B
      nibble conditional mutual information I(H;L|C).

Floats are used throughout: this is offline analysis tooling in the
sense of engineering-spec R-3.6, not the coding path.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for sub in ("bench", "data", "impl"):
    p = str(REPO / sub)
    if p not in sys.path:
        sys.path.insert(0, p)

import baselines  # noqa: E402  (bench/)
import loader  # noqa: E402  (data/)
import metrics  # noqa: E402  (bench/)
import rnr1  # noqa: E402  (impl/)

K_SYNC = 65536  # typical sync spacing (exp-design H10, BUG-003-D "K = 64 KiB")

# ---------------------------------------------------------------------------
# Reference-coder instrumentation (term-resolved Type-I measurement)
# ---------------------------------------------------------------------------


def encode_block_stats(block: bytes, W: int):
    """Encode one sub-block exactly as rnr1.encode_subblock does, while
    accumulating the ideal (pre-quantization) code length.

    Returns (blob, ideal_bits): blob is byte-identical to the reference
    coder's output (same code path, same state machine); ideal_bits is
    sum_i -log2(freq_i[x_i]/total_i), the exact chain-rule cost of the
    class+payload pair before arithmetic-coder quantization, i.e. the
    H(X|Y) proxy realized by this predictor as the side-information
    instrument.
    """
    pred = rnr1.NGramPredictor(W)
    enc = rnr1.ArithmeticEncoder()
    ideal_bits = 0.0
    log2 = math.log2
    for x in block:
        freq, total = pred.dist()
        ideal_bits += log2(total / int(freq[x]))
        bstar, cand, classid, class_freq = rnr1._partition(freq)
        cls = int(classid[x])
        ccum = np.concatenate(([0], np.cumsum(class_freq)))
        enc.encode(int(ccum[cls]), int(ccum[cls + 1]), total)
        if cls != rnr1.CLS_0:
            members = rnr1._class_members(cls, cand, classid)
            mcum = np.concatenate(([0], np.cumsum(freq[members])))
            idx = int(np.nonzero(members == x)[0][0])
            enc.encode(int(mcum[idx]), int(mcum[idx + 1]), int(mcum[-1]))
        pred.update(x)
    return enc.finish(), ideal_bits


def rnr_archive_measure(data: bytes, W: int, K: int = K_SYNC) -> dict:
    """Term-resolved measurement of the RNR1 archive for `data`.

    The archive size is computed analytically from the same quantities
    rnr1.pack assembles (header + index entries + effective blobs); the
    identity len(rnr1.pack(data,W,K)) == result["archive_bytes"] is
    asserted in run_checks.py on a probe input.

    Term accounting (inequality (10.1), per symbol):
      a_ideal_bpb  : ideal predictor cross-entropy (H(X|Y) proxy),
                     coded blocks only, over coded bytes;
      c_rho_bpb    : realized repair rate minus ideal rate (coder
                     residual; arithmetic-coder quantization + flush),
                     coded blocks only;
      d_nu_bpb     : verification bytes = 8 per sub-block hash + 8 (v);
      e_mu_bpb     : header/index metadata minus (d) and minus the
                     8-byte model hash;
      b_model_bpb  : the 8-byte model-description hash (the L(M) charged
                     by this instrument; the counting predictor has no
                     trained parameters).
    """
    n = len(data)
    m = (n + K - 1) // K
    per_block = []
    repair_bytes = 0
    ideal_bits_coded = 0.0
    coded_src_bytes = 0
    raw_blocks = 0
    for j in range(m):
        block = data[j * K: (j + 1) * K]
        blob, ideal_bits = encode_block_stats(block, W)
        mode_raw = len(blob) >= len(block)
        eff = len(block) if mode_raw else len(blob)
        if mode_raw:
            raw_blocks += 1
        else:
            ideal_bits_coded += ideal_bits
            coded_src_bytes += len(block)
        repair_bytes += eff
        per_block.append(
            {
                "index": j,
                "src_bytes": len(block),
                "eff_bytes": eff,
                "mode": "raw" if mode_raw else "coded",
                "ideal_bits": ideal_bits,
                "bpb": 8.0 * eff / len(block),
            }
        )
    archive_bytes = rnr1.HEADER_SIZE + m * rnr1.INDEX_ENTRY_SIZE + repair_bytes
    nu_bytes = 8 * m + 8          # per-sub-block hashes + archive v (sec 10.3)
    b_bytes = 8                   # model-description hash (R-13.3)
    mu_bytes = archive_bytes - repair_bytes - nu_bytes - b_bytes
    out = {
        "W": W,
        "K": K,
        "n": n,
        "m": m,
        "raw_blocks": raw_blocks,
        "coded_fraction": (m - raw_blocks) / m,
        "archive_bytes": archive_bytes,
        "repair_bytes": repair_bytes,
        "total_bpb": 8.0 * archive_bytes / n,
        "repair_bpb": 8.0 * repair_bytes / n,
        "a_ideal_bpb": (ideal_bits_coded / coded_src_bytes)
        if coded_src_bytes
        else None,
        "c_rho_bpb": (
            (8.0 * sum(b["eff_bytes"] for b in per_block if b["mode"] == "coded")
             - ideal_bits_coded) / coded_src_bytes
        )
        if coded_src_bytes
        else None,
        "b_model_bpb": 8.0 * b_bytes / n,
        "d_nu_bpb": 8.0 * nu_bytes / n,
        "e_mu_bpb": 8.0 * mu_bytes / n,
        "per_block": per_block,
    }
    return out


def baseline_block_bpbs(blocks: list, coder: str) -> list:
    """Per-block bits-per-byte of a baseline coder on the SAME blocks.

    This is the matched-granularity comparison: like the RNR1 archive,
    the baseline restarts cold at every K-block boundary, so it is the
    in-process proxy for a seekable container (bgzip / zstd-seekable) at
    sync spacing K.
    """
    return [
        baselines.run_baseline(coder, blk).bits_per_byte for blk in blocks
    ]


def baseline_stream_bpb(data: bytes, coder: str) -> float:
    """Whole-slice (streaming) baseline rate: the harder comparison."""
    return baselines.run_baseline(coder, data).bits_per_byte


# ---------------------------------------------------------------------------
# Entropy / mixing indicators (numpy, plug-in and held-out)
# ---------------------------------------------------------------------------


def _entropy_from_counts(counts: np.ndarray, n: int) -> float:
    c = counts[counts > 0].astype(np.float64)
    p = c / n
    return float(-(p * np.log2(p)).sum())


def _pack_keys(x: np.ndarray, order: int, with_symbol: bool) -> np.ndarray:
    """Pack length-`order` contexts (optionally + next symbol) into uint64."""
    n = x.size
    m = n - order
    if m <= 0:
        raise ValueError("slice shorter than context order")
    keys = np.zeros(m, dtype=np.uint64)
    for j in range(order):
        keys = (keys << np.uint64(8)) | x[j: j + m].astype(np.uint64)
    if with_symbol:
        keys = (keys << np.uint64(8)) | x[order: order + m].astype(np.uint64)
    return keys


def plugin_cond_entropy(x: np.ndarray, order: int) -> dict:
    """Plug-in conditional empirical entropy H(X_i | X_{i-order..i-1}).

    Returns {"h": bits/byte, "contexts": distinct contexts,
    "positions": sample count}. Plug-in estimates are biased LOW when
    contexts are sparse; the coverage fields make that auditable.
    """
    n = x.size
    if order == 0:
        cnt = np.bincount(x, minlength=256)
        return {"h": _entropy_from_counts(cnt, n), "contexts": 1, "positions": n}
    joint = _pack_keys(x, order, with_symbol=True)
    ctx = _pack_keys(x, order, with_symbol=False)
    _, cj = np.unique(joint, return_counts=True)
    uc, cc = np.unique(ctx, return_counts=True)
    m = joint.size
    h = _entropy_from_counts(cj, m) - _entropy_from_counts(cc, m)
    return {"h": h, "contexts": int(uc.size), "positions": int(m)}


def _lookup_counts(uniq: np.ndarray, cnt: np.ndarray, q: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(uniq, q)
    idx_c = np.clip(idx, 0, uniq.size - 1)
    hit = uniq[idx_c] == q
    return np.where(hit, cnt[idx_c], 0)


def holdout_cross_entropy(x: np.ndarray, order: int) -> dict:
    """Held-out conditional cross-entropy: counts from the first half,
    code length evaluated on the second half with Laplace smoothing
    p = (c(ctx,x)+1)/(c(ctx)+256).

    This is a REAL sequential code length, hence an honest upper bound
    on the order-`order` conditional entropy rate (no sparsity
    under-bias, unlike the plug-in estimate). coverage = fraction of
    evaluated positions whose context was seen in training.
    """
    half = x.size // 2
    train, test = x[:half], x[half:]
    if order == 0:
        cnt = np.bincount(train, minlength=256).astype(np.float64)
        p = (cnt[test] + 1.0) / (half + 256.0)
        return {"ce": float(-np.log2(p).mean()), "coverage": 1.0}
    tj = _pack_keys(train, order, with_symbol=True)
    tc = _pack_keys(train, order, with_symbol=False)
    uj, cj = np.unique(tj, return_counts=True)
    uc, cc = np.unique(tc, return_counts=True)
    qj = _pack_keys(test, order, with_symbol=True)
    qc = _pack_keys(test, order, with_symbol=False)
    c_joint = _lookup_counts(uj, cj, qj).astype(np.float64)
    c_ctx = _lookup_counts(uc, cc, qc).astype(np.float64)
    p = (c_joint + 1.0) / (c_ctx + 256.0)
    return {
        "ce": float(-np.log2(p).mean()),
        "coverage": float((c_ctx > 0).mean()),
    }


def thirds_ladder(x: np.ndarray, orders=(0, 1, 2, 3)) -> list:
    """Plug-in ladder per file third: stationarity instrument."""
    n = x.size
    cut = [0, n // 3, 2 * n // 3, n]
    out = []
    for t in range(3):
        seg = x[cut[t]: cut[t + 1]]
        out.append({f"h{k}": plugin_cond_entropy(seg, k)["h"] for k in orders})
    return out


def rel_spread(vals) -> float:
    v = np.asarray(vals, dtype=float)
    mu = v.mean()
    if mu == 0.0:
        return 0.0
    return float((v.max() - v.min()) / mu)


ACF_LAGS = (1, 2, 4, 8, 16, 32, 64, 77, 256, 1024, 4096, 8192)


def byte_acf(x: np.ndarray, lags=ACF_LAGS) -> dict:
    """Autocorrelation of the byte stream (bytes as numeric values).

    Lag 77 targets the base64 line stride (76+newline), 4096 the SQLite
    page size, 8192 the telemetry row stride (2048 float32 = 8192 B).
    """
    xf = x.astype(np.float64)
    mu = xf.mean()
    var = xf.var()
    out = {}
    for lag in lags:
        if lag >= x.size:
            out[lag] = None
            continue
        if var == 0.0:
            out[lag] = 0.0
            continue
        out[lag] = float(((xf[:-lag] - mu) * (xf[lag:] - mu)).mean() / var)
    return out


def nibble_cmi(x: np.ndarray) -> float:
    """Type-III-B instrument: I(H ; L | C) with C = previous byte,
    H/L = high/low nibble of the current byte (exp-design section 5.3,
    plug-in over the empirical order-1 joint law).

    Miller-Madow bias correction with observed supports: the plug-in
    CMI is biased UP by ~(K_CHL - K_CH - K_CL + K_C)/(2 n ln 2); on an
    i.i.d.-uniform 2 MiB control the raw plug-in reads ~0.020 bits
    while the corrected value is ~0 (validated in run_checks.py).
    """
    c = x[:-1].astype(np.int64)
    s = x[1:].astype(np.int64)
    n = s.size
    cnt_c = np.bincount(c, minlength=256)
    cnt_j = np.bincount(c * 256 + s, minlength=65536)
    cnt_ch = np.bincount(c * 16 + (s >> 4), minlength=4096)
    cnt_cl = np.bincount(c * 16 + (s & 15), minlength=4096)
    hc = _entropy_from_counts(cnt_c, n)
    hx_c = _entropy_from_counts(cnt_j, n) - hc
    hh_c = _entropy_from_counts(cnt_ch, n) - hc
    hl_c = _entropy_from_counts(cnt_cl, n) - hc
    k_c = int((cnt_c > 0).sum())
    k_j = int((cnt_j > 0).sum())
    k_ch = int((cnt_ch > 0).sum())
    k_cl = int((cnt_cl > 0).sum())
    bias = (k_j - k_ch - k_cl + k_c) / (2.0 * n * math.log(2.0))
    return float(hh_c + hl_c - hx_c - bias)


# ---------------------------------------------------------------------------
# Statistics (frozen methodology, bench/METHODS.md)
# ---------------------------------------------------------------------------


def bootstrap_p_two_sided(diffs, seed: int, n_resamples: int = metrics.N_RESAMPLES) -> float:
    """Two-sided bootstrap achieved significance for mean(diffs) != 0,
    with the (r+1)/(B+1) finite-resample correction. Labeled SECONDARY:
    the precondition screen is not one of the frozen H1-H18 primaries
    (METHODS.md section 1 allows labeled secondary procedures)."""
    d = np.asarray(diffs, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, d.size, size=(n_resamples, d.size))
    boots = d[idx].mean(axis=1)
    p_lo = (np.sum(boots <= 0.0) + 1.0) / (n_resamples + 1.0)
    p_hi = (np.sum(boots >= 0.0) + 1.0) / (n_resamples + 1.0)
    return float(min(1.0, 2.0 * min(p_lo, p_hi)))


def slice_bytes_of(name: str, want: int) -> bytes:
    """Deterministic dev slice: the first min(want, size) bytes."""
    n = min(want, loader.corpus_bytes(name))
    return loader.read_range(name, 0, n)


def corpus_sha(name: str) -> str:
    return loader.entry(name)["sha256"]
