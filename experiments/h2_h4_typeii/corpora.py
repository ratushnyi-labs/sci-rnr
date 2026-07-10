#!/usr/bin/env python3
"""Record/field decompositions for the H2/H3/H4 Type-II separation experiment.

Every corpus is reduced to a matrix of M fixed-length records of R bytes
(fields = byte positions within the record), the exact setting of
Theorem 5.5: the block is the record, the per-position marginals are the
field marginals, and TC(D_N) = sum_f H(field_f) - H(record).

Two real corpora (from data/MANIFEST.json via data/loader.py):

  sqlite_hdr      first 16 bytes of every 4096-byte page of sqlite_synth
                  (page-header region: page type, freeblock ptr, cell
                  count, content offset, fragment count, first cell
                  pointers -- the positional structure Type-II targets).
  telemetry_expo  bytes (b3, b2) of each little-endian float32 of
                  telemetry_grid (sign+exponent-high, exponent-low+
                  mantissa-high): the cross-field dependence is the
                  exponent structure; mantissa noise bytes are excluded
                  so the plug-in joint entropy is estimable.

Three seeded synthetic corpora (generators documented inline; every
random draw uses a seed derived by bench.metrics.derive_seed with
tag="corpus", so the corpora are reproducible from MANIFEST master seed):

  record_log      fixed-record binary log, R=8 chained fields
                  (src -> level -> code -> len -> flag chain + independent
                  dt + 2 payload fields). All field dependencies are at
                  byte distance <= 3, i.e. learnable by an order-3
                  context coder. Exact model TC computable.
  columnar_int16  columnar int16 sensor array, 8 channels interleaved
                  per sample (R=16). Channels share a latent level m;
                  even channels add {0,1} jitter. High bytes are zero
                  (small sensor values). Exact model TC computable.
  calib_dup       calibration anchor: records (a, a, u, v) with a,u,v
                  independent uniform on 16 values. TC = H(a) = 4 bits
                  per 4-byte record = 1.0 bpb, in closed form.

All generators return records as an (M, R) uint8 array.
"""

from __future__ import annotations

import hashlib
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "data"))
sys.path.insert(0, str(REPO / "bench"))

import loader  # noqa: E402
import metrics  # noqa: E402

EXPERIMENT_ID = "h2h4_typeii"


def _rng(corpus: str, index: int = 0) -> np.random.Generator:
    """Seeded generator per METHODS.md section 4 (tag='corpus')."""
    return np.random.default_rng(
        metrics.derive_seed(f"{EXPERIMENT_ID}:{corpus}", index, tag="corpus")
    )


def records_sha256(records: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(records).tobytes()).hexdigest()


# ---------------------------------------------------------------------------
# Real corpora
# ---------------------------------------------------------------------------

def sqlite_hdr(M: int | None = None) -> dict:
    """First 16 bytes of each 4096-byte page of sqlite_synth."""
    page = 4096
    R = 16
    n_pages = loader.corpus_bytes("sqlite_synth") // page
    take = n_pages if M is None else min(M, n_pages)
    rows = np.empty((take, R), dtype=np.uint8)
    for j, block in enumerate(loader.iter_blocks("sqlite_synth", page,
                                                 drop_last=True)):
        if j >= take:
            break
        rows[j] = np.frombuffer(block[:R], dtype=np.uint8)
    return {
        "name": "sqlite_hdr", "records": rows, "R": R,
        "source": "sqlite_synth",
        "corpus_sha256": loader.entry("sqlite_synth")["sha256"],
        "analytic_tc_bpb": None,
        "field_names": [f"hdr_byte_{i}" for i in range(R)],
    }


def telemetry_expo(M: int = 262144) -> dict:
    """(b3, b2) little-endian bytes of a seeded float32 subsample."""
    raw = np.fromfile(loader.payload_path("telemetry_grid"), dtype=np.uint8)
    quads = raw.reshape(-1, 4)  # little-endian float32 bytes b0..b3
    rng = _rng("telemetry_expo")
    idx = rng.choice(quads.shape[0], size=min(M, quads.shape[0]),
                     replace=False)
    idx.sort()
    rows = quads[idx][:, [3, 2]].copy()  # sign+exp-high, exp-low+mant-high
    return {
        "name": "telemetry_expo", "records": rows, "R": 2,
        "source": "telemetry_grid",
        "corpus_sha256": loader.entry("telemetry_grid")["sha256"],
        "analytic_tc_bpb": None,
        "field_names": ["f32_byte3", "f32_byte2"],
    }


# ---------------------------------------------------------------------------
# Synthetic corpora (seeded, documented generators)
# ---------------------------------------------------------------------------

def _entropy_bits(p: np.ndarray) -> float:
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def _cond_entropy(joint: np.ndarray) -> float:
    """H(col | row) from a joint pmf table."""
    h = 0.0
    for row in joint:
        s = row.sum()
        if s > 0:
            h += s * _entropy_bits(row / s)
    return float(h)


def record_log(M: int = 65536) -> dict:
    """Fixed-record binary log with a learnable dependence chain.

    Fields (one byte each, R=8), with all alphabets small so that the
    order-3 context tables of the reference coder converge:

      f0 src   ~ Zipf-like over 32 sources, p(k) proportional to 1/(k+1)
      f1 level = T1[src] w.p. 0.9, else uniform over 8      (dist 1)
      f2 code  = T2[level] w.p. 0.85, else uniform over 32  (dist 1)
      f3 dt    ~ geometric-like over 8, p(k) prop. 0.6^k    (independent)
      f4 len   = T4[code] w.p. 0.8, else uniform over 16    (dist 2)
      f5 flag  = len & 3 w.p. 0.95, else uniform over 4     (dist 1)
      f6,f7 payload ~ uniform over 4 each                   (independent)

    T1/T2/T4 are seeded random lookup tables. The exact model TC is
    computed by propagating the chain marginals (records are iid draws).
    """
    rng = _rng("record_log")
    n_src, n_lvl, n_code, n_dt, n_len, n_flag, n_pay = 32, 8, 32, 8, 16, 4, 4
    T1 = rng.integers(0, n_lvl, size=n_src)
    T2 = rng.integers(0, n_code, size=n_lvl)
    T4 = rng.integers(0, n_len, size=n_code)

    p_src = 1.0 / (np.arange(n_src) + 1.0)
    p_src /= p_src.sum()
    p_dt = 0.6 ** np.arange(n_dt)
    p_dt /= p_dt.sum()

    src = rng.choice(n_src, size=M, p=p_src)
    lvl = np.where(rng.random(M) < 0.9, T1[src], rng.integers(0, n_lvl, M))
    code = np.where(rng.random(M) < 0.85, T2[lvl], rng.integers(0, n_code, M))
    dt = rng.choice(n_dt, size=M, p=p_dt)
    ln = np.where(rng.random(M) < 0.8, T4[code], rng.integers(0, n_len, M))
    flag = np.where(rng.random(M) < 0.95, ln & 3, rng.integers(0, n_flag, M))
    pay0 = rng.integers(0, n_pay, M)
    pay1 = rng.integers(0, n_pay, M)
    rows = np.stack([src, lvl, code, dt, ln, flag, pay0, pay1],
                    axis=1).astype(np.uint8)

    # Exact model TC: chain factorization  H(joint) = H(f0) + H(f1|f0)
    # + H(f2|f1) + H(f3) + H(f4|f2) + H(f5|f4) + H(f6) + H(f7).
    def cond_table(par_p, n_par, n_ch, table, p_true):
        j = np.zeros((n_par, n_ch))
        for a in range(n_par):
            j[a, :] += par_p[a] * (1 - p_true) / n_ch
            j[a, table[a]] += par_p[a] * p_true
        return j

    j01 = cond_table(p_src, n_src, n_lvl, T1, 0.9)
    p_lvl = j01.sum(axis=0)
    j12 = cond_table(p_lvl, n_lvl, n_code, T2, 0.85)
    p_code = j12.sum(axis=0)
    j24 = cond_table(p_code, n_code, n_len, T4, 0.8)
    p_len = j24.sum(axis=0)
    flag_tab = np.arange(n_len) & 3
    j45 = cond_table(p_len, n_len, n_flag, flag_tab, 0.95)
    p_flag = j45.sum(axis=0)

    h_joint = (_entropy_bits(p_src) + _cond_entropy(j01) + _cond_entropy(j12)
               + _entropy_bits(p_dt) + _cond_entropy(j24) + _cond_entropy(j45)
               + 2 * math.log2(n_pay))
    h_marg = (_entropy_bits(p_src) + _entropy_bits(p_lvl)
              + _entropy_bits(p_code) + _entropy_bits(p_dt)
              + _entropy_bits(p_len) + _entropy_bits(p_flag)
              + 2 * math.log2(n_pay))
    tc_bpb = (h_marg - h_joint) / 8.0
    return {
        "name": "record_log", "records": rows, "R": 8, "source": "generated",
        "corpus_sha256": records_sha256(rows),
        "analytic_tc_bpb": tc_bpb,
        "field_names": ["src", "level", "code", "dt", "len", "flag",
                        "pay0", "pay1"],
    }


def columnar_int16(M: int = 32768) -> dict:
    """Columnar int16 sensor array, sample-major (interleaved) layout.

    8 channels per time sample; channel c reads v_c = m + e_c where the
    latent level m is uniform on {0..63} per sample (iid records) and
    e_c is a {0,1} fair-coin jitter on even channels, 0 on odd channels.
    Values are stored little-endian int16, so high bytes are zero
    (realistic for small-range sensors). R = 16 bytes per record.

    Exact model TC per record: sum_c H(v_c) - H(joint), with
    H(joint) = H(m) + 4 bits (the four jitter bits) and H(v_c) either
    log2 64 (odd channels) or H(uniform64 conv Bernoulli(1/2)) (even).
    """
    rng = _rng("columnar_int16")
    C = 8
    m = rng.integers(0, 64, size=M)
    vals = np.empty((M, C), dtype=np.int64)
    for c in range(C):
        e = rng.integers(0, 2, size=M) if c % 2 == 0 else 0
        vals[:, c] = m + e
    rows = np.zeros((M, 2 * C), dtype=np.uint8)
    rows[:, 0::2] = vals & 0xFF          # low bytes
    rows[:, 1::2] = (vals >> 8) & 0xFF   # high bytes (all zero here)

    p_m = np.full(64, 1 / 64)
    p_me = np.zeros(65)
    p_me[:64] += p_m * 0.5
    p_me[1:65] += p_m * 0.5
    h_odd = math.log2(64)
    h_even = _entropy_bits(p_me)
    h_joint = math.log2(64) + 4.0
    tc_bpb = (4 * h_odd + 4 * h_even - h_joint) / 16.0
    return {
        "name": "columnar_int16", "records": rows, "R": 16,
        "source": "generated",
        "corpus_sha256": records_sha256(rows),
        "analytic_tc_bpb": tc_bpb,
        "field_names": [f"ch{c}_{p}" for c in range(C) for p in ("lo", "hi")],
    }


def calib_dup(M: int = 65536) -> dict:
    """Calibration anchor: (a, a, u, v), a,u,v iid uniform on 16 values.

    TC = H(a) = 4 bits per 4-byte record = 1.0 bpb in closed form; the
    duplication is at byte distance 1, trivially learnable at W >= 1.
    """
    rng = _rng("calib_dup")
    a = rng.integers(0, 16, size=M, dtype=np.uint8)
    u = rng.integers(0, 16, size=M, dtype=np.uint8)
    v = rng.integers(0, 16, size=M, dtype=np.uint8)
    rows = np.stack([a, a, u, v], axis=1)
    return {
        "name": "calib_dup", "records": rows, "R": 4, "source": "generated",
        "corpus_sha256": records_sha256(rows),
        "analytic_tc_bpb": 1.0,
        "field_names": ["a", "a_dup", "u", "v"],
    }


# ---------------------------------------------------------------------------
# Variants: seeded record shuffle (iid framing) and the H4 permutation
# control (destroys cross-field structure, preserves marginals exactly)
# ---------------------------------------------------------------------------

def shuffle_records(records: np.ndarray, corpus: str,
                    seed_index: int) -> np.ndarray:
    """Seeded record shuffle: records become exchangeable, so the
    empirical record distribution is the D_N of Theorem 5.5."""
    rng = np.random.default_rng(
        metrics.derive_seed(f"{EXPERIMENT_ID}:{corpus}:shuffle",
                            seed_index, tag="run"))
    perm = rng.permutation(records.shape[0])
    return records[perm]


def permute_fields(records: np.ndarray, corpus: str,
                   seed_index: int) -> np.ndarray:
    """H4 control: independently permute each field column across records.

    Marginals are preserved exactly; the joint distribution becomes the
    product of marginals, so TC -> 0 (up to estimator bias) and the
    factorized-vs-joint gap must vanish (Remark 5.6 control direction).
    """
    rng = np.random.default_rng(
        metrics.derive_seed(f"{EXPERIMENT_ID}:{corpus}:permute",
                            seed_index, tag="run"))
    out = np.empty_like(records)
    M = records.shape[0]
    for f in range(records.shape[1]):
        out[:, f] = records[rng.permutation(M), f]
    return out


BUILDERS = {
    "sqlite_hdr": sqlite_hdr,
    "telemetry_expo": telemetry_expo,
    "record_log": record_log,
    "columnar_int16": columnar_int16,
    "calib_dup": calib_dup,
}

FULL_M = {
    "sqlite_hdr": None,        # all 3016 pages
    "telemetry_expo": 262144,
    "record_log": 65536,
    "columnar_int16": 32768,
    "calib_dup": 65536,
}

REDUCED_M = {
    "sqlite_hdr": 1024,
    # telemetry's cross-field TC is only ~0.07 bpb; below 2^16 records
    # the structured-vs-permuted CI separation sits at the resolution
    # limit (observed overlap 0.001 bpb at M=32768), so the reduced
    # scale keeps a quarter of the full sample rather than an eighth.
    "telemetry_expo": 65536,
    "record_log": 8192,
    "columnar_int16": 4096,
    "calib_dup": 8192,
}
