#!/usr/bin/env python3
"""H7 -- minimal Type-III-C dictionary coder with greedy marginal-MDL
admission and bits-back accounting (exp-design section 2.6; main paper
Theorems 6.5/6.6, Lemma 6.6a, Theorem 6.6b's o(N) overhead shape).

Coder (real, byte-exact round-trip):
  * Candidate phrases: substrings of length 2..7 counted on the slice.
  * Admission: Theorem 6.5's frozen-rest greedy rule.  For candidate
    tau of length l with realized non-overlapping count f on positions
    still parsed as raw literals,
        Delta_L(tau) = dictionary_cost + f * (code_cost - 8 l),
    dictionary_cost = 8 (l + 1) bits (length byte + raw bytes in the
    header), parse_overhead = 0 (boundaries are implied by token
    lengths in a single entropy-coded symbol stream -- declared, not
    hidden), and code_cost = -log2(f / S') with S' the projected
    symbol-stream length after admission ("all else frozen").  Admit
    iff Delta_L < 0.
  * Format: header (magic, N, K, K entries of [len byte | raw bytes])
    + one adaptive arithmetic-coded stream (impl/rnr1.py coder) over
    the alphabet {256 literals} + {K token ids}, greedy longest-match
    parse.  The decoder mirrors the adaptive counts exactly.

Measurements per ladder point N:
  * rate vs N*h_hat (h_hat = held-out cross-entropy min over orders
    0..3, an honest upper bound on the order-k entropy rate);
  * overhead fraction = dictionary header bits / total bits (o(N)
    check: must decrease along the N-ladder, Theorem 6.6b shape);
  * Lemma 6.6a rate-neutral index check on the token-id substream:
    ideal two-stage (family/variant grid) code length must equal the
    ideal flat length EXACTLY (chain rule), and the REALIZED adaptive
    two-stage bits must match the realized flat bits and the empirical
    entropy within the declared adaptive-redundancy tolerance;
  * bits-back accounting (Theorem 6.6): (a) exact identity (6.4) on a
    micro model where the parse posterior is enumerable; (b) on the
    real coder, the deterministic parse is the q = point-mass case, so
    net = -log2 P(X,z*) = -log2 P_X(X) + [-log2 P(z*|X)]; a DP over
    all parses under the frozen final token law computes the marginal
    and hence the KL term exactly.
"""

from __future__ import annotations

import math
from collections import Counter

import numpy as np

import common
from common import loader, rnr1

MAGIC = b"D3C1"
MAX_PHRASE = 7
MIN_COUNT = 4
CAND_CAP = 3000
K_CAP = 1024


# ---------------------------------------------------------------------------
# Greedy marginal-MDL dictionary builder (Theorem 6.5 rule)
# ---------------------------------------------------------------------------

def _candidates(data: bytes):
    cnt = Counter()
    n = len(data)
    for ell in range(2, MAX_PHRASE + 1):
        for i in range(n - ell + 1):
            cnt[data[i:i + ell]] += 1
    cand = [(t, c) for t, c in cnt.items() if c >= MIN_COUNT]
    cand.sort(key=lambda tc: (-tc[1] * (len(tc[0]) - 1), tc[0]))
    return cand[:CAND_CAP]


def _realized_count(data: bytes, tau: bytes, covered: np.ndarray):
    """Non-overlapping occurrences of tau entirely on uncovered
    positions (frozen-parser realized count f of Theorem 6.5)."""
    ell = len(tau)
    hits = []
    i = data.find(tau)
    while i != -1:
        if not covered[i:i + ell].any():
            hits.append(i)
            i = data.find(tau, i + ell)          # non-overlapping
        else:
            i = data.find(tau, i + 1)
    return hits


def build_dictionary(data: bytes):
    """Greedy admission; returns (dict entries, admission log)."""
    n = len(data)
    covered = np.zeros(n, dtype=bool)
    total_syms = float(n)                        # all-literal start
    admitted, log = [], []
    for tau, _raw_cnt in _candidates(data):
        if len(admitted) >= K_CAP:
            break
        hits = _realized_count(data, tau, covered)
        f = len(hits)
        if f < 2:
            continue
        ell = len(tau)
        s_after = total_syms - f * ell + f
        code_cost = -math.log2(f / s_after)
        dict_cost = 8.0 * (ell + 1)
        delta = dict_cost + f * (code_cost - 8.0 * ell)
        if delta < 0:
            admitted.append(tau)
            log.append({"tau_len": ell, "f": f, "delta_bits": delta,
                        "code_cost": code_cost})
            for h in hits:
                covered[h:h + ell] = True
            total_syms = s_after
    return admitted, log


# ---------------------------------------------------------------------------
# Container: header + one adaptive arithmetic-coded token/literal stream
# ---------------------------------------------------------------------------

def _parse(data: bytes, dictionary):
    """Greedy longest-match parse into symbols (0..255 literals,
    256+k tokens)."""
    by_first = {}
    for k, tau in enumerate(dictionary):
        by_first.setdefault(tau[0], []).append((len(tau), k, tau))
    for lst in by_first.values():
        lst.sort(reverse=True)                   # longest first
    syms = []
    i, n = 0, len(data)
    while i < n:
        hit = None
        for ell, k, tau in by_first.get(data[i], ()):
            if data[i:i + ell] == tau:
                hit = (ell, k)
                break
        if hit:
            syms.append(256 + hit[1])
            i += hit[0]
        else:
            syms.append(data[i])
            i += 1
    return syms


def _code_adaptive(symbols, alphabet, enc=None, dec=None, nsym=None):
    """Adaptive (add-one) arithmetic coding over `alphabet` symbols.
    Encoder mode: symbols given, returns (blob, ideal_bits_adaptive).
    Decoder mode: dec + nsym given, returns list of symbols."""
    counts = np.ones(alphabet, dtype=np.int64)
    if dec is None:
        enc = rnr1.ArithmeticEncoder()
        bits = 0.0
        for s in symbols:
            cum = np.concatenate(([0], np.cumsum(counts)))
            tot = int(cum[-1])
            bits += -math.log2(counts[s] / tot)
            enc.encode(int(cum[s]), int(cum[s + 1]), tot)
            counts[s] += 1
        return enc.finish(), bits
    out = []
    for _ in range(nsym):
        cum = np.concatenate(([0], np.cumsum(counts)))
        tot = int(cum[-1])
        t = dec.target(tot)
        s = int(np.searchsorted(cum, t, side="right")) - 1
        dec.consume(int(cum[s]), int(cum[s + 1]), tot)
        out.append(s)
        counts[s] += 1
    return out


def encode(data: bytes):
    dictionary, adm_log = build_dictionary(data)
    syms = _parse(data, dictionary)
    K = len(dictionary)
    header = bytearray()
    header += MAGIC
    header += len(data).to_bytes(8, "little")
    header += K.to_bytes(2, "little")
    header += len(syms).to_bytes(8, "little")
    for tau in dictionary:
        header.append(len(tau))
        header += tau
    blob, ideal_bits = _code_adaptive(syms, 256 + K)
    archive = bytes(header) + blob
    stats = {
        "n": len(data), "K": K, "n_syms": len(syms),
        "header_bits": 8 * len(header), "dict_bits": 8 * len(header) - 8 * 22,
        "stream_bits": 8 * len(blob), "total_bits": 8 * len(archive),
        "ideal_stream_bits": ideal_bits,
        "admission_log": adm_log,
        "symbols": syms, "dictionary": dictionary,
    }
    return archive, stats


def decode(archive: bytes) -> bytes:
    assert archive[:4] == MAGIC
    n = int.from_bytes(archive[4:12], "little")
    K = int.from_bytes(archive[12:14], "little")
    nsym = int.from_bytes(archive[14:22], "little")
    off = 22
    dictionary = []
    for _ in range(K):
        ell = archive[off]
        dictionary.append(archive[off + 1:off + 1 + ell])
        off += 1 + ell
    dec = rnr1.ArithmeticDecoder(archive[off:])
    syms = _code_adaptive(None, 256 + K, dec=dec, nsym=nsym)
    out = bytearray()
    for s in syms:
        if s < 256:
            out.append(s)
        else:
            out += dictionary[s - 256]
    assert len(out) == n
    return bytes(out)


# ---------------------------------------------------------------------------
# Lemma 6.6a: rate-neutral hierarchical index
# ---------------------------------------------------------------------------

def index_neutrality(symbols, K: int) -> dict:
    """Token-id substream: ideal flat vs ideal two-stage (must be equal
    by the chain rule), realized adaptive flat vs realized adaptive
    two-stage vs empirical entropy."""
    ids = np.array([s - 256 for s in symbols if s >= 256], dtype=np.int64)
    if ids.size == 0 or K == 0:
        return {"n_ids": 0}
    b = int(math.ceil(math.sqrt(K)))
    a = int(math.ceil(K / b))
    fam, var = ids // b, ids % b
    n = ids.size
    cnt = np.bincount(ids, minlength=a * b).astype(float)
    p = cnt / n
    cf = cnt.reshape(a, b).sum(axis=1)
    pf = cf / n
    pv_f = cnt.reshape(a, b) / np.maximum(cf[:, None], 1)
    ideal_flat = float(-(np.log2(p[ids])).sum())
    ideal_two = float(-(np.log2(pf[fam])).sum()
                      - (np.log2(pv_f[fam, var])).sum())
    _, real_flat = _code_adaptive(ids.tolist(), K)
    blob_f, bits_f = _code_adaptive(fam.tolist(), a)
    # per-family adaptive variant coding
    counts = np.ones((a, b), dtype=np.int64)
    bits_v = 0.0
    for fa, va in zip(fam.tolist(), var.tolist()):
        tot = int(counts[fa].sum())
        bits_v += -math.log2(counts[fa, va] / tot)
        counts[fa, va] += 1
    real_two = bits_f + bits_v
    h_emp = float(-(p[p > 0] * np.log2(p[p > 0])).sum()) * n
    return {
        "n_ids": int(n), "grid": [a, b],
        "ideal_flat_bits": ideal_flat, "ideal_two_stage_bits": ideal_two,
        "ideal_identity_gap": abs(ideal_flat - ideal_two),
        "realized_flat_bits": float(real_flat),
        "realized_two_stage_bits": float(real_two),
        "empirical_entropy_bits": h_emp,
    }


# ---------------------------------------------------------------------------
# Bits-back accounting (Theorem 6.6)
# ---------------------------------------------------------------------------

def bitsback_micro_identity(n_len: int = 8) -> dict:
    """Exact verification of identity (6.4) on an enumerable model.

    Model: token alphabet D = {a, b, ab, ba, abab} with weights theta;
    P(X = x, z) = prod theta(tokens of z) / Z_n over parses z of
    strings of length exactly n (Z_n by enumeration).  Latent z = the
    parse.  For q = (i) greedy-longest point mass and (ii) the exact
    posterior, E_q[net] is compared against
    -log2 P_X(x) + KL(q || P(z|x)) for EVERY x in {a,b}^n, and the
    expectation is taken under a data law D != P_X (iid uniform)."""
    toks = [b"a", b"b", b"ab", b"ba", b"abab"]
    theta = np.array([0.3, 0.2, 0.25, 0.15, 0.1])

    def parses(x: bytes):
        memo = {len(x): [([], 1.0)]}

        def rec(i):
            if i in memo:
                return memo[i]
            out = []
            for t, w in zip(toks, theta):
                if x[i:i + len(t)] == t:
                    for tail, tw in rec(i + len(t)):
                        out.append(([t] + tail, w * tw))
            memo[i] = out
            return out
        return rec(0)

    # Z_n: total weight over all strings of length n
    zx = {}
    for m in range(2 ** n_len):
        x = bytes(b"a"[0] if (m >> i) & 1 else b"b"[0]
                  for i in range(n_len))
        ps = parses(x)
        if ps:
            zx[x] = ps
    Z = sum(w for ps in zx.values() for _, w in ps)

    max_gap = 0.0
    exp_net_post, exp_net_greedy = 0.0, 0.0
    exp_rhs_post, exp_rhs_greedy = 0.0, 0.0
    dprob = 1.0 / (2 ** n_len)                     # data law: iid uniform
    covered = 0.0
    for x, ps in zx.items():
        px = sum(w for _, w in ps) / Z
        post = [w / (px * Z) for _, w in ps]
        # (ii) q = exact posterior
        net_post = sum(q * (-math.log2(w / Z) + math.log2(q))
                       for (_, w), q in zip(ps, post))
        rhs_post = -math.log2(px)                  # KL = 0
        # (i) q = point mass on the heaviest parse (greedy stand-in)
        j = max(range(len(ps)), key=lambda i: ps[i][1])
        net_greedy = -math.log2(ps[j][1] / Z)      # + log2 1
        rhs_greedy = -math.log2(px) + (-math.log2(post[j]))
        max_gap = max(max_gap, abs(net_post - rhs_post),
                      abs(net_greedy - rhs_greedy))
        exp_net_post += dprob * net_post
        exp_net_greedy += dprob * net_greedy
        exp_rhs_post += dprob * rhs_post
        exp_rhs_greedy += dprob * rhs_greedy
        covered += dprob
    return {
        "n_len": n_len, "strings_with_parse": len(zx),
        "data_mass_covered": covered,
        "max_pointwise_identity_gap": max_gap,
        "E_net_posterior": exp_net_post, "E_rhs_posterior": exp_rhs_post,
        "E_net_greedy": exp_net_greedy, "E_rhs_greedy": exp_rhs_greedy,
        "bitsback_saving_bits": exp_net_greedy - exp_net_post,
    }


def bitsback_real_coder(data: bytes, symbols, dictionary) -> dict:
    """KL term of (6.4) for the real coder's deterministic parse under
    the FROZEN final empirical token law: DP marginal over all parses
    vs the realized single-parse cost.  net(point mass) - (-log2 P_X)
    = -log2 P(z*|X) >= 0 is the bits-back saving left on the table."""
    K = len(dictionary)
    n = len(data)
    cnt = np.bincount(np.asarray(symbols), minlength=256 + K).astype(float)
    tot = cnt.sum()
    logp = np.full(256 + K, -np.inf)
    nz = cnt > 0
    logp[nz] = np.log2(cnt[nz] / tot)
    by_first = {}
    for k, tau in enumerate(dictionary):
        by_first.setdefault(tau[0], []).append((k, tau))
    # single-parse cost under frozen law
    single = float(-sum(logp[s] for s in symbols))
    # DP marginal: M[i] = log2 sum over parses of data[i:]
    M = np.full(n + 1, -np.inf)
    M[n] = 0.0
    for i in range(n - 1, -1, -1):
        opts = []
        if np.isfinite(logp[data[i]]) and np.isfinite(M[i + 1]):
            opts.append(logp[data[i]] + M[i + 1])
        for k, tau in by_first.get(data[i], ()):
            ell = len(tau)
            if data[i:i + ell] == tau and np.isfinite(M[i + ell]) \
                    and np.isfinite(logp[256 + k]):
                opts.append(logp[256 + k] + M[i + ell])
        if opts:
            m = max(opts)
            M[i] = m + math.log2(sum(2.0 ** (o - m) for o in opts))
    marginal = float(-M[0])
    return {
        "frozen_single_parse_bits": single,
        "frozen_marginal_bits": marginal,
        "kl_term_bits": single - marginal,          # -log2 P(z*|X) >= 0
        "kl_term_bpb": (single - marginal) / n,
    }


# ---------------------------------------------------------------------------
# Campaign driver
# ---------------------------------------------------------------------------

def run_h7(ladder: dict, bitsback_slice: int) -> dict:
    out = {"ladder": {}, "bitsback_micro": bitsback_micro_identity()}
    for corpus, sizes in ladder.items():
        rows = []
        for N in sizes:
            data = loader.read_range(corpus, 0, N)
            archive, st = encode(data)
            ok = decode(archive) == data
            x = np.frombuffer(data, dtype=np.uint8)
            h = common.entropy_rate_estimate(x)
            neut = index_neutrality(st["symbols"], st["K"])
            row = {
                "N": N, "roundtrip_ok": ok,
                "rate_bpb": st["total_bits"] / N,
                "h_hat_bpb": h["h_hat"], "h_hat_order": h["order"],
                "rate_over_Nh": st["total_bits"] / (h["h_hat"] * N),
                "K": st["K"], "n_syms": st["n_syms"],
                "dict_bits": st["dict_bits"],
                "overhead_fraction": st["dict_bits"] / st["total_bits"],
                "stream_bits": st["stream_bits"],
                "ideal_stream_bits": st["ideal_stream_bits"],
                "admitted_all_negative_delta":
                    all(a["delta_bits"] < 0 for a in st["admission_log"]),
                "index_neutrality": neut,
            }
            rows.append(row)
        out["ladder"][corpus] = rows
    # bits-back on the real coder, small slice
    data = loader.read_range("enwik8", 0, bitsback_slice)
    _archive, st = encode(data)
    out["bitsback_real"] = bitsback_real_coder(
        data, st["symbols"], st["dictionary"])
    out["bitsback_real"]["slice_bytes"] = bitsback_slice
    return out
