#!/usr/bin/env python3
"""End-to-end smoke run of the baseline+metrics+results harness.

Generates 1 MiB of synthetic Markov text-like data (order-1 Markov
chain over a Zipf-weighted synthetic vocabulary, deterministic seed),
runs every available baseline coder on it in isolated-child mode,
records the measurements in bench/results/smoke.jsonl, and prints the
Markdown summary table.

This smoke run uses SYNTHETIC data only and 3 runs per coder (the
frozen protocol in METHODS.md prescribes 5 runs on the real corpora;
smoke results are labeled corpus "synthetic-markov-1MiB" and are not
hypothesis evidence). The real-corpus loader hookup happens at
integration.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import baselines
import metrics
import results

SMOKE_CORPUS = "synthetic-markov-1MiB"
SMOKE_RUNS = 3
STORE = Path(__file__).resolve().parent / "results" / "smoke.jsonl"


def markov_text(n_bytes: int = 1 << 20, seed: int = metrics.MASTER_SEED) -> bytes:
    """Synthetic text-like bytes from an order-1 Markov chain over words.

    A seeded vocabulary of 64 lowercase words receives a sparse random
    transition matrix (each word has ~8 plausible successors with
    Dirichlet-like weights), which yields text-like local statistics:
    skewed word frequencies, repeated collocations, spaces, and
    sentence-like punctuation. The small vocabulary makes it more
    compressible than natural text (about 1 bpb under gzip -9 versus
    2-3 bpb for English) but it is structured, non-trivial input --
    adequate for exercising the harness, not for rate claims.
    """
    rng = np.random.default_rng(seed)
    vocab_size = 64
    letters = np.array(list("abcdefghijklmnopqrstuvwxyz"))
    vocab = []
    seen = set()
    while len(vocab) < vocab_size:
        length = int(rng.integers(2, 10))
        word = "".join(rng.choice(letters, size=length))
        if word not in seen:
            seen.add(word)
            vocab.append(word.encode("ascii"))

    # Sparse row-stochastic transition matrix: 8 successors per word.
    successors = np.zeros((vocab_size, 8), dtype=np.int64)
    cumweights = np.zeros((vocab_size, 8), dtype=np.float64)
    for i in range(vocab_size):
        successors[i] = rng.choice(vocab_size, size=8, replace=False)
        w = rng.gamma(shape=0.5, size=8) + 1e-3  # skewed, Zipf-like
        cumweights[i] = np.cumsum(w / w.sum())

    chunks = []
    total = 0
    state = int(rng.integers(vocab_size))
    words_in_sentence = 0
    sentence_len = int(rng.integers(6, 18))
    uniforms = rng.random(2 * n_bytes // 5 + 1024)
    u_idx = 0
    while total < n_bytes:
        u = uniforms[u_idx]
        u_idx += 1
        if u_idx >= len(uniforms):  # pragma: no cover - sized generously
            uniforms = rng.random(len(uniforms))
            u_idx = 0
        state = int(successors[state, np.searchsorted(cumweights[state], u)])
        word = vocab[state]
        words_in_sentence += 1
        if words_in_sentence >= sentence_len:
            chunks.append(word + b".\n")
            words_in_sentence = 0
            sentence_len = int(rng.integers(6, 18))
        else:
            chunks.append(word + b" ")
        total += len(chunks[-1])
    return b"".join(chunks)[:n_bytes]


def main() -> int:
    data = markov_text()
    sha = hashlib.sha256(data).hexdigest()
    versions = baselines.codec_versions()
    coders = baselines.available_codecs()
    print(f"smoke corpus: {SMOKE_CORPUS}, sha256={sha[:16]}..., "
          f"{len(data)} bytes")
    print(f"coders: {', '.join(coders)}")

    for coder in coders:
        for run_idx in range(SMOKE_RUNS):
            seed = metrics.derive_seed(f"smoke-{coder}", run_idx)
            res = baselines.run_baseline(coder, data, isolated=True)
            rec = results.make_record(
                corpus=SMOKE_CORPUS,
                coder=coder,
                config=res.config,
                seed=seed,
                measured={
                    "bpb": metrics.bits_per_byte(
                        res.compressed_size, res.original_size
                    ),
                    "compressed_size": res.compressed_size,
                    "original_size": res.original_size,
                    "enc_s": res.enc_s,
                    "dec_s": res.dec_s,
                    "peak_rss_estimate": res.peak_rss_estimate,
                },
                corpus_sha256=sha,
                rss_mode=res.rss_mode,
                coder_version=versions.get(coder, "unknown"),
                run_index=run_idx,
            )
            results.append(STORE, rec)
            print(
                f"  {coder} run {run_idx}: {rec['measured']['bpb']:.4f} bpb, "
                f"enc {res.enc_s:.3f}s, dec {res.dec_s:.3f}s"
            )

    print()
    print(results.summarize(STORE))
    return 0


if __name__ == "__main__":
    sys.exit(main())
