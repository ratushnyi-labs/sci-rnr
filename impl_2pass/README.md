# impl_2pass — two-pass (train-then-freeze) coding mode, lab prototype

Two-pass variant of the reference Type-I coder `impl/rnr1.py`
(container: `impl/FORMAT.md`): PASS 1 trains the order-W count tables
over the whole input with the reference integer update rules (no sync
resets), prunes them to the profitable contexts, and serializes the
frozen model into the archive; PASS 2 codes every sub-block with the
frozen tables — no adaptation, no per-sync cold start.  Sync semantics
and random access are unchanged (a frozen model has no cross-block
state).  The fallback rule makes the mode lossless in size: if
model + two-pass stream is not strictly smaller than the single-pass
archive, the single-pass archive is emitted bit-identically.

Format, prune rule and honest caveats: `FORMAT2.md` (archive version
1.1.0, header flag bit 1 = mode flag; fallback archives stay 1.0.0).

## Files

| file | purpose |
|------|---------|
| `two_pass.py` | coder: `pack2` / `unpack2` / `open_archive` (+ CLI `pack/unpack/read/info`) |
| `FORMAT2.md` | model-section format, prune-v1 rule, caveats |
| `run_checks.py` | PASS/FAIL gate, tiny inputs (<= 4 MB): round-trip incl. random access, determinism across processes, accounting identity (< 1e-3 bpb coder redundancy), white-noise fallback, 2 MiB size signal |
| `measure_smoke.py` | deferred 16-64 MB campaign driver (JSONL + CSV in the style of `experiments/type1_scale/to_csv.py`) |
| `measure_smoke.sh` | wrapper for the deferred campaign; prints the plan unless `--run` is given — **not** executed by the gate |

## Gate

```
/Users/para/.venvs/rnr/bin/python impl_2pass/run_checks.py [--quick]
```

Gate run 2026-07-11: OVERALL PASS 7/7 (`run_checks_2026-07-11.log`).
First size signal (2 MiB, W=3): two-pass beats single-pass by
1.18 bpb (K=16 KiB) / 0.63 bpb (K=64 KiB) on enwik8 and by
0.88 / 0.40 bpb on scripts-src — the advantage grows as K shrinks
(more per-sync cold starts amortized into one shared model), as
predicted.  Model header share ~0.37-0.40 bpb at 2 MiB.

## Deferred measurement

```
sh impl_2pass/measure_smoke.sh          # plan only
sh impl_2pass/measure_smoke.sh --run    # 16-64 MB campaign (see cost
                                        # and memory caveats inside)
```

Build `impl_fast` first (`sh impl_fast/build.sh`) — the single-pass
side of the campaign uses the byte-identity-gated C engine; without it
the pure reference makes the grid unusably slow.

## Known limitations (honest)

* Encoder-side prune decisions use float `log2` (model selection only;
  decode is exact and the model ships in the archive).  Same-platform
  encoder determinism is gated (T3); cross-platform encoder
  bit-reproducibility would need a fixed-point log.
* The prune saving estimate charges each context's full training usage
  to it even where a kept longer context takes over in PASS 2
  (overestimates shorter-context usage); the whole-archive fallback
  comparison bounds any loss.
* Pure-Python throughput (~0.5 MB/s PASS 2) and PASS-1 memory
  (sparse dicts; expected 2-4 GB at 64 MB text, W=3).
* The repo has no 16-64 MB code/structured payloads; the campaign caps
  those types at the actual payload sizes and records the truth.
