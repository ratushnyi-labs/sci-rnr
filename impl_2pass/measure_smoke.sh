#!/bin/sh
# measure_smoke.sh -- DEFERRED 16-64 MB smoke campaign for the two-pass mode.
#
# NOT executed by run_checks.py and NOT to be started while a timing
# campaign is running on this host: it is a single-core, single-process
# measurement of roughly 30-60 minutes WITH the impl_fast C engine
# built (run `sh impl_fast/build.sh` first), and of DAYS without it
# (the single-pass side then falls back to the pure reference).  The
# 64 MB text tier needs an expected 2-4 GB peak RSS for PASS-1
# training (see the memory caveat in measure_smoke.py).
#
# Grid: {text, code, structured} x tiers {16, 64} decimal MB
#       x K {16 KiB, 64 KiB, 256 KiB} x {rnr1-single, rnr1-2pass}.
# Caveat (recorded per row in the CSV): the repo has no 16-64 MB code
# or structured payload; those types are capped at the actual payload
# size (3.4 MB scripts-src tar, 12.4 MB sqlite db) rather than padded
# or repeated.
#
# Output: impl_2pass/out/smoke_results.jsonl (append-only, resumable)
#         impl_2pass/out/smoke_results.csv   (style of
#         experiments/type1_scale/to_csv.py)
#
# Usage:
#   sh impl_2pass/measure_smoke.sh            # print the plan only
#   sh impl_2pass/measure_smoke.sh --run      # measure (multi-hour)
#   sh impl_2pass/measure_smoke.sh --run --no-decode   # faster variant:
#       verified random-access spot reads instead of full two-pass
#       decodes (weaker round-trip evidence, recorded in the status
#       column)
#
# The script is restart-safe: completed cells are skipped on rerun.

set -eu

REPO="$(cd "$(dirname "$0")/.." && pwd)"
PY="/Users/para/.venvs/rnr/bin/python"

exec "$PY" "$REPO/impl_2pass/measure_smoke.py" "$@"
