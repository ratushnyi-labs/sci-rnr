#!/bin/sh
# Type-I scale benchmark -- single entrypoint.
#
#   ./run_benchmarks.sh                 # 500 MB tier (all types, all coders)
#   ./run_benchmarks.sh 500 1500        # given tiers, in order
#   ./run_benchmarks.sh all             # 500, then 1500, then 5000
#   TYPES=video,images ./run_benchmarks.sh 500     # type filter
#   CODERS=rnr1-fast,zstd-19 ./run_benchmarks.sh 500  # coder filter
#
# Each tier run is checkpointed per cell (out/<tier>/cells.journal): rerunning
# after an interruption resumes where it stopped.  RUNS=N adds repetitions
# (rates are bit-deterministic -- repetitions estimate the timing/RSS
# distribution; every repetition's archive sha256 is recorded as the
# determinism witness).  For mass repetition (e.g. 100000 runs) use
# repeat_stats.py, which batches, resumes, and asserts the witness.  After every tier this
# script re-exports out/results.csv (the flat table for graphics) and
# regenerates the figures from the recorded JSONL.
#
# Prerequisites: data/big masters fetched (data/big/fetch_big.py) and the
# fast engine built (impl_fast/build.sh) -- both are checked below.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
PY="/Users/para/.venvs/rnr/bin/python"

TIERS="${*:-500}"
[ "$TIERS" = "all" ] && TIERS="500 1500 5000"

# fast engine present?
if [ ! -f "$REPO/impl_fast/librnr1fast.dylib" ]; then
    echo "[bench] building fast engine ..."
    (cd "$REPO/impl_fast" && sh build.sh)
fi

# data/big masters present? (fetch_big verifies hashes; cheap when complete)
$PY - <<'EOF'
import json, pathlib, sys
big = pathlib.Path("/Users/para/work/rnr/data/big")
man = json.load(open(big / "MANIFEST-BIG.json"))
missing = [e["file"] for e in man["entries"].values()
           if not (big / e["file"]).exists()]
if missing:
    sys.exit(f"[bench] data/big masters missing: {missing}\n"
             f"        run: {big}/fetch_big.py")
print("[bench] data/big masters present")
EOF

for TIER in $TIERS; do
    echo "== tier $TIER =="
    ARGS="--tier $TIER"
    [ -n "$TYPES" ]  && ARGS="$ARGS --types $TYPES"
    [ -n "$CODERS" ] && ARGS="$ARGS --coders $CODERS"
    [ -n "$RUNS" ]   && ARGS="$ARGS --runs $RUNS"
    $PY -u "$HERE/run_campaign.py" $ARGS
    $PY "$HERE/to_csv.py"
    $PY "$HERE/make_figures.py" || echo "[bench] figures: skipped (partial data ok)"
done

echo "== done: results at $HERE/out/results.csv, figures at $HERE/out/figures/ =="
