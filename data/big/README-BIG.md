# data/big — Type-I scale-benchmark corpus (500 / 1500 / 5000 MB ladder)

Size convention: **decimal MB, 1 MB = 10^6 bytes**.  Ladder tiers are exactly
500 000 000 / 1 500 000 000 / 5 000 000 000 bytes.

Design: **prefix-nested masters**.  Each content type has ONE master payload
(`masters/`); a tier is the byte-prefix of that master, so results at
different scales within a type are measured on nested data.  The manifest
records the sha256 of the full master AND of every tier prefix;
`run_checks.py` recomputes all of them from the payload in a single streaming
pass, which verifies the prefix-nesting property directly.

Content types (see `MANIFEST-BIG.json` for authoritative hashes, sources,
licenses, and per-type status):

| type          | ladder (MB)     | recipe |
|---------------|-----------------|--------|
| white_noise   | 500/1500/5000   | seeded SHAKE-256 XOF, 64 MiB chunks, regenerable from recorded seed |
| text          | 500/1000 (capped)| enwik9 (10^9 bytes, Hutter-prize corpus); text is never looped/repeated — repetition corrupts compressibility, so 1500/5000 tiers do not exist by design |
| video         | 500/1500/5000   | concatenation of distinct open-licensed Blender open-movie files (distinct titles first, then distinct encodes of used titles; no file repeated; tail truncated at the 5000 MB boundary) |
| images        | 500/1500/5000   | video-derived JPEG frames: concatenated MJPEG q:v 2 stream extracted from the downloaded movies (fixed source order, every-Nth-frame, per-source caps) |
| archive_mix   | 500/1500/5000   | deterministic GNU tar of the existing `data/payloads/*` small corpora plus seeded 32–96 MB slices of the other masters |
| precompressed | 500/1500/5000   | zstd -19 -T0 of the text master, then of the video master, concatenated (structured incompressible) |

Payloads (`downloads/`, `masters/`, `logs/`, `tmp/`) are git-ignored; only
scripts and `MANIFEST-BIG.json` are committed.

Workflow:

```
python fetch_big.py --fetch --jobs 4     # resumable downloads + sha256
python fetch_big.py --build all          # build masters (deterministic)
python fetch_big.py --manifest           # (re)write MANIFEST-BIG.json
python run_checks.py                     # full integrity gate (PASS/FAIL)
python run_checks.py --quick             # sizes only, skips full hashing
```

Access from experiment code goes through `loader_big.py`
(`iter_blocks_big(name, tier_mb, K)` etc.), never through raw paths; the
loader refuses tiers the manifest does not mark as materialized.
