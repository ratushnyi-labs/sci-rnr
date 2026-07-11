# Figure and table captions (Type-I scale campaign)

Generated deterministically by make_figures.py from: out/500/results.jsonl, out/runset16/results.jsonl, out/smoke/results.jsonl, out/smoke48/results.jsonl

**fig_rate_vs_size**: Compression rate (bits per byte, lower is better) versus corpus size, one panel per content type; size axis logarithmic. Hue = coder family, linestyle = variant, marker = coder (colorblind-safe; validated worst adjacent CVD deltaE 25).

**fig_enc_time_vs_size**: Encode wall-clock seconds versus corpus size, log-log, one panel per content type. Single-threaded CLI invocations measured around the whole child process.

**fig_dec_time_vs_size**: Decode wall-clock seconds versus corpus size, log-log, one panel per content type; decode output is streamed to a hash, not written to disk (except RNR coders, see protocol 4.3).

**fig_pareto_16.777MB**: Rate-versus-encode-throughput Pareto plot at corpus size 16.777 MB, one panel per content type. Each point is one coder (marker/hue as in the shared legend); the gray staircase joins the Pareto-efficient coders (no coder is both faster to encode and denser). Throughput axis is logarithmic. Censored/failed cells are absent (see the summary CSVs).

**fig_pareto_50.332MB**: Rate-versus-encode-throughput Pareto plot at corpus size 50.332 MB, one panel per content type. Each point is one coder (marker/hue as in the shared legend); the gray staircase joins the Pareto-efficient coders (no coder is both faster to encode and denser). Throughput axis is logarithmic. Censored/failed cells are absent (see the summary CSVs).

**fig_pareto_500MB**: Rate-versus-encode-throughput Pareto plot at corpus size 500 MB, one panel per content type. Each point is one coder (marker/hue as in the shared legend); the gray staircase joins the Pareto-efficient coders (no coder is both faster to encode and denser). Throughput axis is logarithmic. Censored/failed cells are absent (see the summary CSVs).

**summary_<size>MB_<type>.csv**: per-(size, type) table view of every metric -- mean over runs of rate (bpb), compression ratio, encode/decode wall seconds and derived MB/s, peak RSS per phase (/usr/bin/time -l child high-watermark, MB), round-trip verification, and censoring/failure counts. This is the accessible table view backing every figure.

