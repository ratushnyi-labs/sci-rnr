#!/bin/zsh
# Sequential smoke campaign driver (idempotent; resumes via journals).
D=/Users/para/work/rnr/experiments/type1_scale
P=/Users/para/.venvs/rnr/bin/python
cd "$D" || exit 1
echo "DRIVER start $(date) pid=$$"
"$P" -u run_campaign.py --tier smoke --types white_noise,text
echo "STAGE1 rc=$?"
"$P" -u run_campaign.py --tier smoke \
  --types video,images,archive_mix,precompressed \
  --coders zstd-19,zstd-22u,xz-6,xz-9e,brotli-q11,brotli-q9,gzip-9,bzip2-9
echo "STAGE2 rc=$?"
"$P" -u run_campaign.py --tier smoke48 --types white_noise,text \
  --coders zstd-19,xz-6,gzip-9,bzip2-9
echo "STAGE3 rc=$?"
echo "DRIVER done $(date)"
