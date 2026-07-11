//! RNR1 container: pack / unpack / random-access read (impl/FORMAT.md).
//!
//! Layout: 48-byte header, m x 32-byte seek-index entries, repair stream
//! (concatenated per-sub-block arithmetic-coder blobs; byte-aligned).
//! All multi-byte integers little-endian.  Validation rules and their
//! order mirror the normative reference `impl/rnr1.py` (error messages
//! match its ValueError texts).
//!
//! Multi-threaded pack/unpack: sub-blocks are self-contained restart
//! points (fresh predictor and coder state), so per-block work is a pure
//! function of (block bytes, W); assembling results in block order
//! yields output byte-identical to the sequential path regardless of
//! scheduling (verified by run_checks.py).

use crate::coder::{decode_subblock, encode_subblock, h_v, model_hash, AEnc, Pred};

pub const MAGIC: &[u8; 4] = b"RNR1";
pub const VERSION: (u8, u8, u8) = (1, 0, 0);
pub const FLAG_SUBBLOCK_HASHES: u16 = 0x0001;
pub const HEADER_SIZE: usize = 48;
pub const INDEX_ENTRY_SIZE: usize = 32;
pub const MODE_CODED: u32 = 0;
pub const MODE_RAW: u32 = 1;

fn clamp_threads(threads: u32, m: u64) -> u32 {
    let mut t = threads.max(1).min(64);
    if (t as u64) > m {
        t = if m > 0 { m as u32 } else { 1 };
    }
    t
}

// ---------------------------------------------------------------------
// pack
// ---------------------------------------------------------------------

struct BlockResult {
    j: u64,
    blob: Option<Vec<u8>>, // None = raw store mode (use source bytes)
    mode: u32,
    hash: [u8; 8],
}

/// Encode one sub-block and apply the fail-safe store-mode escape
/// (Theorem 7.22 semantics: never pay more than raw for a sub-block).
fn pack_block(block: &[u8], pred: &mut Pred, enc: &mut AEnc, j: u64) -> BlockResult {
    encode_subblock(block, pred, enc);
    if enc.buf.len() >= block.len() {
        BlockResult { j, blob: None, mode: MODE_RAW, hash: h_v(block) }
    } else {
        BlockResult { j, blob: Some(enc.buf.clone()), mode: MODE_CODED, hash: h_v(block) }
    }
}

/// Encode data into an RNR1 archive.  `threads > 1` parallelizes over
/// independent sub-blocks; output is byte-identical to `threads == 1`.
pub fn pack(data: &[u8], w: u32, k: u32, threads: u32) -> Result<Vec<u8>, String> {
    if !(1..=8).contains(&w) {
        return Err("W must be in [1, 8]".into());
    }
    if k < 256 {
        return Err("K must be >= 256".into());
    }
    let n = data.len() as u64;
    let m = if n > 0 { (n + k as u64 - 1) / k as u64 } else { 0 };
    let threads = clamp_threads(threads, m);

    let mut results: Vec<BlockResult> = if threads <= 1 {
        let mut pred = Pred::new(w);
        let mut enc = AEnc::new();
        (0..m)
            .map(|j| {
                let lo = (j * k as u64) as usize;
                let hi = ((j + 1) * k as u64).min(n) as usize;
                pack_block(&data[lo..hi], &mut pred, &mut enc, j)
            })
            .collect()
    } else {
        // Static striping: thread t handles blocks j with j % threads == t.
        std::thread::scope(|s| {
            let mut handles = Vec::new();
            for t in 0..threads as u64 {
                handles.push(s.spawn(move || {
                    let mut pred = Pred::new(w);
                    let mut enc = AEnc::new();
                    let mut out = Vec::new();
                    let mut j = t;
                    while j < m {
                        let lo = (j * k as u64) as usize;
                        let hi = ((j + 1) * k as u64).min(n) as usize;
                        out.push(pack_block(&data[lo..hi], &mut pred, &mut enc, j));
                        j += threads as u64;
                    }
                    out
                }));
            }
            let mut all: Vec<BlockResult> = Vec::with_capacity(m as usize);
            for h in handles {
                all.extend(h.join().expect("worker thread panicked"));
            }
            all.sort_by_key(|r| r.j);
            all
        })
    };

    // Assemble container in block order (identical to the sequential path).
    let mut entries = Vec::with_capacity(m as usize * INDEX_ENTRY_SIZE);
    let mut rep: Vec<u8> = Vec::new();
    let mut repair_bit_offset: u64 = 0;
    for r in results.iter_mut() {
        let lo = (r.j * k as u64) as usize;
        let hi = ((r.j + 1) * k as u64).min(n) as usize;
        let blob: &[u8] = match &r.blob {
            Some(b) => b,
            None => &data[lo..hi],
        };
        entries.extend_from_slice(&(r.j * k as u64).to_le_bytes());
        entries.extend_from_slice(&repair_bit_offset.to_le_bytes());
        entries.extend_from_slice(&0u32.to_le_bytes()); // snapshot_length
        entries.extend_from_slice(&r.mode.to_le_bytes());
        entries.extend_from_slice(&r.hash);
        rep.extend_from_slice(blob);
        repair_bit_offset += 8 * blob.len() as u64;
    }

    let mut out = Vec::with_capacity(HEADER_SIZE + entries.len() + rep.len());
    out.extend_from_slice(MAGIC);
    out.push(VERSION.0);
    out.push(VERSION.1);
    out.push(VERSION.2);
    out.push(0); // reserved
    out.extend_from_slice(&FLAG_SUBBLOCK_HASHES.to_le_bytes());
    out.extend_from_slice(&(w as u16).to_le_bytes());
    out.extend_from_slice(&k.to_le_bytes());
    out.extend_from_slice(&n.to_le_bytes());
    out.extend_from_slice(&model_hash(w));
    out.extend_from_slice(&h_v(data)); // v = H_v(X), Definition 3.1
    out.extend_from_slice(&m.to_le_bytes());
    out.extend_from_slice(&entries);
    out.extend_from_slice(&rep);
    Ok(out)
}

// ---------------------------------------------------------------------
// Archive parsing (header + seek index validation)
// ---------------------------------------------------------------------

pub struct Entry {
    pub byte_offset: u64,
    pub bit_offset: u64,
    pub mode: u32,
    pub hash: [u8; 8],
}

pub struct Archive<'a> {
    pub flags: u16,
    pub w: u32,
    pub k: u32,
    pub n: u64,
    pub model: [u8; 8],
    pub v: [u8; 8],
    pub m: u64,
    pub entries: Vec<Entry>,
    pub repair: &'a [u8],
}

fn rd16(p: &[u8]) -> u16 {
    u16::from_le_bytes([p[0], p[1]])
}
fn rd32(p: &[u8]) -> u32 {
    u32::from_le_bytes([p[0], p[1], p[2], p[3]])
}
fn rd64(p: &[u8]) -> u64 {
    u64::from_le_bytes(p[..8].try_into().unwrap())
}

impl<'a> Archive<'a> {
    /// Parse and validate header + seek index (validation rules and
    /// order per the reference `rnr1.Archive.__init__`).
    pub fn parse(raw: &'a [u8]) -> Result<Self, String> {
        if raw.len() < HEADER_SIZE {
            return Err("truncated archive: no header".into());
        }
        if &raw[0..4] != MAGIC {
            return Err("bad magic: not an RNR1 archive".into());
        }
        // R-13.2: same MAJOR required; refuse newer MINOR.
        if raw[4] != VERSION.0 || raw[5] > VERSION.1 {
            return Err(format!(
                "unsupported archive format version ({}, {}, {})",
                raw[4], raw[5], raw[6]
            ));
        }
        let flags = rd16(&raw[8..]);
        let w = rd16(&raw[10..]) as u32;
        let k = rd32(&raw[12..]);
        let n = rd64(&raw[16..]);
        let mut model = [0u8; 8];
        model.copy_from_slice(&raw[24..32]);
        let mut v = [0u8; 8];
        v.copy_from_slice(&raw[32..40]);
        let m = rd64(&raw[40..]);
        // R-13.3: predictor hash mismatch is a hard failure.
        if w > 8 || model != model_hash(w) {
            return Err("model description hash mismatch".into());
        }
        if k == 0 {
            return Err("unsupported archive parameters".into());
        }
        let expect_m = if n > 0 { (n + k as u64 - 1) / k as u64 } else { 0 };
        if m != expect_m {
            return Err("inconsistent sub-block count".into());
        }
        if m > ((raw.len() - HEADER_SIZE) / INDEX_ENTRY_SIZE) as u64 {
            return Err("truncated archive".into());
        }
        let mut entries = Vec::with_capacity(m as usize);
        let mut off = HEADER_SIZE;
        for _ in 0..m {
            let e = &raw[off..off + INDEX_ENTRY_SIZE];
            let mut hash = [0u8; 8];
            hash.copy_from_slice(&e[24..32]);
            entries.push(Entry {
                byte_offset: rd64(e),
                bit_offset: rd64(&e[8..]),
                mode: rd32(&e[20..]),
                hash,
            });
            off += INDEX_ENTRY_SIZE;
        }
        // R-8.2.1: index sorted by byte_offset ascending (strictly).
        for j in 1..m as usize {
            if entries[j].byte_offset <= entries[j - 1].byte_offset {
                return Err("seek index not sorted".into());
            }
        }
        // Blocks are byte-aligned by the flush.
        for e in &entries {
            if e.bit_offset % 8 != 0 {
                return Err("seek index invalid".into());
            }
        }
        Ok(Archive { flags, w, k, n, model, v, m, entries, repair: &raw[off..] })
    }

    /// Repair-stream byte range of sub-block j (slice-clamped like the
    /// reference's Python slicing so corrupt offsets fail later via the
    /// hash checks, not via a panic).
    fn block_blob(&self, j: usize) -> &'a [u8] {
        let start = (self.entries[j].bit_offset / 8) as usize;
        let end = if j + 1 < self.m as usize {
            (self.entries[j + 1].bit_offset / 8) as usize
        } else {
            self.repair.len()
        };
        let start = start.min(self.repair.len());
        let end = end.min(self.repair.len()).max(start);
        &self.repair[start..end]
    }

    fn block_len(&self, j: usize) -> usize {
        (self.k as u64).min(self.n - j as u64 * self.k as u64) as usize
    }

    /// Decode sub-block j into `dst` (dst.len() = number of positions to
    /// decode, <= block_len).  `pred` is reused across calls.
    fn decode_block_into(
        &self,
        j: usize,
        pred: &mut Pred,
        dst: &mut [u8],
    ) -> Result<(), String> {
        let blob = self.block_blob(j);
        if self.entries[j].mode == MODE_RAW {
            if blob.len() < dst.len() {
                return Err("truncated archive".into());
            }
            dst.copy_from_slice(&blob[..dst.len()]);
            Ok(())
        } else {
            decode_subblock(blob, dst.len(), pred, dst)
        }
    }

    fn check_block_hash(&self, j: usize, out: &[u8]) -> Result<(), String> {
        if self.flags & FLAG_SUBBLOCK_HASHES != 0 && h_v(out) != self.entries[j].hash {
            return Err(format!("sub-block {} hash mismatch", j));
        }
        Ok(())
    }

    /// Full sequential/parallel decode; checks H_v(X) = v (Def 3.1).
    pub fn unpack(&self, verify: bool, threads: u32) -> Result<Vec<u8>, String> {
        let threads = clamp_threads(threads, self.m);
        let mut out = vec![0u8; self.n as usize];
        if threads <= 1 {
            let mut pred = Pred::new(self.w);
            for j in 0..self.m as usize {
                let lo = j * self.k as usize;
                let hi = (lo + self.block_len(j)).min(self.n as usize);
                self.decode_block_into(j, &mut pred, &mut out[lo..hi])?;
                if verify {
                    self.check_block_hash(j, &out[lo..hi])?;
                }
            }
        } else {
            // Static striping over disjoint output chunks (safe, no copy).
            let chunks: Vec<(usize, &mut [u8])> =
                out.chunks_mut(self.k as usize).enumerate().collect();
            let mut buckets: Vec<Vec<(usize, &mut [u8])>> =
                (0..threads).map(|_| Vec::new()).collect();
            for (j, ch) in chunks {
                buckets[j % threads as usize].push((j, ch));
            }
            let errs: Vec<Result<(), String>> = std::thread::scope(|s| {
                let handles: Vec<_> = buckets
                    .into_iter()
                    .map(|bucket| {
                        s.spawn(move || {
                            let mut pred = Pred::new(self.w);
                            for (j, dst) in bucket {
                                self.decode_block_into(j, &mut pred, dst)?;
                                if verify {
                                    self.check_block_hash(j, dst)?;
                                }
                            }
                            Ok(())
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .map(|h| h.join().expect("worker thread panicked"))
                    .collect()
            });
            for e in errs {
                e?;
            }
        }
        if verify && h_v(&out) != self.v {
            return Err("archive verification failed: H_v mismatch".into());
        }
        Ok(out)
    }

    /// Binary search for the latest sync point <= p (R-8.2.2).
    fn find_sync(&self, p: u64) -> usize {
        let (mut lo, mut hi) = (0usize, self.m as usize - 1);
        while lo < hi {
            let mid = (lo + hi + 1) / 2;
            if self.entries[mid].byte_offset <= p {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        lo
    }

    /// Random access: source bytes [p, p+k') with k' = min(len, n-p),
    /// decoding only the covering sub-block(s) (Theorem 10.3 pattern).
    /// Returns (data, warmup_bytes, decoded_bytes).
    pub fn read(
        &self,
        p: u64,
        len: u64,
        verify: bool,
    ) -> Result<(Vec<u8>, u64, u64), String> {
        if p > self.n {
            return Err("position out of range".into());
        }
        let kprime = len.min(self.n - p);
        if kprime == 0 {
            return Ok((Vec::new(), 0, 0));
        }
        let j0 = self.find_sync(p);
        let j1 = self.find_sync(p + kprime - 1);
        let mut pieces: Vec<u8> = Vec::with_capacity(kprime as usize);
        let mut decoded: u64 = 0;
        let mut pred = Pred::new(self.w);
        for j in j0..=j1 {
            let base = self.entries[j].byte_offset;
            let blk_len = self.block_len(j) as u64;
            let take = if j < j1 || verify {
                blk_len // full block
            } else {
                (p + kprime - base).min(blk_len)
            };
            let mut blk = vec![0u8; take as usize];
            self.decode_block_into(j, &mut pred, &mut blk)?;
            if verify && take == blk_len {
                self.check_block_hash(j, &blk)?;
            }
            decoded += take;
            let lo = (p.max(base) - base) as usize;
            let hi = ((p + kprime).min(base + take) - base) as usize;
            pieces.extend_from_slice(&blk[lo..hi]);
        }
        Ok((pieces, p - self.entries[j0].byte_offset, decoded))
    }
}
