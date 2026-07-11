//! Core RNR1 Type-I coding engine.
//!
//! Direct port of the validated integer state machine of the normative
//! reference `impl/rnr1.py` (see also the C port `impl_fast/rnr1_fast.c`,
//! whose derived-cache design this file reuses).  Every component is
//! bit-exact with respect to the reference by construction and verified
//! by `rnr-cli/run_checks.py` against the full corpus matrix:
//!
//!  * Code12 E12-A systematic (12,8) map (Section 4.2 parity equations).
//!    E12-A is linear with minimum distance 3, so for the byte-aligned
//!    argmax predictor the error mask between two valid codewords has
//!    weight 0 or >= 3: classes C1/C2 are structurally empty and the
//!    class partition below exploits exactly that (checked at init).
//!
//!  * Order-W adaptive byte-context predictor: per-order exact-key
//!    context tables (open addressing, full 64-bit key compare), u16
//!    counts with the +1 floor applied at read time
//!    (freq = 1 + 32*count, total = 256 + 32*ctx_total), deterministic
//!    longest-context-first backoff, cap-halving at ctx_total >= 2^16,
//!    full reset to Q_init at every sub-block boundary (generation tags).
//!
//!  * 32-bit binary arithmetic coder, Witten-Neal-Cleary renormalization
//!    with pending-bit carry handling, MSB-first bit packing, identical
//!    termination (2 bits + carries, zero-padded to a byte boundary).
//!    All intermediates fit in u64: range <= 2^32, totals < 2^24.
//!
//!  * Type-I class/payload partition (Section 4.3): top-16 ordering
//!    (count desc, byte asc) identical to the reference's stable argsort,
//!    exact partition sums, candidate-rank payload for CAND,
//!    ascending-byte payload for ESC.
//!
//! Integer determinism: the whole coding path is integer-only, no
//! floating point, no randomness, fixed evaluation order.

use sha2::{Digest, Sha256};

pub const COUNT_SCALE: u64 = 32;
pub const COUNT_CAP: u32 = 65536;
pub const NUM_CAND: u64 = 15;

// ---------------------------------------------------------------------
// Code12 E12-A table (Section 4.2)
// ---------------------------------------------------------------------

/// E12-A lookup table: bit i of the entry holds codeword position i+1.
pub fn build_e12() -> [u16; 256] {
    let mut table = [0u16; 256];
    for b in 0..256usize {
        let bits: Vec<u16> = (0..8).map(|j| ((b >> j) & 1) as u16).collect();
        let mut c = [0u16; 13]; // 1-based positions
        c[3] = bits[0];
        c[5] = bits[1];
        c[6] = bits[2];
        c[7] = bits[3];
        c[9] = bits[4];
        c[10] = bits[5];
        c[11] = bits[6];
        c[12] = bits[7];
        c[1] = c[3] ^ c[5] ^ c[7] ^ c[9] ^ c[11];
        c[2] = c[3] ^ c[6] ^ c[7] ^ c[10] ^ c[11];
        c[4] = c[5] ^ c[6] ^ c[7] ^ c[12];
        c[8] = c[9] ^ c[10] ^ c[11] ^ c[12];
        let mut w = 0u16;
        for i in 1..=12 {
            w |= c[i] << (i - 1);
        }
        table[b] = w;
    }
    // Self-check: minimum pairwise mask weight 3 (linear, d_min = 3).
    // This is what makes classes C1/C2 structurally empty; the class
    // partition in encode/decode relies on it.
    for a in 0..256usize {
        for b in (a + 1)..256usize {
            let w = ((table[a] ^ table[b]) & 0x0fff).count_ones();
            assert!(w >= 3, "E12-A self-check failed: d_min < 3");
        }
    }
    table
}

/// Truncated SHA-256 of the canonical model description (R-13.3);
/// byte-identical to `rnr1.model_description_hash(W)`.
pub fn model_hash(w: u32) -> [u8; 8] {
    let e12 = build_e12();
    // E12 serialized as 256 little-endian int64 values (numpy "<i8").
    let mut e12b = [0u8; 256 * 8];
    for (i, v) in e12.iter().enumerate() {
        e12b[8 * i..8 * i + 8].copy_from_slice(&(*v as i64).to_le_bytes());
    }
    let d = Sha256::digest(e12b);
    let hex: String = d[..8].iter().map(|b| format!("{:02x}", b)).collect();
    let desc = format!(
        "rnr-type1-code12a;predictor=ngram;W={};scale={};cap={};cand={};\
         floor=1;classes=c0,c1,c2,cand,esc;e12={}",
        w, COUNT_SCALE, COUNT_CAP, NUM_CAND, hex
    );
    let d = Sha256::digest(desc.as_bytes());
    let mut out = [0u8; 8];
    out.copy_from_slice(&d[..8]);
    out
}

/// Verification value H_v (Definition 3.1): truncated SHA-256.
pub fn h_v(data: &[u8]) -> [u8; 8] {
    let d = Sha256::digest(data);
    let mut out = [0u8; 8];
    out.copy_from_slice(&d[..8]);
    out
}

// ---------------------------------------------------------------------
// Binary arithmetic coder (32-bit registers, WNC renormalization)
// ---------------------------------------------------------------------

const AC_MASK: u64 = 0xffff_ffff;
const AC_HALF: u64 = 0x8000_0000;
const AC_Q1: u64 = 0x4000_0000;
const AC_Q3: u64 = 0xc000_0000;

pub struct AEnc {
    low: u64,
    high: u64,
    pending: u64,
    pub buf: Vec<u8>,
    cur: u32,
    ncur: u32,
}

impl AEnc {
    pub fn new() -> Self {
        AEnc { low: 0, high: AC_MASK, pending: 0, buf: Vec::new(), cur: 0, ncur: 0 }
    }

    pub fn reset(&mut self) {
        self.low = 0;
        self.high = AC_MASK;
        self.pending = 0;
        self.buf.clear();
        self.cur = 0;
        self.ncur = 0;
    }

    #[inline]
    fn put(&mut self, bit: u32) {
        self.cur = (self.cur << 1) | bit;
        self.ncur += 1;
        if self.ncur == 8 {
            self.buf.push(self.cur as u8);
            self.cur = 0;
            self.ncur = 0;
        }
    }

    #[inline]
    fn emit(&mut self, bit: u32) {
        self.put(bit);
        let inv = 1 - bit;
        while self.pending > 0 {
            self.put(inv);
            self.pending -= 1;
        }
    }

    /// Encode a symbol occupying [cum_lo, cum_hi) of [0, total).
    #[inline]
    pub fn encode(&mut self, cum_lo: u64, cum_hi: u64, total: u64) {
        let rng = self.high - self.low + 1; // <= 2^32
        // (rng*total)/total == rng and (rng*0)/total == 0 exactly.
        self.high = self.low
            + (if cum_hi == total { rng } else { rng * cum_hi / total })
            - 1;
        self.low += if cum_lo == 0 { 0 } else { rng * cum_lo / total };
        loop {
            if self.high < AC_HALF {
                self.emit(0);
            } else if self.low >= AC_HALF {
                self.emit(1);
                self.low -= AC_HALF;
                self.high -= AC_HALF;
            } else if self.low >= AC_Q1 && self.high < AC_Q3 {
                self.pending += 1;
                self.low -= AC_Q1;
                self.high -= AC_Q1;
            } else {
                break;
            }
            self.low <<= 1;
            self.high = (self.high << 1) | 1;
        }
    }

    /// Flush (Definition 10.3(b)): terminate and byte-align the block.
    pub fn finish(&mut self) {
        self.pending += 1;
        let bit = if self.low < AC_Q1 { 0 } else { 1 };
        self.emit(bit);
        while self.ncur != 0 {
            self.put(0);
        }
    }
}

pub struct ADec<'a> {
    low: u64,
    high: u64,
    value: u64,
    data: &'a [u8],
    bitpos: usize,
}

impl<'a> ADec<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        let mut d = ADec { low: 0, high: AC_MASK, value: 0, data, bitpos: 0 };
        for _ in 0..32 {
            let bit = d.get() as u64;
            d.value = (d.value << 1) | bit;
        }
        d
    }

    /// Read the next bit; reads 0 past the end of the block.
    #[inline]
    fn get(&mut self) -> u32 {
        let i = self.bitpos >> 3;
        let bit = if i < self.data.len() {
            ((self.data[i] >> (7 - (self.bitpos & 7))) & 1) as u32
        } else {
            0
        };
        self.bitpos += 1;
        bit
    }

    /// Scaled position of the coded point inside the current interval.
    #[inline]
    pub fn target(&self, total: u64) -> u64 {
        let rng = self.high - self.low + 1;
        ((self.value - self.low + 1) * total - 1) / rng
    }

    /// Advance past a symbol occupying [cum_lo, cum_hi) of [0, total).
    #[inline]
    pub fn consume(&mut self, cum_lo: u64, cum_hi: u64, total: u64) {
        let rng = self.high - self.low + 1;
        self.high = self.low
            + (if cum_hi == total { rng } else { rng * cum_hi / total })
            - 1;
        self.low += if cum_lo == 0 { 0 } else { rng * cum_lo / total };
        loop {
            if self.high < AC_HALF {
                // nothing
            } else if self.low >= AC_HALF {
                self.low -= AC_HALF;
                self.high -= AC_HALF;
                self.value -= AC_HALF;
            } else if self.low >= AC_Q1 && self.high < AC_Q3 {
                self.low -= AC_Q1;
                self.high -= AC_Q1;
                self.value -= AC_Q1;
            } else {
                break;
            }
            self.low <<= 1;
            self.high = (self.high << 1) | 1;
            let bit = self.get() as u64;
            self.value = (self.value << 1) | bit;
        }
    }
}

// ---------------------------------------------------------------------
// Context entry (counts + derived top-16 list and bucket sums)
// ---------------------------------------------------------------------

/// Context entry.  `counts`/`total` are the reference's adaptive state;
/// the top-16 list (count desc, byte asc -- identical to the reference's
/// stable argsort prefix) and the 16 bucket sums are derived caches
/// maintained exactly (see `bump`).
#[derive(Clone)]
pub struct Ctx {
    total: u32,
    c16: [u16; 16],
    t16: [u8; 16],
    bucket: [u16; 16],
    counts: [u16; 256],
}

impl Ctx {
    /// Q_init image: zero counts; t16 = 0..15 mirrors the stable argsort
    /// of the all-zero count vector (byte-ascending tie-break).
    fn zeroed() -> Self {
        let mut t16 = [0u8; 16];
        for (i, t) in t16.iter_mut().enumerate() {
            *t = i as u8;
        }
        Ctx { total: 0, c16: [0; 16], t16, bucket: [0; 16], counts: [0; 256] }
    }

    /// Full stable top-16 selection over all 256 counts: order by
    /// (count desc, byte asc), identical to np.argsort(-freq,
    /// kind="stable")[:16].  Used at cap-halving only.
    fn rebuild_top16(&mut self) {
        let mut c16 = [0u16; 16];
        let mut t16 = [0u8; 16];
        let mut nt = 0usize;
        for b in 0..256usize {
            let v = self.counts[b];
            if nt == 16 {
                if v <= c16[15] {
                    continue; // tie: smaller byte stays
                }
                let mut i = 15;
                while i > 0 && v > c16[i - 1] {
                    c16[i] = c16[i - 1];
                    t16[i] = t16[i - 1];
                    i -= 1;
                }
                c16[i] = v;
                t16[i] = b as u8;
            } else {
                let mut i = nt;
                while i > 0 && v > c16[i - 1] {
                    c16[i] = c16[i - 1];
                    t16[i] = t16[i - 1];
                    i -= 1;
                }
                c16[i] = v;
                t16[i] = b as u8;
                nt += 1;
            }
        }
        self.c16 = c16;
        self.t16 = t16;
    }

    /// Cap-halving (reference: counts >>= 1; total = counts.sum()).
    /// Equal counts can collapse under halving, which reorders byte-asc
    /// ties, so the top-16 list and buckets are rebuilt from scratch.
    /// `wrapped` repairs the one transient u16 overflow case (true
    /// count 2^16, all mass on one byte at the cap).
    fn halve(&mut self, x: u8, wrapped: bool) {
        if wrapped {
            self.counts[x as usize] = 0; // placeholder; fixed below
        }
        for b in 0..256usize {
            self.counts[b] >>= 1;
        }
        if wrapped {
            self.counts[x as usize] = 32768; // 65536 >> 1
        }
        let mut s = 0u32;
        for k in 0..16usize {
            let mut bs = 0u32;
            for j in 0..16usize {
                bs += self.counts[16 * k + j] as u32;
            }
            self.bucket[k] = bs as u16;
            s += bs;
        }
        self.total = s;
        self.rebuild_top16();
    }

    /// The reference's per-context update: counts[x] += 1, total += 1,
    /// halve at total >= 2^16; maintain the derived caches exactly.
    fn bump(&mut self, x: u8) {
        let xi = x as usize;
        let old = self.counts[xi];
        self.counts[xi] = old.wrapping_add(1);
        self.total += 1;
        self.bucket[xi >> 4] = self.bucket[xi >> 4].wrapping_add(1);
        if self.total >= COUNT_CAP {
            // rare; rebuilds caches
            self.halve(x, old == 0xffff);
            return;
        }
        let v = old + 1;
        // Early outs.  Any member of the top-16 has cached count >=
        // c16[15], so x (old count v-1) can be a member only when
        // v > c16[15]; and an outside byte enters iff it now precedes
        // the 16th element in the total order (count desc, byte asc).
        let mut pos: usize;
        if v < self.c16[15] {
            return;
        }
        if v == self.c16[15] {
            // x is outside the list; tie-break entry
            if x > self.t16[15] {
                return;
            }
            pos = 15;
        } else {
            match find16(&self.t16, x) {
                Some(p) => pos = p,
                None => pos = 15, // v > c16[15]: always enters
            }
        }
        // Re-insert (v, x) into the sorted prefix: shift down while x
        // precedes its predecessor.
        while pos > 0
            && (v > self.c16[pos - 1]
                || (v == self.c16[pos - 1] && x < self.t16[pos - 1]))
        {
            self.c16[pos] = self.c16[pos - 1];
            self.t16[pos] = self.t16[pos - 1];
            pos -= 1;
        }
        self.c16[pos] = v;
        self.t16[pos] = x;
    }
}

/// Position of byte x in the 16-byte list, or None.
#[inline]
fn find16(t16: &[u8; 16], x: u8) -> Option<usize> {
    t16.iter().position(|&t| t == x)
}

/// freq[b] = 1 + 32 * count[b] (the +1 probability floor).
#[inline]
fn fr(count: u16) -> u64 {
    1 + COUNT_SCALE * count as u64
}

/// Sum of the 15 candidate counts c16[1..16].
#[inline]
fn cand_count_sum(c: &Ctx) -> u32 {
    c.c16[1..].iter().map(|&v| v as u32).sum()
}

/// Sum of freq over ESC bytes strictly below x (ascending-byte prefix).
fn esc_prefix(c: &Ctx, x: u8) -> u64 {
    let xi = x as usize;
    let kb = xi >> 4;
    let mut cnt = 0u32;
    for k in 0..kb {
        cnt += c.bucket[k] as u32;
    }
    for b in (kb << 4)..xi {
        cnt += c.counts[b] as u32;
    }
    let mut pref = xi as u64 + COUNT_SCALE * cnt as u64;
    for i in 0..16 {
        if c.t16[i] < x {
            pref -= fr(c.c16[i]);
        }
    }
    pref
}

// ---------------------------------------------------------------------
// Order-W adaptive predictor (exact-key open-addressing context tables)
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Default)]
struct Slot {
    key: u64,
    gen: u32,
    idx: u32,
}

const KMASK: [u64; 9] = [
    0x0,
    0xff,
    0xffff,
    0xff_ffff,
    0xffff_ffff,
    0xff_ffff_ffff,
    0xffff_ffff_ffff,
    0xff_ffff_ffff_ffff,
    0xffff_ffff_ffff_ffff,
];

pub struct Pred {
    w: u32,
    slots: Vec<Vec<Slot>>, // per order 0..=W; mask = len-1
    used: [u32; 9],        // live entries this generation
    pool: Vec<Ctx>,
    pool_head: u32,
    gen: u32,
    roll: u64,    // last <=8 bytes, newest in LSB
    histlen: u32, // min(bytes since sync, W)
}

#[inline]
fn key_hash(key: u64) -> u32 {
    let h = key.wrapping_mul(0x9e37_79b9_7f4a_7c15);
    ((h >> 32) as u32) ^ (h as u32)
}

impl Pred {
    pub fn new(w: u32) -> Self {
        let mut slots = Vec::with_capacity(w as usize + 1);
        for o in 0..=w {
            let nslots = if o == 0 { 16 } else if o == 1 { 1024 } else { 8192 };
            slots.push(vec![Slot::default(); nslots]);
        }
        Pred {
            w,
            slots,
            used: [0; 9],
            pool: Vec::with_capacity(4096),
            pool_head: 0,
            gen: 0,
            roll: 0,
            histlen: 0,
        }
    }

    /// Reset to Q_init at a sync point (Definition 10.3(a)): O(1) via
    /// generation tags; slots whose gen != current gen are empty.
    pub fn reset(&mut self) {
        self.gen = self.gen.wrapping_add(1);
        if self.gen == 0 {
            // wrapped: hard-clear and restart
            for tab in self.slots.iter_mut() {
                for s in tab.iter_mut() {
                    *s = Slot::default();
                }
            }
            self.gen = 1;
        }
        self.used = [0; 9];
        self.pool_head = 0;
        self.roll = 0;
        self.histlen = 0;
    }

    /// Double a table, re-inserting only current-generation entries.
    /// Table layout never influences coded output.
    fn grow(&mut self, o: usize) {
        let old = std::mem::take(&mut self.slots[o]);
        let new_n = old.len() * 2;
        let mut ns = vec![Slot::default(); new_n];
        let nmask = new_n - 1;
        for s in old {
            if s.gen != self.gen {
                continue;
            }
            let mut j = key_hash(s.key) as usize & nmask;
            while ns[j].gen == self.gen {
                j = (j + 1) & nmask;
            }
            ns[j] = s;
        }
        self.slots[o] = ns;
    }

    /// Exact-key get-or-create; returns pool index.  Entries created
    /// here start with total == 0 and are therefore invisible to the
    /// backoff rule (which requires total > 0), exactly like absent
    /// dict entries in the reference.
    fn insert(&mut self, o: usize, key: u64) -> u32 {
        if (self.used[o] as usize + 1) * 2 > self.slots[o].len() {
            self.grow(o);
        }
        let mask = self.slots[o].len() - 1;
        let mut i = key_hash(key) as usize & mask;
        loop {
            let s = self.slots[o][i];
            if s.gen != self.gen {
                let idx = self.pool_head;
                if (idx as usize) < self.pool.len() {
                    self.pool[idx as usize] = Ctx::zeroed();
                } else {
                    self.pool.push(Ctx::zeroed());
                }
                self.pool_head += 1;
                self.used[o] += 1;
                self.slots[o][i] = Slot { key, gen: self.gen, idx };
                return idx;
            }
            if s.key == key {
                return s.idx;
            }
            i = (i + 1) & mask;
        }
    }

    /// One shared table pass per position: get-or-create the context
    /// entry of every order 0..avail, then pick the longest entry with
    /// observations (backoff; equals the reference's dist-then-update
    /// semantics).  Returns (avail, chosen pool index or None=Q_init).
    #[inline]
    pub fn lookup(&mut self, idx: &mut [u32; 9]) -> (u32, Option<u32>) {
        let avail = self.histlen.min(self.w);
        for o in 0..=avail as usize {
            idx[o] = self.insert(o, self.roll & KMASK[o]);
        }
        for o in (0..=avail as usize).rev() {
            if self.pool[idx[o] as usize].total > 0 {
                return (avail, Some(idx[o]));
            }
        }
        (avail, None)
    }

    #[inline]
    pub fn ctx_ref<'s>(&'s self, found: &Option<u32>, zero: &'s Ctx) -> &'s Ctx {
        match found {
            Some(i) => &self.pool[*i as usize],
            None => zero,
        }
    }

    /// Deterministic integer count update (Definition 10.1 semantics).
    #[inline]
    pub fn update(&mut self, x: u8, avail: u32, idx: &[u32; 9]) {
        for o in 0..=avail as usize {
            self.pool[idx[o] as usize].bump(x);
        }
        self.roll = (self.roll << 8) | x as u64;
        if self.histlen < self.w {
            self.histlen += 1;
        }
    }
}

// ---------------------------------------------------------------------
// Sub-block encode / decode (Type-I class/payload, Section 4.3)
//
// Classes C1/C2 are structurally empty (E12-A min distance 3 and both
// codewords valid), so the exact partition sums reduce to
// ccum = [0, c0, c0, c0, c0+ccand, total], as the reference computes
// them via np.add.at over the classid map.
// ---------------------------------------------------------------------

pub fn encode_subblock(block: &[u8], pred: &mut Pred, enc: &mut AEnc) {
    pred.reset();
    enc.reset();
    let zero = Ctx::zeroed();
    let mut idx = [0u32; 9];
    for &x in block {
        let (avail, found) = pred.lookup(&mut idx);
        let c = pred.ctx_ref(&found, &zero);
        let total = 256 + COUNT_SCALE * c.total as u64;
        let c0 = fr(c.c16[0]);
        let ccand = NUM_CAND + COUNT_SCALE * cand_count_sum(c) as u64;
        match find16(&c.t16, x) {
            Some(0) => {
                enc.encode(0, c0, total);
            }
            Some(rank) => {
                enc.encode(c0, c0 + ccand, total);
                let mut mlo = 0u64;
                for i in 1..rank {
                    mlo += fr(c.c16[i]);
                }
                enc.encode(mlo, mlo + fr(c.c16[rank]), ccand);
            }
            None => {
                let cesc = total - c0 - ccand;
                enc.encode(c0 + ccand, total, total);
                let pref = esc_prefix(c, x);
                enc.encode(pref, pref + fr(c.counts[x as usize]), cesc);
            }
        }
        pred.update(x, avail, &idx);
    }
    enc.finish();
}

pub fn decode_subblock(
    blob: &[u8],
    nbytes: usize,
    pred: &mut Pred,
    out: &mut [u8],
) -> Result<(), String> {
    pred.reset();
    let mut dec = ADec::new(blob);
    let zero = Ctx::zeroed();
    let mut idx = [0u32; 9];
    for pos in 0..nbytes {
        let (avail, found) = pred.lookup(&mut idx);
        let c = pred.ctx_ref(&found, &zero);
        let total = 256 + COUNT_SCALE * c.total as u64;
        let c0 = fr(c.c16[0]);
        let ccand = NUM_CAND + COUNT_SCALE * cand_count_sum(c) as u64;
        let t = dec.target(total);
        let x: u8;
        if t < c0 {
            dec.consume(0, c0, total);
            x = c.t16[0];
        } else if t < c0 + ccand {
            dec.consume(c0, c0 + ccand, total);
            let t2 = dec.target(ccand);
            let mut cum = 0u64;
            let mut rank = 1usize;
            while rank < 16 {
                let f = fr(c.c16[rank]);
                if cum + f > t2 {
                    break;
                }
                cum += f;
                rank += 1;
            }
            if rank == 16 {
                return Err("internal error".into());
            }
            dec.consume(cum, cum + fr(c.c16[rank]), ccand);
            x = c.t16[rank];
        } else {
            let cesc = total - c0 - ccand;
            dec.consume(c0 + ccand, total, total);
            let t2 = dec.target(cesc);
            // Bitmap of the 16 special bytes + per-bucket corrections.
            let mut spec = [0u64; 4];
            let mut bcor = [0u32; 16];
            for i in 0..16 {
                let b = c.t16[i] as usize;
                spec[b >> 6] |= 1u64 << (b & 63);
                bcor[b >> 4] += fr(c.c16[i]) as u32;
            }
            let mut cum = 0u64;
            let mut k = 0usize;
            while k < 16 {
                let besc = 16 + COUNT_SCALE * c.bucket[k] as u64 - bcor[k] as u64;
                if cum + besc > t2 {
                    break;
                }
                cum += besc;
                k += 1;
            }
            if k == 16 {
                return Err("internal error".into());
            }
            let mut b = k << 4;
            let bend = (k + 1) << 4;
            let mut fx = 0u64;
            while b < bend {
                if (spec[b >> 6] >> (b & 63)) & 1 == 1 {
                    b += 1;
                    continue;
                }
                fx = fr(c.counts[b]);
                if cum + fx > t2 {
                    break;
                }
                cum += fx;
                b += 1;
            }
            if b == bend {
                return Err("internal error".into());
            }
            dec.consume(cum, cum + fx, cesc);
            x = b as u8;
        }
        out[pos] = x;
        pred.update(x, avail, &idx);
    }
    Ok(())
}
