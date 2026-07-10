/* rnr1_fast.c -- fast C implementation of the RNR1 Type-I coder.
 *
 * This is a performance port of the normative Python reference
 * impl/rnr1.py (format: impl/FORMAT.md).  The reference stays normative:
 * this implementation is required to produce archives BYTE-IDENTICAL to
 * rnr1.pack() and to decode exactly the archives rnr1.py accepts
 * (verified by impl_fast/run_checks.py against the full data/ corpus).
 *
 * Ported components, each bit-exact by construction:
 *
 *  1. Code12 E12-A systematic (12,8) map (Section 4.2 parity equations).
 *     At init we verify the min nonzero pairwise mask weight is 3
 *     (linear code, min distance 3), which proves classes C1/C2 are
 *     structurally empty for the byte-aligned argmax predictor; the
 *     class partition below exploits exactly that (and only that).
 *
 *  2. Order-W adaptive byte-context predictor: per-order EXACT-key
 *     context tables (open addressing, full 64-bit key compare -- the
 *     reference uses exact dict keys, so no lossy hashing is allowed),
 *     uint16 counts with the +1 floor applied at read time
 *     (freq = 1 + 32*count, total = 256 + 32*ctx_total), deterministic
 *     backoff longest-context-first, cap-halving at ctx_total >= 2^16.
 *     Adaptive state is reset to Q_init at every sub-block boundary
 *     (generation tags, O(1) reset).
 *
 *  3. 32-bit binary arithmetic coder, Witten-Neal-Cleary renormalization
 *     with pending-bit carry handling, MSB-first bit packing, identical
 *     termination (2 bits + carries, zero-padded to byte boundary).
 *     All intermediates fit in uint64: range <= 2^32, totals < 2^24.
 *
 *  4. Type-I class/payload partition (Section 4.3): top-16 ordering
 *     (freq desc, byte value asc on ties) identical to the reference's
 *     stable argsort, exact partition sums, candidate-rank payload for
 *     CAND, ascending-byte payload for ESC.
 *
 *  5. Container per impl/FORMAT.md: 48-byte header, 32-byte seek-index
 *     entries, per-sub-block truncated SHA-256 (CommonCrypto), model
 *     description hash, fail-safe raw store mode when the coded blob is
 *     not smaller than the raw sub-block.
 *
 * Performance structure (bit-identity preserving):
 *   - Each context entry carries its top-16 list (count desc, byte asc)
 *     maintained INCREMENTALLY under +1 count updates.  Exactness: only
 *     one key changes per update, so re-inserting that key into the
 *     sorted prefix reproduces the full stable sort's prefix; an outside
 *     byte enters iff it now precedes the current 16th element in the
 *     total order.  On cap-halving (which can reorder equal-count ties)
 *     the list is rebuilt by a full stable selection.
 *   - Each entry carries 16 bucket sums (16 bytes per bucket) so the
 *     ESC-class prefix sums need not scan all 256 counts.
 *   - dist() and update() share one exact-key table pass per position;
 *     entries created before coding have total == 0 and are therefore
 *     invisible to the backoff rule (which requires total > 0), exactly
 *     like absent dict entries in the reference.
 *   - uint16 counts: any per-context count is <= its total < 2^16 at
 *     rest; the single transient 2^16 case (all mass on one byte at the
 *     cap) wraps to 0 and is repaired to 32768 inside the halving step,
 *     which fires in the same update.
 *
 * Integer determinism: the whole coding path is integer-only (uint16/
 * uint32/uint64), no floating point, no randomness, fixed evaluation
 * order.
 *
 * Build (macOS/arm64):
 *   clang -O3 -std=c11 -Wall -Wextra -fPIC -shared \
 *         -o librnr1fast.dylib rnr1_fast.c
 */

#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <pthread.h>

#include <CommonCrypto/CommonDigest.h>

#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

/* ------------------------------------------------------------------ */
/* Format constants (impl/FORMAT.md)                                   */
/* ------------------------------------------------------------------ */

#define VER_MAJOR 1
#define VER_MINOR 0
#define VER_PATCH 0
#define FLAG_SUBBLOCK_HASHES 0x0001u

#define COUNT_SCALE 32u
#define COUNT_CAP 65536u
#define NUM_CAND 15

#define HEADER_SIZE 48u
#define INDEX_ENTRY_SIZE 32u
#define MODE_CODED 0u
#define MODE_RAW 1u

/* Error codes (negative), mapped to messages by rnr1_strerror. */
#define RNR1_OK 0
#define RNR1_EBADW (-1)
#define RNR1_EBADK (-2)
#define RNR1_ENOMEM (-3)
#define RNR1_ETRUNC (-4)
#define RNR1_EMAGIC (-5)
#define RNR1_EVERSION (-6)
#define RNR1_EMODEL (-7)
#define RNR1_ECOUNT (-8)
#define RNR1_EINDEX (-9)
#define RNR1_EBLOCKHASH (-10)
#define RNR1_EARCHHASH (-11)
#define RNR1_EUNSUP (-12)
#define RNR1_EINTERNAL (-13)

/* ------------------------------------------------------------------ */
/* Little-endian scalar I/O                                            */
/* ------------------------------------------------------------------ */

static inline void wr16(uint8_t *p, uint16_t v) {
    p[0] = (uint8_t)v; p[1] = (uint8_t)(v >> 8);
}
static inline void wr32(uint8_t *p, uint32_t v) {
    p[0] = (uint8_t)v; p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16); p[3] = (uint8_t)(v >> 24);
}
static inline void wr64(uint8_t *p, uint64_t v) {
    for (int i = 0; i < 8; i++) p[i] = (uint8_t)(v >> (8 * i));
}
static inline uint16_t rd16(const uint8_t *p) {
    return (uint16_t)(p[0] | ((uint16_t)p[1] << 8));
}
static inline uint32_t rd32(const uint8_t *p) {
    return p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16)
         | ((uint32_t)p[3] << 24);
}
static inline uint64_t rd64(const uint8_t *p) {
    uint64_t v = 0;
    for (int i = 7; i >= 0; i--) v = (v << 8) | p[i];
    return v;
}

/* ------------------------------------------------------------------ */
/* SHA-256 helpers (truncated to 8 bytes, Definition 3.1 H_v)          */
/* ------------------------------------------------------------------ */

static void sha256_trunc8(const uint8_t *data, uint64_t len, uint8_t out[8]) {
    CC_SHA256_CTX c;
    CC_SHA256_Init(&c);
    while (len > 0) {
        CC_LONG chunk = len > (1u << 30) ? (1u << 30) : (CC_LONG)len;
        CC_SHA256_Update(&c, data, chunk);
        data += chunk;
        len -= chunk;
    }
    uint8_t d[CC_SHA256_DIGEST_LENGTH];
    CC_SHA256_Final(d, &c);
    memcpy(out, d, 8);
}

/* ------------------------------------------------------------------ */
/* Code12 E12-A table (Section 4.2) + one-time self-check              */
/* ------------------------------------------------------------------ */

static uint16_t E12[256];
static int g_init_done = 0;
static int g_init_rc = RNR1_EINTERNAL;

static uint8_t g_model_hash_cache[9][8]; /* per W in 0..8 */
static uint8_t g_model_hash_valid[9];

/* Context entry.  counts/total are the reference's adaptive state;
 * top-16 list and bucket sums are derived caches (see header comment). */
typedef struct {
    uint32_t total;         /* sum of counts, < 2^16 at rest            */
    uint16_t c16[16];       /* top-16 counts, order (count desc, b asc) */
    uint8_t t16[16];        /* top-16 byte values                       */
    uint16_t bucket[16];    /* bucket[k] = sum counts[16k .. 16k+15]    */
    uint8_t pad_[44];       /* header = 128B when pool is 64B-aligned   */
    uint16_t counts[256];
} Ctx;                      /* 640 bytes                                */

static Ctx g_zero_ctx;      /* Q_init image: zero counts, t16 = 0..15   */

static int popcount12(uint32_t v) { return __builtin_popcount(v & 0xfffu); }

static int rnr1_init(void) {
    if (g_init_done) return g_init_rc;
    for (int b = 0; b < 256; b++) {
        int bits[8], c[13];
        for (int j = 0; j < 8; j++) bits[j] = (b >> j) & 1;
        memset(c, 0, sizeof c);
        c[3] = bits[0]; c[5] = bits[1]; c[6] = bits[2]; c[7] = bits[3];
        c[9] = bits[4]; c[10] = bits[5]; c[11] = bits[6]; c[12] = bits[7];
        c[1] = c[3] ^ c[5] ^ c[7] ^ c[9] ^ c[11];
        c[2] = c[3] ^ c[6] ^ c[7] ^ c[10] ^ c[11];
        c[4] = c[5] ^ c[6] ^ c[7] ^ c[12];
        c[8] = c[9] ^ c[10] ^ c[11] ^ c[12];
        uint16_t w = 0;
        for (int i = 1; i <= 12; i++) w |= (uint16_t)(c[i] << (i - 1));
        E12[b] = w;
    }
    /* Self-check: E12-A is linear with min distance 3, so every nonzero
     * pairwise mask has weight >= 3.  This is what makes classes C1/C2
     * structurally empty; the class partition below relies on it. */
    for (int a = 0; a < 256; a++)
        for (int b = a + 1; b < 256; b++) {
            int w = popcount12((uint32_t)(E12[a] ^ E12[b]));
            if (w < 3) { g_init_done = 1; return g_init_rc; }
        }
    memset(&g_zero_ctx, 0, sizeof g_zero_ctx);
    for (int i = 0; i < 16; i++) g_zero_ctx.t16[i] = (uint8_t)i;
    g_init_rc = RNR1_OK;
    g_init_done = 1;
    return g_init_rc;
}

/* Truncated SHA-256 of the canonical model description (R-13.3);
 * byte-identical to rnr1.model_description_hash(W). */
static void model_hash(uint32_t W, uint8_t out[8]) {
    if (W <= 8 && g_model_hash_valid[W]) {
        memcpy(out, g_model_hash_cache[W], 8);
        return;
    }
    /* E12 serialized as 256 little-endian int64 values (numpy "<i8"). */
    uint8_t e12b[256 * 8];
    memset(e12b, 0, sizeof e12b);
    for (int i = 0; i < 256; i++) wr64(e12b + 8 * i, E12[i]);
    uint8_t d[CC_SHA256_DIGEST_LENGTH];
    CC_SHA256(e12b, sizeof e12b, d);
    char hex[17];
    for (int i = 0; i < 8; i++)
        snprintf(hex + 2 * i, 3, "%02x", d[i]);
    char desc[192];
    int len = snprintf(desc, sizeof desc,
        "rnr-type1-code12a;predictor=ngram;W=%u;scale=%u;cap=%u;cand=%d;"
        "floor=1;classes=c0,c1,c2,cand,esc;e12=%s",
        W, COUNT_SCALE, COUNT_CAP, NUM_CAND, hex);
    CC_SHA256(desc, (CC_LONG)len, d);
    memcpy(out, d, 8);
    if (W <= 8) {
        memcpy(g_model_hash_cache[W], d, 8);
        g_model_hash_valid[W] = 1;
    }
}

/* ------------------------------------------------------------------ */
/* Binary arithmetic coder (32-bit registers, WNC renormalization)     */
/* ------------------------------------------------------------------ */

#define AC_MASK 0xffffffffull
#define AC_HALF 0x80000000ull
#define AC_Q1 0x40000000ull
#define AC_Q3 0xc0000000ull

typedef struct {
    uint64_t low, high;     /* invariants: low <= high <= AC_MASK */
    uint64_t pending;
    uint8_t *buf;
    size_t len, cap;
    uint32_t cur;
    int ncur;
    int oom;
} AEnc;

static void aenc_reset(AEnc *e) {
    e->low = 0; e->high = AC_MASK; e->pending = 0;
    e->len = 0; e->cur = 0; e->ncur = 0; e->oom = 0;
}

static void aenc_grow(AEnc *e) {
    size_t ncap = e->cap ? e->cap * 2 : 65536;
    uint8_t *nb = (uint8_t *)realloc(e->buf, ncap);
    if (!nb) { e->oom = 1; return; }
    e->buf = nb; e->cap = ncap;
}

static inline void aenc_put(AEnc *e, int bit) {
    e->cur = (e->cur << 1) | (uint32_t)bit;
    if (++e->ncur == 8) {
        if (e->len == e->cap) { aenc_grow(e); if (e->oom) return; }
        e->buf[e->len++] = (uint8_t)e->cur;
        e->cur = 0; e->ncur = 0;
    }
}

static inline void aenc_emit(AEnc *e, int bit) {
    aenc_put(e, bit);
    int inv = 1 - bit;
    while (e->pending) { aenc_put(e, inv); e->pending--; }
}

static inline void aenc_encode(AEnc *e, uint64_t cum_lo, uint64_t cum_hi,
                               uint64_t total) {
    uint64_t rng = e->high - e->low + 1;      /* <= 2^32 */
    /* (rng*total)/total == rng and (rng*0)/total == 0 exactly; skipping
     * the division in those cases changes nothing but the cycle count. */
    e->high = e->low + (cum_hi == total ? rng : (rng * cum_hi) / total) - 1;
    e->low = e->low + (cum_lo == 0 ? 0 : (rng * cum_lo) / total);
    for (;;) {
        if (e->high < AC_HALF) {
            aenc_emit(e, 0);
        } else if (e->low >= AC_HALF) {
            aenc_emit(e, 1);
            e->low -= AC_HALF; e->high -= AC_HALF;
        } else if (e->low >= AC_Q1 && e->high < AC_Q3) {
            e->pending++;
            e->low -= AC_Q1; e->high -= AC_Q1;
        } else {
            break;
        }
        e->low <<= 1;
        e->high = (e->high << 1) | 1;
    }
}

static void aenc_finish(AEnc *e) {
    e->pending++;
    aenc_emit(e, e->low < AC_Q1 ? 0 : 1);
    while (e->ncur) aenc_put(e, 0);
}

typedef struct {
    uint64_t low, high, value;
    const uint8_t *data;
    size_t len;
    size_t bitpos;
} ADec;

static inline int adec_get(ADec *d) {
    size_t i = d->bitpos >> 3;
    int bit = 0;
    if (i < d->len) bit = (d->data[i] >> (7 - (d->bitpos & 7))) & 1;
    d->bitpos++;
    return bit;
}

static void adec_init(ADec *d, const uint8_t *data, size_t len) {
    d->low = 0; d->high = AC_MASK; d->value = 0;
    d->data = data; d->len = len; d->bitpos = 0;
    for (int i = 0; i < 32; i++) d->value = (d->value << 1) | (uint64_t)adec_get(d);
}

static inline uint64_t adec_target(ADec *d, uint64_t total) {
    uint64_t rng = d->high - d->low + 1;
    return ((d->value - d->low + 1) * total - 1) / rng;
}

static inline void adec_consume(ADec *d, uint64_t cum_lo, uint64_t cum_hi,
                                uint64_t total) {
    uint64_t rng = d->high - d->low + 1;
    d->high = d->low + (cum_hi == total ? rng : (rng * cum_hi) / total) - 1;
    d->low = d->low + (cum_lo == 0 ? 0 : (rng * cum_lo) / total);
    for (;;) {
        if (d->high < AC_HALF) {
            /* nothing */
        } else if (d->low >= AC_HALF) {
            d->low -= AC_HALF; d->high -= AC_HALF; d->value -= AC_HALF;
        } else if (d->low >= AC_Q1 && d->high < AC_Q3) {
            d->low -= AC_Q1; d->high -= AC_Q1; d->value -= AC_Q1;
        } else {
            break;
        }
        d->low <<= 1;
        d->high = (d->high << 1) | 1;
        d->value = (d->value << 1) | (uint64_t)adec_get(d);
    }
}

/* ------------------------------------------------------------------ */
/* Context entry maintenance (top-16 list + bucket sums)               */
/* ------------------------------------------------------------------ */

/* Full stable top-16 selection over all 256 counts: order by
 * (count desc, byte asc), identical to np.argsort(-freq, kind="stable")
 * [:16].  Used at cap-halving only; the hot path is incremental. */
static void ctx_rebuild_top16(Ctx *c) {
    uint16_t c16[16];
    uint8_t t16[16];
    int nt = 0;
    for (int b = 0; b < 256; b++) {
        uint16_t v = c->counts[b];
        if (nt == 16) {
            if (v <= c16[15]) continue;   /* tie: smaller byte stays */
            int i = 15;
            while (i > 0 && v > c16[i - 1]) {
                c16[i] = c16[i - 1]; t16[i] = t16[i - 1]; i--;
            }
            c16[i] = v; t16[i] = (uint8_t)b;
        } else {
            int i = nt;
            while (i > 0 && v > c16[i - 1]) {
                c16[i] = c16[i - 1]; t16[i] = t16[i - 1]; i--;
            }
            c16[i] = v; t16[i] = (uint8_t)b; nt++;
        }
    }
    memcpy(c->c16, c16, sizeof c16);
    memcpy(c->t16, t16, sizeof t16);
}

/* Cap-halving (reference: counts >>= 1; total = counts.sum()).  Equal
 * counts can collapse under halving, which reorders byte-asc ties, so
 * the top-16 list and buckets are rebuilt from scratch.  `wrapped`
 * repairs the one transient uint16 overflow case (true count 2^16). */
static void ctx_halve(Ctx *c, uint8_t x, int wrapped) {
    if (wrapped) c->counts[x] = 0;      /* placeholder; fixed below */
    uint32_t s = 0;
    for (int b = 0; b < 256; b++) {
        uint16_t hv = (uint16_t)(c->counts[b] >> 1);
        c->counts[b] = hv;
    }
    if (wrapped) c->counts[x] = 32768;  /* 65536 >> 1 */
    for (int k = 0; k < 16; k++) {
        uint32_t bs = 0;
        for (int j = 0; j < 16; j++) bs += c->counts[16 * k + j];
        c->bucket[k] = (uint16_t)bs;
        s += bs;
    }
    c->total = s;
    ctx_rebuild_top16(c);
}

/* Position of byte x in a 16-byte list, or -1. */
static inline int find16(const uint8_t *t16, uint8_t x) {
#if defined(__ARM_NEON)
    uint8x16_t eq = vceqq_u8(vld1q_u8(t16), vdupq_n_u8(x));
    uint8x8_t nib = vshrn_n_u16(vreinterpretq_u16_u8(eq), 4);
    uint64_t mask = vget_lane_u64(vreinterpret_u64_u8(nib), 0);
    if (mask == 0) return -1;
    return __builtin_ctzll(mask) >> 2;
#else
    for (int i = 0; i < 16; i++)
        if (t16[i] == x) return i;
    return -1;
#endif
}

/* Apply the reference's per-context update: counts[x] += 1, total += 1,
 * halve at total >= 2^16; maintain the derived caches exactly. */
static inline void ctx_bump(Ctx *c, uint8_t x) {
    uint16_t old = c->counts[x];
    c->counts[x] = (uint16_t)(old + 1);
    c->total++;
    c->bucket[x >> 4]++;
    if (c->total >= COUNT_CAP) {        /* rare; rebuilds caches */
        ctx_halve(c, x, old == 0xffff);
        return;
    }
    uint16_t v = (uint16_t)(old + 1);
    uint16_t *c16 = c->c16;
    uint8_t *t16 = c->t16;
    /* Early outs.  Any member of the top-16 has cached count >= c16[15],
     * so x (old count v-1) can be a member only when v > c16[15]; and an
     * outside byte enters iff it now precedes the 16th element in the
     * total order (count desc, byte asc). */
    int pos;
    if (v < c16[15]) return;
    if (v == c16[15]) {                 /* x is outside; tie-break entry */
        if (x > t16[15]) return;
        pos = 15;
    } else {
        pos = find16(t16, x);
        if (pos < 0) pos = 15;          /* v > c16[15]: always enters */
    }
    /* Re-insert (v, x) into the sorted prefix: shift down while x
     * precedes its predecessor. */
    int i = pos;
    while (i > 0 && (v > c16[i - 1] || (v == c16[i - 1] && x < t16[i - 1]))) {
        c16[i] = c16[i - 1]; t16[i] = t16[i - 1]; i--;
    }
    c16[i] = v; t16[i] = x;
}

/* ------------------------------------------------------------------ */
/* Order-W adaptive byte-context predictor (exact-key context tables)  */
/* ------------------------------------------------------------------ */

typedef struct {
    uint64_t key;
    uint32_t gen;
    uint32_t idx;
} Slot;

typedef struct {
    uint32_t W;
    Slot *slots[9];
    uint32_t mask[9];                   /* slot-count-1 per order */
    uint32_t used[9];                   /* live entries this generation */
    Ctx *pool;
    uint32_t pool_head, pool_cap;
    uint32_t gen;
    uint64_t roll;                      /* last <=8 bytes, newest in LSB */
    uint32_t histlen;                   /* min(bytes since sync, W) */
    int oom;
} Pred;

static const uint64_t KMASK[9] = {
    0x0ull, 0xffull, 0xffffull, 0xffffffull, 0xffffffffull,
    0xffffffffffull, 0xffffffffffffull, 0xffffffffffffffull,
    0xffffffffffffffffull,
};

static int pred_create(Pred *p, uint32_t W) {
    memset(p, 0, sizeof *p);
    p->W = W;
    for (uint32_t o = 0; o <= W; o++) {
        uint32_t nslots = o == 0 ? 16 : (o == 1 ? 1024 : 8192);
        p->slots[o] = (Slot *)calloc(nslots, sizeof(Slot));
        if (!p->slots[o]) return RNR1_ENOMEM;
        p->mask[o] = nslots - 1;
    }
    p->pool_cap = 4096;
    if (posix_memalign((void **)&p->pool, 64,
                       (size_t)p->pool_cap * sizeof(Ctx)) != 0)
        return RNR1_ENOMEM;
    p->gen = 0;
    return RNR1_OK;
}

static void pred_destroy(Pred *p) {
    for (int o = 0; o < 9; o++) free(p->slots[o]);
    free(p->pool);
    memset(p, 0, sizeof *p);
}

/* Reset to Q_init at a sync point (Definition 10.3(a)): O(1) via
 * generation tags; slots whose gen != current gen are empty. */
static void pred_reset(Pred *p) {
    p->gen++;
    if (p->gen == 0) {              /* wrapped: hard-clear and restart */
        for (uint32_t o = 0; o <= p->W; o++)
            memset(p->slots[o], 0, ((size_t)p->mask[o] + 1) * sizeof(Slot));
        p->gen = 1;
    }
    memset(p->used, 0, sizeof p->used);
    p->pool_head = 0;
    p->roll = 0;
    p->histlen = 0;
}

static inline uint32_t key_hash(uint64_t key) {
    uint64_t h = key * 0x9e3779b97f4a7c15ull;
    return (uint32_t)(h >> 32) ^ (uint32_t)h;
}

/* Double a table, re-inserting only current-generation entries.  Table
 * layout never influences coded output, so growth points are irrelevant
 * to bit-identity; growth is itself deterministic anyway. */
static int tab_grow(Pred *p, uint32_t o) {
    uint32_t old_n = p->mask[o] + 1;
    uint64_t new_n = (uint64_t)old_n * 2;
    if (new_n > (1u << 31)) return RNR1_ENOMEM;
    Slot *ns = (Slot *)calloc((size_t)new_n, sizeof(Slot));
    if (!ns) return RNR1_ENOMEM;
    uint32_t nmask = (uint32_t)new_n - 1;
    for (uint32_t i = 0; i < old_n; i++) {
        Slot *s = &p->slots[o][i];
        if (s->gen != p->gen) continue;
        uint32_t j = key_hash(s->key) & nmask;
        while (ns[j].gen == p->gen) j = (j + 1) & nmask;
        ns[j] = *s;
    }
    free(p->slots[o]);
    p->slots[o] = ns;
    p->mask[o] = nmask;
    return RNR1_OK;
}

static int pool_grow(Pred *p) {
    uint64_t ncap = (uint64_t)p->pool_cap * 2;
    if (ncap > (1u << 31)) return RNR1_ENOMEM;
    Ctx *np;
    if (posix_memalign((void **)&np, 64, (size_t)ncap * sizeof(Ctx)) != 0)
        return RNR1_ENOMEM;
    memcpy(np, p->pool, (size_t)p->pool_head * sizeof(Ctx));
    free(p->pool);
    p->pool = np;
    p->pool_cap = (uint32_t)ncap;
    return RNR1_OK;
}

/* Exact-key get-or-create; returns pool INDEX (the pool may move when
 * it grows, so callers resolve pointers only after all inserts). */
static inline uint32_t tab_insert(Pred *p, uint32_t o, uint64_t key) {
    if ((uint64_t)(p->used[o] + 1) * 2 > p->mask[o] + 1) {
        if (tab_grow(p, o) != RNR1_OK) { p->oom = 1; return 0; }
    }
    Slot *slots = p->slots[o];
    uint32_t mask = p->mask[o];
    uint32_t i = key_hash(key) & mask;
    for (;;) {
        Slot *s = &slots[i];
        if (s->gen != p->gen) {
            if (p->pool_head >= p->pool_cap) {
                if (pool_grow(p) != RNR1_OK) { p->oom = 1; return 0; }
                slots = p->slots[o];    /* unchanged, but keep tidy */
            }
            s->gen = p->gen;
            s->key = key;
            s->idx = p->pool_head++;
            p->used[o]++;
            memcpy(&p->pool[s->idx], &g_zero_ctx, sizeof(Ctx));
            return s->idx;
        }
        if (s->key == key) return s->idx;
        i = (i + 1) & mask;
    }
}

/* One shared table pass per position: get-or-create the context entry
 * of every order 0..avail (creation is invisible to the backoff rule,
 * which requires total > 0, so this equals the reference's dist-then-
 * update semantics), then pick the longest entry with observations. */
static inline const Ctx *pred_lookup(Pred *p, uint32_t *avail_out,
                                     uint32_t idx[9]) {
    uint32_t avail = p->histlen < p->W ? p->histlen : p->W;
    *avail_out = avail;
    /* Warm the slot cache lines of every order before probing. */
    for (uint32_t o = 0; o <= avail; o++) {
        uint64_t key = p->roll & KMASK[o];
        __builtin_prefetch(&p->slots[o][key_hash(key) & p->mask[o]], 0, 1);
    }
    for (uint32_t o = 0; o <= avail; o++) {
        idx[o] = tab_insert(p, o, p->roll & KMASK[o]);
        /* Warm the entry for the coding read and the coming update. */
        __builtin_prefetch(&p->pool[idx[o]], 1, 2);
    }
    if (p->oom) return NULL;
    for (int o = (int)avail; o >= 0; o--) {
        Ctx *c = &p->pool[idx[o]];
        if (c->total > 0) return c;
    }
    return &g_zero_ctx;
}

/* Deterministic integer count update (Definition 10.1 semantics). */
static inline void pred_update(Pred *p, uint8_t x, uint32_t avail,
                               const uint32_t idx[9]) {
    for (uint32_t o = 0; o <= avail; o++)
        ctx_bump(&p->pool[idx[o]], x);
    p->roll = (p->roll << 8) | x;
    if (p->histlen < p->W) p->histlen++;
}

/* ------------------------------------------------------------------ */
/* Type-I class/payload coding (Section 4.3)                           */
/*                                                                     */
/* Classes C1/C2 are structurally empty (init self-check: E12-A min    */
/* distance 3 and both codewords valid), so the exact partition sums   */
/* reduce to ccum = [0, c0, c0, c0, c0+ccand, total], as the reference */
/* computes them via np.add.at over the classid map.                   */
/* ------------------------------------------------------------------ */

/* freq[b] = 1 + 32*counts[b]; class frequencies from the cached top16. */
static inline uint64_t fr(uint32_t count) {
    return 1 + (uint64_t)COUNT_SCALE * count;
}

/* Sum of the 15 candidate counts c16[1..15]. */
static inline uint32_t cand_count_sum(const Ctx *c) {
#if defined(__ARM_NEON)
    uint32_t s = vaddlvq_u16(vld1q_u16(c->c16))
               + vaddlvq_u16(vld1q_u16(c->c16 + 8));
    return s - c->c16[0];
#else
    uint32_t s = 0;
    for (int i = 1; i < 16; i++) s += c->c16[i];
    return s;
#endif
}

/* Sum of freq over ESC bytes below x (ascending-byte prefix). */
static inline uint64_t esc_prefix(const Ctx *c, uint8_t x) {
    uint32_t cnt = 0;
    int kb = x >> 4;
    for (int k = 0; k < kb; k++) cnt += c->bucket[k];
    for (int b = kb << 4; b < x; b++) cnt += c->counts[b];
    uint64_t pref = (uint64_t)x + (uint64_t)COUNT_SCALE * cnt;
    for (int i = 0; i < 16; i++)
        if (c->t16[i] < x) pref -= fr(c->c16[i]);
    return pref;
}

/* ------------------------------------------------------------------ */
/* Sub-block encode / decode                                           */
/* ------------------------------------------------------------------ */

static int encode_subblock(const uint8_t *block, uint64_t blen,
                           Pred *p, AEnc *e) {
    pred_reset(p);
    aenc_reset(e);
    uint32_t idx[9];
    for (uint64_t pos = 0; pos < blen; pos++) {
        uint8_t x = block[pos];
        uint32_t avail;
        const Ctx *c = pred_lookup(p, &avail, idx);
        if (!c) return RNR1_ENOMEM;
        uint64_t total = 256 + (uint64_t)COUNT_SCALE * c->total;
        uint64_t c0 = fr(c->c16[0]);
        uint64_t ccand = NUM_CAND
            + (uint64_t)COUNT_SCALE * cand_count_sum(c);
        int rank = find16(c->t16, x);
        if (rank == 0) {
            aenc_encode(e, 0, c0, total);
        } else if (rank > 0) {
            aenc_encode(e, c0, c0 + ccand, total);
            uint64_t mlo = 0;
            for (int i = 1; i < rank; i++) mlo += fr(c->c16[i]);
            aenc_encode(e, mlo, mlo + fr(c->c16[rank]), ccand);
        } else {
            uint64_t cesc = total - c0 - ccand;
            aenc_encode(e, c0 + ccand, total, total);
            uint64_t pref = esc_prefix(c, x);
            aenc_encode(e, pref, pref + fr(c->counts[x]), cesc);
        }
        pred_update(p, x, avail, idx);
        if (e->oom || p->oom) return RNR1_ENOMEM;
    }
    aenc_finish(e);
    if (e->oom) return RNR1_ENOMEM;
    return RNR1_OK;
}

static int decode_subblock(const uint8_t *blob, uint64_t blob_len,
                           uint64_t nbytes, Pred *p, uint8_t *out) {
    pred_reset(p);
    ADec d;
    adec_init(&d, blob, (size_t)blob_len);
    uint32_t idx[9];
    for (uint64_t pos = 0; pos < nbytes; pos++) {
        uint32_t avail;
        const Ctx *c = pred_lookup(p, &avail, idx);
        if (!c) return RNR1_ENOMEM;
        uint64_t total = 256 + (uint64_t)COUNT_SCALE * c->total;
        uint64_t c0 = fr(c->c16[0]);
        uint64_t ccand = NUM_CAND
            + (uint64_t)COUNT_SCALE * cand_count_sum(c);
        uint64_t t = adec_target(&d, total);
        uint8_t x;
        if (t < c0) {
            adec_consume(&d, 0, c0, total);
            x = c->t16[0];
        } else if (t < c0 + ccand) {
            adec_consume(&d, c0, c0 + ccand, total);
            uint64_t t2 = adec_target(&d, ccand);
            uint64_t cum = 0;
            int rank = 1;
            for (; rank < 16; rank++) {
                uint64_t f = fr(c->c16[rank]);
                if (cum + f > t2) break;
                cum += f;
            }
            if (rank == 16) return RNR1_EINTERNAL;
            adec_consume(&d, cum, cum + fr(c->c16[rank]), ccand);
            x = c->t16[rank];
        } else {
            uint64_t cesc = total - c0 - ccand;
            adec_consume(&d, c0 + ccand, total, total);
            uint64_t t2 = adec_target(&d, cesc);
            /* Bitmap of the 16 special bytes + per-bucket corrections. */
            uint64_t spec[4] = {0, 0, 0, 0};
            uint32_t bcor[16];
            memset(bcor, 0, sizeof bcor);
            for (int i = 0; i < 16; i++) {
                uint8_t b = c->t16[i];
                spec[b >> 6] |= 1ull << (b & 63);
                bcor[b >> 4] += (uint32_t)fr(c->c16[i]);
            }
            uint64_t cum = 0;
            int k = 0;
            for (; k < 16; k++) {
                uint64_t besc = 16 + (uint64_t)COUNT_SCALE * c->bucket[k]
                              - bcor[k];
                if (cum + besc > t2) break;
                cum += besc;
            }
            if (k == 16) return RNR1_EINTERNAL;
            int b = k << 4;
            int bend = (k + 1) << 4;
            uint64_t fx = 0;
            for (; b < bend; b++) {
                if ((spec[b >> 6] >> (b & 63)) & 1) continue;
                fx = fr(c->counts[b]);
                if (cum + fx > t2) break;
                cum += fx;
            }
            if (b == bend) return RNR1_EINTERNAL;
            adec_consume(&d, cum, cum + fx, cesc);
            x = (uint8_t)b;
        }
        out[pos] = x;
        pred_update(p, x, avail, idx);
        if (p->oom) return RNR1_ENOMEM;
    }
    return RNR1_OK;
}

/* ------------------------------------------------------------------ */
/* Container: pack / unpack (impl/FORMAT.md)                           */
/* ------------------------------------------------------------------ */

const char *rnr1_strerror(int64_t rc) {
    switch (rc) {
    case RNR1_OK: return "ok";
    case RNR1_EBADW: return "W must be in [1, 8]";
    case RNR1_EBADK: return "K must be >= 256";
    case RNR1_ENOMEM: return "out of memory";
    case RNR1_ETRUNC: return "truncated archive";
    case RNR1_EMAGIC: return "bad magic: not an RNR1 archive";
    case RNR1_EVERSION: return "unsupported archive format version";
    case RNR1_EMODEL: return "model description hash mismatch";
    case RNR1_ECOUNT: return "inconsistent sub-block count";
    case RNR1_EINDEX: return "seek index invalid";
    case RNR1_EBLOCKHASH: return "sub-block hash mismatch";
    case RNR1_EARCHHASH: return "archive verification failed: H_v mismatch";
    case RNR1_EUNSUP: return "unsupported archive parameters";
    default: return "internal error";
    }
}

void rnr1_free(uint8_t *ptr) { free(ptr); }

int64_t rnr1_pack(const uint8_t *data, uint64_t n, uint32_t W, uint32_t K,
                  uint8_t **out, uint64_t *out_len) {
    int rc = rnr1_init();
    if (rc != RNR1_OK) return rc;
    if (W < 1 || W > 8) return RNR1_EBADW;
    if (K < 256) return RNR1_EBADK;

    uint64_t m = n ? (n + K - 1) / K : 0;

    Pred pred;
    rc = pred_create(&pred, W);
    if (rc != RNR1_OK) { pred_destroy(&pred); return rc; }
    AEnc enc;
    memset(&enc, 0, sizeof enc);

    uint8_t *entries = NULL;
    if (m) {
        entries = (uint8_t *)malloc((size_t)m * INDEX_ENTRY_SIZE);
        if (!entries) { pred_destroy(&pred); return RNR1_ENOMEM; }
    }
    /* Repair stream accumulator. */
    size_t rep_cap = (size_t)(n / 2 + 4096);
    uint8_t *rep = (uint8_t *)malloc(rep_cap);
    size_t rep_len = 0;
    if (!rep) { free(entries); pred_destroy(&pred); return RNR1_ENOMEM; }

    uint64_t repair_bit_offset = 0;
    for (uint64_t j = 0; j < m; j++) {
        const uint8_t *block = data + j * K;
        uint64_t blen = (j + 1) * K <= n ? K : n - j * K;
        rc = encode_subblock(block, blen, &pred, &enc);
        if (rc != RNR1_OK) goto fail;
        const uint8_t *blob = enc.buf;
        uint64_t blob_len = enc.len;
        uint32_t mode = MODE_CODED;
        if (blob_len >= blen) {         /* fail-safe store mode */
            blob = block;
            blob_len = blen;
            mode = MODE_RAW;
        }
        uint8_t *ent = entries + (size_t)j * INDEX_ENTRY_SIZE;
        wr64(ent, j * (uint64_t)K);     /* byte_offset */
        wr64(ent + 8, repair_bit_offset);
        wr32(ent + 16, 0);              /* snapshot_length */
        wr32(ent + 20, mode);
        sha256_trunc8(block, blen, ent + 24);
        while (rep_len + blob_len > rep_cap) {
            rep_cap = rep_cap * 2 + 65536;
            uint8_t *nr = (uint8_t *)realloc(rep, rep_cap);
            if (!nr) { rc = RNR1_ENOMEM; goto fail; }
            rep = nr;
        }
        memcpy(rep + rep_len, blob, blob_len);
        rep_len += blob_len;
        repair_bit_offset += 8 * blob_len;
    }

    {
        uint64_t total_len = HEADER_SIZE + m * INDEX_ENTRY_SIZE + rep_len;
        uint8_t *arch = (uint8_t *)malloc((size_t)total_len);
        if (!arch) { rc = RNR1_ENOMEM; goto fail; }
        uint8_t *h = arch;
        memcpy(h, "RNR1", 4);
        h[4] = VER_MAJOR; h[5] = VER_MINOR; h[6] = VER_PATCH;
        h[7] = 0;                       /* reserved */
        wr16(h + 8, FLAG_SUBBLOCK_HASHES);
        wr16(h + 10, (uint16_t)W);
        wr32(h + 12, K);
        wr64(h + 16, n);
        model_hash(W, h + 24);
        sha256_trunc8(data, n, h + 32); /* v = H_v(X) */
        wr64(h + 40, m);
        if (m) memcpy(arch + HEADER_SIZE, entries, (size_t)m * INDEX_ENTRY_SIZE);
        memcpy(arch + HEADER_SIZE + m * INDEX_ENTRY_SIZE, rep, rep_len);
        *out = arch;
        *out_len = total_len;
        rc = RNR1_OK;
    }
fail:
    free(rep);
    free(entries);
    free(enc.buf);
    pred_destroy(&pred);
    return rc;
}

/* ------------------------------------------------------------------ */
/* Multi-threaded pack/unpack over independent sub-blocks.             */
/*                                                                     */
/* Sub-blocks are self-contained restart points (fresh predictor and   */
/* coder state), so per-block work is a pure function of (block bytes, */
/* W); assembling results in block order yields output byte-identical  */
/* to the sequential path regardless of scheduling (verified by        */
/* run_checks F3b/F4b).                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    const uint8_t *data;
    uint64_t n, m;
    uint32_t W, K;
    _Atomic uint64_t next;
    _Atomic int err;
    uint8_t **blobs;            /* per block: malloc'd coded blob or NULL */
    uint64_t *blob_lens;
    uint8_t *entries;           /* 32B per block; bit offsets patched later */
} PackMt;

static void *pack_worker(void *arg) {
    PackMt *w = (PackMt *)arg;
    Pred pred;
    if (pred_create(&pred, w->W) != RNR1_OK) {
        atomic_store(&w->err, RNR1_ENOMEM);
        pred_destroy(&pred);
        return NULL;
    }
    AEnc enc;
    memset(&enc, 0, sizeof enc);
    for (;;) {
        uint64_t j = atomic_fetch_add(&w->next, 1);
        if (j >= w->m || atomic_load(&w->err) != RNR1_OK) break;
        const uint8_t *block = w->data + j * w->K;
        uint64_t blen = (j + 1) * (uint64_t)w->K <= w->n ? w->K
                                                         : w->n - j * w->K;
        int rc = encode_subblock(block, blen, &pred, &enc);
        if (rc != RNR1_OK) { atomic_store(&w->err, rc); break; }
        uint32_t mode = MODE_CODED;
        uint64_t blob_len = enc.len;
        if (blob_len >= blen) {         /* fail-safe store mode */
            mode = MODE_RAW;
            blob_len = blen;
            w->blobs[j] = NULL;         /* raw: use source directly */
        } else {
            uint8_t *b = (uint8_t *)malloc(blob_len ? blob_len : 1);
            if (!b) { atomic_store(&w->err, RNR1_ENOMEM); break; }
            memcpy(b, enc.buf, blob_len);
            w->blobs[j] = b;
        }
        w->blob_lens[j] = blob_len;
        uint8_t *ent = w->entries + (size_t)j * INDEX_ENTRY_SIZE;
        wr64(ent, j * (uint64_t)w->K);
        wr32(ent + 16, 0);
        wr32(ent + 20, mode);
        sha256_trunc8(block, blen, ent + 24);
    }
    free(enc.buf);
    pred_destroy(&pred);
    return NULL;
}

int64_t rnr1_pack_mt(const uint8_t *data, uint64_t n, uint32_t W, uint32_t K,
                     uint32_t nthreads, uint8_t **out, uint64_t *out_len);

int64_t rnr1_pack_mt(const uint8_t *data, uint64_t n, uint32_t W, uint32_t K,
                     uint32_t nthreads, uint8_t **out, uint64_t *out_len) {
    int rc = rnr1_init();
    if (rc != RNR1_OK) return rc;
    if (W < 1 || W > 8) return RNR1_EBADW;
    if (K < 256) return RNR1_EBADK;
    uint64_t m = n ? (n + K - 1) / K : 0;
    if (nthreads < 1) nthreads = 1;
    if (nthreads > 64) nthreads = 64;
    if (nthreads > m) nthreads = m ? (uint32_t)m : 1;
    if (nthreads <= 1) return rnr1_pack(data, n, W, K, out, out_len);

    PackMt w;
    memset(&w, 0, sizeof w);
    w.data = data; w.n = n; w.m = m; w.W = W; w.K = K;
    atomic_store(&w.next, 0);
    atomic_store(&w.err, RNR1_OK);
    w.blobs = (uint8_t **)calloc((size_t)m, sizeof(uint8_t *));
    w.blob_lens = (uint64_t *)calloc((size_t)m, sizeof(uint64_t));
    w.entries = (uint8_t *)malloc((size_t)m * INDEX_ENTRY_SIZE);
    if (!w.blobs || !w.blob_lens || !w.entries) { rc = RNR1_ENOMEM; goto done; }

    {
        pthread_t tid[64];
        uint32_t started = 0;
        for (uint32_t t = 0; t < nthreads; t++) {
            if (pthread_create(&tid[t], NULL, pack_worker, &w) != 0) break;
            started++;
        }
        if (started == 0) { rc = RNR1_ENOMEM; goto done; }
        for (uint32_t t = 0; t < started; t++) pthread_join(tid[t], NULL);
    }
    rc = atomic_load(&w.err);
    if (rc != RNR1_OK) goto done;

    {
        uint64_t rep_len = 0;
        for (uint64_t j = 0; j < m; j++) rep_len += w.blob_lens[j];
        uint64_t total_len = HEADER_SIZE + m * INDEX_ENTRY_SIZE + rep_len;
        uint8_t *arch = (uint8_t *)malloc((size_t)total_len);
        if (!arch) { rc = RNR1_ENOMEM; goto done; }
        uint8_t *h = arch;
        memcpy(h, "RNR1", 4);
        h[4] = VER_MAJOR; h[5] = VER_MINOR; h[6] = VER_PATCH;
        h[7] = 0;
        wr16(h + 8, FLAG_SUBBLOCK_HASHES);
        wr16(h + 10, (uint16_t)W);
        wr32(h + 12, K);
        wr64(h + 16, n);
        model_hash(W, h + 24);
        sha256_trunc8(data, n, h + 32);
        wr64(h + 40, m);
        uint8_t *rep = arch + HEADER_SIZE + m * INDEX_ENTRY_SIZE;
        uint64_t bit_off = 0;
        for (uint64_t j = 0; j < m; j++) {
            uint8_t *ent = w.entries + (size_t)j * INDEX_ENTRY_SIZE;
            wr64(ent + 8, bit_off);
            const uint8_t *src = w.blobs[j] ? w.blobs[j] : data + j * K;
            memcpy(rep, src, (size_t)w.blob_lens[j]);
            rep += w.blob_lens[j];
            bit_off += 8 * w.blob_lens[j];
        }
        memcpy(arch + HEADER_SIZE, w.entries, (size_t)m * INDEX_ENTRY_SIZE);
        *out = arch;
        *out_len = total_len;
        rc = RNR1_OK;
    }
done:
    if (w.blobs)
        for (uint64_t j = 0; j < m; j++) free(w.blobs[j]);
    free(w.blobs);
    free(w.blob_lens);
    free(w.entries);
    return rc;
}

typedef struct {
    const uint8_t *idx, *rep;
    uint64_t rep_len, n, m;
    uint32_t W, K;
    uint16_t flags;
    int verify;
    uint8_t *data;
    _Atomic uint64_t next;
    _Atomic int err;
} UnpackMt;

static void *unpack_worker(void *arg) {
    UnpackMt *w = (UnpackMt *)arg;
    Pred pred;
    if (pred_create(&pred, w->W) != RNR1_OK) {
        atomic_store(&w->err, RNR1_ENOMEM);
        pred_destroy(&pred);
        return NULL;
    }
    for (;;) {
        uint64_t j = atomic_fetch_add(&w->next, 1);
        if (j >= w->m || atomic_load(&w->err) != RNR1_OK) break;
        const uint8_t *ent = w->idx + j * INDEX_ENTRY_SIZE;
        uint64_t start = rd64(ent + 8) / 8;
        uint64_t end = (j + 1 < w->m)
            ? rd64(w->idx + INDEX_ENTRY_SIZE * (j + 1) + 8) / 8 : w->rep_len;
        if (start > w->rep_len) start = w->rep_len;
        if (end > w->rep_len) end = w->rep_len;
        if (end < start) end = start;
        uint64_t blk_len = (j + 1) * (uint64_t)w->K <= w->n
                         ? w->K : w->n - j * w->K;
        uint32_t mode = rd32(ent + 20);
        uint8_t *dst = w->data + j * w->K;
        int rc = RNR1_OK;
        if (mode == MODE_RAW) {
            if (end - start < blk_len) rc = RNR1_ETRUNC;
            else memcpy(dst, w->rep + start, (size_t)blk_len);
        } else {
            rc = decode_subblock(w->rep + start, end - start, blk_len,
                                 &pred, dst);
        }
        if (rc == RNR1_OK && w->verify && (w->flags & FLAG_SUBBLOCK_HASHES)) {
            uint8_t hv[8];
            sha256_trunc8(dst, blk_len, hv);
            if (memcmp(hv, ent + 24, 8) != 0) rc = RNR1_EBLOCKHASH;
        }
        if (rc != RNR1_OK) { atomic_store(&w->err, rc); break; }
    }
    pred_destroy(&pred);
    return NULL;
}

int64_t rnr1_unpack_mt(const uint8_t *raw, uint64_t raw_len, int verify,
                       uint32_t nthreads, uint8_t **out, uint64_t *out_len);

int64_t rnr1_unpack(const uint8_t *raw, uint64_t raw_len, int verify,
                    uint8_t **out, uint64_t *out_len) {
    int rc = rnr1_init();
    if (rc != RNR1_OK) return rc;
    if (raw_len < HEADER_SIZE) return RNR1_ETRUNC;
    if (memcmp(raw, "RNR1", 4) != 0) return RNR1_EMAGIC;
    if (raw[4] != VER_MAJOR || raw[5] > VER_MINOR) return RNR1_EVERSION;
    uint16_t flags = rd16(raw + 8);
    uint32_t W = rd16(raw + 10);
    uint32_t K = rd32(raw + 12);
    uint64_t n = rd64(raw + 16);
    uint64_t m = rd64(raw + 40);
    uint8_t mh[8];
    model_hash(W, mh);
    if (memcmp(mh, raw + 24, 8) != 0) return RNR1_EMODEL;
    /* The reference accepts any W whose model hash matches; our coder
     * supports the range the reference writer can produce. */
    if (W > 8) return RNR1_EUNSUP;
    if (K == 0) return RNR1_EUNSUP;
    uint64_t expect_m = n ? (n + K - 1) / K : 0;
    if (m != expect_m) return RNR1_ECOUNT;
    if (m > (raw_len - HEADER_SIZE) / INDEX_ENTRY_SIZE) return RNR1_ETRUNC;

    const uint8_t *idx = raw + HEADER_SIZE;
    const uint8_t *rep = idx + m * INDEX_ENTRY_SIZE;
    uint64_t rep_len = raw_len - HEADER_SIZE - m * INDEX_ENTRY_SIZE;

    for (uint64_t j = 1; j < m; j++) {
        if (rd64(idx + j * INDEX_ENTRY_SIZE) <=
            rd64(idx + (j - 1) * INDEX_ENTRY_SIZE))
            return RNR1_EINDEX;
    }
    for (uint64_t j = 0; j < m; j++) {
        if (rd64(idx + j * INDEX_ENTRY_SIZE + 8) % 8 != 0)
            return RNR1_EINDEX;
    }

    uint8_t *data = (uint8_t *)malloc(n ? (size_t)n : 1);
    if (!data) return RNR1_ENOMEM;

    Pred pred;
    rc = pred_create(&pred, W);
    if (rc != RNR1_OK) { pred_destroy(&pred); free(data); return rc; }

    for (uint64_t j = 0; j < m; j++) {
        const uint8_t *ent = idx + j * INDEX_ENTRY_SIZE;
        uint64_t bit_off = rd64(ent + 8);
        uint64_t start = bit_off / 8;
        uint64_t end = (j + 1 < m) ? rd64(idx + INDEX_ENTRY_SIZE * (j + 1) + 8) / 8
                                   : rep_len;
        if (start > rep_len) start = rep_len;   /* mimic slice clamping */
        if (end > rep_len) end = rep_len;
        if (end < start) end = start;
        uint64_t blk_len = (j + 1) * K <= n ? K : n - j * K;
        uint32_t mode = rd32(ent + 20);
        uint8_t *dst = data + j * K;
        if (mode == MODE_RAW) {
            if (end - start < blk_len) { rc = RNR1_ETRUNC; goto fail; }
            memcpy(dst, rep + start, (size_t)blk_len);
        } else {
            rc = decode_subblock(rep + start, end - start, blk_len, &pred, dst);
            if (rc != RNR1_OK) goto fail;
        }
        if (verify && (flags & FLAG_SUBBLOCK_HASHES)) {
            uint8_t hv[8];
            sha256_trunc8(dst, blk_len, hv);
            if (memcmp(hv, ent + 24, 8) != 0) { rc = RNR1_EBLOCKHASH; goto fail; }
        }
    }
    if (verify) {
        uint8_t hv[8];
        sha256_trunc8(data, n, hv);
        if (memcmp(hv, raw + 32, 8) != 0) { rc = RNR1_EARCHHASH; goto fail; }
    }
    *out = data;
    *out_len = n;
    pred_destroy(&pred);
    return RNR1_OK;
fail:
    free(data);
    pred_destroy(&pred);
    return rc;
}

int64_t rnr1_unpack_mt(const uint8_t *raw, uint64_t raw_len, int verify,
                       uint32_t nthreads, uint8_t **out, uint64_t *out_len) {
    int rc = rnr1_init();
    if (rc != RNR1_OK) return rc;
    if (raw_len < HEADER_SIZE) return RNR1_ETRUNC;
    if (memcmp(raw, "RNR1", 4) != 0) return RNR1_EMAGIC;
    if (raw[4] != VER_MAJOR || raw[5] > VER_MINOR) return RNR1_EVERSION;
    uint16_t flags = rd16(raw + 8);
    uint32_t W = rd16(raw + 10);
    uint32_t K = rd32(raw + 12);
    uint64_t n = rd64(raw + 16);
    uint64_t m = rd64(raw + 40);
    uint8_t mh[8];
    model_hash(W, mh);
    if (memcmp(mh, raw + 24, 8) != 0) return RNR1_EMODEL;
    if (W > 8) return RNR1_EUNSUP;
    if (K == 0) return RNR1_EUNSUP;
    uint64_t expect_m = n ? (n + K - 1) / K : 0;
    if (m != expect_m) return RNR1_ECOUNT;
    if (m > (raw_len - HEADER_SIZE) / INDEX_ENTRY_SIZE) return RNR1_ETRUNC;
    if (nthreads < 1) nthreads = 1;
    if (nthreads > 64) nthreads = 64;
    if (nthreads > m) nthreads = m ? (uint32_t)m : 1;
    if (nthreads <= 1) return rnr1_unpack(raw, raw_len, verify, out, out_len);

    UnpackMt w;
    memset(&w, 0, sizeof w);
    w.idx = raw + HEADER_SIZE;
    w.rep = w.idx + m * INDEX_ENTRY_SIZE;
    w.rep_len = raw_len - HEADER_SIZE - m * INDEX_ENTRY_SIZE;
    w.n = n; w.m = m; w.W = W; w.K = K;
    w.flags = flags; w.verify = verify;
    atomic_store(&w.next, 0);
    atomic_store(&w.err, RNR1_OK);

    for (uint64_t j = 1; j < m; j++) {
        if (rd64(w.idx + j * INDEX_ENTRY_SIZE) <=
            rd64(w.idx + (j - 1) * INDEX_ENTRY_SIZE))
            return RNR1_EINDEX;
    }
    for (uint64_t j = 0; j < m; j++) {
        if (rd64(w.idx + j * INDEX_ENTRY_SIZE + 8) % 8 != 0)
            return RNR1_EINDEX;
    }

    w.data = (uint8_t *)malloc(n ? (size_t)n : 1);
    if (!w.data) return RNR1_ENOMEM;

    {
        pthread_t tid[64];
        uint32_t started = 0;
        for (uint32_t t = 0; t < nthreads; t++) {
            if (pthread_create(&tid[t], NULL, unpack_worker, &w) != 0) break;
            started++;
        }
        if (started == 0) { free(w.data); return RNR1_ENOMEM; }
        for (uint32_t t = 0; t < started; t++) pthread_join(tid[t], NULL);
    }
    rc = atomic_load(&w.err);
    if (rc != RNR1_OK) { free(w.data); return rc; }
    if (verify) {
        uint8_t hv[8];
        sha256_trunc8(w.data, n, hv);
        if (memcmp(hv, raw + 32, 8) != 0) { free(w.data); return RNR1_EARCHHASH; }
    }
    *out = w.data;
    *out_len = n;
    return RNR1_OK;
}
