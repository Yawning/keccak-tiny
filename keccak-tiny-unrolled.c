/** libkeccak-tiny
 *
 * A single-file implementation of SHA-3 and SHAKE.
 *
 * Implementor: David Leon Gil
 * License: CC0, attribution kindly requested. Blame taken too,
 * but not liability.
 */
#include "keccak-tiny.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/******** The Keccak-f[1600] permutation ********/

/*** Constants. ***/
static const uint8_t rho[24] = \
  { 1,  3,   6, 10, 15, 21,
    28, 36, 45, 55,  2, 14,
    27, 41, 56,  8, 25, 43,
    62, 18, 39, 61, 20, 44};
static const uint8_t pi[24] = \
  {10,  7, 11, 17, 18, 3,
    5, 16,  8, 21, 24, 4,
   15, 23, 19, 13, 12, 2,
   20, 14, 22,  9, 6,  1};
static const uint64_t RC[24] = \
  {1ULL, 0x8082ULL, 0x800000000000808aULL, 0x8000000080008000ULL,
   0x808bULL, 0x80000001ULL, 0x8000000080008081ULL, 0x8000000000008009ULL,
   0x8aULL, 0x88ULL, 0x80008009ULL, 0x8000000aULL,
   0x8000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL, 0x8000000000008003ULL,
   0x8000000000008002ULL, 0x8000000000000080ULL, 0x800aULL, 0x800000008000000aULL,
   0x8000000080008081ULL, 0x8000000000008080ULL, 0x80000001ULL, 0x8000000080008008ULL};

/*** Helper macros to unroll the permutation. ***/
#define rol(x, s) (((x) << s) | ((x) >> (64 - s)))
#define REPEAT6(e) e e e e e e
#define REPEAT24(e) REPEAT6(e e e e)
#define REPEAT5(e) e e e e e
#define FOR5(v, s, e) \
  v = 0;            \
  REPEAT5(e; v += s;)

/*** Keccak-f[1600] ***/
static inline void keccakf(void* state) {
  uint64_t* a = (uint64_t*)state;
  uint64_t b[5] = {0};
  uint64_t t = 0;
  uint8_t x, y, i = 0;

  REPEAT24(
      // Theta
      FOR5(x, 1,
           b[x] = 0;
           FOR5(y, 5,
                b[x] ^= a[x + y]; ))
      FOR5(x, 1,
           FOR5(y, 5,
                a[y + x] ^= b[(x + 4) % 5] ^ rol(b[(x + 1) % 5], 1); ))
      // Rho and pi
      t = a[1];
      x = 0;
      REPEAT24(b[0] = a[pi[x]];
               a[pi[x]] = rol(t, rho[x]);
               t = b[0];
               x++; )
      // Chi
      FOR5(y,
         5,
         FOR5(x, 1,
              b[x] = a[y + x];)
         FOR5(x, 1,
              a[y + x] = b[x] ^ ((~b[(x + 1) % 5]) & b[(x + 2) % 5]); ))
      // Iota
      a[0] ^= RC[i];
      i++; )
}

/******** The FIPS202-defined functions. ********/

/*** Some helper macros. ***/

#define _(S) do { S } while (0)
#define FOR(i, ST, L, S) \
  _(for (size_t i = 0; i < L; i += ST) { S; })
#define mkapply_ds(NAME, S)                                          \
  static inline void NAME(uint8_t* dst,                              \
                          const uint8_t* src,                        \
                          size_t len) {                              \
    FOR(i, 1, len, S);                                               \
  }
#define mkapply_sd(NAME, S)                                          \
  static inline void NAME(const uint8_t* src,                        \
                          uint8_t* dst,                              \
                          size_t len) {                              \
    FOR(i, 1, len, S);                                               \
  }

mkapply_ds(xorin, dst[i] ^= src[i])  // xorin
mkapply_sd(setout, dst[i] = src[i])  // setout

#define P keccakf
#define Plen 200

// Fold P*F over the full blocks of an input.
#define foldP(I, L, F) \
  while (L >= s->rate) {  \
    F(s->a, I, s->rate);  \
    P(s->a);              \
    I += s->rate;         \
    L -= s->rate;         \
  }

typedef struct keccak_state_s {
  uint8_t a[Plen];
  size_t rate;
  uint8_t delim;

  uint8_t block[Plen];
  size_t offset;
  uint8_t finalized : 1;
} keccak_state;

static inline void keccak_init(keccak_state *s,
                               size_t rate, uint8_t delim)
{
  memset_s(s, sizeof(*s), 0, sizeof(*s));
  s->rate = rate;
  s->delim = delim;
}

static inline void keccak_absorb_blocks(keccak_state *s,
                                        const uint8_t *in, size_t nr_blocks)
{
  size_t inlen = nr_blocks * s->rate;
  foldP(in, inlen, xorin);
}

static int keccak_update(keccak_state *s,
                         const uint8_t *in, size_t inlen) {
  size_t remaining;

  if ((s->finalized) || ((in == NULL) && inlen != 0))
    return -1;

  for (remaining = inlen; remaining > 0; ) {
    size_t buf_avail;
    size_t buf_bytes;

    /* Process full blocks directly if possible. */
    if (s->offset == 0) {
      size_t blocks = remaining / s->rate;
      size_t direct_bytes = blocks * s->rate;
      if (direct_bytes > 0) {
        keccak_absorb_blocks(s, in, blocks);
        remaining -= direct_bytes;
        in += direct_bytes;
      }
    }

    /* Buffer up to 1 block worth of data... */
    buf_avail = s->rate - s->offset;
    buf_bytes = (buf_avail > remaining) ? remaining : buf_avail;
    if (buf_bytes > 0) {
      memcpy(&s->block[s->offset], in, buf_bytes);
      s->offset += buf_bytes;
      remaining -= buf_bytes;
      in += buf_bytes;
    }
    if (s->offset == s->rate) { /* ... and process it. */
      keccak_absorb_blocks(s, s->block, 1);
      s->offset = 0;
    }
  }
  return 0;
}

static int keccak_finalize(keccak_state *s) {
  if (s->finalized)
    return -1;

  // Xor in the DS and pad frame.
  s->a[s->offset] ^= s->delim;
  s->a[s->rate - 1] ^= 0x80;
  // Xor in the last block.
  xorin(s->a, s->block, s->offset);

  // Update bookkeeping to prepare for squeezing.
  s->finalized = 1;
  s->offset = s->rate;
  return 0;
}

static inline void keccak_squeeze_blocks(keccak_state *s,
                                         uint8_t *out, size_t nr_blocks)
{
  size_t i, rate;
  for (i = 0, rate = s->rate; i < nr_blocks; i++) {
    keccakf(s->a);
    setout(s->a, out, rate);
    out += rate;
  }
}

static int keccak_squeeze(keccak_state *s,
                          uint8_t *out, size_t outlen)
{
  size_t remaining;
  if (!s->finalized)
    return -1;

  for (remaining = outlen; remaining > 0; ) {
    if (s->offset == s->rate) {
      /* Process full blocks directly if possible */
      size_t blocks = remaining / s->rate;
      size_t direct_bytes = blocks * s->rate;
      if (blocks > 0) {
        keccak_squeeze_blocks(s, out, blocks);
        out += direct_bytes;
        remaining -= direct_bytes;
      }

      /* Squeeze out another block into the internal buffer. */
      if (remaining > 0) {
        keccak_squeeze_blocks(s, s->block, 1);
        s->offset = 0;
      }
    }

    /* If there's a (partial) buffered block, drain it. */
    size_t buf_bytes = s->rate - s->offset;
    size_t indirect_bytes = (buf_bytes > remaining) ? remaining : buf_bytes;
    if (indirect_bytes > 0) {
      memcpy(out, &s->block[s->offset], indirect_bytes);
      out += indirect_bytes;
      s->offset += indirect_bytes;
      remaining -= indirect_bytes;
    }
  }
  return 0;
}

/** The sponge-based hash construction. **/
static inline int hash(uint8_t* out, size_t outlen,
                       const uint8_t* in, size_t inlen,
                       size_t rate, uint8_t delim) {
  keccak_state s;
  int ret = 0;

  if ((out == NULL) || ((in == NULL) && inlen != 0) || (rate >= Plen)) {
    return -1;
  }

  keccak_init(&s, rate, delim);
  ret |= keccak_update(&s, in, inlen);
  ret |= keccak_finalize(&s);
  ret |= keccak_squeeze(&s, out, outlen);
  memset_s(&s, sizeof(s), 0, sizeof(s));

  return ret;
}


/*** Helper macros to define SHA3 and SHAKE instances. ***/
#define defshake(bits)                                            \
  int shake##bits(uint8_t* out, size_t outlen,                    \
                  const uint8_t* in, size_t inlen) {              \
    return hash(out, outlen, in, inlen, 200 - (bits / 4), 0x1f);  \
  }                                                               \
  shake_##bits##_ctx *shake_##bits##_init(void) {                 \
    shake_##bits##_ctx *ctx = malloc(sizeof(shake_##bits##_ctx)); \
    keccak_init(ctx, 200 - (bits / 4), 0x1f);                     \
    return ctx;                                                   \
  }                                                               \
  int shake_##bits##_absorb(shake_##bits##_ctx *s,                \
                            const uint8_t *in, size_t inlen) {    \
    return keccak_update(s, in, inlen);                           \
  }                                                               \
  int shake_##bits##_squeeze(shake_##bits##_ctx *s,               \
                             uint8_t *out, size_t outlen) {       \
    int ret = 0;                                                  \
    if (!s->finalized)                                            \
      ret |= keccak_finalize(s);                                  \
    ret |= keccak_squeeze(s, out, outlen);                        \
    return ret;                                                   \
  }                                                               \
  void shake_##bits##_free(shake_##bits##_ctx *s) {               \
    memset_s(s, sizeof(*s), 0, sizeof(*s));                       \
    free(s);                                                      \
  }

#define defsha3(bits) \
  int sha3_##bits(uint8_t* out, size_t outlen,                    \
                  const uint8_t* in, size_t inlen) {              \
    if (outlen > (bits/8)) {                                      \
      return -1;                                                  \
    }                                                             \
    return hash(out, outlen, in, inlen, 200 - (bits / 4), 0x06);  \
  }                                                               \
  sha3_##bits##_ctx *sha3_##bits##_init(void) {                   \
    sha3_##bits##_ctx *ctx = malloc(sizeof(sha3_##bits##_ctx));   \
    keccak_init(ctx, 200 - (bits / 4), 0x06);                     \
    return ctx;                                                   \
  }                                                               \
  int sha3_##bits##_update(sha3_##bits##_ctx *s,                  \
                         const uint8_t *in, size_t inlen) {       \
    return keccak_update(s, in, inlen);                           \
  }                                                               \
  int sha3_##bits##_sum(const sha3_##bits##_ctx *s,               \
                       uint8_t *out, size_t outlen) {             \
    sha3_##bits##_ctx tmp;                                        \
    int ret = 0;                                                  \
    if (outlen > (bits / 8))                                      \
      return -1;                                                  \
    memcpy(&tmp, s, sizeof(tmp));                                 \
    ret |= keccak_finalize(&tmp);                                 \
    ret |= keccak_squeeze(&tmp, out, outlen);                     \
    memset_s(&tmp, sizeof(tmp), 0, sizeof(tmp));                  \
    return ret;                                                   \
  }                                                               \
  void sha3_##bits##_free(sha3_##bits##_ctx *s) {                 \
    memset_s(s, sizeof(*s), 0, sizeof(*s));                       \
    free(s);                                                      \
  }

/*** FIPS202 SHAKE VOFs ***/
defshake(128)
defshake(256)

/*** FIPS202 SHA3 FOFs ***/
defsha3(224)
defsha3(256)
defsha3(384)
defsha3(512)
