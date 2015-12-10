#ifndef KECCAK_FIPS202_H
#define KECCAK_FIPS202_H
#define __STDC_WANT_LIB_EXT1__ 1
#include <stdint.h>
#include <stdlib.h>

struct keccak_state_s;

#define decshake(bits) \
  typedef struct keccak_state_s shake_##bits##_ctx; \
  shake_##bits##_ctx *shake_##bits##_init(void); \
  int shake_##bits##_absorb(shake_##bits##_ctx *, const uint8_t *, size_t); \
  int shake_##bits##_squeeze(shake_##bits##_ctx *, uint8_t*, size_t); \
  void shake_##bits##_free(shake_##bits##_ctx *); \
  int shake##bits(uint8_t*, size_t, const uint8_t*, size_t);

#define decsha3(bits) \
  typedef struct keccak_state_s sha3_##bits##_ctx; \
  sha3_##bits##_ctx *sha3_##bits##_init(void); \
  int sha3_##bits##_update(sha3_##bits##_ctx *, const uint8_t *, size_t); \
  int sha3_##bits##_sum(const sha3_##bits##_ctx *, uint8_t*, size_t); \
  void sha3_##bits##_free(sha3_##bits##_ctx *); \
  int sha3_##bits(uint8_t*, size_t, const uint8_t*, size_t);

decshake(128)
decshake(256)
decsha3(224)
decsha3(256)
decsha3(384)
decsha3(512)
#endif
