extern "C" {
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <execinfo.h>
#include <time.h>
#include <openssl/rand.h>
}

#include <stxxl/map>

#define __cplusplus

#define DATA_NODE_BLOCK_SIZE (4096)
#define DATA_LEAF_BLOCK_SIZE (4096)

#define SUB_BLOCK_SIZE (8192)
#define SUB_BLOCKS_PER_BLOCK (256)

#define SEL_CODE_LEN (56)
#define SEL_CODE_BYTES (SEL_CODE_LEN >> 3)
#define MAX_SELECTOR 6

#ifndef AQF_EXT_H
#define AQF_EXT_H

typedef struct ext_t {
  uint64_t bits; // store bits of extension
  int len;       // store which bits are valid; if len=0, then ext is empty
} Ext;

// Store 0.875 extension bits per element
#define EXT_CODE_LEN 56
#define EXT_CODE_BYTES (EXT_CODE_LEN >> 3)

#endif //AQF_EXT_H

/*
 *  Bit rank and select, taken from CQF (Pandey et al.)
 */

// Returns the trailing zero count of val
int tzcnt(uint64_t val) {
  asm("tzcnt %[val], %[val]"
  : [val] "+r" (val)
  :
  : "cc");
  return val;
}

// Returns the number of 1s in val
int popcnt(uint64_t val) {
  asm("popcnt %[val], %[val]"
      : [val] "+r" (val)
      :
      : "cc");
  return val;
}

// Returns the number of 1s up to (including) the pos'th bit, indexing from 0
uint64_t bitrank(uint64_t val, uint64_t pos) {
  assert(pos >= 0 && pos < 64 && "pos should be in [0, 63]");
  val &= ((2ull << pos) - 1); // 2^(pos+1) - 1
  asm("popcnt %[val], %[val]"
      : [val] "+r" (val)
      :
      : "cc");
  return val;
}

/**
 * Returns the position of the k-th 1 in the 64-bit word x.
 * k is 0-based, so k=0 returns the position of the first 1.
 *
 * Uses the broadword selection algorithm by Vigna [1], improved by Gog
 * and Petri [2] and Vigna [3].
 *
 * [1] Sebastiano Vigna. Broadword Implementation of Rank/Select
 *    Queries. WEA, 2008
 *
 * [2] Simon Gog, Matthias Petri. Optimized succinct data
 * structures for massive data. Softw. Pract. Exper., 2014
 *
 * [3] Sebastiano Vigna. MG4J 5.2.1. http://mg4j.di.unimi.it/
 * The following code is taken from
 * https://github.com/facebook/folly/blob/b28186247104f8b90cfbe094d289c91f9e413317/folly/experimental/Select64.h
 */
const uint8_t kSelectInByteAQF[2048] = {
                                     8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0,
                                     1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0,
                                     2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
                                     1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0,
                                     3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0,
                                     1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0,
                                     2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0,
                                     1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
                                     4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0,
                                     1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8, 8, 8, 1,
                                     8, 2, 2, 1, 8, 3, 3, 1, 3, 2, 2, 1, 8, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
                                     2, 1, 8, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1,
                                     4, 3, 3, 1, 3, 2, 2, 1, 8, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4,
                                     4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1,
                                     3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 7, 7, 1, 7, 2,
                                     2, 1, 7, 3, 3, 1, 3, 2, 2, 1, 7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
                                     7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3,
                                     3, 1, 3, 2, 2, 1, 7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1,
                                     4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2,
                                     2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 8, 8, 8, 8, 8, 8, 2,
                                     8, 8, 8, 3, 8, 3, 3, 2, 8, 8, 8, 4, 8, 4, 4, 2, 8, 4, 4, 3, 4, 3, 3, 2, 8, 8,
                                     8, 5, 8, 5, 5, 2, 8, 5, 5, 3, 5, 3, 3, 2, 8, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3,
                                     4, 3, 3, 2, 8, 8, 8, 6, 8, 6, 6, 2, 8, 6, 6, 3, 6, 3, 3, 2, 8, 6, 6, 4, 6, 4,
                                     4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 8, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2,
                                     6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 7, 8, 7, 7, 2, 8, 7,
                                     7, 3, 7, 3, 3, 2, 8, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, 8, 7, 7, 5,
                                     7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3,
                                     3, 2, 8, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6, 4, 4, 2,
                                     6, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6, 5,
                                     5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 3, 8, 8, 8, 8, 8, 8,
                                     8, 5, 8, 8, 8, 5, 8, 5, 5, 3, 8, 8, 8, 5, 8, 5, 5, 4, 8, 5, 5, 4, 5, 4, 4, 3,
                                     8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 3, 8, 8, 8, 6, 8, 6, 6, 4, 8, 6,
                                     6, 4, 6, 4, 4, 3, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5, 5, 3, 8, 6, 6, 5,
                                     6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7,
                                     7, 3, 8, 8, 8, 7, 8, 7, 7, 4, 8, 7, 7, 4, 7, 4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 5,
                                     8, 7, 7, 5, 7, 5, 5, 3, 8, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3, 8, 8,
                                     8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 3, 8, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4,
                                     6, 4, 4, 3, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5,
                                     5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 4, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6,
                                     6, 4, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 6, 8, 6, 6, 5,
                                     8, 6, 6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8,
                                     8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7,
                                     8, 7, 7, 5, 8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 4, 8, 8, 8, 8, 8, 8,
                                     8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 4,
                                     8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6,
                                     6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6,
                                     8, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
                                     8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 8,
                                     8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6,
                                     6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7
};

static inline uint64_t select64(uint64_t x, int k) {
  if (k >= popcnt(x)) { return 64; }

  const uint64_t kOnesStep4  = 0x1111111111111111ull;
  const uint64_t kOnesStep8  = 0x0101010101010101ull;
  const uint64_t kMSBsStep8  = 0x80ull * kOnesStep8;

  uint64_t s = x;
  s = s - ((s & 0xAu * kOnesStep4) >> 1u);
  s = (s & 0x3 * kOnesStep4) + ((s >> 2u) & 0x3 * kOnesStep4);
  s = (s + (s >> 4u)) & 0xF * kOnesStep8;
  uint64_t byteSums = s * kOnesStep8;

  uint64_t kStep8 = k * kOnesStep8;
  uint64_t geqKStep8 = (((kStep8 | kMSBsStep8) - byteSums) & kMSBsStep8);
  uint64_t place = popcnt(geqKStep8) * 8;
  uint64_t byteRank = k - (((byteSums << 8u) >> place) & (uint64_t)(0xFF));
  return place + kSelectInByteAQF[((x >> place) & 0xFFu) | (byteRank << 8u)];
}

// Returns the position of the rank'th 1.  (rank = 0 returns the 1st 1)
// Returns 64 if there are fewer than rank+1 1s.
uint64_t bitselect(uint64_t val, uint64_t rank) {
#ifdef __SSE4_2_
  uint64_t i = 1ull << rank; // i = 2^rank
  asm("pdep %[val], %[mask], %[val]"
      : [val] "+r" (val)
      : [mask] "r" (i));
  asm("tzcnt %[bit], %[index]"
      : [index] "=r" (i)
      : [bit] "g" (val)
      : "cc");
  return i;
#endif
  return select64(val, rank);
}

#ifndef MACROS_H
#define MACROS_H

#ifdef __cplusplus
extern "C" {
#endif

/** Set bit at i, indexing from LSB */
#define ONE(i) (1ULL << (i))

/** Get i ones from bit 0 to i-1, equivalent to MASK_CLOSED(0, i-1) */
#define ONES(i) (ONE(i)-1)

/** Bits in the closed interval [a,b] set, all others unset; a <= b */
#define MASK_CLOSED(a,b) (ONES((b)-(a)+1) << (a))

/** Bits in the half-open interval [a,b) set, all others unset; a <= b */
#define MASK_HALF_OPEN(a,b) (ONES((b)-(a)) << (a))

/* Miscellaneous */
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

/*
   Access metadata
   - occupieds, runends index from *LSB*
   - GET returns 0 or a nonzero value, not 0 or 1,
     but this works just fine for conditionals
*/
#define GET(bitarr, i) (ONE(i) & (bitarr))
#define SET(bitarr, i) ((bitarr) |= ONE(i))
#define UNSET(bitarr, i) ((bitarr) &= ~ONE(i))

// Get the block containing absolute index x
#define block_containing(filter, x) ((filter)->blocks[(x)/64])

/* Occupieds (i is an absolute index) */
#define get_occupied(filter, i)                         \
  (GET((filter)->blocks[(i)/64].occupieds, (i)%64))
#define set_occupied_to(filter, i, x)                               \
  ((x) ?                                                            \
   SET((filter)->blocks[(i)/64].occupieds, (i)%64) :                \
   UNSET((filter)->blocks[(i)/64].occupieds, (i)%64))
#define set_occupied(filter, i) SET((filter)->blocks[(i)/64].occupieds, (i)%64)
#define unset_occupied(filter, i) UNSET((filter)->blocks[(i)/64].occupieds, (i)%64)

/* Runends (i is an absolute index) */
#define get_runend(filter, i)                           \
  (GET((filter)->blocks[(i)/64].runends, (i)%64))
#define set_runend_to(filter, i, x)                                 \
  ((x) ?                                                            \
   SET((filter)->blocks[(i)/64].runends, (i)%64) :                  \
   UNSET((filter)->blocks[(i)/64].runends, (i)%64))
#define set_runend(filter, i) SET((filter)->blocks[(i)/64].runends, (i)%64)
#define unset_runend(filter, i) UNSET((filter)->blocks[(i)/64].runends, (i)%64)

/* Remainders */
/* Shorthand to get i-th remainder */
#define remainder(filter, i) ((filter)->blocks[(i)/64].remainders[(i)%64])

// Round v to nearest power of 2
// Pre: v >= 0
// https://graphics.stanford.edu/~seander/bithacks.html
static inline size_t nearest_pow_of_2(size_t v) {
  v--;
  v |= v >> 1UL;
  v |= v >> 2UL;
  v |= v >> 4UL;
  v |= v >> 8UL;
  v |= v >> 16UL;
  v |= v >> 32UL;
  v++;
  return v;
}

#ifdef __cplusplus
}
#endif

#endif

#ifndef _MURMURHASH3_H_
#define _MURMURHASH3_H_

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------

void MurmurHash3_x86_32 (const void *key, int len, uint32_t seed, void *out);

void MurmurHash3_x86_128(const void *key, int len, uint32_t seed, void *out);

void MurmurHash3_x64_128(const void *key, int len, uint32_t seed, void *out);

//-----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#endif // _MURMURHASH3_H_
//-----------------------------------------------------------------------------
// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.
// C implementation by Peter Scott:
// https://github.com/PeterScott/murmur3/blob/master/murmur3.h

// Note - The x86 and x64 versions do _not_ produce the same results, as the
// algorithms are optimized for their respective platforms. You can still
// compile and run any of them on any platform, but your performance with the
// non-native version will be less than optimal.


//-----------------------------------------------------------------------------
// Platform-specific functions and macros

#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif

static FORCE_INLINE uint32_t rotl32 ( uint32_t x, int8_t r )
{
  return (x << r) | (x >> (32 - r));
}

static FORCE_INLINE uint64_t rotl64 ( uint64_t x, int8_t r )
{
  return (x << r) | (x >> (64 - r));
}

#define ROTL32(x,y)  rotl32(x,y)
#define ROTL64(x,y) rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

#define getblock(p, i) (p[i])

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

static FORCE_INLINE uint32_t fmix32 ( uint32_t h )
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

//----------

static FORCE_INLINE uint64_t fmix64 ( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}

//-----------------------------------------------------------------------------

void MurmurHash3_x86_32 ( const void * key, int len,
                          uint32_t seed, void * out )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 4;
  int i;

  uint32_t h1 = seed;

  uint32_t c1 = 0xcc9e2d51;
  uint32_t c2 = 0x1b873593;

  //----------
  // body

  const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);

  for(i = -nblocks; i; i++)
  {
    uint32_t k1 = getblock(blocks,i);

    k1 *= c1;
    k1 = ROTL32(k1,15);
    k1 *= c2;

    h1 ^= k1;
    h1 = ROTL32(h1,13);
    h1 = h1*5+0xe6546b64;
  }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

  uint32_t k1 = 0;

  switch(len & 3)
  {
  case 3: k1 ^= tail[2] << 16;
  case 2: k1 ^= tail[1] << 8;
  case 1: k1 ^= tail[0];
          k1 *= c1; k1 = ROTL32(k1,15); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len;

  h1 = fmix32(h1);

  *(uint32_t*)out = h1;
}

//-----------------------------------------------------------------------------

void MurmurHash3_x86_128 ( const void * key, const int len,
                           uint32_t seed, void * out )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 16;
  int i;

  uint32_t h1 = seed;
  uint32_t h2 = seed;
  uint32_t h3 = seed;
  uint32_t h4 = seed;

  uint32_t c1 = 0x239b961b;
  uint32_t c2 = 0xab0e9789;
  uint32_t c3 = 0x38b34ae5;
  uint32_t c4 = 0xa1e38b93;

  //----------
  // body

  const uint32_t * blocks = (const uint32_t *)(data + nblocks*16);

  for(i = -nblocks; i; i++)
  {
    uint32_t k1 = getblock(blocks,i*4+0);
    uint32_t k2 = getblock(blocks,i*4+1);
    uint32_t k3 = getblock(blocks,i*4+2);
    uint32_t k4 = getblock(blocks,i*4+3);

    k1 *= c1; k1  = ROTL32(k1,15); k1 *= c2; h1 ^= k1;

    h1 = ROTL32(h1,19); h1 += h2; h1 = h1*5+0x561ccd1b;

    k2 *= c2; k2  = ROTL32(k2,16); k2 *= c3; h2 ^= k2;

    h2 = ROTL32(h2,17); h2 += h3; h2 = h2*5+0x0bcaa747;

    k3 *= c3; k3  = ROTL32(k3,17); k3 *= c4; h3 ^= k3;

    h3 = ROTL32(h3,15); h3 += h4; h3 = h3*5+0x96cd1c35;

    k4 *= c4; k4  = ROTL32(k4,18); k4 *= c1; h4 ^= k4;

    h4 = ROTL32(h4,13); h4 += h1; h4 = h4*5+0x32ac3b17;
  }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

  uint32_t k1 = 0;
  uint32_t k2 = 0;
  uint32_t k3 = 0;
  uint32_t k4 = 0;

  switch(len & 15)
  {
  case 15: k4 ^= tail[14] << 16;
  case 14: k4 ^= tail[13] << 8;
  case 13: k4 ^= tail[12] << 0;
           k4 *= c4; k4  = ROTL32(k4,18); k4 *= c1; h4 ^= k4;

  case 12: k3 ^= tail[11] << 24;
  case 11: k3 ^= tail[10] << 16;
  case 10: k3 ^= tail[ 9] << 8;
  case  9: k3 ^= tail[ 8] << 0;
           k3 *= c3; k3  = ROTL32(k3,17); k3 *= c4; h3 ^= k3;

  case  8: k2 ^= tail[ 7] << 24;
  case  7: k2 ^= tail[ 6] << 16;
  case  6: k2 ^= tail[ 5] << 8;
  case  5: k2 ^= tail[ 4] << 0;
           k2 *= c2; k2  = ROTL32(k2,16); k2 *= c3; h2 ^= k2;

  case  4: k1 ^= tail[ 3] << 24;
  case  3: k1 ^= tail[ 2] << 16;
  case  2: k1 ^= tail[ 1] << 8;
  case  1: k1 ^= tail[ 0] << 0;
           k1 *= c1; k1  = ROTL32(k1,15); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len; h2 ^= len; h3 ^= len; h4 ^= len;

  h1 += h2; h1 += h3; h1 += h4;
  h2 += h1; h3 += h1; h4 += h1;

  h1 = fmix32(h1);
  h2 = fmix32(h2);
  h3 = fmix32(h3);
  h4 = fmix32(h4);

  h1 += h2; h1 += h3; h1 += h4;
  h2 += h1; h3 += h1; h4 += h1;

  ((uint32_t*)out)[0] = h1;
  ((uint32_t*)out)[1] = h2;
  ((uint32_t*)out)[2] = h3;
  ((uint32_t*)out)[3] = h4;
}

//-----------------------------------------------------------------------------

void MurmurHash3_x64_128 ( const void * key, const int len,
                           const uint32_t seed, void * out )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 16;
  int i;

  uint64_t h1 = seed;
  uint64_t h2 = seed;

  uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
  uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

  //----------
  // body

  const uint64_t * blocks = (const uint64_t *)(data);

  for(i = 0; i < nblocks; i++)
  {
    uint64_t k1 = getblock(blocks,i*2+0);
    uint64_t k2 = getblock(blocks,i*2+1);

    k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

    h1 = ROTL64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

    k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

    h2 = ROTL64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
  }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

  uint64_t k1 = 0;
  uint64_t k2 = 0;

  switch(len & 15)
  {
  case 15: k2 ^= (uint64_t)(tail[14]) << 48;
  case 14: k2 ^= (uint64_t)(tail[13]) << 40;
  case 13: k2 ^= (uint64_t)(tail[12]) << 32;
  case 12: k2 ^= (uint64_t)(tail[11]) << 24;
  case 11: k2 ^= (uint64_t)(tail[10]) << 16;
  case 10: k2 ^= (uint64_t)(tail[ 9]) << 8;
  case  9: k2 ^= (uint64_t)(tail[ 8]) << 0;
           k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

  case  8: k1 ^= (uint64_t)(tail[ 7]) << 56;
  case  7: k1 ^= (uint64_t)(tail[ 6]) << 48;
  case  6: k1 ^= (uint64_t)(tail[ 5]) << 40;
  case  5: k1 ^= (uint64_t)(tail[ 4]) << 32;
  case  4: k1 ^= (uint64_t)(tail[ 3]) << 24;
  case  3: k1 ^= (uint64_t)(tail[ 2]) << 16;
  case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;
  case  1: k1 ^= (uint64_t)(tail[ 0]) << 0;
           k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len; h2 ^= len;

  h1 += h2;
  h2 += h1;

  h1 = fmix64(h1);
  h2 = fmix64(h2);

  h1 += h2;
  h2 += h1;

  ((uint64_t*)out)[0] = h1;
  ((uint64_t*)out)[1] = h2;
}

//-----------------------------------------------------------------------------

#define HIGH (~0ull >> (64 - SEL_CODE_LEN))

int encode_ext(const Ext exts[64], uint64_t* code) {
  uint64_t low = 0;
  uint64_t high = HIGH;

  for (int i=0; i<64; i++) {
    Ext ext = exts[i];
    uint64_t range = high - low;
    // Multiply range by ~0.90624 (Pr[ext is empty])
    uint64_t gap = (range >> 1) + (range >> 2) + (range >> 3) + (range >> 5);
    /* printf("range: %f\n", ((double)range)/HIGH); */
    /* printf("big gap: %f\n", ((double)gap)/HIGH); */
    if (ext.len == 0) {
      // If extension is empty, lower top of range
      high = low + gap;
      /* gap = (range >> 5) + (range >> 6); */
      /* /\* gap = range - ((range >> 1) + (range >> 2) + (range >> 3) + (range >> 5)); *\/ */
      /* printf("small gap: %f\n", ((double)gap)/HIGH); */

    } else {
      // If extension is nonempty, raise bottom of range
      low = low + gap;
      // Set gap to range * ~0.04687
      gap = (range >> 5) + (range >> 6);
      /* gap = range - ((range >> 1) + (range >> 2) + (range >> 3) + (range >> 5)); */
      /* printf("small gap: %f\n", ((double)gap)/HIGH); */
      /* exit(1); */
       /* gap = high - low; */
      // Account for probability of extension length:
      // extension length k>0 has probability 2^{-k}
      for (int j=1; j<ext.len; j++) {
        low += gap;
        gap >>= 1;
      }
      // Account for probability of a particular extension of length k:
      // all equally likely -> 1/(2^k)
      gap >>= ext.len; // Divide gap by 2^k
      low += ext.bits * gap; // Take bits-th 1/(2^k)-bit piece
      high = low + gap;
    }
    if (high - low < 2) {
      return -1;
    }
  }
  *code = low;
  return 0;
}

void decode_ext(uint64_t code, Ext exts[64]) {
  uint64_t low = 0;
  uint64_t high = HIGH;
  for (int i=0; i<64; i++) {
    uint64_t range = high - low;
    // Multiply range by ~0.90624 (Pr[ext is empty])
    uint64_t gap = (range >> 1) + (range >> 2) + (range >> 3) + (range >> 5);
    if (low + gap > code) {
      high = low + gap;
      exts[i].len = 0;
      exts[i].bits = 0;
    } else {
      low = low + gap;
      // Multiply range by ~0.04687
      gap = (range >> 5) + (range >> 6);
      /* gap = range - ((range >> 1) + (range >> 2) + (range >> 3) + (range >> 5)); */
      /* gap = range * 0.04701; */

      // Compute k, the length of the extension, by
      // iteratively shrinking the gap in proportion to
      // the probability of each length k=1, 2, ...
      int len = 1;
      while (low + gap <= code) {
        low += gap;
        gap >>= 1;
        len += 1;
      }
      // Get the bits given k, the length of the extension,
      // by dividing the interval into 2^k even pieces and
      // determining which piece the input belongs to
      gap >>= len;
      uint64_t bits = (code - low)/gap;
      low += bits * gap;
      high = low + gap;

      exts[i].bits = bits;
      exts[i].len = len;
    }
  }
}

int encode_sel(const int sels[64], uint64_t *code) {
  uint64_t low = 0;
  uint64_t high = HIGH;

  for(int i=0; i < 64; i++){
    int letter = sels[i];
    if (letter > MAX_SELECTOR) letter %= MAX_SELECTOR;

    uint64_t range = high - low;
    //takes advantage of "fall through" switch behavior
    switch (letter) {
      default:
        return 0;
      case 6:
        low += (range >> 19) + (range >> 20) + (range >> 23);
      case 5:
        low += (range >> 14) + (range >> 16);
      case 4:
        low += (range >> 10) + (range >> 11);
      case 3:
        low += (range >> 6) + (range >> 8);
      case 2:
        low += (range >> 3) + (range >> 4) + (range >> 7) + (range >> 9);
      case 1:
        low += (range >> 1) + (range >> 2) + (range >> 5);
      case 0: ;
    }

    switch (letter) {
      default:
        //high  = low + (range >> 22);
        return 0;
      case 6:
        high  = low + (range >> 24) + (range >> 25) + (range >> 26);
        break;
      case 5:
        high  = low + (range >> 19) + (range >> 20) + (range >> 23);
        break;
      case 4:
        high  = low + (range >> 14) + (range >> 16);
        break;
      case 3:
        high  = low + (range >> 10) + (range >> 11);
        break;
      case 2:
        high  = low + (range >> 6) + (range >> 8);
        break;
      case 1:
        high  = low + (range >> 3) + (range >> 4) + (range >> 7) + (range >> 9);
        break;
      case 0:
        high  = low + (range >> 1) + (range >> 2) + (range >> 5);
        break;
    }
    // Check if ran out of bits
    if(high - low < 2) return -1;
  }
  *code = low;
  return 0;
}

void decode_sel(uint64_t code, int out[64]) {
  uint64_t low = 0;
  uint64_t high = HIGH;

  for (int i=0; i < 64; i++){
    uint64_t range = high - low;
    uint64_t gap = (range >> 1) + (range >> 2) + (range >> 5);

    if(low + gap > code) {
      high = low + gap;
      out[i] = 0;
      continue;
    }
    low += gap;
    gap = (range >> 3) + (range >> 4) + (range >> 7) + (range >> 9);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 1;
      continue;
    }
    low += gap;
    gap = (range >> 6) + (range >> 8);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 2;
      continue;
    }
    low += gap;
    gap = (range >> 10) + (range >> 11);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 3;
      continue;
    }
    low += gap;
    gap = (range >> 14) + (range >> 16);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 4;
      continue;
    }
    low += gap;
    gap = (range >> 19) + (range >> 20) + (range >> 23);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 5;
      continue;
    }
    low += gap;
    gap = (range >> 24) + (range >> 25) + (range >> 26);
    if(low + gap > code) {
      high = low + gap;
      out[i] = 6;
      continue;
    }
    out[i] = 7;
  }
}

/* Tests */
//#define TEST_ARCD 1
#ifdef TEST_ARCD

#define assert_eq(a, b) assert((a) == (b))

#define test_assert_eq(a, b, msg, ...)   \
  if ((a) != (b)) {                 \
    do {                            \
      fprintf(stderr, "Assertion failed: %s != %s: ", #a, #b); \
      fprintf(stderr, msg"\n", __VA_ARGS__); \
      assert_eq(a, b);              \
    } while (0);                    \
  }

#define panic(...) \
  do {             \
    fprintf(stderr, __VA_ARGS__);  \
    exit(1);       \
  } while (0)

/**
 * Convert a string into an ext
 */
void str_to_ext(char* str, Ext* ext) {
  int len = (int)strlen(str);
  if (len == 0) {
    ext->bits = 0;
  } else {
    ext->bits = strtol(str, NULL, 2);
  }
  ext->len = len;
}

/**
 * Convert an array of strs into an array of exts
 */
void strs_to_exts(char* strs[64], Ext exts[64]) {
  for (int i=0; i<64; i++) {
    str_to_ext(strs[i], &exts[i]);
  }
}

int ext_arr_eq(Ext a[64], Ext b[64]) {
  for (int i=0; i<64; i++) {
    if (a[i].len != b[i].len ||
        (a[i].bits & ONES(a[i].len)) != (b[i].bits & ONES(b[i].len))) {
      return 0;
    }
  }
  return 1;
}

void test_strs_to_exts() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char* strs[64] = {
          "", "", "", "", "", "", "", "",
          "0", "0", "0", "0", "0", "0", "0", "0",
          "1", "1", "1", "1", "1", "1", "1", "1",
          "01", "01", "01", "01", "01", "01", "01", "01",
          "10", "10", "10", "10", "10", "10", "10", "10",
          "11", "11", "11", "11", "11", "11", "11", "11",
          "100", "100", "100", "100", "100", "100", "100", "100",
          "1000", "1000", "1000", "1000", "1000", "1000", "1000", "1000",
  };
  strs_to_exts(strs, exts);
  for (int i=0; i<64; i++) {
    switch (i / 8) {
      case 0:
        assert_eq(exts[i].bits, 0);
        assert_eq(exts[i].len, 0);
        break;
      case 1:
        assert_eq(exts[i].bits, 0);
        assert_eq(exts[i].len, 1);
        break;
      case 2:
        assert_eq(exts[i].bits, 1);
        assert_eq(exts[i].len, 1);
        break;
      case 3:
        assert_eq(exts[i].bits, 1);
        assert_eq(exts[i].len, 2);
        break;
      case 4:
        assert_eq(exts[i].bits, 0b10);
        assert_eq(exts[i].len, 2);
        break;
      case 5:
        assert_eq(exts[i].bits, 0b11);
        assert_eq(exts[i].len, 2);
        break;
      case 6:
        assert_eq(exts[i].bits, 0b100);
        assert_eq(exts[i].len, 3);
        break;
      case 7:
        assert_eq(exts[i].bits, 0b1000);
        assert_eq(exts[i].len, 4);
        break;
    }
  }
  printf("passed.\n");
}

/*
 * Check that exts = decode(encode(exts))
 */
void test_encode_decode_w_input(Ext exts[64]) {
  uint64_t code;
  Ext decoded[64];
  int out;
  out = encode_ext(exts, &code);
  if (out < 0) {
    panic("Encoding failed\n");
  } else {
    decode_ext(code, decoded);
    assert(ext_arr_eq(exts, decoded));
  }
}

void test_encode_decode_w_input_expect_fail(Ext exts[64]) {
  uint64_t code;
  int out;
  out = encode_ext(exts, &code);
  assert_eq(out, -1);
}

void test_encode_decode_empty() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_one() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "0",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_few() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "000", "10", "1", "0",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_many() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "0", "0", "0", "0", "0", "0", "0", "0",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_many_rev() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "0", "0", "0", "0", "0", "0", "0", "0",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

void test_encode_decode_long() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "1111111111111111111", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input(exts);
  printf("passed.\n");
}

int will_overflow(Ext exts[64]) {
  uint64_t code;
  Ext decoded[64];
  if (encode_ext(exts, &code) < 0) {
    return 1;
  } else {
    decode_ext(code, decoded);
    assert(ext_arr_eq(exts, decoded));
    return 0;
  }
}

void print_exts(Ext exts[64]) {
  for (int k=0; k<64; k++) {
    printf(" %d", exts[k].len);
  }
  printf("\n");
}


void test_encode_decode_capacity() {
  printf("Testing %s...\n", __FUNCTION__);
  int limit = 20;
  Ext exts[64];
  // Try making ext arrays with different configs to stress test arcd
  // i: length of ext
  // j: number of i-length exts
  // k: iterate over exts
  for (int len=1; len < limit; len++) {
    for (int i=0; i < 64; i++) {
      exts[i].bits = 0;
      exts[i].len = 0;
    }
    int n;
    for (n=1; n < 64; n++) {
      for (int i=0; i < 64; i++) {
        exts[i].len = i < n ? len : 0;
      }
      if (will_overflow(exts)) {
        n--;
        break;
      }
    }
    printf("Can hold %d %d-length exts\n", n, len);
  }
}

void test_encode_decode_too_many() {
  printf("Testing %s...", __FUNCTION__);
  Ext exts[64];
  char *strs[64] = {
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "",
          "", "", "", "", "", "", "", "0",
          "0", "0", "0", "0", "0", "0", "0", "0",
  };
  strs_to_exts(strs, exts);
  test_encode_decode_w_input_expect_fail(exts);
  printf("passed.\n");
}

int main() {
  test_strs_to_exts();
  test_encode_decode_empty();
  test_encode_decode_one();
  test_encode_decode_few();
  test_encode_decode_many();
  test_encode_decode_many_rev();
  test_encode_decode_long();
  test_encode_decode_capacity();
  test_encode_decode_too_many();
}

#endif // TEST_ARCD

#ifndef EXAF_CONSTANTS_H
#define EXAF_CONSTANTS_H

/** Size of a remainder in the filter
 * Must be power of 2 since mod is taken using & */
#define REM_SIZE 7

#endif //EXAF_CONSTANTS_H

#ifndef AQF_REMAINDERS_H
#define AQF_REMAINDERS_H

#ifdef __cplusplus
extern "C" {
#endif

#if REM_SIZE <= 8
typedef uint8_t rem_t;
#elif REM_SIZE <= 16
typedef uint16_t rem_t;
#elif REM_SIZE <= 32
typedef uint32_t rem_t;
#else
typedef uint64_t rem_t;
#endif

#ifdef __cplusplus
}
#endif

#endif //AQF_REMAINDERS_H

#ifndef HASH_SET_H
#define HASH_SET_H

#ifdef __cplusplus
extern "C" {
#endif

#define HASH_SET_SEED 26571997

//linked list node for hashing with chaining
typedef struct setnode_t {
  struct setnode_t* next;
  char* value;
  int sources;
} Setnode;

int set_insert(char* word, int length, int source, Setnode* set, int set_size);

int set_lookup(char* word, int length, Setnode* set, int set_size);

/**
 * Returns a char* array of all values stored in the set.
 */
char** get_values(Setnode *set, int set_size);

void set_deallocate(Setnode* set, int set_size);

void print_set(Setnode* set, int set_size);

int run_set_tests();


#ifdef __cplusplus
}
#endif

#endif

//insert the string word into set
//return 1 if inserted in the set; 0 if already found (and not inserted)
int set_insert(char* word, int length, int source, Setnode* set, int set_size) {
  uint32_t hash[4];
  MurmurHash3_x64_128(word, length, HASH_SET_SEED, hash);
  int index = hash[0] % set_size;

  Setnode* currentNode = &set[index];

  if(currentNode->value == NULL) {
      currentNode->value = strdup(word);
      //currentNode->sources = (1 << source);
      currentNode->sources = source + 1;
  } else {
    while(currentNode->next != NULL){
      if(strcmp(currentNode->value, word) == 0) {
        //currentNode->sources |= (1 << source);
        currentNode->sources = source + 1;
        return 0;
      }
      currentNode = currentNode->next;
    }
    if(strcmp(currentNode->value, word) == 0) {
      //currentNode->sources |= (1 << source);
      currentNode->sources = source + 1;
      return 0;
    }
    Setnode* newNode = (Setnode*)malloc(sizeof(*currentNode));
    currentNode->next = newNode;

    newNode->next = NULL;
    newNode->value = strdup(word);
    //newNode->sources = (1 << source);
    currentNode->sources = source + 1;
  }
  return 1;
}

//look up word in the set
//returns sources if it is found; returns 0 otherwise
int set_lookup(char* word, int length, Setnode* set, int set_size) {
  uint32_t hash[4];
  MurmurHash3_x64_128(word, length, HASH_SET_SEED, hash);
  int index = hash[0] % set_size;

  Setnode* currentNode = &set[index];
  if(currentNode->value == NULL) {
    return 0;
  } else {
    while(currentNode->next != NULL){
      if(strcmp(currentNode->value, word) == 0) {
        return currentNode->sources;
      }
      currentNode = currentNode->next;
    }
    if(strcmp(currentNode->value, word)){
      return 0;
    } else{
      return currentNode->sources;
    }
  }
}

void set_deallocate(Setnode* set, int set_size) {
  for(int i = 0; i < set_size; i++) {
    Setnode* current = &set[i];
    Setnode* next = current->next;
    while(next != NULL) {
      current = next;
      next = current->next;
      free(current->value);
      free(current);
    }
    free(set[i].value);
  }
  free(set);
}

char** get_values(Setnode *set, int set_size) {
  char** values = (char**)calloc(set_size, sizeof(char*));
  int j=0;
  for (int i=0; i<set_size; i++) {
    Setnode* node = &set[i];
    if (node->value != NULL) {
      do {
        values[j] = node->value;
        node = node->next;
        j++;
      } while (node != NULL);
    }
  }
  return values;
}

/* Printing */

void print_set(Setnode* set, int set_size) {
  printf("SET (size=%d):\n", set_size);
  for(int i = 0; i < set_size; i++) {
    if (set[i].value != NULL) {
      printf(" %d: [%s]", i, set[i].value);
      Setnode* next = set[i].next;
      while (next != NULL) {
        printf("-> [%s]", next->value);
        next = next->next;
      }
      printf("\n");
    }
  }
}

/* Tests */

/**
  * @return 1 if the set contains the word, 0 otherwise
 */
static int contains(char* set[], int set_size, char* word) {
  for (int i=0; i<set_size; i++)
    if (strcmp(set[i], word) == 0) return 1;
  return 0;
}

static int test_set_values() {
  printf("Starting get_values test...\n");
  char* inputs[] = {
          "premiere", "partie", "combray", "i", "longtemps", "me", "suis", "couche", "de", "bonne", "heure", "parfois",
          "a", "peine", "ma", "bougie", "eteinte", "mes", "yeux", "se", "fermaient", "si", "vite", "que", "n", "pas",
          "le", "temps", "dire", "m", "et", "une", "demi", "apres", "des", "epoques", "vecues", "par", "eux",
          "distantes", "entre", "lesquelles", "tant", "jours", "sont", "venus", "placer", "dans",
  };
  int num_words = sizeof(inputs) / sizeof(char*);
  int set_size = (int)(num_words * 1.5);
  Setnode* set = new Setnode[set_size];
  for (int i=0; i < num_words; i++) {
    set_insert(inputs[i], (int)strlen(inputs[i]), 0, set, set_size);
  }
  char** values = get_values(set, set_size);

  // Check that values = inputs by checking subsets both ways
  for (int i=0; i < num_words; i++) {
    if (!contains(values, num_words, inputs[i])) {
      printf("get_values test failed: input %s not found in set.\n", inputs[i]);
      return 0;
    }
    if (!contains(inputs, num_words, values[i])) {
      printf("get_values test failed: %s found in set but not in inputs.\n", values[i]);
      return 0;
    }
  }
  printf("get_values test passed.\n");
  return 1;
}

int run_set_tests() {
  return test_set_values();
}

#ifndef TAF_H
#define TAF_H

#ifdef __cplusplus
extern "C" {
#endif

#define TAF_MODE_NORMAL 0
#define TAF_MODE_ARCD_OVERWRITE 1

typedef struct taf_block_t {
  rem_t remainders[64];
  uint64_t occupieds;
  uint64_t runends;
  size_t offset;
  uint8_t sel_code[SEL_CODE_BYTES];
} TAFBlock;

typedef uint64_t elt_t;

typedef struct remote_elt_t {
  uint64_t elt;
  uint64_t hash;
} Remote_elt ;

typedef struct taf_t {
  size_t p;                     /* fingerprint prefix size = log2(n/E) to get false-pos rate E */
  size_t q;                     /* length of quotient */
  size_t r;                     /* length of remainder */
  size_t nslots;                /* number of slots available (2^q) */
  size_t nblocks;               /* nslots/64 */
  size_t nelts;                 /* number of elements stored  */
  int seed;                     /* seed for Murmurhash */
  TAFBlock* blocks;           /* blocks of 64 remainders with metadata  */
  Remote_elt* remote;           /* array of inserted elements (up to 64 bits) */

  // Extra modes
  int mode;            // mode flag: handle non-adaptive case
} TAF;

void taf_init(TAF *filter, size_t n, int seed);
void taf_destroy(TAF* filter);
int taf_lookup(TAF *filter, elt_t elt);
void taf_insert(TAF *filter, elt_t elt);
void taf_clear(TAF* filter);

// Printing
double taf_load(TAF *filter);
void print_taf(TAF* filter);
void print_taf_metadata(TAF* filter);
void print_taf_block(TAF* filter, size_t block_index);

#ifdef __cplusplus
}
#endif

#endif //TAF_H

//! [hash]
struct HashFunctor {
        size_t operator () (uint64_t key) const {
               return (size_t)(key * 2654435761u);
        }
};
//! [hash]

//! [comparator]
struct CompareGreater {
        bool operator () (const uint64_t& a, const uint64_t& b) const {
                return a > b;
        }
        static uint64_t max_value() {
                return 0ULL;
        }
};
//! [comparator]

typedef stxxl::map<uint64_t, uint64_t, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> ordered_map_t;
typedef std::pair<uint64_t, uint64_t> pair_t;

/**
 * Generate a hash for the input word.
 * Returns the full 128 bit murmurhash.
 * Expects a two-slot array of uint64_t's.
 */
static uint64_t taf_hash(const TAF *filter, elt_t elt) {
  uint64_t buf[2];
  MurmurHash3_x64_128(&elt, 8, filter->seed, buf);
  return buf[0];
}

/**
 * Returns the quotient for a 64-bit fingerprint hash.
 */
static size_t calc_quot(const TAF* filter, uint64_t hash) {
  return hash & ONES(filter->q);
}

/**
 * Returns the k-th remainder for h
 */
static rem_t calc_rem(const TAF* filter, uint64_t hash, int k) {
  int n_rems = (64 - (int)filter->q)/(int)filter->r;
  if (k >= n_rems) k %= n_rems;
  return (hash >> (filter->q + k * filter->r)) & ONES(filter->r);
}

/* TAF Helpers */

/**
 * @return The selector code at block `block_i`, padded with zeros as a uint64_t.
 */
static uint64_t get_sel_code(const TAF* filter, size_t block_i) {
  uint64_t code = 0;
  memcpy(&code, filter->blocks[block_i].sel_code, SEL_CODE_BYTES);
  return code;
}

/**
 * Set the selector code for the `block_i`-th block to the first `CODE_BYTES`
 * bits of `code`.
 */
static void set_sel_code(TAF* filter, size_t block_i, uint64_t code) {
  memcpy(filter->blocks[block_i].sel_code, &code, SEL_CODE_BYTES);
}

/**
 * Print an array of selectors.
 */
static void print_sels(const int sels[64]) {
  for (int i=0; i<8; i++) {
    printf("   ");
    for (int j=0; j<8; j++) {
      int sel = sels[i*8 + j];
      if (sel == 0) {
        printf(" _");
      } else {
        printf(" %d", sel);
      }
    }
    printf("\n");
  }
}

/**
 * Returns the absolute index of the `rank`-th 1 bit in Q.runends past the start of
 * the block at `block_index`. `rank` indexes from 0.
 *
 * Returns -1 if result is invalid (out of bounds).
 */
static int select_runend(const TAF* filter, size_t block_index, size_t rank) {
  assert(block_index < filter->nblocks && "block_index out of bounds");

  size_t step;
  size_t loc = block_index * 64;
  while (1) {
    TAFBlock* b = &filter->blocks[loc / 64];
    step = bitselect(b->runends, rank >= 64 ? 63 : (int)rank);
    loc += step;
    if (step != 64 || loc >= filter->nslots) {
      break;
    }
    rank -= popcnt(b->runends);
  }
  if (loc >= filter->nslots) {
    return -1;
  } else {
    return (int)loc;
  }
}

#define RANK_SELECT_EMPTY (-1)
#define RANK_SELECT_OVERFLOW (-2)
/** Performs the blocked equivalent of the unblocked operation
 *    y = select(Q.runends, rank(Q.occupieds, x)).
 *  Note: x indexes from 0.
 *
 *  Return behavior:
 *  - If y <= x, returns Empty
 * - If y > x, returns Full(y)
 * - If y runs off the edge, returns Overflow
 */
static int rank_select(const TAF* filter, size_t x) {
  // Exit early if x obviously out of range
  if (x >= filter->nslots) {
    return RANK_SELECT_OVERFLOW;
  }
  size_t block_i = x/64;
  size_t slot_i = x%64;
  TAFBlock *b = &filter->blocks[block_i];

  // Compute i + O_i where i = x - (x mod 64)
  if (!GET(b->occupieds, 0) && b->offset == 0 && !GET(b->runends, 0)) {
    // b[0] unoccupied, b.offset = 0, b[0] not a runend =>
    // negative offset
    if (slot_i == 0) {
      return RANK_SELECT_EMPTY;
    }
  } else {
    // non-negative offset
    if (slot_i == 0) {
      return (int)(block_i * 64 + b->offset);
    } else {
      block_i += b->offset/64;
    }
  }

  // Handle case where offset runs off the edge
  if (block_i >= filter->nblocks) {
    return RANK_SELECT_OVERFLOW;
  }

  // Count the number of occupied quotients between i+1 (b.start + i) and j (x)
  uint64_t d = bitrank(b->occupieds, slot_i) - GET(b->occupieds, 0);

  // Advance offset to relevant value for the block that b.offset points to
  size_t offset = b->offset % 64;
  b = &filter->blocks[block_i];

  // Account for the runends in [0, offset] of the new block
  d += bitrank(b->runends, offset);

  // If rank(Q.occupieds, x) == 0, then there's nothing to see here
  if (d == 0) {
    return RANK_SELECT_EMPTY;
  } else {
    // (rank-1) accounts for select's indexing from 0
    int loc = select_runend(filter, block_i, d-1);
    if (loc == -1) {
      return RANK_SELECT_OVERFLOW;
    } else if (loc < x) {
      return RANK_SELECT_EMPTY;
    } else {
      return loc;
    }
  }
}

#define NO_UNUSED (-1)
/**
 * Finds the first unused slot at or after absolute location x.
 */
static int first_unused(const TAF* filter, size_t x) {
  while (1) {
    int loc = rank_select(filter, x);
    switch (loc) {
      case RANK_SELECT_EMPTY: return x;
      case RANK_SELECT_OVERFLOW: return NO_UNUSED;
      default:
        if (x <= loc) {
          x = loc + 1;
        } else {
          return x;
        }
    }
  }
}

/**
 * Shift the remainders and runends in [a, b] forward by 1 into [a+1, b+1]
 */
static void shift_rems_and_runends(TAF* filter, int a, int b) {
  if (a > b) return;
  for (int i=b; i>=a; i--) {
    remainder(filter, i+1) = remainder(filter, i);
    set_runend_to(filter, i+1, get_runend(filter, i));
  }
  set_runend_to(filter, a, 0);
}

/**
 * Shift the remote elements in [a,b] forward by 1
 */
static void shift_remote_elts(TAF* filter, int a, int b) {
  if (a > b) return;
  for (int i=b; i>=a; i--) {
    filter->remote[i+1] = filter->remote[i];
  }
  filter->remote[a].elt = 0;
  filter->remote[a].hash = 0;
}

static void inline swap_ptrs(int **a, int **b) {
  int *tmp = *a;
  *a = *b;
  *b = tmp;
}

/**
 * Helper for `shift_sels`.  Shifts sels in `[0, b]` a single block.
 */
static void shift_block_sels(TAF *filter, int block_i, int sels[64], const int prev_sels[64], int b) {
  uint64_t code;
  for (int i=b; i > 0; i--) {
    sels[i] = sels[i-1];
  }
  sels[0] = prev_sels[63];
  if (encode_sel(sels, &code) == -1) {
    code = 0;
  }
  set_sel_code(filter, block_i, code);
}

/**
 * Shift the hash selectors in [a,b] forward by 1
 */
static void shift_sels(TAF* filter, int a, int b) {
  if (a > b) return;
  uint64_t code;
  if (a/64 == (b+1)/64) {
    // a and b+1 in the same block
    int sels[64];
    decode_sel(get_sel_code(filter, a/64), sels);
    for (int i = (b+1)%64; i > a%64; i--) {
      sels[i] = sels[i-1];
    }
    sels[a%64] = 0;
    if (encode_sel(sels, &code) == -1) {
      code = 0;
    }
    set_sel_code(filter, a/64, code);
  } else {
    // a and b+1 in different blocks
    int* sels = new int[64];
    int* prev_sels = new int[64];
    // (1) last block
    int block_i = (b+1)/64;
    decode_sel(get_sel_code(filter, block_i), sels);
    decode_sel(get_sel_code(filter, block_i - 1), prev_sels);
    shift_block_sels(filter, block_i, sels, prev_sels, (b + 1) % 64);
    swap_ptrs(&sels, &prev_sels);
    // (2) middle blocks
    for (block_i--; block_i > a/64; block_i--) {
      decode_sel(get_sel_code(filter, block_i - 1), prev_sels);
      shift_block_sels(filter, block_i, sels, prev_sels, 63);
      swap_ptrs(&sels, &prev_sels);
    }
    // (3) first block
    for (int i=63; i>a%64; i--) {
      sels[i] = sels[i-1];
    }
    sels[a%64] = 0;
    if (encode_sel(sels, &code) == -1) {
      code = 0;
    }
    set_sel_code(filter, a/64, code);
    free(sels);
    free(prev_sels);
  }
}

/**
 * Increment all non-negative offsets with targets in [a,b]
 */
static void inc_offsets(TAF* filter, size_t a, size_t b) {
  assert(a < filter->nslots && b < filter->nslots);
  // Exit early if invalid range
  if (a > b) {
    return;
  }
  // Start i at the first block after b, clamping it so it doesn't go off the end, and work backwards
  size_t start = min(b/64 + 1, filter->nblocks - 1);
  for (int i = start; i>=0; i--) {
    TAFBlock *block = &filter->blocks[i];
    size_t block_start = i * 64;
    // Skip this block if it has a negative offset
    if (!GET(block->occupieds, 0) &&
        block->offset == 0 &&
        !GET(block->runends, 0)) {
      continue;
    }
    // Exit if the target for b.offset is before the interval;
    // if it's within the interval, increment offset
    size_t target = block_start + block->offset;
    if (target < a) {
      break;
    } else if (target <= b) {
      block->offset++;
    }
  }
}

/**
 * Increment non-negative offsets to accommodate insertion of a new run
 * for `quot` at `loc`.
 *
 * Concretely, this function increments unowned offsets in blocks whose
 * first slot `s` is not after `quot`: `s >= quot`.
 */
static void inc_offsets_for_new_run(TAF* filter, size_t quot, size_t loc) {
  assert(loc < filter->nslots);
  // Start i at the first block after loc,
  // clamping it so it doesn't go off the end
  size_t start = min(loc/64 + 1, filter->nblocks - 1);
  for (int i=start; i>=0; i--) {
    TAFBlock *b = &filter->blocks[i];
    size_t b_start = i*64;
    // Skip this block if it has a negative offset
    if (!GET(b->occupieds, 0) && b->offset == 0 && !GET(b->runends, 0)) {
      continue;
    }
    // Exit if the target for b.offset is before the interval;
    // if the target is within the interval, increment b.offset
    size_t target = b_start + b->offset;
    if (target < loc) {
      break;
    } else if (target == loc && !GET(b->occupieds, 0) && quot <= b_start) {
      b->offset++;
    }
  }
}

static void add_block(TAF *filter) {
  // Add block to new_blocks
  TAFBlock *new_blocks = (TAFBlock*)realloc(filter->blocks, (filter->nblocks + 1) * sizeof(TAFBlock));
  if (new_blocks == NULL) {
    printf("add_block failed to realloc new blocks\n");
    exit(1);
  }
  filter->blocks = new_blocks;
  memset(filter->blocks + filter->nblocks, 0, sizeof(TAFBlock));

  // Reallocate remote rep
  Remote_elt *new_remote = (Remote_elt*)realloc(filter->remote,(filter->nslots + 64) * sizeof(Remote_elt));
  if (new_remote == NULL) {
    printf("add_block failed to realloc new remote rep\n");
    exit(1);
  }
  filter->remote = new_remote;
  memset(filter->remote + filter->nslots, 0, 64 * sizeof(Remote_elt));

  // Update counters
  filter->nblocks += 1;
  filter->nslots += 64;
}

/**
 * Adapt a fingerprint at a particular location by incrementing the selector and
 * updating the remainder.
 */
static void adapt_loc(TAF *filter, size_t loc, int sels[64]) {
  // Increment selector at loc%64
  int old_sel = sels[loc%64];
  int new_sel = (old_sel + 1) % MAX_SELECTOR;
  sels[loc%64] = new_sel;
  // Write encoding to block
  uint64_t code;
  if (encode_sel(sels, &code) == -1) {
    // Encoding failed: rebuild
    // Reset all remainders and selectors in block
    memset(sels, 0, 64 * sizeof(sels[0]));
    TAFBlock *b = &filter->blocks[loc/64];
    uint64_t b_start = loc - (loc % 64);
    for (int i=0; i<64; i++) {
      b->remainders[i] = calc_rem(filter, filter->remote[b_start + i].hash, 0);
    }
    // Set sel to new_sel and attempt encode
    sels[loc % 64] = new_sel;
    if (encode_sel(sels, &code) == -1) {
      fprintf(stderr, "Encoding (sel=%d) failed after rebuild!\n", new_sel);
      sels[loc % 64] = 0;
      new_sel = 0;
      code = 0;
    }
  }
  // Encoding succeeded: update sel_code and remainder
  rem_t new_rem = calc_rem(filter, filter->remote[loc].hash, new_sel);
  switch (filter->mode) {
    case TAF_MODE_NORMAL:
      remainder(filter, loc) = new_rem;
      set_sel_code(filter, loc/64, code);
      break;
    case TAF_MODE_ARCD_OVERWRITE:
      remainder(filter, loc) = remainder(filter, loc);
      set_sel_code(filter, loc/64, 0);
      break;
  }
}

/**
 * Adapt on a query element that collided with a stored fingerprint at loc.
 *
 * Go through the rest of the run and fix any other remaining collisions.
 */
static void adapt(TAF *filter, elt_t query, int loc, size_t quot, uint64_t hash, int sels[64]) {
  assert(quot <= loc && loc < filter->nslots);
  // Make sure the query elt isn't mapped to an earlier index in the sequence
  for (int i=loc; i>=(int)quot && (i == loc || !get_runend(filter, i)); i--) {
    if (filter->remote[i].elt == query) {
      return;
    }
  }
  // Adapt on all collisions in the run
  for (int i=loc; i>=(int)quot && (i == loc || !get_runend(filter, i)); i--) {
    // Re-decode if at a new block
    if (i != loc && i % 64 == 63) {
      decode_sel(get_sel_code(filter, i/64), sels);
    }
    // Check collision
    int sel = sels[i % 64];
    if (remainder(filter, i) == calc_rem(filter, hash, sel)) {
      adapt_loc(filter, i, sels);
    }
  }
}

/* TAF */

void taf_init(TAF *filter, size_t n, int seed) {
  filter->seed = seed;
  filter->nelts = 0;
  filter->nblocks = max(1, nearest_pow_of_2(n)/64);
  filter->nslots = filter->nblocks * 64;
  filter->q = (size_t)log2((double)filter->nslots); // nslots = 2^q
  filter->r = REM_SIZE;
  filter->p = filter->q + filter->r;
  filter->blocks = new TAFBlock[filter->nblocks];
  filter->remote = new Remote_elt[filter->nslots];
  filter->mode = TAF_MODE_NORMAL;
}

void taf_destroy(TAF* filter) {
  free(filter->blocks);
  free(filter->remote);
  free(filter);
}

void taf_clear(TAF* filter) {
  filter->nelts = 0;
  free(filter->blocks);
  free(filter->remote);
  filter->blocks = new TAFBlock[filter->nblocks];
  filter->remote = new Remote_elt[filter->nslots];
}

int extra_blocks = 0;
static void raw_insert(TAF* filter, elt_t elt, uint64_t hash) {
  size_t quot = calc_quot(filter, hash);
  rem_t rem = calc_rem(filter, hash, 0);
  filter->nelts++;

  // Find the appropriate runend
  int r = rank_select(filter, quot);
  switch (r) {
    case RANK_SELECT_EMPTY: {
      set_occupied(filter, quot);
      set_runend(filter, quot);
      remainder(filter, quot) = rem;
      filter->remote[quot].elt = elt;
      filter->remote[quot].hash = hash;
      break;
    }
    case RANK_SELECT_OVERFLOW: {
      printf("TAF failed to find runend (nslots=%lu, quot=(block=%lu, slot=%lu))\n",
             filter->nslots, quot/64, quot%64);
      exit(1);
    }
    default: {
      // Find u, the first open slot after r, and
      // shift everything in [r+1, u-1] forward by 1 into [r+2, u],
      // leaving r+1 writable
      size_t u = first_unused(filter, r+1);
      if (u == NO_UNUSED) {
        // Extend filter by one block and use the first empty index
        add_block(filter);
        u = filter->nslots - 64;
	extra_blocks++;
      }
      inc_offsets(filter, r+1, u-1);
      shift_rems_and_runends(filter, r + 1, (int)u - 1);
      shift_remote_elts(filter, r + 1, (int)u - 1);
      shift_sels(filter, r + 1, (int)u - 1);

      // Start a new run or extend an existing one
      if (get_occupied(filter, quot)) {
        // quot occupied: extend an existing run
        inc_offsets(filter, r, r);
        unset_runend(filter, r);
      } else {
        // quot unoccupied: start a new run
        inc_offsets_for_new_run(filter, quot, r);
        set_occupied(filter, quot);
      }
      set_runend(filter, r+1);
      remainder(filter, r+1) = rem;
      filter->remote[r+1].elt = elt;
      filter->remote[r+1].hash = hash;
    }
  }
}

static int raw_lookup(TAF* filter, elt_t elt, uint64_t hash) {
  size_t quot = calc_quot(filter, hash);

  if (get_occupied(filter, quot)) {
    int loc = rank_select(filter, quot);
    if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
      return 0;
    }
    // Cache decoded selectors
    int decoded[64];
    int decoded_i = -1;
    do {
      // Refresh cached code
      if (decoded_i != loc/64) {
        decoded_i = loc/64;
        uint64_t code = get_sel_code(filter, loc/64);
        decode_sel(code, decoded);
      }
      int sel = decoded[loc%64];
      rem_t rem = calc_rem(filter, hash, sel);
      if (remainder(filter, loc) == rem) {
        // Check remote
        if (elt != filter->remote[loc].elt) {
          adapt(filter, elt, loc, quot, hash, decoded);
        }
        return loc;
      }
      loc--;
    } while (loc >= (int)quot && !get_runend(filter, loc));
  }
  return 0;
}

/**
 * Return 1 if word is in the filter.
 *
 * Exits with 0 immediately if quot(word) is unoccupied.
 * Otherwise, linear probes through the run to see if the
 * run contains rem(word).
 */
int taf_lookup(TAF *filter, elt_t elt) {
  uint64_t hash = taf_hash(filter, elt);
  return raw_lookup(filter, elt, hash);
}

void taf_insert(TAF *filter, elt_t elt) {
  uint64_t hash = taf_hash(filter, elt);
  raw_insert(filter, elt, hash);
}

double taf_load(TAF *filter) {
  return (double)filter->nelts/(double)filter->nslots;
}

/* Printing */

void print_taf_metadata(TAF* filter) {
  printf("FILTER METADATA:\n");
  printf("  p=%ld, q=%ld, r=%ld\n",
         filter->p, filter->q, filter->r);
  printf("  nslots=%ld, nblocks=%ld, blocksize=%ld, nelts=%ld\n",
         filter->nslots, filter->nslots/64, sizeof(TAFBlock), filter->nelts);
  printf("  seed=%d\n", filter->seed);
  printf("  load factor=%f\n", taf_load(filter));
}

void print_taf_block(TAF* filter, size_t block_index) {
  assert(0 <= block_index && block_index < filter->nslots/64);
  TAFBlock block = filter->blocks[block_index];
  printf("BLOCK %lu:\n", block_index);
  printf("  occupieds=0x%lx\n", block.occupieds);
  printf("  runends=0x%lx\n", block.runends);
  printf("  offset=%ld\n", block.offset);
  printf("  remainders=\n");
  // Print out 8x8
  for (int i=0; i<8; i++) {
    printf("   ");
    for (int j=0; j<8; j++) {
      printf(get_occupied(filter, block_index*64 + i*8 + j) ? "o" : " ");
      printf(get_runend(filter, block_index*64 + i*8 + j) ? "r" : " ");
      printf(" 0x%-*x", (int)(filter->r / 8 + 3), block.remainders[i*8+j]);
    }
    printf("\n");
  }
  printf("  selector code=0x%lx\n", get_sel_code(filter, block_index));
  printf("  selectors=\n");
  int sels [64];
  decode_sel(get_sel_code(filter, block_index), sels);
  print_sels(sels);
  printf("  remote elts=\n");
  for (int i=0; i<8; i++) {
    printf("   ");
    for (int j=0; j<8; j++) {
      printf(" 0x%-*lx", 8, filter->remote[block_index * 64 + i*8 + j].elt);
    }
    printf("\n");
  }
}

void print_taf(TAF* filter) {
  print_taf_metadata(filter);
  for (int i=0; i<filter->nblocks; i++) {
    print_taf_block(filter, i);
  }
}

void print_taf_stats(TAF* filter) {
  printf("TAF stats:\n");
  // Hash selector counts
  int sel_counts[MAX_SELECTOR];
  for (int i=0; i<MAX_SELECTOR; i++) {
    sel_counts[i] = 0;
  }
  int sels[64];
  for (int i=0; i<filter->nslots; i++) {
    if (i%64 == 0) {
      decode_sel(get_sel_code(filter, i/64), sels);
    }
    sel_counts[sels[i%64]]++;
  }
  printf("Hash selector counts:\n");
  for (int i=0; i<MAX_SELECTOR; i++) {
    printf(" %d: %d (%f%%)\n", i, sel_counts[i],
           100 * (double)sel_counts[i]/(double)filter->nslots);
  }
}

// Tests
#define TEST_TAF 1

void print_backtrace() {
  void* callstack[128];
  int i, frames = backtrace(callstack, 128);
  char** strs = backtrace_symbols(callstack, frames);
  printf("\n");
  for (i = 0; i < frames; ++i) {
    printf("%s\n", strs[i]);
  }
  free(strs);
}

#define assert_eq(a, b) assert((a) == (b))

#define test_assert_eq(a, b, msg, ...)   \
  if ((a) != (b)) {                 \
    do {                            \
      fprintf(stderr, "Assertion failed: %s != %s: ", #a, #b); \
      fprintf(stderr, msg"\n", __VA_ARGS__); \
      assert_eq(a, b);              \
    } while (0);                    \
  }

#define TAF_SEED 32776517

TAF *new_taf(size_t n) {
  TAF *filter = new TAF();
  taf_init(filter, n, TAF_SEED);
  return filter;
}

void test_calc_rem() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(128);

  // q = 7, r = 8
  for (int i=0; i<16; i++) {
    assert_eq(calc_rem(filter, 0, i), 0);
    assert_eq(calc_rem(filter, 0b1111111, i), 0);
  }
  assert_eq(calc_rem(filter, 0b10000000, 0), 1);
  for (int i=1; i<7; i++) {
    assert_eq(calc_rem(filter, 0b10000000, i), 0);
  }
  assert_eq(calc_rem(filter, 0b111000001110000000, 0), 0b111);
  assert_eq(calc_rem(filter, 0b111000001110000000, 1), 0b111);
  for (int i=2; i<7; i++) {
    assert_eq(calc_rem(filter, 0b111000001110000000, i), 0);
  }
  assert_eq(calc_rem(filter, 0b111000001110000000, 7), 0b111);
  assert_eq(calc_rem(filter, 0b111000001110000000, 8), 0b111);

  taf_destroy(filter);
  printf("passed.\n");
}

void test_add_block() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(64 * 2);
  assert_eq(filter->nslots, 128);
  assert_eq(filter->nblocks, 2);
  add_block(filter);
  // Check metadata
  assert_eq(filter->nslots, 192);
  assert_eq(filter->nblocks, 3);
  // Check new block
  TAFBlock b = filter->blocks[2];
  assert_eq(b.occupieds, 0);
  assert_eq(b.runends, 0);
  assert_eq(b.offset, 0);
  for (int i=0; i<64; i++) {
    assert_eq(b.remainders[i], 0);
  }
  // Check remote rep
  for (int i=0; i<64; i++) {
    assert_eq(filter->remote[128 + i].elt, 0);
    assert_eq(filter->remote[128 + i].hash, 0);
  }
  taf_destroy(filter);
  printf("passed.\n");
}

/// Check that adding a block doesn't overwrite existing data
void test_add_block_no_clobber() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(128);

  // Setup
  for (int i=0; i<filter->nslots; i++) {
    set_occupied(filter, i);
    set_runend(filter, i);
    remainder(filter, i) = i%16;
    filter->remote[i].elt = i;
    filter->remote[i].hash = i;
  }
  add_block(filter);
  // Check that data in first 2 blocks is preserved
  for (int i=0; i<128; i++) {
    assert(get_occupied(filter, i));
    assert(get_runend(filter, i));
    assert_eq(remainder(filter, i), i%16);
    assert_eq(filter->remote[i].elt, i);
    assert_eq(filter->remote[i].hash, i);
  }
  // Check that 3rd block is empty
  for (int i=128; i<filter->nslots; i++) {
    assert(!get_occupied(filter, i));
    assert(!get_runend(filter, i));
    assert_eq(remainder(filter, i), 0);
    assert_eq(filter->remote[i].elt, 0);
    assert_eq(filter->remote[i].hash, 0);
  }
  // Check filter metadata
  assert_eq(filter->nslots, 192);
  assert_eq(filter->nblocks, 3);

  taf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_loc_1() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(128);
  // q=7, r=8
  uint64_t code;
  int sels[64];
  code = get_sel_code(filter, 0);
  assert_eq(code, 0);
  decode_sel(code, sels);
  for (int i=0; i<64; i++) {
    assert_eq(sels[i], 0);
  }
  adapt_loc(filter, 0, sels);
  decode_sel(get_sel_code(filter, 0), sels);
  assert_eq(sels[0], 1);
  for (int i=1; i<64; i++) {
    assert_eq(sels[i], 0);
  }
  taf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_loc_2() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(128);
  // q=7, r=8
  int sels[64];
  for (int i=0; i<64; i++) {
    sels[i] = 0;
  }
  // Push encoding to limit (15 1's)
  int limit = 15;
  for (int i=0; i<limit; i++) {
    adapt_loc(filter, i, sels);
  }
  decode_sel(get_sel_code(filter, 0), sels);
  for (int i=0; i<limit; i++) {
    assert_eq(sels[i], 1);
  }
  for (int i=limit; i<64; i++) {
    assert_eq(sels[i], 0);
  }
  adapt_loc(filter, limit, sels);
  for (int i=0; i<64; i++) {
    test_assert_eq(sels[i], i == limit, "i=%d", i);
  }
  taf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_loc_3() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(128);
  // q=7, r=8
  int sels[64];
  for (int i=0; i<64; i++) {
    sels[i] = 0;
  }
  // Push encoding to limit (five 2's, one 1)
  int limit = 5;
  for (int i=0; i<5; i++) {
    adapt_loc(filter, i, sels);
    adapt_loc(filter, i, sels);
  }
  decode_sel(get_sel_code(filter, 0), sels);
  for (int i=0; i<limit; i++) {
    assert_eq(sels[i], 2);
  }
  for (int i=limit; i<64; i++) {
    assert_eq(sels[i], 0);
  }
  adapt_loc(filter, limit, sels);
  for (int i=0; i<64; i++) {
    test_assert_eq(sels[i], i < limit ? 2 : (i == limit ? 1 : 0), "i=%d", i);
  }
  adapt_loc(filter, limit, sels);
  for (int i=0; i<64; i++) {
    test_assert_eq(sels[i], i == limit ? 2 : 0, "i=%d", i);
  }
  taf_destroy(filter);
  printf("passed.\n");
}

void test_adapt_1() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(64);

  printf("[Unimplemented] ");

  taf_destroy(filter);
  printf("passed.\n");
}

void test_raw_lookup_1() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(128);

  printf(" [Unimplemented] ");

  taf_destroy(filter);
  printf("passed.\n");
}

void test_raw_insert_1() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(128);

  printf(" [Unimplemented] ");

  taf_destroy(filter);
  printf("passed.\n");
}

void test_shift_remote_elts() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(128);
  for (int i=0; i<filter->nslots; i++) {
    assert_eq(filter->remote[i].elt, 0);
    assert_eq(filter->remote[i].hash, 0);
  }
  for (int i=0; i<filter->nslots; i++) {
    filter->remote[i].elt = i;
    filter->remote[i].hash = i;
  }
  // Shift elts in [32, 64+32] to [33, 64+33]
  shift_remote_elts(filter, 32, 64+32);
  for (int i=0; i<=31; i++) {
    assert_eq(filter->remote[i].elt, i);
    assert_eq(filter->remote[i].hash, i);
  }
  assert_eq(filter->remote[32].elt, 0);
  assert_eq(filter->remote[32].hash, 0);
  for (int i=33; i<=64+33; i++) {
    assert_eq(filter->remote[i].elt, i-1);
    assert_eq(filter->remote[i].hash, i-1);
  }
  for (int i=64+34; i<filter->nslots; i++) {
    assert_eq(filter->remote[i].elt, i);
    assert_eq(filter->remote[i].hash, i);
  }
  taf_destroy(filter);
  printf("passed.\n");
}

double rand_zipfian(double s, double max) {
        double p = (double)rand() / RAND_MAX;

        double pD = p * (12 * (pow(max, -s + 1) - 1) / (1 - s) + 6 + 6 * pow(max, -s) + s - s * pow(max, -s + 1));
        double x = max / 2;
        while (1) {
                double m = pow(x, -s - 2);
                double mx = m * x;
                double mxx = mx * x;
                double mxxx = mxx * x;

                double b = 12 * (mxxx - 1) / (1 - s) + 6 + 6 * mxx + s - (s * mx) - pD;
                double c = 12 * mxx - (6 * s * mx) + (m * s * (s + 1));
                double newx = x - b / c > 1 ? x - b / c : 1;
                if (abs(newx - x) <= 0.01) { // this is the tolerance for approximation
                        return newx;
                }
                x = newx;
        }
}

uint64_t hash_str(char *str) {
        uint64_t hash = 5381;
        int c;
        while ((c = *str++)) {
                hash = ((hash << 5) + hash) + c;
        }
        return hash;
}

void csv_get(char* buffer, int col) {
        int i, j;
        for (i = 0; buffer[i] != '\0' && col > 0; i++) {
                if (buffer[i] == ',') col--;
        }
        for (j = 0; buffer[i + j] != '\0' && buffer[i + j] != ','; j++) {
                buffer[j] = buffer[i + j];
        }
        buffer[j] = '\0';
}

double avg_insert_time = 0;
double avg_query_time = 0;
double avg_fp_rate = 0;

/// General integration test: insert and query elts, ensuring that there
/// are no false negatives
void test_insert_and_query() {
  printf("Testing %s...", __FUNCTION__);
  size_t a = 1 << 20;
  double a_s = 100.0; // a/s
  double load = 0.9;
  size_t s = nearest_pow_of_2((size_t)((double)a / a_s));
  s = (size_t)((double)s * load);
  TAF* filter = new_taf(a);
  s = a * load;

  // Generate query set
  srandom(time(NULL));
  int nset = (int)(1.5*(double)s);
  Setnode* set = new Setnode[nset];

  elt_t *elts_to_insert = new elt_t[s];
  RAND_bytes((unsigned char*)elts_to_insert, s * sizeof(elt_t));

  char str[64];
	double measurement_interval = 0.05f;
        double current_measurement_point = measurement_interval;
        uint64_t next_measurement = current_measurement_point * a;
        clock_t end_time, start_time = clock();
        for (int i = 0; i < s; i++) {
		sprintf(str, "%lu", elts_to_insert[i]);
                set_insert(str, (int)strlen(str), 0, set, nset);
		taf_insert(filter, elts_to_insert[i]);
                if (i > next_measurement) {
                        printf("%f:\t%f\n", current_measurement_point, 1000000.0f * (measurement_interval * a) / (clock() - start_time));
                        current_measurement_point += measurement_interval;
                        next_measurement = current_measurement_point * a;
                        start_time = clock();
                }
        }
        printf("%f:\t%f\n", current_measurement_point, 1000000.0f * (measurement_interval * a) / (clock() - start_time));

  /*char str[64];
  clock_t start_time = clock();
  for (int i=0; i<s; i++) {
    elt_t elt = random();
    sprintf(str, "%lu", elt);
    set_insert(str, (int)strlen(str), 0, set, nset);
    taf_insert(filter, elt);
    //assert(set_lookup(str, (int)strlen(str), set, nset));
    //assert(taf_lookup(filter, elt));
  }
  clock_t end_time = clock();
  avg_insert_time += (double)(end_time - start_time) / s;*/

  int num_queries = 10000000;
  // Query [0, a] and ensure that all items in the set return true
  int fps = 0;
  int fns = 0;

  elt_t *query_set = new elt_t[num_queries];
  /*for (int i = 0; i < num_queries; i++) {
	  query_set[i] = rand_zipfian(1.5f, 1lu << 30);
  }*/
  RAND_bytes((unsigned char*)query_set, num_queries * sizeof(elt_t));

  start_time = clock();
  for (int i=0; i<num_queries; i++) {
	  elt_t elt;
    //elt = i;
    //elt = rand_zipfian(1.5f, 1lu << 30);
    //elt = random();
    elt = query_set[i];
    sprintf(str, "%lu", elt);
    int in_set = set_lookup(str, (int)strlen(str), set, nset);
    int in_taf = taf_lookup(filter, elt);
    if (in_set && !in_taf) {
      fns++;
      uint64_t hash = taf_hash(filter, elt);
      size_t quot = calc_quot(filter, hash);
      if (get_occupied(filter, quot)) {
        int loc = rank_select(filter, quot);
        if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu)"
                 " was occupied but didn't have an associated runend\n",
                 elt, elt, quot, quot/64, quot%64);
          print_taf_block(filter, quot/64);
          exit(1);
        } else {
          int sels[64];
          decode_sel(get_sel_code(filter, loc/64), sels);
          int sel = sels[loc%64];
          rem_t query_rem = calc_rem(filter, hash, sel);
          rem_t stored_rem = remainder(filter, loc);
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu),"
                 "loc=%d (block=%d, slot=%d); stored rem=0x%hhx doesn't match query rem=0x%hhx\n",
                 elt, elt, quot, quot/64, quot%64, loc, loc/64, loc%64, stored_rem, query_rem);
          print_taf_metadata(filter);
          print_taf_block(filter, loc/64);
          exit(1);
        }
      } else {
        printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu) wasn't occupied\n",
               elt, elt, quot, quot/64, quot%64);
        exit(1);
      }
    } else if (!in_set && in_taf) {
      //fps += in_taf;
      fps++;
    }
    /*if (i % 100 == 0) {
	    printf("%d,%f\n", i, (double)fps / i);
    }*/
  }
  end_time = clock();
  avg_query_time += (double)(end_time - start_time) / num_queries;
  avg_fp_rate += (double)(fps) / num_queries;

  printf("passed. ");
  printf("FPs: %d (%f%%), FNs: %d (%f%%)\n",
         fps, (double)fps/(double)num_queries, fns, (double)fns/(double)a * 100);
  print_taf_metadata(filter);
  taf_destroy(filter);
}

void test_insert_and_query_w_repeats(elt_t *query_set, int query_set_size, int n_queries, int step_size, double *fprates, int iter) {
  printf("Testing %s...\n", __FUNCTION__);
  int nslots = 1 << 20;
  double load = 0.95;
  double a_s = 100;
  int queries_per_elt = 10;
  printf("nslots=%d, load=%f, a/s=%f, queries_per_elt = %d\n", nslots, load, a_s, queries_per_elt);

  int s = (int)((double)nearest_pow_of_2(nslots) * load);
  int a = (int)((double)s * a_s);
  printf("%d\n", a);
  //int n_queries = a * queries_per_elt;
  //n_queries = 1000000;

  int fps = 0;  // false positives
  int rfps = 0; // repeated false positives
  int fns = 0;  // false negatives
  int tot_queries = n_queries * queries_per_elt;

  TAF *filter = new_taf(s);
  int nset = (int)(s * 1.5);
  Setnode *set = new Setnode[nset];

  srandom(time(NULL));
  char str[64];
  int len;
  fprintf(stderr, "Initializing membership set and filter...\n");
  for (int i = 0; i < s; i++) {
    elt_t elt = random();
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    set_insert(str, len, 0, set, nset);
    taf_insert(filter, elt);
  }
  /*fprintf(stderr, "Initializing query set...\n");
  FILE *caida = fopen("../../../aqf/AdaptiveQF/data/20140619-140100.csv", "r");
  char buffer[1024];
  fgets(buffer, sizeof(buffer), caida);
  elt_t *query_set = calloc(a, sizeof(elt_t));
  for (int i=0; i<a; i++) {
    //query_set[i] = random();
    //query_set[i] = rand_zipfian(1.5f, 1lu << 30);
    fgets(buffer, sizeof(buffer), caida);
    csv_get(buffer, 3);
    query_set[i] = hash_str(buffer);
  }
  fclose(caida);*/
  fprintf(stderr, "Querying set and filter...\n");
  int nseen = (int)(s * 1.5);
  Setnode *seen = new Setnode[nseen];
  for (int i=0; i<n_queries; i++) {
    elt_t elt = query_set[random() % query_set_size];
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    int in_filter = taf_lookup(filter, elt);
    int in_set = set_lookup(str, len, set, nset);
    if (in_filter && !in_set) {
      fps++;
      if (set_lookup(str, len, seen, nseen)) {
        rfps++;
      } else {
        set_insert(str, len, 0, seen, nseen);
      }
    } else if (!in_filter && in_set) {
      fns++;
      uint64_t hash = taf_hash(filter, elt);
      size_t quot = calc_quot(filter, hash);
      if (get_occupied(filter, quot)) {
        int loc = rank_select(filter, quot);
        if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu)"
                 " was occupied but didn't have an associated runend\n",
                 elt, elt, quot, quot/64, quot%64);
          print_taf_block(filter, quot/64);
          exit(1);
        } else {
          int sels[64];
          decode_sel(get_sel_code(filter, loc/64), sels);
          int sel = sels[loc%64];
          rem_t query_rem = calc_rem(filter, hash, sel);
          rem_t stored_rem = remainder(filter, loc);
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu),"
                 "loc=%d (block=%d, slot=%d); stored rem=0x%hhx doesn't match query rem=0x%hhx\n",
                 elt, elt, quot, quot/64, quot%64, loc, loc/64, loc%64, stored_rem, query_rem);
          print_taf_metadata(filter);
          print_taf_block(filter, loc/64);
          exit(1);
        }
      } else {
        printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu) wasn't occupied\n",
               elt, elt, quot, quot/64, quot%64);
        exit(1);
      }
    }
    if (i % 100 == 0) {
	    fprates[i / 100] += (double)fps / i;
    }
  }
  printf("Test results:\n");
  printf("FPs: %d (%f%%), RFPs: %d (%f%%)\n",
         fps, (double)fps/tot_queries, rfps, (double)rfps/tot_queries * 100);
  printf("FNs: %d (%f%%)\n", fns, (double)fns/tot_queries * 100);

  taf_destroy(filter);
  printf("Done testing %s.\n", __FUNCTION__);
}

void test_dataset_evolution(char* input_file_name, char* output_file_name, int num_trials, int query_space_size, int num_queries, int step_size) {
	FILE* outfp = fopen(output_file_name, "w");
	fclose(outfp);
	outfp = fopen(output_file_name, "a");
	char buffer[1024];
	//FILE* infp = fopen(input_file_name, "r");
	//fgets(buffer, sizeof(buffer), infp);

	elt_t *query_set = new elt_t[query_space_size];
	for (int i = 0; i < query_space_size; i++) {
		/*fgets(buffer, sizeof(buffer), infp);
		csv_get(buffer, 3);
		query_set[i] = hash_str(buffer);*/
		query_set[i] = rand();
	}

	double *fprates = new double[num_queries / 100];

	for (int i = 0; i < num_trials; i++) {
		test_insert_and_query_w_repeats(query_set, query_space_size, num_queries, step_size, fprates, i);
	}

	for (int i = 0; i < num_queries / step_size; i++) {
		sprintf(buffer, "%d,%f\n", i * step_size, fprates[i] / num_trials);
		fputs(buffer, outfp);
	}
	
	//fclose(infp);
	fclose(outfp);
	free(query_set);
	free(fprates);
}

void test_hash_accesses(int qbits, int rbits, double load, uint64_t num_queries, uint64_t seed) {
	if (seed == -1) seed = time(NULL);
	printf("testing hash accesses on seed %lu\n", seed);
	srand(seed);

	uint64_t nslots = 1ull << qbits;
	//uint64_t xnslots = nslots + 10 * sqrt(nslots);

	uint64_t num_inserts = nslots * load;

	int fps = 0;  // false positives
	int rfps = 0; // repeated false positives
	int fns = 0;
	int tps = 0;
	int tns = 0;
	int hash_accesses = 0;
	int negatives = 0;

	TAF *filter = new_taf(nslots);
	int nset = 1.3 * num_inserts;
	Setnode *set = new Setnode[nset];
	elt_t *inserts = new elt_t[num_inserts];

	char str[64];
	int len;
	printf("starting %lu inserts\n", num_inserts);
	for (int i = 0; i < num_inserts; i++) {
		elt_t elt = rand();
		sprintf(str, "%lu", elt);
		len = (int)strlen(str);

		set_insert(str, len, i, set, nset);
		taf_insert(filter, elt);
		inserts[i] = elt;
	}

	FILE *fp = fopen("target/hash_accesses.txt", "w");
	fprintf(fp, "queries\taccesses\tfps\trfps\ttps\tnegatives\tfns\ttns\n");
	fclose(fp);
	fp = fopen("target/hash_accesses.txt", "a");

	int nseen = 1.3 * num_inserts;
	Setnode *seen = new Setnode[nseen];
	printf("starting %lu queries\n", num_queries);
	for (int i = 0; i < num_queries; i++) {
		elt_t elt = rand_zipfian(1.5f, 1ull << 30);
		sprintf(str, "%lu", elt);
		len = (int)strlen(str);

		int in_filter = taf_lookup(filter, elt);
		int in_set = set_lookup(str, len, set, nset);

		if (in_filter) {
			hash_accesses++;
			if (in_set) {
				tps++;
				if (inserts[in_set - 1] != elt) {
					printf("original insert was %lu but set was triggered by %lu\n", inserts[in_set - 1], elt);
				}
			}
		}
		else {
			negatives++;
			if (!in_set) tns++;
		}
		if (in_filter && !in_set) {
			fps++;
			if (set_lookup(str, len, seen, nseen)) {
				rfps++;
			} else {
				set_insert(str, len, 0, seen, nseen);
			}
		} else if (!in_filter && in_set) {
			fns++;
			uint64_t hash = taf_hash(filter, elt);
			size_t quot = calc_quot(filter, hash);
			if (get_occupied(filter, quot)) {
				int loc = rank_select(filter, quot);
				if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
					printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu)"
							" was occupied but didn't have an associated runend\n",
							elt, elt, quot, quot/64, quot%64);
					print_taf_block(filter, quot/64);
					exit(1);
				} else {
					int sels[64];
					decode_sel(get_sel_code(filter, loc/64), sels);
					int sel = sels[loc%64];
					rem_t query_rem = calc_rem(filter, hash, sel);
					rem_t stored_rem = remainder(filter, loc);
					printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu),"
							"loc=%d (block=%d, slot=%d); stored rem=0x%hhx doesn't match query rem=0x%hhx\n",
							elt, elt, quot, quot/64, quot%64, loc, loc/64, loc%64, stored_rem, query_rem);
					print_taf_metadata(filter);
					print_taf_block(filter, loc/64);
					exit(1);
				}
			} else {
				printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu) wasn't occupied\n",
						elt, elt, quot, quot/64, quot%64);
				exit(1);
			}
		}

		if (i % 100000 == 0) {
			fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, hash_accesses, fps, rfps, tps, negatives, fns, tns);
		}
	}

	fclose(fp);

	printf("Test results:\n");
	printf("FPs: %d (%f%%), RFPs: %d (%f%%)\n",
			fps, (double)fps/num_queries, rfps, (double)rfps/num_queries * 100);
	printf("FNs: %d (%f%%)\n", fns, (double)fns/num_queries * 100);

	taf_destroy(filter);
	printf("Done testing %s.\n", __FUNCTION__);
}


void test_mixed_insert_and_query_w_repeats() {
  printf("Testing %s...\n", __FUNCTION__);
  int nslots = 1 << 14;
  double load = 0.95;
  double a_s = 100;
  int queries_per_elt = 10;
  printf("nslots=%d, load=%f, a/s=%f, queries_per_elt = %d\n", nslots, load, a_s, queries_per_elt);

  int s = (int)((double)nearest_pow_of_2(nslots) * load);
  int a = (int)((double)s * a_s);
  int n_queries = a * queries_per_elt;

  int fps = 0;  // false positives
  int rfps = 0; // repeated false positives
  int fns = 0;  // false negatives
  int tot_queries = n_queries * queries_per_elt;

  TAF *filter = new_taf(s);
  int nset = (int)(s * 1.5);
  Setnode *set = new Setnode[nset];

  srandom(TAF_SEED);
  char str[64];
  int len;
  fprintf(stderr, "Initializing query set...\n");
  elt_t *query_set = new elt_t[a];
  for (int i=0; i<a; i++) {
    query_set[i] = random();
  }
  fprintf(stderr, "Initializing membership set and filter...\n");
  for (int i=0; i<s/2; i++) {
    elt_t elt = random();
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    set_insert(str, len, 0, set, nset);
    taf_insert(filter, elt);
  }
  fprintf(stderr, "Performing interleaved queries...\n");
  for (int i=0; i<n_queries; i++) {
    elt_t elt = query_set[random() % a];
    taf_lookup(filter, elt);
  }
  fprintf(stderr, "Finishing initialization of membership set and filter...\n");
  for (int i=s/2; i<s; i++) {
    elt_t elt = random();
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    set_insert(str, len, 0, set, nset);
    taf_insert(filter, elt);
  }
  fprintf(stderr, "Querying set and filter...\n");
  int nseen = (int)(s * 1.5);
  Setnode *seen = new Setnode[nseen];
  for (int i=0; i<n_queries; i++) {
    elt_t elt = query_set[random() % a];
    sprintf(str, "%lu", elt);
    len = (int)strlen(str);
    int in_filter = taf_lookup(filter, elt);
    int in_set = set_lookup(str, len, set, nset);
    if (in_filter && !in_set) {
      fps++;
      if (set_lookup(str, len, seen, nseen)) {
        rfps++;
      } else {
        set_insert(str, len, 0, seen, nseen);
      }
    } else if (!in_filter && in_set) {
      fns++;
      uint64_t hash = taf_hash(filter, elt);
      size_t quot = calc_quot(filter, hash);
      if (get_occupied(filter, quot)) {
        int loc = rank_select(filter, quot);
        if (loc == RANK_SELECT_EMPTY || loc == RANK_SELECT_OVERFLOW) {
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu)"
                 " was occupied but didn't have an associated runend\n",
                 elt, elt, quot, quot/64, quot%64);
          print_taf_block(filter, quot/64);
          exit(1);
        } else {
          int sels[64];
          decode_sel(get_sel_code(filter, loc/64), sels);
          int sel = sels[loc%64];
          rem_t query_rem = calc_rem(filter, hash, sel);
          rem_t stored_rem = remainder(filter, loc);
          printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu),"
                 "loc=%d (block=%d, slot=%d); stored rem=0x%hhx doesn't match query rem=0x%hhx\n",
                 elt, elt, quot, quot/64, quot%64, loc, loc/64, loc%64, stored_rem, query_rem);
          print_taf_metadata(filter);
          print_taf_block(filter, loc/64);
          exit(1);
        }
      } else {
        printf("False negative (elt=%lu, 0x%lx): quot=%lu (block=%lu, slot=%lu) wasn't occupied\n",
               elt, elt, quot, quot/64, quot%64);
        exit(1);
      }
    }
  }
  printf("Test results:\n");
  printf("FPs: %d (%f%%), RFPs: %d (%f%%)\n",
         fps, (double)fps/tot_queries, rfps, (double)rfps/tot_queries * 100);
  printf("FNs: %d (%f%%)\n", fns, (double)fns/tot_queries * 100);
  print_taf_stats(filter);
  taf_destroy(filter);
  printf("Done testing %s.\n", __FUNCTION__);
}

// a, b+1 in the same block
void test_shift_sels_single_block() {
  printf("Testing %s...", __FUNCTION__);

  // Setup
  TAF *filter = new_taf(64 * 3);
  int sels[64];
  for (int i=0; i<64; i++) {
    sels[i] = i % 8 == 0;
  }
  uint64_t code;
  encode_sel(sels, &code);
  set_sel_code(filter, 0, code);

  // Shift sels in [0, 62] -> [1, 63] and check
  for (int j=1; j <= 64; j++) {
    shift_sels(filter, 0, 62);
    decode_sel(get_sel_code(filter, 0), sels);
    for (int i=0; i<64; i++) {
      test_assert_eq(sels[i],
                     (i - j) % 8 == 0 && i > j-1,
                     "j=%d, i=%d", j, i);
    }
  }
  taf_destroy(filter);
  printf("passed.\n");
}

void test_swap_sels() {
  printf("Testing %s...", __FUNCTION__);
  int* xs = new int[64];
  int* ys = new int[64];
  for (int i=0; i<64; i++) {
    xs[i] = i;
    ys[i] = 64 + i;
  }
  swap_ptrs(&xs, &ys);
  for (int i=0; i<64; i++) {
    assert_eq(ys[i], i);
    assert_eq(xs[i], 64+i);
  }
  printf("passed.\n");
}

TAF* sel_setup() {
  TAF* filter = new_taf(64 * 4);
  int sels[64];
  uint64_t code;
  for (int i=0; i<filter->nblocks; i++) {
    for (int j=0; j<64; j++) {
      sels[j] = j % 8 == 0;
    }
    encode_sel(sels, &code);
    set_sel_code(filter, i, code);
  }
  return filter;
}

// a, b+1 in different blocks
void test_shift_sels_multi_block() {
  printf("Testing %s...", __FUNCTION__);
  int sels[64];
  TAF *filter = sel_setup();
  // (1) Shift sels in [0, 127] -> [1,128]
  shift_sels(filter, 0, 127);
  for (int i=0; i<filter->nblocks; i++) {
    decode_sel(get_sel_code(filter, i), sels);
    for (int j=0; j<64; j++) {
      test_assert_eq(sels[j],
                     (i < 2) ? ((j-1)%8 == 0) : (i == 2 && j == 0 ? 0 : j%8 == 0),
                     "i=%d, j=%d", i, j);
    }
  }
  taf_destroy(filter);
  filter = sel_setup();
  // (2) Shift sels in [32, 64+32] -> [33, 64+33]
  shift_sels(filter, 32, 64+32);
  for (int i=0; i<filter->nslots; i++) {
    if (i%64 == 0) {
      decode_sel(get_sel_code(filter, i/64), sels);
    }
    test_assert_eq(sels[i%64],
                   i < 32 ? i%8 == 0 :
                   (i == 32 ? 0 :
                    (i <= 64 + 33) ? (i-1)%8 == 0 : i%8 == 0),
                   "i=%d", i);
  }
  taf_destroy(filter);
  printf("passed.\n");
}

void test_template() {
  printf("Testing %s...", __FUNCTION__);
  TAF *filter = new_taf(64 * 3);

  printf("[unimplemented]");

  taf_destroy(filter);
  printf("passed.\n");
}

int main(int argc, char **argv) {
	TAF *filter = new_taf(18);
	ordered_map_t backing_map((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
	ordered_map_t database((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);

	uint64_t num_slots = 1ULL << 18;
	uint64_t num_inserts = num_slots * 0.95;
	uint64_t *insert_set = new uint64_t[num_inserts];
	RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));

	double measure_interval = 0.05f, current_interval = measure_interval;
	uint64_t next_point = measure_interval * num_slots, last_point = 0;
	clock_t start_time = clock(), interval_time = start_time, end_time;
	for (size_t i = 0; i < num_inserts; i++) {
		taf_insert(filter, insert_set[i]);
		backing_map.insert(pair_t(insert_set[i], i));
		database.insert(pair_t(insert_set[i], i));
		if (i > next_point) {
			printf("interval %f - %f ops/sec\n", current_interval, (double)CLOCKS_PER_SEC * (next_point - last_point) / (clock() - interval_time));
			current_interval += measure_interval;
			last_point = next_point;
			next_point = current_interval * num_slots;
			interval_time = clock();
		}
	}
	printf("interval %f - %f ops/sec\n", current_interval, (double)CLOCKS_PER_SEC * (num_inserts - last_point) / (clock() - interval_time));

	taf_destroy(filter);
}
