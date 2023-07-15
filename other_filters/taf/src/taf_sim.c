#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <execinfo.h>

#include <time.h>
#include <sys/time.h>
#include <openssl/rand.h>

#include "murmur3.h"
#include "macros.h"
#include "arcd.h"
#include "taf.h"
#include "bit_util.h"
#include "set.h"

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
    int* sels = malloc(64 * sizeof(int));
    int* prev_sels = malloc(64 * sizeof(int));
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
  TAFBlock *new_blocks = realloc(filter->blocks, (filter->nblocks + 1) * sizeof(TAFBlock));
  if (new_blocks == NULL) {
    printf("add_block failed to realloc new blocks\n");
    exit(1);
  }
  filter->blocks = new_blocks;
  memset(filter->blocks + filter->nblocks, 0, sizeof(TAFBlock));

  // Reallocate remote rep
  Remote_elt *new_remote = realloc(filter->remote,(filter->nslots + 64) * sizeof(Remote_elt));
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
  filter->blocks = calloc(filter->nblocks, sizeof(TAFBlock));
  filter->remote = calloc(filter->nslots, sizeof(Remote_elt));
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
  filter->blocks = calloc(filter->nblocks, sizeof(TAFBlock));
  filter->remote = calloc(filter->nslots, sizeof(Remote_elt));
}

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
	  return -1;
        }
        return 1;
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
//#define TEST_TAF_SIM 1
#ifdef TEST_TAF_SIM

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
  TAF *filter = malloc(sizeof(TAF));
  taf_init(filter, n, TAF_SEED);
  return filter;
}

double rand_zipfian(double s, double max, uint64_t source) {
        double p = (double)source / (-1ULL);

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

void test_splinter_throughput() {
        printf("Testing %s...\n", __FUNCTION__);
        size_t num_slots = 1ull << 26;
        double load = 0.9;

        size_t xnum_slots = num_slots + (10 * sqrt(num_slots));

        TAF *filter = new_taf(xnum_slots);

        printf("Generating inserts\n");
        uint64_t num_inserts = num_slots * load;
        uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
        RAND_bytes((void*)insert_set, num_inserts * sizeof(uint64_t));

        printf("Performing inserts\n");
        uint64_t num_segments = 100;
        uint64_t measure_point = num_slots / num_segments - 1;
        uint64_t cur_segment = 0;

        clock_t start_clock = clock(), end_clock;
        struct timeval timecheck;
        gettimeofday(&timecheck, NULL);
        uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;

        for (uint64_t i = 0; i < num_inserts; i++) {
                taf_insert(filter, insert_set[i]);

                if (i >= measure_point) {
                        fprintf(stderr, "\rload factor: %d%%", (int)((i + 1) * 100 / num_slots));

                        cur_segment++;
                        measure_point = (cur_segment + 1) * num_slots / num_segments - 1;
                }
        }

        gettimeofday(&timecheck, NULL);
        end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
        end_clock = clock();

        printf("\n");
        printf("time for inserts:     %f sec\n", (double)(end_time - start_time) / 1000000);
        printf("insert throughput:    %f ops/sec\n", (double)num_inserts * 1000000 / (end_time - start_time));
        printf("cpu time for inserts: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	printf("cpu insert through:   %f ops/sec\n", (double)num_inserts * CLOCKS_PER_SEC / (end_clock - start_clock));

        printf("Generating queries\n");
        uint64_t num_queries = 100000000ull;

        uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
        RAND_bytes((void*)query_set, num_queries * sizeof(uint64_t));
	/*int murmur_seed = rand();
	for (uint64_t i = 0; i < num_queries; i++) {
		query_set[i] = (uint64_t)rand_zipfian(1.5f, 10000000ull, query_set[i]);
		query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
	}*/

        printf("Performing queries\n");
        measure_point = num_queries / num_segments - 1;
        cur_segment = 0;

        uint64_t fp_count = 0;

        start_clock = clock();
        gettimeofday(&timecheck, NULL);
	start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
        for (uint64_t i = 0; i < num_queries; i++) {
                int ret = taf_lookup(filter, query_set[i]);
                if (ret == -1) {
                        fp_count++;
                }

                if (i >= measure_point) {
                        fprintf(stderr, "\rqueries done: %d%%", (int)((i + 1) * 100 / num_queries));

                        cur_segment++;
                        measure_point = (cur_segment + 1) * num_queries / num_segments - 1;
                }
        }

        gettimeofday(&timecheck, NULL);
        end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
        end_clock = clock();

        printf("\n");
        printf("time for queries:     %f sec\n", (double)(end_time - start_time) / 1000000);
        printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
        printf("cpu time for queries: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	printf("cpu query through:    %f ops/sec\n", (double)num_queries * CLOCKS_PER_SEC / (end_clock - start_clock));

        printf("false positive rate:  %f%%\n", (double)fp_count / num_queries * 100);

        free(insert_set);
        free(query_set);

        taf_destroy(filter);
}

int main() {
	test_splinter_throughput();
}
#endif // TEST_TAF_SIM
