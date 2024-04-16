#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <execinfo.h>

#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <openssl/rand.h>

#include <splinterdb/splinterdb.h>
#include <splinterdb/data.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/default_data_config.h>

#include "murmur3.h"
#include "macros.h"
#include "arcd.h"
#include "splinter_taf.h"
#include "bit_util.h"
#include "set.h"


#define TEST_DB_NAME "db"

#define Kilo (1024UL)
#define Mega (1024UL * Kilo)
#define Giga (1024UL * Mega)

#define MAX_KEY_SIZE 16
#define MAX_VAL_SIZE 16

void csv_get_col(char* buffer, int col) {
	int i, j;
	for (i = 0; buffer[i] != '\0' && col > 0; i++) {
		if (buffer[i] == ',') col--;
	}
	for (j = 0; buffer[i + j] != '\0' && buffer[i + j] != ','; j++) {
		buffer[j] = buffer[i + j];
	}
	buffer[j] = '\0';
}

uint64_t hash_str(char *str) {
	uint64_t hash = 5381;
	int c;
	while ((c = *str++)) {
		hash = ((hash << 5) + hash) + c;
	}
	return hash;
}

uint64_t db_insert_count = 0, db_update_count = 0, db_query_count = 0;

struct backing_data {
	uint64_t elt;
	uint64_t hash;
} typedef backing_data;

void pad_data(void *dest, const void *src, const size_t dest_len, const size_t src_len) {
	assert(dest_len >= src_len);
	bzero(dest, dest_len);
	memcpy(dest, src, src_len);
}

slice padded_slice(const void *data, const size_t dest_len, const size_t src_len, void *buffer) {
	pad_data(buffer, data, dest_len, src_len);

	return slice_create(dest_len, buffer);
}

void db_insert(TAF *filter, const uint64_t loc, const uint64_t key, const uint64_t hash) {
	char key_padded[MAX_KEY_SIZE];
	char val_padded[MAX_VAL_SIZE];
	bzero(key_padded, MAX_KEY_SIZE);
	bzero(val_padded, MAX_VAL_SIZE);
	memcpy(key_padded, &loc, sizeof(loc));
	memcpy(val_padded, &key, sizeof(key));
	memcpy(val_padded + sizeof(key), &hash, sizeof(hash));

	slice key_slice = slice_create(MAX_KEY_SIZE, key_padded);
	slice val_slice = slice_create(MAX_VAL_SIZE, val_padded);

	splinterdb_insert(filter->database, key_slice, val_slice);
	db_insert_count++;
}

void db_update(TAF *filter, const uint64_t loc, const uint64_t key, const uint64_t hash) {
	char buffer[MAX_KEY_SIZE];
	slice loc_slice = padded_slice(&loc, MAX_KEY_SIZE, sizeof(loc), buffer);

	splinterdb_delete(filter->database, loc_slice);

	char val_padded[MAX_VAL_SIZE];
	bzero(val_padded, MAX_VAL_SIZE);
	memcpy(val_padded, &key, sizeof(key));
	memcpy(val_padded + sizeof(key), &hash, sizeof(hash));

	slice val_slice = slice_create(MAX_VAL_SIZE, val_padded);

	splinterdb_insert(filter->database, loc_slice, val_slice);
	db_update_count++;
}

backing_data db_query(TAF *filter, const uint64_t loc) {
	char buffer[MAX_KEY_SIZE];
	slice loc_slice = padded_slice(&loc, MAX_KEY_SIZE, sizeof(loc), buffer);
	splinterdb_lookup(filter->database, loc_slice, &filter->db_result);
	db_query_count++;

	backing_data data;
	if (!splinterdb_lookup_found(&filter->db_result)) {
		data.elt = data.hash = 0;
		return data;
	}

	slice result_val;
	splinterdb_lookup_result_value(&filter->db_result, &result_val);

	memcpy(&data.elt, slice_data(result_val), sizeof(uint64_t));
	memcpy(&data.hash, slice_data(result_val) + sizeof(uint64_t), sizeof(uint64_t));
	return data;
}


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
	//return hash % (1ull << filter->q);
}

/**
 * Returns the k-th remainder for h
 */
static rem_t calc_rem(const TAF* filter, uint64_t hash, int k) {
	int n_rems = (64 - (int)filter->q)/(int)filter->r;
	if (k >= n_rems) k %= n_rems;
	return (hash >> (filter->q + k * filter->r)) & ONES(filter->r);
	//return (hash >> (filter->q + k * filter->r)) % (1ull << filter->r);
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
		backing_data data = db_query(filter, i);
		db_update(filter, i + 1, data.elt, data.hash);
	}
	db_update(filter, a, 0, 0);
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
			//backing_data orig_data = db_query(filter, b_start + i);
			b->remainders[i] = calc_rem(filter, db_query(filter, b_start + i).hash, 0);
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
	rem_t new_rem = calc_rem(filter, db_query(filter, loc).hash, new_sel);
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
		if (db_query(filter, i).elt == query) {
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
	filter->nblocks += 10 * sqrt(filter->nslots) / 64;
	filter->q = (size_t)log2((double)filter->nslots); // nslots = 2^q
	filter->r = REM_SIZE;
	filter->p = filter->q + filter->r;
	filter->blocks = calloc(filter->nblocks, sizeof(TAFBlock));
	filter->mode = TAF_MODE_NORMAL;

        default_data_config_init(MAX_KEY_SIZE, &filter->data_cfg);
	filter->splinterdb_cfg = (splinterdb_config){
		.filename   = TEST_DB_NAME,
		.cache_size = 32 * Mega,
		.disk_size  = 20 * Giga,
		.data_cfg   = &filter->data_cfg,
		.io_flags   = O_RDWR | O_CREAT | O_DIRECT
	};
	splinterdb_create(&filter->splinterdb_cfg, &filter->database);
	splinterdb_lookup_result_init(filter->database, &filter->db_result, 0, NULL);
}

void taf_destroy(TAF* filter) {
	free(filter->blocks);
	splinterdb_lookup_result_deinit(&filter->db_result);
	splinterdb_close(&filter->database);
	free(filter);
}

void taf_clear(TAF* filter) {
	filter->nelts = 0;
	free(filter->blocks);
	splinterdb_close(&filter->database);
	filter->blocks = calloc(filter->nblocks, sizeof(TAFBlock));
	splinterdb_create(&filter->splinterdb_cfg, &filter->database);
	splinterdb_lookup_result_init(filter->database, &filter->db_result, 0, NULL);
}

static void raw_insert(TAF* filter, elt_t elt, uint64_t hash) {
	size_t quot = calc_quot(filter, hash);
	rem_t rem = calc_rem(filter, hash, 0);
	filter->nelts++;

	// Find the appropriate runend
	int r = rank_select(filter, quot);
	switch (r) {
		case RANK_SELECT_EMPTY: { // empty slot
						set_occupied(filter, quot);
						set_runend(filter, quot);
						remainder(filter, quot) = rem;
						db_insert(filter, quot, elt, hash);
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
				 db_insert(filter, r + 1, elt, hash);
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
				if (elt != db_query(filter, loc).elt) {
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
			printf(" 0x%-*lx", 8, db_query(filter, block_index * 64 + i*8 + j).elt);
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
//#define TEST_SPLINTER_TAF 1
#ifdef TEST_SPLINTER_TAF

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
		backing_data data = db_query(filter, 128 + i);
		assert_eq(data.elt, 0);
		assert_eq(data.hash, 0);
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
		db_insert(filter, i, i, i);
	}
	add_block(filter);
	// Check that data in first 2 blocks is preserved
	for (int i=0; i<128; i++) {
		assert(get_occupied(filter, i));
		assert(get_runend(filter, i));
		assert_eq(remainder(filter, i), i%16);
		backing_data data = db_query(filter, i);
		assert_eq(data.elt, i);
		assert_eq(data.hash, i);
	}
	// Check that 3rd block is empty
	for (int i=128; i<filter->nslots; i++) {
		assert(!get_occupied(filter, i));
		assert(!get_runend(filter, i));
		assert_eq(remainder(filter, i), 0);
		backing_data data = db_query(filter, i);
		assert_eq(data.elt, 0);
		assert_eq(data.hash, 0);
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
		assert_eq(db_query(filter, i).elt, 0);
		assert_eq(db_query(filter, i).hash, 0);
	}
	for (int i=0; i<filter->nslots; i++) {
		db_insert(filter, i, i, i);
	}
	// Shift elts in [32, 64+32] to [33, 64+33]
	shift_remote_elts(filter, 32, 64+32);
	for (int i=0; i<=31; i++) {
		assert_eq(db_query(filter, i).elt, i);
		assert_eq(db_query(filter, i).hash, i);
	}
	assert_eq(db_query(filter, 32).elt, 0);
	assert_eq(db_query(filter, 32).hash, 0);
	for (int i=33; i<=64+33; i++) {
		assert_eq(db_query(filter, i).elt, i-1);
		assert_eq(db_query(filter, i).hash, i-1);
	}
	for (int i=64+34; i<filter->nslots; i++) {
		assert_eq(db_query(filter, i).elt, i);
		assert_eq(db_query(filter, i).hash, i);
	}
	taf_destroy(filter);
	printf("passed.\n");
}

/// General integration test: insert and query elts, ensuring that there
/// are no false negatives
void test_insert_and_query() {
	printf("Testing %s...", __FUNCTION__);
	size_t a = 1 << 20;
	double a_s = 100.0; // a/s
	double load = 0.95;
	size_t s = nearest_pow_of_2((size_t)((double)a / a_s));
	s = (size_t)((double)s * load);
	TAF* filter = new_taf(s);

	// Generate query set
	srandom(TAF_SEED);
	int nset = (int)(1.5*(double)s);
	Setnode* set = calloc(nset, sizeof(set[0]));
	char str[64];
	for (int i=0; i<s; i++) {
		elt_t elt = random() % a;
		sprintf(str, "%lu", elt);
		set_insert(str, (int)strlen(str), 0, set, nset);
		taf_insert(filter, elt);
		assert(set_lookup(str, (int)strlen(str), set, nset));
		assert(taf_lookup(filter, elt));
	}
	// Query [0, a] and ensure that all items in the set return true
	int fps = 0;
	int fns = 0;
	for (int i=0; i<a; i++) {
		elt_t elt = i;
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
			fps += in_taf;
		}
	}
	printf("passed. ");
	printf("FPs: %d (%f%%), FNs: %d (%f%%)\n",
			fps, (double)fps/(double)a * 100, fns, (double)fns/(double)a * 100);
	print_taf_metadata(filter);
	taf_destroy(filter);
}

void test_insert_and_query_w_repeats() {
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
	Setnode *set = calloc(nset, sizeof(set[0]));

	srandom(TAF_SEED);
	char str[64];
	int len;
	fprintf(stderr, "Initializing membership set and filter...\n");
	for (int i=0; i<s; i++) {
		elt_t elt = random();
		sprintf(str, "%lu", elt);
		len = (int)strlen(str);
		set_insert(str, len, 0, set, nset);
		taf_insert(filter, elt);
	}
	fprintf(stderr, "Initializing query set...\n");
	elt_t *query_set = calloc(a, sizeof(elt_t));
	for (int i=0; i<a; i++) {
		query_set[i] = random();
	}
	fprintf(stderr, "Querying set and filter...\n");
	int nseen = (int)(s * 1.5);
	Setnode *seen = calloc(nseen, sizeof(seen[0]));
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
	Setnode *set = calloc(nset, sizeof(set[0]));

	srandom(TAF_SEED);
	char str[64];
	int len;
	fprintf(stderr, "Initializing query set...\n");
	elt_t *query_set = calloc(a, sizeof(elt_t));
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
	Setnode *seen = calloc(nseen, sizeof(seen[0]));
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
	int* xs = malloc(64 * sizeof(int));
	int* ys = malloc(64 * sizeof(int));
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

void test_splinter_ops() {
	data_config data_cfg;
	default_data_config_init(8, &data_cfg);
	splinterdb_config splinterdb_cfg = (splinterdb_config){
		.filename   = TEST_DB_NAME,
		.cache_size = 64 * Mega,
		.disk_size  = 127 * Mega,
		.data_cfg   = &data_cfg
	};

	splinterdb *database;
	int rc = splinterdb_create(&splinterdb_cfg, &database);
	assert(rc == 0);

	char key_val_data[16];
	RAND_bytes((void*)key_val_data, 16);

	slice key = slice_create(8, key_val_data);
	slice val = slice_create(8, key_val_data + 8);

	rc = splinterdb_insert(database, key, val);
	assert(rc == 0);

	splinterdb_lookup_result result;
	splinterdb_lookup_result_init(database, &result, 0, NULL);

	rc = splinterdb_lookup(database, key, &result);
	assert(rc == 0);
	assert(splinterdb_lookup_found(&result));
	slice result_val;
	splinterdb_lookup_result_value(&result, &result_val);
	assert(slice_length(result_val) == 8);
	assert(strncmp(slice_data(result_val), key_val_data + 8, 8) == 0);

	rc = splinterdb_delete(database, key);
	assert(rc == 0);

	rc = splinterdb_lookup(database, key, &result);
	assert(rc == 0);
	assert(!splinterdb_lookup_found(&result));

	splinterdb_lookup_result_deinit(&result);
	splinterdb_close(&database);

	printf("Test completed\n");
}

void test_splinter_veracity() {
	printf("Testing %s...\n", __FUNCTION__);
	size_t num_slots = 1ull << 20;
	double load = 0.85;

	//size_t xnum_slots = num_slots + (10 * sqrt(num_slots));

	TAF *filter = new_taf(num_slots);

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

	printf("Verifying inserts\n");
	for (uint64_t i = 0; i < num_inserts; i++) {
		assert(taf_lookup(filter, insert_set[i]) != 0);
	}

	printf("Generating queries\n");
	uint64_t num_queries = 1000000;

	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	RAND_bytes((void*)query_set, num_queries * sizeof(uint64_t));
	int murmur_seed = rand();
	for (uint64_t i = 0; i < num_queries; i++) {
		query_set[i] = query_set[i] % (1ull << 24);
		query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
	}

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

	printf("false positive rate:  %f%%\n", (double)fp_count / num_queries * 100);

	free(insert_set);
	free(query_set);

	taf_destroy(filter);
}

void test_splinter_throughput(int qbits, uint64_t num_queries, uint64_t num_inc_queries) {
	printf("Testing %s...\n", __FUNCTION__);
	size_t num_slots = 1ull << qbits;
	double load = 0.9;

	//size_t xnum_slots = num_slots + (10 * sqrt(num_slots));

	TAF *filter = new_taf(num_slots);

	printf("Generating inserts\n");
	uint64_t num_inserts = num_slots * load;
	uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
	RAND_bytes((void*)insert_set, num_inserts * sizeof(uint64_t));

	printf("Performing inserts\n");
	FILE *inserts_fp = fopen("stats_splinter_inserts.csv", "w");
	fprintf(inserts_fp, "fill through\n");
	fclose(inserts_fp);

	uint64_t num_segments = 100;
	uint64_t measure_point = num_slots / num_segments - 1;
	uint64_t cur_segment = 0;

	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time, interval_time = start_time;

	for (uint64_t i = 0; i < num_inserts; i++) {
		taf_insert(filter, insert_set[i]);

		if (i >= measure_point) {
			gettimeofday(&timecheck, NULL);
			inserts_fp = fopen("stats_splinter_inserts.csv", "a");
			fprintf(inserts_fp, "%f %f\n", 100. * i / num_slots, (double)num_slots / num_segments * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time));
			fclose(inserts_fp);
			fprintf(stderr, "\rload factor: %d%%", (int)((i + 1) * 100 / num_slots));

			cur_segment++;
			measure_point = (cur_segment + 1) * num_slots / num_segments - 1;
			
			gettimeofday(&timecheck, NULL);
			interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		}
	}

	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("\n");
	printf("time for inserts:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughput:    %f ops/sec\n", (double)num_inserts * 1000000 / (end_time - start_time));
	printf("cpu time for inserts: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	printf("database inserts:     %lu\n", db_insert_count);
	printf("database updates:     %lu\n", db_update_count);
	printf("database queries:     %lu\n", db_query_count);

	printf("Generating queries\n");
	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	RAND_bytes((void*)query_set, num_queries * sizeof(uint64_t));
	//int murmur_seed = rand();
	for (uint64_t i = 0; i < num_queries; i++) {
		//query_set[i] %= (1ull << 24);
		//query_set[i] = query_set[i] % (1ull << 24);
		//query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
	}

	char dataset_buffer[256];
	FILE *shalla = fopen("../../../data/shalla.txt", "r");
	FILE *caida = fopen("../../../data/20140619-140100.csv", "r");
	fgets(dataset_buffer, sizeof(dataset_buffer), caida);
	for (int q = 0; q < num_queries; q++) {
		fgets(dataset_buffer, sizeof(dataset_buffer), shalla);
		query_set[q] = hash_str(dataset_buffer);
		//fgets(dataset_buffer, sizeof(dataset_buffer), caida);
		//csv_get_col(dataset_buffer, 3);
		//query_set[q] = hash_str(dataset_buffer);
	}
	fclose(shalla);
	fclose(caida);

	printf("Performing queries\n");
	FILE *queries_fp = fopen("stats_splinter_queries.csv", "w");
	fprintf(queries_fp, "queries through\n");
	fclose(queries_fp);
	FILE *fprates_fp = fopen("stats_splinter_fprates.csv", "w");
	fprintf(fprates_fp, "queries fprate\n");
	fclose(fprates_fp);

	measure_point = num_queries / num_segments - 1;
	cur_segment = 0;

	uint64_t fp_count = 0;
	uint64_t attempted_adapts = 0, successful_adapts = 0;

	start_clock = clock();
	gettimeofday(&timecheck, NULL);
	interval_time = start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

	for (uint64_t i = 0; i < num_queries; i++) {
		int ret = taf_lookup(filter, query_set[i]);
		if (ret == -1) {
			fp_count++;
			attempted_adapts++;
			ret = taf_lookup(filter, query_set[i]);
			if (ret != -1) {
				successful_adapts++;
			}
		}

		if (i >= measure_point) {
			gettimeofday(&timecheck, NULL);

			queries_fp = fopen("stats_splinter_queries.csv", "a");
			fprintf(queries_fp, "%lu %f\n", i, (double)num_queries / num_segments * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time));
			fclose(queries_fp);
			
			fprates_fp = fopen("stats_splinter_fprates.csv", "a");
			fprintf(fprates_fp, "%lu %f\n", i, 100. * fp_count / i);
			fclose(fprates_fp);

			fprintf(stderr, "\rqueries done: %d%%", (int)((i + 1) * 100 / num_queries));

			cur_segment++;
			measure_point = (cur_segment + 1) * num_queries / num_segments - 1;

			gettimeofday(&timecheck, NULL);
			interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		}
	}

	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("\n");
	printf("time for queries:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
	printf("cpu time for queries: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	printf("attempted adapts:     %lu\n", attempted_adapts);
	printf("successful adapts:    %lu\n", successful_adapts);

	printf("false positive rate:  %f%%\n", (double)fp_count / num_queries * 100);

	free(insert_set);
	free(query_set);

	taf_destroy(filter);
}

void test_splinter_adversarial(int qbits, uint64_t num_queries, char **argv, int argc) {
	printf("Testing %s...\n", __FUNCTION__);
	size_t num_slots = 1ull << qbits;
	double load = 0.9;

	//size_t xnum_slots = num_slots + (10 * sqrt(num_slots));

	TAF *filter = new_taf(num_slots);

	printf("Generating inserts\n");
	uint64_t num_inserts = num_slots * load;
	uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
	RAND_bytes((void*)insert_set, num_inserts * sizeof(uint64_t));

	printf("Performing inserts\n");
	FILE *inserts_fp = fopen("stats_splinter_adv_inserts.csv", "w");
	fprintf(inserts_fp, "fill through\n");
	fclose(inserts_fp);

	uint64_t num_segments = 100;
	uint64_t measure_point = num_slots / num_segments - 1;
	uint64_t cur_segment = 0;

	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time, interval_time = start_time;

	for (uint64_t i = 0; i < num_inserts; i++) {
		taf_insert(filter, insert_set[i]);

		if (i >= measure_point) {
			gettimeofday(&timecheck, NULL);
			inserts_fp = fopen("stats_splinter_adv_inserts.csv", "a");
			fprintf(inserts_fp, "%f %f\n", 100. * i / num_slots, (double)num_slots / num_segments * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time));
			fclose(inserts_fp);
			fprintf(stderr, "\rload factor: %d%%", (int)((i + 1) * 100 / num_slots));

			cur_segment++;
			measure_point = (cur_segment + 1) * num_slots / num_segments - 1;
			
			gettimeofday(&timecheck, NULL);
			interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		}
	}

	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("\n");
	printf("time for inserts:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughput:    %f ops/sec\n", (double)num_inserts * 1000000 / (end_time - start_time));
	printf("cpu time for inserts: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));

	printf("Performing queries\n");
	char filename[100];

	uint64_t fp_count = 0;

	uint64_t max_adv_set_len = 1000000ull;
	uint64_t *adv_set = malloc(max_adv_set_len * sizeof(uint64_t));
	uint64_t adv_set_len = 0;
	uint64_t adv_index = 0;
	uint64_t adv_attempted = 0, adv_successful = 0, adv_failed = 0, adv_found = 0;

	for (int trial = 0; trial < argc; trial++) {
		RAND_bytes((void*)query_set, num_queries * sizeof(uint64_t));

		uint64_t adv_freq = strtoull(argv[trial], NULL, 10);

		sprintf(filename, "%d-%lu-%lu-adversarial.csv", qbits, num_queries, adv_freq);
		FILE *adv_fp = fopen(filename, "w");
		fprintf(adv_fp, "queries through fprate\n");

		adv_set_len = 0;
		adv_index = 0;
		adv_attempted = adv_successful = adv_failed = adv_found = 0;
		fp_count = 0;

		measure_point = num_queries / num_segments - 1;
		cur_segment = 0;

		start_clock = clock();
		gettimeofday(&timecheck, NULL);
		interval_time = start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

		for (uint64_t i = 0; i < num_queries; i++) {
			if (i % adv_freq == 0 && adv_set_len > 0) {
				adv_attempted++;
				if (adv_index >= adv_set_len) adv_index = 0;

				int ret = taf_lookup(filter, adv_set[adv_index]);
				if (ret == -1) {
					fp_count++;
					adv_successful++;
				}
				else {
					adv_set[adv_index] = adv_set[--adv_set_len];
					adv_failed++;
				}

				adv_index++;
			}
			else {
				int ret = taf_lookup(filter, query_set[i]);
				if (ret == -1) {
					fp_count++;
					adv_set[adv_set_len++] = query_set[i];
					adv_found++;
				}
			}

			if (i >= measure_point) {
				gettimeofday(&timecheck, NULL);
				fprintf(adv_fp, "%lu %f %f\n", i, (double)num_queries / num_segments * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time), (double)fp_count / i);

				fprintf(stderr, "\rqueries done: %d%%", (int)((i + 1) * 100 / num_queries));

				cur_segment++;
				measure_point = (cur_segment + 1) * num_queries / num_segments - 1;

				gettimeofday(&timecheck, NULL);
				interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
			}
		}

		gettimeofday(&timecheck, NULL);
		end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		end_clock = clock();

		fclose(adv_fp);

		printf("\n");
		printf("time for queries:     %f sec\n", (double)(end_time - start_time) / 1000000);
		printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
		printf("cpu time for queries: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

		printf("adv found:            %lu\n", adv_found);
		printf("adv attempted:        %lu\n", adv_attempted);
		printf("adv successful:       %lu\n", adv_successful);
		printf("adv failed:           %lu\n", adv_failed);

		printf("false positive rate:  %f%%\n", (double)fp_count / num_queries * 100);
	}

	free(insert_set);
	free(query_set);
	free(adv_set);

	taf_destroy(filter);
}

int main(int argc, char **argv) {
	if (argc == 1) {
		printf("need more arguments\n");
		exit(0);
	}
	//test_splinter_ops();
	//test_splinter_veracity();
	test_splinter_throughput(atoi(argv[1]), strtoull(argv[2], NULL, 10), strtoull(argv[3], NULL, 10)); // ./splinter_taf [log of nslots] [num queries] [num inc queries]
	//test_splinter_adversarial(atoi(argv[1]), strtoull(argv[2], NULL, 10), &argv[3], argc - 3); // ./splinter_taf [log of nslots] [num queries] [adv freq]
}
#endif // TEST_SPLINTER_TAF
