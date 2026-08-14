#include "aqf.h"
#include "hashutil.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/******************************************************************
 * Code for managing the metadata bits and slots w/o interpreting *
 * the content of the slots.
 ******************************************************************/

#define METADATA_INC_MODE 1

#define MAX_VALUE(nbits) ((1ULL << (nbits)) - 1)
#define BITMASK(nbits)                                    \
  ((nbits) == 64 ? 0xffffffffffffffff : MAX_VALUE(nbits))
#define CLUSTER_SIZE (1ULL<<14)
#define METADATA_WORD(qf,field,slot_index)                              \
  (get_block((qf), (slot_index) /                                       \
             QF_SLOTS_PER_BLOCK)->field[((slot_index)  % QF_SLOTS_PER_BLOCK) / 64])

#define GET_NO_LOCK(flag) (flag & QF_NO_LOCK)
#define GET_KEY_HASH(flag) (flag & QF_KEY_IS_HASH)
#define GET_WAIT_FOR_LOCK(flag) (flag & QF_WAIT_FOR_LOCK)
#define QF_WAIT_FOR_LOCK (0x04)


static void modify_metadata(pc_t *metadata, int cnt)
{
	pc_add(metadata, cnt);
	return;
}

static inline bool qf_spin_lock(volatile int *lock, uint8_t flag)
{
	if (GET_WAIT_FOR_LOCK(flag) != QF_WAIT_FOR_LOCK) {
		return !__sync_lock_test_and_set(lock, 1);
	} else {
		while (__sync_lock_test_and_set(lock, 1))
			while (*lock);
		return true;
	}

	return false;
}

static inline void qf_spin_unlock(volatile int *lock)
{
	__sync_lock_release(lock);
	return;
}

static inline bool qf_lock(QF *qf, uint64_t hash_bucket_index, bool small, uint8_t runtime_lock)
{
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (small) {
		if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], runtime_lock))
			return false;
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1], runtime_lock)) {
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
	} else {
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1], runtime_lock))
				return false;
		}
		if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], runtime_lock)) {
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <= CLUSTER_SIZE)
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
		if (!qf_spin_lock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1], runtime_lock)) {
			qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <= CLUSTER_SIZE)
				qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
	}
	return true;
}

static inline void qf_unlock(QF *qf, uint64_t hash_bucket_index, bool small)
{
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (small) {
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		}
		qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
	} else {
		qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <= CLUSTER_SIZE)
			qf_spin_unlock(&qf->runtimedata->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
	}
}

static inline int popcnt(uint64_t val)
{
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

static inline int64_t bitscanreverse(uint64_t val)
{
	if (val == 0) {
		return -1;
	} else {
		asm("bsr %[val], %[val]"
				: [val] "+r" (val)
				:
				: "cc");
		return val;
	}
}

static inline int popcntv(const uint64_t val, int ignore)
{
	if (ignore % 64)
		return popcnt (val & ~BITMASK(ignore % 64));
	else
		return popcnt(val);
}

static inline int bitrank(uint64_t val, int pos) {
	val = val & ((2ULL << pos) - 1);
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
static const uint8_t kSelectInByte[2048] = {
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

static inline uint64_t _select64(uint64_t x, int k)
{
	if (k >= popcnt(x)) { return 64; }

	const uint64_t kOnesStep4  = 0x1111111111111111ULL;
	const uint64_t kOnesStep8  = 0x0101010101010101ULL;
	const uint64_t kMSBsStep8  = 0x80ULL * kOnesStep8;

	uint64_t s = x;
	s = s - ((s & 0xA * kOnesStep4) >> 1);
	s = (s & 0x3 * kOnesStep4) + ((s >> 2) & 0x3 * kOnesStep4);
	s = (s + (s >> 4)) & 0xF * kOnesStep8;
	uint64_t byteSums = s * kOnesStep8;

	uint64_t kStep8 = k * kOnesStep8;
	uint64_t geqKStep8 = (((kStep8 | kMSBsStep8) - byteSums) & kMSBsStep8);
	uint64_t place = popcnt(geqKStep8) * 8;
	uint64_t byteRank = k - (((byteSums << 8) >> place) & (uint64_t)(0xFF));
	return place + kSelectInByte[((x >> place) & 0xFF) | (byteRank << 8)];
}

// Returns the position of the rank'th 1.  (rank = 0 returns the 1st 1)
// Returns 64 if there are fewer than rank+1 1s.
static inline uint64_t bitselect(uint64_t val, int rank) {
#ifdef __SSE4_2_
	uint64_t i = 1ULL << rank;
	asm("pdep %[val], %[mask], %[val]" : [val] "+r" (val) : [mask] "r" (i));
	asm("tzcnt %[bit], %[index]" : [index] "=r" (i) : [bit] "g" (val) : "cc");
	return i;
#endif
	return _select64(val, rank);
}

static inline uint64_t bitselectv(const uint64_t val, int ignore, int rank)
{
	return bitselect(val & ~BITMASK(ignore % 64), rank);
}

// metadata checks
static inline int is_runend(const QF *qf, uint64_t index)
{
	return ~(METADATA_WORD(qf, extensions, index) >> ((index % QF_SLOTS_PER_BLOCK) % 64)) &
	(METADATA_WORD(qf, runends, index) >> ((index % QF_SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

static inline int is_occupied(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, occupieds, index) >> ((index % QF_SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

static inline int is_extension(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, extensions, index) >> ((index % QF_SLOTS_PER_BLOCK) % 64)) &
	~(METADATA_WORD(qf, runends, index) >> ((index % QF_SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

static inline int is_counter(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, extensions, index) >> ((index % QF_SLOTS_PER_BLOCK) % 64)) &
	(METADATA_WORD(qf, runends, index) >> ((index % QF_SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

static inline int is_runend_or_counter(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, runends, index) >> ((index % QF_SLOTS_PER_BLOCK) & 0b111111)) & 1ULL;
}

static inline int is_extension_or_counter(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, extensions, index) >> ((index % QF_SLOTS_PER_BLOCK) & 0b111111)) & 1ULL;
}

// slot accessors
#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32 || QF_BITS_PER_SLOT == 64

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	return get_block(qf, index / QF_SLOTS_PER_BLOCK)->slots[index % QF_SLOTS_PER_BLOCK];
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	get_block(qf, index / QF_SLOTS_PER_BLOCK)->slots[index % QF_SLOTS_PER_BLOCK] =
		value & BITMASK(qf->metadata->bits_per_slot);
}

#elif QF_BITS_PER_SLOT > 0

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	//printf("index%lu\n", index);
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 QF_SLOTS_PER_BLOCK)->slots[(index %
																																QF_SLOTS_PER_BLOCK)
																			 * QF_BITS_PER_SLOT / 8];
	return (uint64_t)(((*p) >> (((index % QF_SLOTS_PER_BLOCK) * QF_BITS_PER_SLOT) %
															8)) & BITMASK(QF_BITS_PER_SLOT));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 QF_SLOTS_PER_BLOCK)->slots[(index %
																																QF_SLOTS_PER_BLOCK)
																			 * QF_BITS_PER_SLOT / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(QF_BITS_PER_SLOT);
	uint64_t v = value;
	int shift = ((index % QF_SLOTS_PER_BLOCK) * QF_BITS_PER_SLOT) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#else

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	if (index > qf->metadata->xnslots) {
		return QF_NO_SPACE;
	}
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 QF_SLOTS_PER_BLOCK)->slots[(index %
																																QF_SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	return (uint64_t)(((*p) >> (((index % QF_SLOTS_PER_BLOCK) *
															 qf->metadata->bits_per_slot) % 8)) &
										BITMASK(qf->metadata->bits_per_slot));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 QF_SLOTS_PER_BLOCK)->slots[(index %
																																QF_SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(qf->metadata->bits_per_slot);
	uint64_t v = value;
	int shift = ((index % QF_SLOTS_PER_BLOCK) * qf->metadata->bits_per_slot) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#endif

// run lookup
static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index);

static inline uint64_t block_offset(const QF *qf, uint64_t blockidx)
{
	if (sizeof(qf->blocks[0].offset) > 1 ||
			get_block(qf, blockidx)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
		return get_block(qf, blockidx)->offset;

	return run_end(qf, QF_SLOTS_PER_BLOCK * blockidx - 1) - QF_SLOTS_PER_BLOCK *
		blockidx + 1;
}

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index)
{
	uint64_t bucket_block_index       = hash_bucket_index / QF_SLOTS_PER_BLOCK;
	uint64_t bucket_intrablock_offset = hash_bucket_index % QF_SLOTS_PER_BLOCK;
	uint64_t bucket_blocks_offset = block_offset(qf, bucket_block_index);

	uint64_t bucket_intrablock_rank   = bitrank(get_block(qf, bucket_block_index)->occupieds[0], bucket_intrablock_offset);

	if (bucket_intrablock_rank == 0) {
		if (bucket_blocks_offset <= bucket_intrablock_offset)
			return hash_bucket_index;
		else
			return QF_SLOTS_PER_BLOCK * bucket_block_index + bucket_blocks_offset - 1;
	}

	uint64_t runend_block_index  = bucket_block_index + bucket_blocks_offset / QF_SLOTS_PER_BLOCK;
	uint64_t runend_ignore_bits  = bucket_blocks_offset % QF_SLOTS_PER_BLOCK;
	uint64_t runend_rank         = bucket_intrablock_rank - 1;
	uint64_t runend_block_offset = bitselectv(get_block(qf, runend_block_index)->runends[0] & (~get_block(qf, runend_block_index)->extensions[0]), runend_ignore_bits, runend_rank);
	if (runend_block_offset == QF_SLOTS_PER_BLOCK) {
		if (bucket_blocks_offset == 0 && bucket_intrablock_rank == 0) {
			/* The block begins in empty space, and this bucket is in that region of
			 * empty space */
			return hash_bucket_index;
		} else {
			do {
				runend_rank        -= popcntv(get_block(qf, runend_block_index)->runends[0] & (~get_block(qf, runend_block_index)->extensions[0]), runend_ignore_bits);
				runend_block_index++;
				runend_ignore_bits  = 0;
				runend_block_offset = bitselectv(get_block(qf, runend_block_index)->runends[0] & (~get_block(qf, runend_block_index)->extensions[0]), runend_ignore_bits, runend_rank);
			} while (runend_block_offset == QF_SLOTS_PER_BLOCK);
		}
	}

	uint64_t runend_index = QF_SLOTS_PER_BLOCK * runend_block_index + runend_block_offset;
	while (is_extension_or_counter(qf, runend_index + 1)) runend_index++;
	if (runend_index < hash_bucket_index)
		return hash_bucket_index;
	else {
		return runend_index;
	}
}

// slot state checks
static inline int offset_lower_bound(const QF *qf, uint64_t slot_index)
{
	const qfblock *b = get_block(qf, slot_index / QF_SLOTS_PER_BLOCK);
	const uint64_t slot_offset = slot_index % QF_SLOTS_PER_BLOCK;
	const uint64_t boffset = b->offset;
	const uint64_t occupieds = b->occupieds[0] & BITMASK(slot_offset + 1);
	assert(QF_SLOTS_PER_BLOCK == 64);
	if (boffset <= slot_offset) {
		const uint64_t runends = (b->runends[0] & ~(b->extensions[0]) & BITMASK(slot_offset)) >> boffset;
		return popcnt(occupieds) - popcnt(runends);
	}
	return boffset - slot_offset + popcnt(occupieds);
}

static inline int is_empty(const QF *qf, uint64_t slot_index)
{
	return offset_lower_bound(qf, slot_index) == 0;
}

static inline int might_be_empty(const QF *qf, uint64_t slot_index)
{
	return !is_occupied(qf, slot_index)
		&& !is_runend_or_counter(qf, slot_index)
		&& !is_extension_or_counter(qf, slot_index);
}

static inline uint64_t find_first_empty_slot(QF *qf, uint64_t from)
{
	do {
		while (is_extension_or_counter(qf, from)) from++;
		int t = offset_lower_bound(qf, from);
		if (t < 0) {
			return -1;
		}
		assert(t>=0);
		if (t == 0)
			break;
		from = from + t;
	} while(1);
	return from;
}

// REMAINDER_WORD macro needed by shift_remainders
#define REMAINDER_WORD(qf, i) ((uint64_t *)&(get_block(qf, (i)/qf->metadata->bits_per_slot)->slots[8 * ((i) % qf->metadata->bits_per_slot)]))

// shift helpers
static inline uint64_t shift_into_b(const uint64_t a, const uint64_t b, const int bstart, const int bend, const int amount)
{
	const uint64_t a_component = bstart == 0 ? (a >> (64 - amount)) : 0;
	const uint64_t b_shifted_mask = BITMASK(bend - bstart) << bstart;
	const uint64_t b_shifted = ((b_shifted_mask & b) << amount) & b_shifted_mask;
	const uint64_t b_mask = ~b_shifted_mask;
	return a_component | b_shifted | (b & b_mask);
}

#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32 || QF_BITS_PER_SLOT == 64

static inline void shift_remainders(QF *qf, uint64_t start_index, uint64_t empty_index)
{
	uint64_t start_block  = start_index / QF_SLOTS_PER_BLOCK;
	uint64_t start_offset = start_index % QF_SLOTS_PER_BLOCK;
	uint64_t empty_block  = empty_index / QF_SLOTS_PER_BLOCK;
	uint64_t empty_offset = empty_index % QF_SLOTS_PER_BLOCK;

	assert (start_index <= empty_index && empty_index < qf->metadata->xnslots);

	while (start_block < empty_block) {
		memmove(&get_block(qf, empty_block)->slots[1], 
						&get_block(qf, empty_block)->slots[0],
						empty_offset * sizeof(qf->blocks[0].slots[0]));
		get_block(qf, empty_block)->slots[0] = get_block(qf,
																			empty_block-1)->slots[QF_SLOTS_PER_BLOCK-1];
		empty_block--;
		empty_offset = QF_SLOTS_PER_BLOCK-1;
	}

	memmove(&get_block(qf, empty_block)->slots[start_offset+1], 
					&get_block(qf, empty_block)->slots[start_offset],
					(empty_offset - start_offset) * sizeof(qf->blocks[0].slots[0]));
}

#else

#define REMAINDER_WORD(qf, i) ((uint64_t *)&(get_block(qf, (i)/qf->metadata->bits_per_slot)->slots[8 * ((i) % qf->metadata->bits_per_slot)]))

static inline void shift_remainders(QF *qf, const uint64_t start_index, const uint64_t empty_index)
{
	uint64_t last_word = (empty_index + 1) * qf->metadata->bits_per_slot / 64;
	const uint64_t first_word = start_index * qf->metadata->bits_per_slot / 64;
	int bend = ((empty_index + 1) * qf->metadata->bits_per_slot) % 64;
	const int bstart = (start_index * qf->metadata->bits_per_slot) % 64;

	while (last_word != first_word) {
		*REMAINDER_WORD(qf, last_word) = shift_into_b(*REMAINDER_WORD(qf, last_word-1), *REMAINDER_WORD(qf, last_word), 0, bend, qf->metadata->bits_per_slot);
		last_word--;
		bend = 64;
	}
	*REMAINDER_WORD(qf, last_word) = shift_into_b(0, *REMAINDER_WORD(qf, last_word), bstart, bend, qf->metadata->bits_per_slot);
}

static inline void shift_runends(QF *qf, int64_t first, uint64_t last, uint64_t distance)
{ // also shifts extensions
	assert(last < qf->metadata->xnslots && distance < 64);
	uint64_t first_word = first / 64;
	uint64_t bstart = first % 64;
	uint64_t last_word = (last + distance + 1) / 64;
	uint64_t bend = (last + distance + 1) % 64;

	if (last_word != first_word) {
		METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)), METADATA_WORD(qf, runends, 64*last_word), 0, bend, distance);
		METADATA_WORD(qf, extensions, 64*last_word) = shift_into_b(METADATA_WORD(qf, extensions, 64*(last_word-1)), METADATA_WORD(qf, extensions, 64*last_word), 0, bend, distance);
		bend = 64;
		last_word--;
		while (last_word != first_word) {
			METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)), METADATA_WORD(qf, runends, 64*last_word), 0, bend, distance);
			METADATA_WORD(qf, extensions, 64*last_word) = shift_into_b(METADATA_WORD(qf, extensions, 64*(last_word-1)), METADATA_WORD(qf, extensions, 64*last_word), 0, bend, distance);
			last_word--;
		}
	}
	METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(0, METADATA_WORD(qf, runends, 64*last_word), bstart, bend, distance);
	METADATA_WORD(qf, extensions, 64*last_word) = shift_into_b(0, METADATA_WORD(qf, extensions, 64*last_word), bstart, bend, distance);
}

// insert — external dependency for count > 1 (not used by AQF, which always passes count=1)
extern int insert_and_extend(QF *qf, uint64_t index, uint64_t key, uint64_t count,
                              uint64_t other_key, uint64_t *ret_hash,
                              uint64_t *ret_other_hash, uint8_t flags);

// target_index is the index of the item's target bucket (this is used to update offset)
// insert_index is the index to actually place the item
// value is the data (remainder/extension) to go in the slot
static inline int insert_one_slot(QF *qf, uint64_t target_index, uint64_t insert_index, uint64_t value) {
	uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);

	if (empty_slot_index >= qf->metadata->xnslots) {
		return QF_NO_SPACE;
	}
	shift_remainders(qf, insert_index, empty_slot_index);

	set_slot(qf, insert_index, value);

	shift_runends(qf, insert_index, empty_slot_index - 1, 1);

	uint64_t i;
	for (i = target_index / QF_SLOTS_PER_BLOCK + 1; i <= empty_slot_index / QF_SLOTS_PER_BLOCK; i++) {
		if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset))) {
			get_block(qf, i)->offset++;
		}
		assert(get_block(qf, i)->offset != 0);
	}

	return 1;
}

// insert
static inline int insert(QF *qf, qf_insert_result *result, uint64_t count, uint8_t runtime_lock)
{
	uint64_t hash_remainder = result->hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index = (result->hash & BITMASK(qf->metadata->quotient_bits + qf->metadata->bits_per_slot)) >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % QF_SLOTS_PER_BLOCK;

	if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
		if (!qf_lock(qf, hash_bucket_index, true, runtime_lock))
			return QF_COULDNT_LOCK;
	}

	if (might_be_empty(qf, hash_bucket_index) && run_end(qf, hash_bucket_index) == hash_bucket_index) {
		set_slot(qf, hash_bucket_index, hash_remainder);
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL << hash_bucket_block_offset;
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << hash_bucket_block_offset;

#if METADATA_INC_MODE == 1
		qf->metadata->ndistinct_elts++;
		qf->metadata->noccupied_slots++;
		qf->metadata->nelts++;
#elif METADATA_INC_MODE == 2
		modify_metadata(&qf->runtimedata->pc_ndistinct_elts, 1);
		modify_metadata(&qf->runtimedata->pc_noccupied_slots, 1);
		modify_metadata(&qf->runtimedata->pc_nelts, 1);
#endif

		if (count > 1) {
			insert_and_extend(qf, hash_bucket_index, result->hash, count - 1, result->hash, result->hash, result->hash, QF_KEY_IS_HASH | QF_NO_LOCK);
		}
	} else {
		int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1;

		if (!is_occupied(qf, hash_bucket_index)) { /* Empty bucket, but its slot is taken. */
			insert_one_slot(qf, hash_bucket_index, runstart_index, hash_remainder);

			METADATA_WORD(qf, runends, runstart_index) |= 1ULL << (runstart_index % 64);
			METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << hash_bucket_block_offset;

#if METADATA_INC_MODE == 1
			qf->metadata->ndistinct_elts++;
			qf->metadata->noccupied_slots++;
			qf->metadata->nelts++;
#elif METADATA_INC_MODE == 2
			modify_metadata(&qf->runtimedata->pc_ndistinct_elts, 1);
			modify_metadata(&qf->runtimedata->pc_noccupied_slots, 1);
			modify_metadata(&qf->runtimedata->pc_nelts, count);
#endif

			if (count > 1) insert_and_extend(qf, hash_bucket_index, result->hash, count - 1, result->hash, result->hash, result->hash, QF_KEY_IS_HASH | QF_NO_LOCK);
		} else {
			uint64_t current_index = runstart_index, temp_index, temp_remainder, old_runend_index;
			uint64_t current_remainder;

			uint64_t count_info;
			int count_slots;
			assert(!is_extension_or_counter(qf, current_index));
			int was_runend = 0;
			uint64_t insert_index;
			do {
				current_remainder = get_slot(qf, current_index);
				if (current_remainder >= hash_remainder) {
					if (current_remainder == hash_remainder) {
						result->minirun_existed = 1;
						result->minirun_rank = 0;
						temp_index = current_index;
						temp_remainder = current_remainder;
						while (temp_remainder == hash_remainder) {
							result->minirun_rank++;
							if (is_runend(qf, temp_index)) {
								old_runend_index = temp_index;
								was_runend = 1;
								while (is_extension_or_counter(qf, ++temp_index));
								break;
							}
							while (is_extension_or_counter(qf, ++temp_index));
							temp_remainder = get_slot(qf, temp_index);
						}
						insert_index = temp_index;
						break;
					}
					insert_index = current_index;
					break;
				}
				else if (is_runend(qf, current_index)) {
					was_runend = 1;
					old_runend_index = current_index;
					insert_index = current_index + 1;
					while (is_extension_or_counter(qf, insert_index)) insert_index++;
					break;
				}
				else {
					current_index++;
					while (is_extension_or_counter(qf, current_index)) current_index++;
				}
			} while (current_index < qf->metadata->xnslots);
			if (current_index >= qf->metadata->xnslots) {
				if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
					qf_unlock(qf, hash_bucket_index, true);
				}
				return QF_NO_SPACE;
			}

			insert_one_slot(qf, hash_bucket_index, insert_index, hash_remainder);
			if (was_runend) {
				METADATA_WORD(qf, runends, insert_index) |= (1ULL << (insert_index % QF_SLOTS_PER_BLOCK));
				METADATA_WORD(qf, runends, old_runend_index) ^= (1ULL << (old_runend_index % QF_SLOTS_PER_BLOCK));
			}

#if METADATA_INC_MODE == 1
			qf->metadata->ndistinct_elts++;
			qf->metadata->noccupied_slots++;
			qf->metadata->nelts += count;
#elif METADATA_INC_MODE == 2
			modify_metadata(&qf->runtimedata->pc_ndistinct_elts, 1);
			modify_metadata(&qf->runtimedata->pc_noccupied_slots, 1);
			modify_metadata(&qf->runtimedata->pc_nelts, count);
#endif

			if (count > 1) insert_and_extend(qf, insert_index, result->hash, count - 1, result->hash, result->hash, result->hash, QF_KEY_IS_HASH | QF_NO_LOCK);
		}
	}

	if (GET_NO_LOCK(runtime_lock) != QF_NO_LOCK) {
		qf_unlock(qf, hash_bucket_index, true);
	}

	return 0;
}

static uint64_t aqf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t value_bits,
                         enum qf_hashmode hash, uint32_t seed, void* buffer, uint64_t
                         buffer_len)
{
	uint64_t num_slots, xnslots, nblocks, qbits;
	uint64_t key_remainder_bits, bits_per_slot;
	uint64_t size;
	uint64_t total_num_bytes;

	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */
	qbits = num_slots = nslots;
	xnslots = nslots + 10*sqrt((double)nslots);
	nblocks = (xnslots + QF_SLOTS_PER_BLOCK - 1) / QF_SLOTS_PER_BLOCK;
	key_remainder_bits = key_bits;
	while (nslots > 1 && key_remainder_bits > 0) {
		key_remainder_bits--;
		nslots >>= 1;
	}
	assert(key_remainder_bits >= 2);

	bits_per_slot = key_remainder_bits + value_bits;
	assert (QF_BITS_PER_SLOT == 0 || QF_BITS_PER_SLOT == qf->metadata->bits_per_slot);
	assert(bits_per_slot > 1);
#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32 || QF_BITS_PER_SLOT == 64
	size = nblocks * sizeof(qfblock);
#else
	size = nblocks * (sizeof(qfblock) + QF_SLOTS_PER_BLOCK * bits_per_slot / 8);
#endif

	total_num_bytes = sizeof(qfmetadata) + size;
	if (buffer == NULL || total_num_bytes > buffer_len)
		return total_num_bytes;

	qf->metadata = (qfmetadata *)(buffer);
	qf->blocks = (qfblock *)(qf->metadata + 1);

	qf->metadata->magic_endian_number = MAGIC_NUMBER;
	qf->metadata->reserved = 0;
	qf->metadata->hash_mode = hash;
	qf->metadata->total_size_in_bytes = size;
	qf->metadata->seed = seed;
	qf->metadata->nslots = num_slots;
	qf->metadata->xnslots = xnslots;
	qf->metadata->key_bits = key_bits;
	qf->metadata->value_bits = value_bits;
	qf->metadata->key_remainder_bits = key_remainder_bits;
	qf->metadata->bits_per_slot = bits_per_slot;
	qf->metadata->quotient_bits = 0;
	while (qbits > 1) {
		qbits >>= 1;
		qf->metadata->quotient_bits++;
	}
	qf->metadata->quotient_remainder_bits = qf->metadata->quotient_bits + qf->metadata->bits_per_slot;

	qf->metadata->range = qf->metadata->nslots;
	qf->metadata->range <<= qf->metadata->key_remainder_bits;
	qf->metadata->nblocks = (qf->metadata->xnslots + QF_SLOTS_PER_BLOCK - 1) /
		QF_SLOTS_PER_BLOCK;
	qf->metadata->nelts = 0;
	qf->metadata->ndistinct_elts = 0;
	qf->metadata->noccupied_slots = 0;

	qf->runtimedata->num_locks = (qf->metadata->xnslots / NUM_SLOTS_TO_LOCK) + 2;

	pc_init(&qf->runtimedata->pc_nelts, (int64_t*)&qf->metadata->nelts, 8, 100);
	pc_init(&qf->runtimedata->pc_ndistinct_elts, (int64_t*)&qf->metadata->ndistinct_elts, 8, 100);
	pc_init(&qf->runtimedata->pc_noccupied_slots, (int64_t*)&qf->metadata->noccupied_slots, 8, 100);
	qf->runtimedata->auto_resize = 0;
	qf->runtimedata->container_resize = NULL;
	qf->runtimedata->metadata_lock = 0;
	qf->runtimedata->locks = (volatile int *)calloc(qf->runtimedata->num_locks,
																									sizeof(volatile int));
	if (qf->runtimedata->locks == NULL) {
		perror("Couldn't allocate memory for runtime locks.");
		exit(EXIT_FAILURE);
	}

	return total_num_bytes;
}

static void *aqf_destroy(QF *qf)
{
	assert(qf->runtimedata != NULL);
	if (qf->runtimedata->locks != NULL)
		free((void*)qf->runtimedata->locks);
	if (qf->runtimedata->wait_times != NULL)
		free(qf->runtimedata->wait_times);
	free(qf->runtimedata);

	return (void*)qf->metadata;
}

bool aqf_malloc(QF *qf, uint64_t nslots, uint64_t key_bits,
								uint64_t value_bits, enum qf_hashmode hash, uint32_t seed)
{
	uint64_t total_num_bytes = aqf_init(qf, nslots, key_bits, value_bits,
																			hash, seed, NULL, 0);

	void *buffer = calloc(1, total_num_bytes);
	if (buffer == NULL) {
		perror("Couldn't allocate memory for the AQF.");
		exit(EXIT_FAILURE);
	}

	qf->runtimedata = (qfruntime *)calloc(sizeof(qfruntime), 1);
	if (qf->runtimedata == NULL) {
		perror("Couldn't allocate memory for runtime data.");
		exit(EXIT_FAILURE);
	}

	uint64_t init_size = aqf_init(qf, nslots, key_bits, value_bits, hash, seed,
																buffer, total_num_bytes);

	if (init_size == total_num_bytes)
		return true;
	else
		return false;
}

void aqf_free(QF *qf)
{
	assert(qf->metadata != NULL);
	void *buffer = aqf_destroy(qf);
	if (buffer != NULL)
		free(buffer);
}

int aqf_insert(QF *qf, uint64_t key, uint64_t count,
               qf_insert_result *result, uint8_t flags)
{
	result->hash = key;
	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			result->hash = MurmurHash64A(((void *)&key), sizeof(key), qf->metadata->seed);
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			result->hash = hash_64(key, -1ULL);
	}

	result->minirun_existed = 0;
	result->minirun_rank = 0;
	result->minirun_id = result->hash & BITMASK(qf->metadata->quotient_bits + qf->metadata->bits_per_slot);
	return insert(qf, result, count, flags);
}

static inline int get_slot_info(const QF *qf, uint64_t index, uint64_t *ext, int *ext_slots, uint64_t *count, int *count_slots) {
	assert(!is_extension_or_counter(qf, index));

	int curr = ++index;
	if (ext != NULL) {
		*ext = 0;
		while (is_extension(qf, curr)) {
			*ext |= get_slot(qf, curr) << ((curr - index) * qf->metadata->bits_per_slot);
			curr++;
		}
		if (ext_slots != NULL) *ext_slots = curr - index;
	}
	else while (is_extension(qf, index++));

	index = curr;
	if (count != NULL) {
		*count = 0;
		if (!is_counter(qf, curr)) {
			*count = 1;
			if (count_slots != NULL) *count_slots = 0;
		}
		while (is_counter(qf, curr)) {
			*count |= get_slot(qf, curr) << ((curr - index) * qf->metadata->bits_per_slot);
			curr++;
		}
		if (count_slots != NULL) *count_slots = curr - index;
	}
	return 1;
}

static inline int query(const QF *qf, uint64_t key, uint64_t *ret_hash, uint8_t *ret_minirun_rank, uint8_t flags) {
	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			*ret_hash = MurmurHash64A(((void *)&key), sizeof(key), qf->metadata->seed);
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE)
			*ret_hash = hash_64(key, -1ULL);
		else
			*ret_hash = key;
	}
	else {
		*ret_hash = key;
	}
	uint64_t hash_remainder   = *ret_hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index = (*ret_hash >> qf->metadata->bits_per_slot) & BITMASK(qf->metadata->quotient_bits);

	if (!is_occupied(qf, hash_bucket_index))
		return 0;

	uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	uint64_t current_index = runstart_index;
	*ret_minirun_rank = 0;
	do {
		if (get_slot(qf, current_index) == hash_remainder) {
			uint64_t ext, count;
			int ext_len, count_len;
			get_slot_info(qf, current_index, &ext, &ext_len, &count, &count_len);
			if (((*ret_hash >> (qf->metadata->quotient_bits + qf->metadata->bits_per_slot)) & BITMASK(qf->metadata->bits_per_slot * ext_len)) == ext) {
				return count;
			}
			(*ret_minirun_rank)++;
			if (is_runend(qf, current_index++)) break;
			current_index += ext_len + count_len;
		}
		else {
			if (is_runend(qf, current_index++)) break;
			while (is_extension_or_counter(qf, current_index)) current_index++;
		}
	} while (current_index < qf->metadata->xnslots);

	return 0;
}

int aqf_query(const QF *qf, uint64_t key, uint64_t *ret_hash,
              uint8_t *ret_minirun_rank, uint8_t flags)
{
	return query(qf, key, ret_hash, ret_minirun_rank, flags);
}

static inline int adapt(QF *qf, uint64_t index, uint64_t hash_bucket_index, uint64_t hash, uint64_t other_hash, uint64_t *ret_hash) {
	uint64_t ext, count;
	int ext_len, count_len;
	if (!get_slot_info(qf, index, &ext, &ext_len, &count, &count_len)) return 0;
	assert((hash & BITMASK(qf->metadata->quotient_bits + qf->metadata->bits_per_slot)) == (get_slot(qf, index) | (hash_bucket_index << qf->metadata->bits_per_slot)));
	int ext_bits = qf->metadata->bits_per_slot * ext_len;
	*ret_hash = hash & BITMASK(ext_bits + qf->metadata->quotient_bits + qf->metadata->bits_per_slot);
	int slots_used = ext_len + 1;
	hash >>= ext_bits + qf->metadata->quotient_bits;
	other_hash >>= ext_bits + qf->metadata->quotient_bits;

	do {
		hash >>= qf->metadata->bits_per_slot;
		other_hash >>= qf->metadata->bits_per_slot;

		uint64_t empty_slot_index = find_first_empty_slot(qf, index + slots_used);
		if (empty_slot_index >= qf->metadata->xnslots) {
			return QF_NO_SPACE;
		}

		shift_remainders(qf, index + slots_used, empty_slot_index);

		set_slot(qf, index + slots_used, hash & BITMASK(qf->metadata->bits_per_slot));
		*ret_hash |= (hash & BITMASK(qf->metadata->bits_per_slot)) << (ext_bits + qf->metadata->quotient_bits + qf->metadata->bits_per_slot);

		shift_runends(qf, index + slots_used, empty_slot_index - 1, 1);

		uint64_t i;
		for (i = hash_bucket_index / QF_SLOTS_PER_BLOCK + 1; i <= empty_slot_index / QF_SLOTS_PER_BLOCK; i++) {
			if (get_block(qf, i)->offset < BITMASK(8 * sizeof(qf->blocks[0].offset))) get_block(qf, i)->offset++;
		}

		METADATA_WORD(qf, extensions, index + slots_used) |= 1ULL << ((index + slots_used) % 64);
		qf->metadata->noccupied_slots++;
		slots_used++;
		ext_bits += qf->metadata->bits_per_slot;
	} while (((hash & BITMASK(qf->metadata->bits_per_slot)) == (other_hash & BITMASK(qf->metadata->bits_per_slot))) && (ext_bits + qf->metadata->quotient_bits + qf->metadata->bits_per_slot < 64));

	return ext_bits + qf->metadata->quotient_bits + qf->metadata->bits_per_slot;
}

int aqf_adapt(QF *qf, uint64_t orig_key, uint64_t fp_key,
              uint64_t minirun_rank, uint8_t flags)
{
	uint64_t hash;
	if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT) {
			hash = MurmurHash64A(((void *)&fp_key), sizeof(fp_key), qf->metadata->seed);
		}
		else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE) {
			hash = hash_64(fp_key, -1ULL);
		}
	}
	else hash = fp_key;

	uint64_t hash_remainder = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index = (hash >> qf->metadata->bits_per_slot) & BITMASK(qf->metadata->quotient_bits);
	if (!is_occupied(qf, hash_bucket_index)) return 0;

	uint64_t current_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1;
	int curr_minirun_rank = 0;

	uint64_t count_info, hash_info;
	int count_slots, hash_slots;
	do {
		if (get_slot(qf, current_index) < hash_remainder) {
			while (is_extension_or_counter(qf, ++current_index));
		}
		else if (get_slot(qf, current_index) == hash_remainder) {
			get_slot_info(qf, current_index, &hash_info, &hash_slots, &count_info, &count_slots);

			if (curr_minirun_rank == minirun_rank) {
				uint64_t key_hash = orig_key, fp_hash = fp_key;
				if (GET_KEY_HASH(flags) != QF_KEY_IS_HASH) {
					if (qf->metadata->hash_mode == QF_HASH_DEFAULT) {
						key_hash = MurmurHash64A(((void *)&orig_key), sizeof(orig_key), qf->metadata->seed);
						fp_hash = MurmurHash64A(((void *)&fp_key), sizeof(fp_key), qf->metadata->seed);
					}
					else if (qf->metadata->hash_mode == QF_HASH_INVERTIBLE) {
						key_hash = hash_64(orig_key, -1ULL);
						fp_hash = hash_64(fp_key, -1ULL);
					}
				}
				if (key_hash == fp_hash) return 0;
				adapt(qf, current_index, (key_hash % qf->metadata->range) >> qf->metadata->bits_per_slot, key_hash, fp_hash, &hash_info);
				return 1;
			}

			if (is_runend(qf, current_index)) break;
			current_index += count_slots + hash_slots + 1;
			curr_minirun_rank++;
		}
		else break;
	} while (current_index < qf->metadata->xnslots);

	return 0;
}

uint64_t aqf_get_noccupied_slots(const QF *qf)
{
	return qf->metadata->noccupied_slots;
}

uint64_t aqf_get_nslots(const QF *qf)
{
	return qf->metadata->nslots;
}

uint64_t aqf_get_total_size_in_bytes(const QF *qf)
{
	return qf->metadata->total_size_in_bytes;
}

int aqf_bulk_insert(QF *qf, uint64_t *keys, uint64_t numKeys, uint8_t flags)
{
	uint64_t *fingerprints = malloc(numKeys * sizeof(uint64_t));
	if (!fingerprints) return -1;

	for (uint64_t i = 0; i < numKeys; i++) {
		if (qf->metadata->hash_mode == QF_HASH_DEFAULT)
			fingerprints[i] = MurmurHash64A(&keys[i], sizeof(keys[i]),
																			 qf->metadata->seed);
		else
			fingerprints[i] = keys[i];
	}

	int ret = 0;
	for (uint64_t i = 0; i < numKeys; i++) {
		qf_insert_result result;
		ret = aqf_insert(qf, fingerprints[i], 1, &result,
								 flags | QF_KEY_IS_HASH);
		if (ret < 0) break;
	}

	free(fingerprints);
	return ret;
}

static inline void aqf_dump_block(const QF *qf, uint64_t i)
{
	uint64_t j;

	printf("%-192d", get_block(qf, i)->offset);
	printf("\n");

	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf("%02lx ", j);
	printf("\n");

	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->occupieds[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->runends[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");

#if QF_BITS_PER_SLOT == 8 || QF_BITS_PER_SLOT == 16 || QF_BITS_PER_SLOT == 32
	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#elif QF_BITS_PER_SLOT == 64
	for (j = 0; j < QF_SLOTS_PER_BLOCK; j++)
		printf("%02lx ", get_block(qf, i)->slots[j]);
#else
	for (j = 0; j < QF_SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#endif

	printf("\n");

	printf("\n");
}

void aqf_dump_metadata(const QF *qf) {
	printf("Slots: %lu Occupied: %lu Elements: %lu Distinct: %lu\n",
				 qf->metadata->nslots,
				 qf->metadata->noccupied_slots,
				 qf->metadata->nelts,
				 qf->metadata->ndistinct_elts);
	printf("Key_bits: %lu Value_bits: %lu Remainder_bits: %lu Bits_per_slot: %lu\n",
				 qf->metadata->key_bits,
				 qf->metadata->value_bits,
				 qf->metadata->key_remainder_bits,
				 qf->metadata->bits_per_slot);
}

void aqf_dump(const QF *qf)
{
	uint64_t i;

	printf("%lu %lu %lu\n",
				 qf->metadata->nblocks,
				 qf->metadata->ndistinct_elts,
				 qf->metadata->nelts);

	for (i = 0; i < qf->metadata->nblocks; i++) {
		aqf_dump_block(qf, i);
	}

}

int aqf_read_quotient_slots(const QF *qf, uint64_t start_quotient,
                             int num_quotients, qf_slot_debug *slots,
                             uint64_t max_slots, uint64_t *num_slots)
{
	uint64_t total = 0;
	uint64_t end_quotient = start_quotient + num_quotients;
	if (end_quotient > qf->metadata->nslots) end_quotient = qf->metadata->nslots;

	for (uint64_t q = start_quotient; q < end_quotient; q++) {
		if (!is_occupied(qf, q)) continue;

		uint64_t run_start = (q == 0) ? 0 : run_end(qf, q - 1) + 1;
		uint64_t run_end_idx = run_end(qf, q);

		uint64_t i = run_start;
		while (i <= run_end_idx) {
			if (total < max_slots) {
				slots[total].slot = get_slot(qf, i);
				slots[total].is_runend = is_runend(qf, i);
				slots[total].is_extension = is_extension(qf, i);
				slots[total].quotient = q;
			}
			total++;
			i++;
		}
	}

	*num_slots = total;
	return 0;
}
