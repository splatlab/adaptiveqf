#ifndef AQF_INT_H
#define AQF_INT_H

#include <inttypes.h>
#include <stdbool.h>

#include "partitioned_counter.h"

#ifdef __cplusplus
extern "C" {
#endif

struct quotient_filter;

/* Hash modes (from gqf.h) */
enum qf_hashmode {
  QF_HASH_DEFAULT,
  QF_HASH_INVERTIBLE,
  QF_HASH_NONE
};

/* Lock and flag defines (from gqf.h) */
#define QF_NO_LOCK (0x01)
#define QF_KEY_IS_HASH (0x08)

/* Return codes (from gqf.h) */
#define QF_NO_SPACE (-1)
#define QF_COULDNT_LOCK (-2)
#define QF_DOESNT_EXIST (-3)

/* Insert result (from gqf.h) */
struct qf_insert_result_t {
  uint64_t hash;
  uint64_t minirun_id;   // base fingerprint (quotient + remainder), not position in block
  int minirun_existed;
  int minirun_rank;      // position within the minirun of matching remainders
} typedef qf_insert_result;

/* Block layout (from gqf_int.h) */
#define MAGIC_NUMBER 1018874902021329732
#define NUM_SLOTS_TO_LOCK (1ULL<<16)
#define MAX_VALUE(nbits) ((1ULL << (nbits)) - 1)
#define BITMASK(nbits)                                    \
  ((nbits) == 64 ? 0xffffffffffffffff : MAX_VALUE(nbits))
#define QF_BITS_PER_SLOT 0
#define QF_BLOCK_OFFSET_BITS (6)
#define QF_SLOTS_PER_BLOCK (1ULL << QF_BLOCK_OFFSET_BITS)
#define QF_METADATA_WORDS_PER_BLOCK ((QF_SLOTS_PER_BLOCK + 63) / 64)

typedef struct __attribute__ ((__packed__)) qfblock {
  uint8_t offset;
  uint64_t occupieds[QF_METADATA_WORDS_PER_BLOCK];
  uint64_t runends[QF_METADATA_WORDS_PER_BLOCK];
  uint64_t extensions[QF_METADATA_WORDS_PER_BLOCK];
#if QF_BITS_PER_SLOT == 8
  uint8_t  slots[QF_SLOTS_PER_BLOCK];
#elif QF_BITS_PER_SLOT == 16
  uint16_t  slots[QF_SLOTS_PER_BLOCK];
#elif QF_BITS_PER_SLOT == 32
  uint32_t  slots[QF_SLOTS_PER_BLOCK];
#elif QF_BITS_PER_SLOT == 64
  uint64_t  slots[QF_SLOTS_PER_BLOCK];
#elif QF_BITS_PER_SLOT != 0
  uint8_t   slots[QF_SLOTS_PER_BLOCK * QF_BITS_PER_SLOT / 8];
#else
  uint8_t   slots[];
#endif
} qfblock;

typedef struct qfblock qfblock;

typedef struct file_info {
  int fd;
  char *filepath;
} file_info;

typedef struct {
  uint64_t total_time_single;
  uint64_t total_time_spinning;
  uint64_t locks_taken;
  uint64_t locks_acquired_single_attempt;
} wait_time_data;

typedef struct quotient_filter_runtime_data {
  file_info f_info;
  uint32_t auto_resize;
  int64_t (*container_resize)(struct quotient_filter *qf, uint64_t nslots);
  pc_t pc_nelts;
  pc_t pc_ndistinct_elts;
  pc_t pc_noccupied_slots;
  uint64_t num_locks;
  volatile int metadata_lock;
  volatile int *locks;
  wait_time_data *wait_times;
} quotient_filter_runtime_data;

typedef quotient_filter_runtime_data qfruntime;

typedef struct quotient_filter_metadata {
  uint64_t magic_endian_number;
  enum qf_hashmode hash_mode;
  uint32_t reserved;
  uint64_t total_size_in_bytes;
  uint32_t seed;
  uint64_t nslots;
  uint64_t xnslots;
  uint64_t key_bits;
  uint64_t value_bits;
  uint64_t key_remainder_bits;
  uint64_t bits_per_slot;
  uint64_t quotient_bits;
  uint64_t quotient_remainder_bits;
  __uint128_t range;
  uint64_t nblocks;
  uint64_t nelts;
  uint64_t ndistinct_elts;
  uint64_t noccupied_slots;
} quotient_filter_metadata;

typedef quotient_filter_metadata qfmetadata;

typedef struct quotient_filter {
  qfruntime *runtimedata;
  qfmetadata *metadata;
  qfblock *blocks;
} quotient_filter;

typedef quotient_filter QF;

#if QF_BITS_PER_SLOT > 0
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
  return &qf->blocks[block_index];
}
#else
static inline qfblock * get_block(const QF *qf, uint64_t block_index)
{
  return (qfblock *)(((char *)qf->blocks)
                     + block_index * (sizeof(qfblock) + QF_SLOTS_PER_BLOCK *
                                      qf->metadata->bits_per_slot / 8));
}
#endif

#ifdef __cplusplus
}

/* External declarations — implemented in gqf.c, linked at build time */
extern "C" {
#endif

typedef struct {
  uint64_t slot;
  int is_runend;
  int is_extension;
  uint64_t quotient;
} qf_slot_debug;

int aqf_read_quotient_slots(const QF *qf, uint64_t start_quotient,
                             int num_quotients, qf_slot_debug *slots,
                             uint64_t max_slots, uint64_t *num_slots);
int qf_insert_using_ll_table(QF *qf, uint64_t key, uint64_t count,
                             qf_insert_result *result, uint8_t flags);
int qf_get_count_using_ll_table(const QF *qf, uint64_t key,
                                uint64_t *ret_hash, uint8_t *ret_minirun_rank,
                                uint8_t flags);
int qf_adapt_using_ll_table(QF *qf, uint64_t orig_key, uint64_t fp_key,
                            uint64_t minirun_rank, uint8_t flags);

#ifdef __cplusplus
}
#endif

#endif /* AQF_INT_H */