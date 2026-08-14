#ifndef AQF_H
#define AQF_H

#include <inttypes.h>
#include <stdbool.h>

#include "aqf_int.h"
#include "hashutil.h"

// [TODO] File-backed QF operations (serialize/deserialize, init from file)
// were removed. These need to be re-added to the AQF API.
// See original gqf_file.h for the interface.

#ifdef __cplusplus
extern "C" {
#endif

// Lifecycle
bool aqf_malloc(QF *qf, uint64_t nslots, uint64_t key_bits,
                uint64_t value_bits, enum qf_hashmode hash, uint32_t seed);
void aqf_free(QF *qf);

// Core operations
int aqf_insert(QF *qf, uint64_t key, uint64_t count,
               qf_insert_result *result, uint8_t flags);
int aqf_query(const QF *qf, uint64_t key, uint64_t *ret_hash,
              uint8_t *ret_minirun_rank, uint8_t flags);
int aqf_adapt(QF *qf, uint64_t orig_key, uint64_t fp_key,
              uint64_t minirun_rank, uint8_t flags);

// Metadata accessors
uint64_t aqf_get_noccupied_slots(const QF *qf);
uint64_t aqf_get_nslots(const QF *qf);
uint64_t aqf_get_total_size_in_bytes(const QF *qf);

// Bulk insert — handles hashing internally based on hash_mode
int aqf_bulk_insert(QF *qf, uint64_t *keys, uint64_t numKeys, uint8_t flags);

// Debug
void aqf_dump_metadata(const QF *qf);
void aqf_dump(const QF *qf);

#ifdef __cplusplus
}
#endif

#endif