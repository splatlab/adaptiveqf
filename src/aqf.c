#include "aqf.h"
#include "hashutil.h"

bool aqf_malloc(QF *qf, uint64_t nslots, uint64_t key_bits,
                uint64_t value_bits, enum qf_hashmode hash, uint32_t seed)
{
  return qf_malloc(qf, nslots, key_bits, value_bits, hash, seed);
}

void aqf_free(QF *qf)
{
  qf_free(qf);
}

int aqf_insert(QF *qf, uint64_t key, uint64_t count,
               qf_insert_result *result, uint8_t flags)
{
  return qf_insert_using_ll_table(qf, key, count, result, flags);
}

int aqf_query(const QF *qf, uint64_t key, uint64_t *ret_hash,
              uint8_t *ret_minirun_rank, uint8_t flags)
{
  return qf_get_count_using_ll_table(qf, key, ret_hash, ret_minirun_rank, flags);
}

int aqf_adapt(QF *qf, uint64_t orig_key, uint64_t fp_key,
              uint64_t minirun_rank, uint8_t flags)
{
  return qf_adapt_using_ll_table(qf, orig_key, fp_key, minirun_rank, flags);
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
      fingerprints[i] = MurmurHash64A(&keys[i], sizeof(keys[i]), qf->metadata->seed);
    else
      fingerprints[i] = keys[i];
  }

  int ret = 0;
  for (uint64_t i = 0; i < numKeys; i++) {
    qf_insert_result result;
    ret = qf_insert_using_ll_table(qf, fingerprints[i], 1, &result,
                                   flags | QF_KEY_IS_HASH);
    if (ret < 0) break;
  }

  free(fingerprints);
  return ret;
}