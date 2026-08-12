#ifndef ADAPTIVEQF_H
#define ADAPTIVEQF_H

#include "qf_filter.hpp"
#include <cstddef>
#include <cstdint>

extern "C" {
#include <gqf.h>
#include <gqf_int.h>
#include <hashutil.h>
}

template <typename ReverseMap>
class AdaptiveQF {
public:
  int init(QFilterConfig cfg) {
    config = cfg;
    size_t num_slots = 1ull << config.qbits;
    if (!qf_malloc(&qf, num_slots, config.qbits + config.rbits, 0,
                   (enum qf_hashmode)config.hash_mode, 0)) {
      return -1;
    }
    reverseMap.init(config.reverse_map, config.qbits + config.rbits);
    fullPoint = config.max_load_factor * (1ULL << config.qbits);
    return 0;
  }

  int bulkInsert(uint64_t *keys, uint64_t numKeys) {
    uint64_t *fingerprints = new uint64_t[numKeys];
    for (uint64_t i = 0; i < numKeys; i++) {
      if (config.hash_mode == QF_HASH_DEFAULT) {
        fingerprints[i] = MurmurHash64A(&keys[i], sizeof(keys[i]),
                                         qf.metadata->seed);
      } else {
        fingerprints[i] = keys[i];  // INVERTIBLE or NONE
      }
    }
    for (uint64_t i = 0; i < numKeys; i++) {
      qf_insert_result result;
      int ret = qf_insert_using_ll_table(&qf, fingerprints[i], 1, &result,
                                          QF_NO_LOCK | QF_KEY_IS_HASH);
      if (ret < 0) {
        delete[] fingerprints;
        return ret;
      }
    }
    reverseMap.bulkInsert(keys, fingerprints, numKeys);
    delete[] fingerprints;
    return 0;
  }

  int queryFilter(uint64_t queryKey, QFilterQueryResult *result) {
    uint8_t minirun_rank;
    uint64_t hash;
    uint64_t minirun_count = qf_get_count_using_ll_table(
        &qf, queryKey, &hash, &minirun_rank, QF_KEY_IS_HASH);
    if (minirun_count > 0) {
      result->key_present = 1;
      result->minirun_count = minirun_count;
    } else {
      result->key_present = 0;
    }
    result->hash = hash;
    result->minirun_rank = minirun_rank;
    result->hash_index = 0;
    return 0;
  }

  int adapt(uint64_t queryKey, QFilterQueryResult *filterResult) {
    if (qf.metadata->noccupied_slots >= fullPoint) {
      return -1;
    }
    uint64_t origKey;
    uint64_t fingerprint = filterResult->hash;
    int ret = reverseMap.getFingerprint(fingerprint, filterResult->minirun_rank,
                                         &origKey);
    if (ret) return -1;
    ret = qf_adapt_using_ll_table(&qf, origKey, queryKey,
                                   filterResult->minirun_rank, QF_KEY_IS_HASH);
    return 1;
  }

  double loadFactor() {
    return (double)qf.metadata->noccupied_slots / qf.metadata->nslots;
  }

  int close() {
    reverseMap.close();
    return 0;
  }

  uint64_t sizeInBytes() {
    return qf.metadata->total_size_in_bytes;
  }

protected:
  QF qf;
  ReverseMap reverseMap;
  QFilterConfig config;
  size_t fullPoint;
};

#endif
