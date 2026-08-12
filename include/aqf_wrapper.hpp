#ifndef AQF_WRAPPER_HPP
#define AQF_WRAPPER_HPP

#include <cstddef>
#include <cstdint>

#include "aqf.h"
#include "reverse_map_config.hpp"

struct QFilterConfig {
  size_t qbits;
  size_t rbits;
  double max_load_factor;
  int hash_mode;
  ReverseMapConfig reverse_map;
};

struct QFilterQueryResult {
  int key_present;
  int minirun_rank;
  int minirun_count;
  uint64_t hash;
  uint64_t hash_index;
};

template <typename ReverseMap>
class AdaptiveQF {
public:
  int init(QFilterConfig cfg) {
    config = cfg;
    size_t num_slots = 1ull << config.qbits;
    if (!aqf_malloc(&qf, num_slots, config.qbits + config.rbits, 0,
                    (enum qf_hashmode)config.hash_mode, 0)) {
      return -1;
    }
    reverseMap.init(config.reverse_map, config.qbits + config.rbits);
    fullPoint = config.max_load_factor * (1ULL << config.qbits);
    return 0;
  }

  int bulkInsert(uint64_t *keys, uint64_t numKeys) {
    return aqf_bulk_insert(&qf, keys, numKeys, QF_NO_LOCK);
  }

  int queryFilter(uint64_t queryKey, QFilterQueryResult *result) {
    uint8_t minirun_rank;
    uint64_t hash;
    uint64_t minirun_count = aqf_query(&qf, queryKey, &hash, &minirun_rank,
                                        QF_KEY_IS_HASH);
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
    if (aqf_get_noccupied_slots(&qf) >= fullPoint) {
      return -1;
    }
    uint64_t origKey;
    uint64_t fingerprint = filterResult->hash;
    int ret = reverseMap.getFingerprint(fingerprint, filterResult->minirun_rank,
                                         &origKey);
    if (ret) return -1;
    return aqf_adapt(&qf, origKey, queryKey, filterResult->minirun_rank,
                     QF_KEY_IS_HASH);
  }

  double loadFactor() {
    return (double)aqf_get_noccupied_slots(&qf) / aqf_get_nslots(&qf);
  }

  int close() {
    aqf_free(&qf);
    reverseMap.close();
    return 0;
  }

  uint64_t sizeInBytes() {
    return aqf_get_total_size_in_bytes(&qf);
  }

protected:
  QF qf;
  ReverseMap reverseMap;
  QFilterConfig config;
  size_t fullPoint;
};

#endif