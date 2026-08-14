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
    uint8_t flags = QF_NO_LOCK;
    if (config.hash_mode == QF_HASH_NONE) flags |= QF_KEY_IS_HASH;

    uint64_t *fingerprints = new uint64_t[numKeys];
    uint64_t *ranks = new uint64_t[numKeys];
    for (uint64_t i = 0; i < numKeys; i++) {
      qf_insert_result result;
      int ret = aqf_insert(&qf, keys[i], 1, &result, flags);
      if (ret < 0) { delete[] fingerprints; delete[] ranks; return ret; }
      fingerprints[i] = result.minirun_id;
      ranks[i] = result.minirun_rank;
    }
    reverseMap.bulkInsert(keys, fingerprints, ranks, numKeys);
    delete[] fingerprints;
    delete[] ranks;
    return 0;
  }

  int queryFilter(uint64_t queryKey, QFilterQueryResult *result) {
    uint8_t flags = (config.hash_mode == QF_HASH_NONE) ? QF_KEY_IS_HASH : 0;
    uint8_t minirun_rank;
    uint64_t hash;
    uint64_t minirun_count = aqf_query(&qf, queryKey, &hash, &minirun_rank,
                                        flags);
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
    uint8_t flags = (config.hash_mode == QF_HASH_NONE) ? QF_KEY_IS_HASH : 0;
    uint64_t origKey;
    uint64_t fingerprint = filterResult->hash & BITMASK(config.qbits + config.rbits);
    int ret = reverseMap.getFingerprint(fingerprint, filterResult->minirun_rank,
                                         &origKey);
    if (ret) return -1;
    return aqf_adapt(&qf, origKey, queryKey, filterResult->minirun_rank,
                     flags);
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

  void dumpMetadata() { aqf_dump_metadata(&qf); }
  void dump() { aqf_dump(&qf); }

  void dumpBlock(uint64_t blockIndex) {
    uint64_t start_q = blockIndex * QF_SLOTS_PER_BLOCK;
    qf_slot_debug buf[128];
    uint64_t n;
    int ret = aqf_read_quotient_slots(&qf, start_q, QF_SLOTS_PER_BLOCK,
                                       buf, 128, &n);
    if (ret < 0) { return; }

    uint64_t bits_per_slot = qf.metadata->bits_per_slot;
    uint64_t qbits = qf.metadata->quotient_bits;

    fprintf(stdout, "Block %lu (offset=%u):\n", blockIndex,
            get_block(&qf, blockIndex)->offset);

    for (uint64_t i = 0; i < n; i++) {
      if (buf[i].is_extension) continue;

      uint64_t quotient = buf[i].quotient;
      uint64_t remainder = buf[i].slot;

      int rank = 0;
      for (uint64_t k = 0; k < i; k++) {
        if (!buf[k].is_extension && buf[k].quotient == quotient &&
            buf[k].slot == remainder) {
          rank++;
        }
      }

      uint64_t minirun_id = (quotient << bits_per_slot) | remainder;
      fprintf(stdout, "  q=%lu rem=%lu rank=%d", quotient, remainder, rank);

      fprintf(stdout, " exts:");
      for (uint64_t k = i + 1; k < n && buf[k].is_extension; k++) {
        fprintf(stdout, " %lu", buf[k].slot);
      }

      uint64_t key;
      if (reverseMap.getFingerprint(minirun_id, rank, &key) == 0) {
        fprintf(stdout, " key=%lx", key);
      }
      fprintf(stdout, "\n");
    }
  }

  QF *getQF() { return &qf; }
  ReverseMap *getReverseMap() { return &reverseMap; }

protected:
  QF qf;
  ReverseMap reverseMap;
  QFilterConfig config;
  size_t fullPoint;
};

#endif