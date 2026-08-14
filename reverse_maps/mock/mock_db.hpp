#ifndef MOCK_DB_H
#define MOCK_DB_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include "reverse_map_config.hpp"

class MockDB {
public:
  int init(ReverseMapConfig config, int base_fingerprint_bits) {
    return 0;
  }

  int insert(uint64_t fingerprint, uint64_t rank, uint64_t key) {
    storage[fingerprint].push_back(key);
    return 0;
  }

  void bulkInsert(uint64_t *keys, uint64_t *fingerprints, uint64_t *ranks, uint64_t numKeys) {
    for (uint64_t i = 0; i < numKeys; i++) {
      insert(fingerprints[i], ranks[i], keys[i]);
    }
  }

  int getFingerprint(uint64_t fingerprint, int rank, uint64_t *value) {
    auto it = storage.find(fingerprint);
    if (it == storage.end()) return -1;
    if (rank >= (int)it->second.size()) return -1;
    *value = it->second[rank];
    return 0;
  }

  int close() { return 0; }

private:
  std::unordered_map<uint64_t, std::vector<uint64_t>> storage;
};

#endif
