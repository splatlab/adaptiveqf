#ifndef MOCK_DB_H
#define MOCK_DB_H

#include <cstdint>
#include <string>
#include "reverse_map_config.hpp"

class MockDB {
public:
  // [TODO] base_fingerprint_bits is an internal detail — find a better way
  // to allocate rank without requiring the reverse map to know about it.
  int init(ReverseMapConfig config, int base_fingerprint_bits) {
    this->base_fingerprint_bits = base_fingerprint_bits;
    return 0;
  }

  int insertFingerprint(uint64_t fingerprint, uint64_t rank, uint64_t key) {
    return 0;
  }

  int insertAndCommitFingerprint(uint64_t fingerprint, uint64_t rank,
                                  uint64_t key) {
    return 0;
  }

  int getFingerprint(uint64_t fingerprint, int rank, uint64_t *value) {
    uint64_t mask = -1;
    mask = mask << base_fingerprint_bits;
    *value = fingerprint | mask;
    return 0;
  }

  int insertKV(uint64_t key, uint64_t value, int isUpdate) { return 0; }

  void bulkInsert(uint64_t *keys, uint64_t *fingerprints, uint64_t numKeys) {
    for (uint64_t i = 0; i < numKeys; i++) {
      insertAndCommitFingerprint(fingerprints[i], 0, keys[i]);
    }
  }

  void commitFingerprints() {}

  int lookupKey(uint64_t key) { return 0; }

  int close() { return 0; }

private:
  int base_fingerprint_bits;
};

#endif