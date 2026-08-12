#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <unordered_set>
#include <vector>

#include "aqf_wrapper.hpp"
#include "mock_db.hpp"

int main() {
  QFilterConfig cfg;
  cfg.qbits = 12;
  cfg.rbits = 8;
  cfg.max_load_factor = 0.75;
  cfg.hash_mode = QF_HASH_INVERTIBLE;
  cfg.reverse_map.storage_path = "reverseMap";
  cfg.reverse_map.cache_size_bytes = 64ULL * 1024 * 1024;
  cfg.reverse_map.overwrite = true;

  AdaptiveQF<MockDB> filter;
  if (filter.init(cfg)) {
    fprintf(stderr, "Failed to init filter\n");
    return 1;
  }

  srand(42);
  std::vector<uint64_t> inserted(1000);
  for (int i = 0; i < 1000; i++) {
    inserted[i] = ((uint64_t)rand() << 32) | rand();
  }

  if (int ret = filter.bulkInsert(inserted.data(), 1000)) {
    fprintf(stderr, "Failed bulk insert: %d\n", ret);
    return 1;
  }

  std::unordered_set<uint64_t> adapted_keys;
  int total_fps = 0;
  int adapted_count = 0;
  int i;

  for (i = 0; i < 100000; i++) {
    uint64_t q = ((uint64_t)rand() << 32) | rand();
    QFilterQueryResult result;
    filter.queryFilter(q, &result);

    if (result.key_present) {
      total_fps++;

      int ret = filter.adapt(q, &result);
      if (ret > 0) {
        adapted_keys.insert(q);
        adapted_count++;

        QFilterQueryResult recheck;
        filter.queryFilter(q, &recheck);
        if (recheck.key_present) {
          fprintf(stderr, "FAIL: adapted key %lx still reports as present\n", q);
          return 1;
        }
      }
    }

    if (filter.loadFactor() >= 0.95) {
      printf("Stopping early at i=%d: load factor reached 0.95\n", i);
      break;
    }
  }

  printf("Inserted: 1000\n");
  printf("Queries: %d\n", i);
  printf("Total false positives: %d\n", total_fps);
  printf("Adapted: %d\n", adapted_count);
  printf("Unique adapted keys: %zu\n", adapted_keys.size());
  printf("Final load factor: %f\n", filter.loadFactor());
  printf("False positive rate: %f\n", (double)total_fps / i);

  filter.close();
  return 0;
}