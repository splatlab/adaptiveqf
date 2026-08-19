#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include "aqf_wrapper.hpp"
#include "mock_db.hpp"

int main() {
  QFilterConfig cfg;
  cfg.qbits = 12;
  cfg.rbits = 8;
  cfg.max_load_factor = 0.75;
  cfg.hash_mode = QF_HASH_NONE;
  cfg.reverse_map = {};

  AdaptiveQF<MockDB> filter;
  if (filter.init(cfg)) return 1;

  srand(42);
  uint64_t keys[1000];
  for (int i = 0; i < 1000; i++) keys[i] = ((uint64_t)rand() << 32) | rand();
  filter.bulkInsert(keys, 1000);
  printf("Inserted 1000 keys\n");

  for (int i = 0; i < 100000; i++) {
    uint64_t q = ((uint64_t)rand() << 32) | rand();
    QFilterQueryResult r;
    filter.queryFilter(q, &r);
    if (!r.key_present) continue;

    printf("False positive at query %d: key=%lx  rank=%d\n",
           i, q, r.minirun_rank);

    do {
      filter.queryFilter(q, &r);
      if (!r.key_present) break;
      filter.adapt(q, &r);
    } while (r.key_present);

    printf("  Adapted, key no longer a false positive\n");
    filter.queryFilter(q, &r);
    printf("  Re-query: present=%d\n", r.key_present);
    break;
  }

  filter.close();
  return 0;
}