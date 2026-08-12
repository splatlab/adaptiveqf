#ifndef QF_FILTER_H
#define QF_FILTER_H

#include <cstdint>
#include "reverse_map_config.hpp"

struct QFilterConfig {
  size_t qbits;
  size_t rbits;
  double max_load_factor;
  int hash_mode;           // QF_HASH_DEFAULT, QF_HASH_INVERTIBLE, or QF_HASH_NONE
  ReverseMapConfig reverse_map;
};

struct QFilterQueryResult {
  int key_present;
  int minirun_rank;
  int minirun_count;
  uint64_t hash;
  uint64_t hash_index;
};

#endif