# Getting Started

- [Building from source](building.md) — prerequisites, cloning, compilation

## Using the library

### Headers

```cpp
#include "aqf_wrapper.h"
```

### Building the filter

> TODO: For adapt to work, the AQF needs a reverse map

The AQF needs a **reverse map** — a separate data structure that stores the
original keys. When a false positive occurs, the AQF looks up the original key
in the reverse map and extends its fingerprint. The reverse map can be large
(it holds all original keys), and we provide several implementations.

```cpp
  QFilterConfig cfg;
  cfg.qbits = 12;
  cfg.rbits = 8;
  cfg.max_load_factor = 0.75;
  cfg.hash_mode = QF_HASH_NONE;

  AdaptiveQF<MockDB> filter;
  filter.init(cfg);

  uint64_t keys[1000] = { /* your data */ };
  filter.bulkInsert(keys, 1000);
```

> TODO: Add an API to construct from file

### Querying and Adapting

> TODO: Describe the query and adapt worklow

```cpp
  QFilterQueryResult r;
  filter.queryFilter(some_key, &r);
  if (r.key_present) {
    filter.adapt(some_key, &r);  // correct the false positive
  }

  filter.queryFilter(some_key, &r);
  assert(!r.key_present);

```

### Updating the filter

> TODO: Add the insert/delete API

### Resizing Policy


## API reference

- [API Reference](api.md) — full C and C++ API documentation
- [Quotient filter](qf.md) — Explain the quotient filter design
- [Reverse maps](reverse-maps.md) — why they're needed, size tradeoffs,
  provided implementations

