# API Reference

## C API (`aqf.h`)

```c
#include "aqf.h"
```

### Lifecycle

```c
bool aqf_malloc(QF *qf, uint64_t nslots, uint64_t key_bits,
                uint64_t value_bits, enum qf_hashmode hash, uint32_t seed);
void aqf_free(QF *qf);
```

### Core operations

```c
int aqf_insert(QF *qf, uint64_t key, uint64_t count,
               qf_insert_result *result, uint8_t flags);
```

Inserts `key` (or increments count). `result` returns the minirun ID and rank
for the reverse map. Returns 0 on success, negative on error.

```c
int aqf_query(const QF *qf, uint64_t key, uint64_t *ret_hash,
              uint8_t *ret_minirun_rank, uint8_t flags);
```

Returns the count if found, 0 if not found (may be a false positive).

```c
int aqf_adapt(QF *qf, uint64_t orig_key, uint64_t fp_key,
              uint64_t minirun_rank, uint8_t flags);
```

Extends the fingerprint of `orig_key` so `fp_key` is no longer a false positive.

### Metadata

```c
uint64_t aqf_get_noccupied_slots(const QF *qf);
uint64_t aqf_get_nslots(const QF *qf);
uint64_t aqf_get_total_size_in_bytes(const QF *qf);
void     aqf_dump_metadata(const QF *qf);
void     aqf_dump(const QF *qf);
```

### Flags

| Flag | Value | Meaning |
|------|-------|---------|
| `QF_NO_LOCK` | `0x01` | Skip locking (caller must synchronize) |
| `QF_KEY_IS_HASH` | `0x08` | `key` is already a hash (skip hashing) |

### Hash modes

| Mode | Description |
|------|-------------|
| `QF_HASH_DEFAULT` | MurmurHash64A |
| `QF_HASH_INVERTIBLE` | Invertible hash |
| `QF_HASH_NONE` | Key is pre-hashed |

---

## C++ Wrapper (`aqf_wrapper.hpp`)

```cpp
#include "aqf_wrapper.hpp"
```

Templated on a reverse map implementation.

```cpp
template <typename ReverseMap>
class AdaptiveQF {
public:
  int init(QFilterConfig cfg);
  int bulkInsert(uint64_t *keys, uint64_t numKeys);
  int queryFilter(uint64_t key, QFilterQueryResult *result);
  int adapt(uint64_t key, QFilterQueryResult *result);
  int close();

  // Accessors
  double loadFactor();
  uint64_t sizeInBytes();
  void dumpMetadata();
  void dump();
  QF *getQF();
  ReverseMap *getReverseMap();
};
```

### `QFilterConfig`

```cpp
struct QFilterConfig {
  size_t qbits;            // log2 of number of slots
  size_t rbits;            // remainder bits
  double max_load_factor;  // adapt returns -1 when exceeded
  int hash_mode;           // QF_HASH_DEFAULT, _INVERTIBLE, or _NONE
  ReverseMapConfig reverse_map;
};
```

### `ReverseMapConfig`

```cpp
struct ReverseMapConfig {
  std::string storage_path;       // file path for persistent stores
  uint64_t cache_size_bytes;      // memory budget for caching
  bool overwrite;                 // delete existing data on init
};
```

---

## Reverse Maps

The `ReverseMap` template parameter must implement:

```cpp
struct ReverseMap {
  int  init(ReverseMapConfig config, int base_fingerprint_bits);
  int  insert(uint64_t fingerprint, uint64_t rank, uint64_t key);
  void bulkInsert(uint64_t *keys, uint64_t *fingerprints,
                  uint64_t *ranks, uint64_t numKeys);
  int  getFingerprint(uint64_t fingerprint, int rank, uint64_t *value);
  int  close();
};
```

| Implementation | File | Storage | Persistence |
|---|---|---|---|
| `MockDB` | `reverse_maps/mock/mock_db.hpp` | `std::unordered_map` | None (volatile) |
| `SplinterDB` | `reverse_maps/splinterdb/splinter_reverse_map.hpp` | SplinterDB | Disk-backed |
| `ll_table` | `reverse_maps/ll_table/` | In-memory hash table + linked lists | None (legacy C) |