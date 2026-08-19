# Reverse Maps

> TODO: This doc needs a review

The AQF decouples the filter from the reverse map. The filter stores compact
fingerprints; the reverse map stores the original keys needed during
adaptation.

## Why a separate structure?

When two keys collide on a fingerprint, the filter cannot distinguish them.
To adapt, the AQF needs to find the original key that caused the collision —
that lookup is the reverse map's job.

## Storage cost

The reverse map stores one entry per inserted key. Entry size depends on the
key size and the implementation overhead. For a dataset of `N` keys of
`K` bytes, expect roughly `N × (K + metadata)` bytes.

## Provided implementations

| Implementation | Location | Storage | Persistence |
|---|---|---|---|
| `MockDB` | `reverse_maps/mock/mock_db.hpp` | `std::unordered_map` | None (volatile) |
| `SplinterDB` | `reverse_maps/splinterdb/splinter_reverse_map.hpp` | SplinterDB | Disk-backed |
| `ll_table` | `reverse_maps/ll_table/` | In-memory hash table + linked lists | None (legacy C) |

**MockDB** is the simplest — an in-memory `std::unordered_map` mapping
fingerprints to vectors of keys. Use it for testing and development.

**SplinterDB** is a persistent key-value store that writes to disk. Use it
when the reverse map is too large for RAM or needs to survive restarts.

**ll_table** is a legacy C implementation using sglib-style linked lists.
It exists for backward compatibility during the refactor.

## Implementing your own

Any class implementing the `ReverseMap` interface can be used:

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

See the [API Reference](api.md) for details.

