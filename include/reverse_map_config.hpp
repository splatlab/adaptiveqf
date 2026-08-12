#ifndef REVERSE_MAP_CONFIG_H
#define REVERSE_MAP_CONFIG_H

#include <string>
#include <cstdint>

// Configuration for initializing a reverse map implementation.
//
// Interpretation of fields varies by implementation:
//   - MockDB: ignores storage_path, cache_size_bytes, overwrite
//   - SplinterDB: storage_path is a file path
//   - WiredTiger/RocksDB: storage_path is a directory path
//
// Implementations should create parent directories if needed.

struct ReverseMapConfig {
  std::string storage_path;            // filesystem path for persistent storage
  uint64_t cache_size_bytes;           // memory budget for caching, in bytes.
                                       // This is a recommendation — implementations
                                       // may not adhere to it (e.g., MockDB ignores it).
  bool overwrite;                      // if true, delete existing data at
                                       // storage_path before initializing
#ifdef BENCH
  int collect_stats;
  uint64_t sleep_us;                   // delay per getFingerprint call to simulate
                                       // storage-level latency in benchmarks
#endif
};

#endif