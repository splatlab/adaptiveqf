# Refactoring TODO — adaptiveqf/

This file lists every code-level work item for refactoring `adaptiveqf/`.
Use it as a checklist. Items are tagged by phase: `[TODO-P1]` through
`[TODO-P5]`. Dead-code items are tagged `[DEAD]`. Vendored items are
tagged `[VENDORED]`.

---

## PHASE 1 — Document & Map

### File Inventory

#### Core Implementation
- [ ] [TODO-P1] `src/gqf.c` (3081 lines) — Legacy CQF implementation.
  Being replaced by `src/aqf.c`. See Phase 4 for the migration plan.
- [ ] [TODO-P1] `src/aqf.c` (~1000 lines) — AQF implementation.
  Self-contained: insert, query, adapt, metadata, debug, locking.
  One file, same structure as gqf.c but simplified (no counters, no
  extensions, no merge/iterate/resize yet).
- [ ] [TODO-P1] `include/aqf.h` — Public C API for the AQF.
  Currently a subset of gqf.h. Needs to grow to cover all features
  that callers depend on.
- [ ] [TODO-P1] `include/aqf_int.h` — Internal struct definitions shared
  between aqf.c and the C++ wrapper.
- [ ] [TODO-P1] [DEAD] `src/gqf_backup.c` (3352 lines) — Near-duplicate
  of gqf.c. Run `diff src/gqf.c src/gqf_backup.c` and document the
  actual differences. Leave in place until Phase 4.
- [ ] [TODO-P1] `src/gqf_ski_rental.cpp` (3581 lines) — C++ version of
  gqf.c with SkiRental adaptation. Reference for correct adapt logic.
  Not compiled into libqf.a.

#### `taf.cc` (Telescoping? Unknown Acronym)

- [ ] [TODO-P1] `src/taf.cc` (3029 lines) — C++ file using STXXL (disk-backed
      STL). Purpose unclear — likely an earlier adaptive filter prototype.
      Investigate: is this dead code? Does any variant or test depend on it?
      If dead, mark [DEAD] and remove in Phase 4.

#### Vendored / Third-Party
- [ ] [TODO-P1] [VENDORED] `include/sglib.h` (1953 lines) — SGLIB data
  structure library. Audit which macros/functions are actually used.
  Candidate for removal once C++ STL replaces it (Phase 4).
- [ ] [TODO-P1] [DEAD] `other_filters/` — Snapshots of ACF, ACFSIM, CF,
  CQF, TAF. The CQF copy has compiled binaries committed. Either remove
  or vendor explicitly.

#### Hash Table Utility
- [ ] [TODO-P1] `src/hash_table.c` / `include/hash_table.h` — simple
  MurmurHash-based hash table. Is this used? Check callers.

#### Test Files — Deleted

All legacy test files (`src/test_*.c`, `src/test_*.cc`, `test/test_*.c`,
`src/test_driver.c`, `include/test_driver.h`) have been deleted.
See `CHANGELOG.md` for a record of what functionality each test covered.

Tests will be re-added from scratch once the refactored structure is stable.

---

## PHASE 2 — C++ Filter Layer

### `aqf_wrapper.hpp` (exists, needs work)
- [X] [TODO-P2] `include/aqf_wrapper.hpp` — C++ template wrapping AQF C API.
  Currently wraps `aqf_insert`, `aqf_query`, `aqf_adapt`.
- [ ] [TODO-P2] Add `AdaptiveQF::loadFromFile("dump.qf", cfg)` that calls
      `aqf_use()` to recreate a filter from a saved file.
- [ ] [TODO-P2] `base_fingerprint_bits` is derived from `qbits + rbits` and passed
      to `reverseMap.init()`. It's an internal detail the fingerprint uses to
      encode rank. Find a way to handle this inside the reverse map so the
      filter doesn't need to pass it explicitly.
- [ ] [TODO-P2] `fullPoint` (max_load_factor * nslots) is computed in the filter
      wrapper but logically belongs inside aqf.c. The filter internals should
      track load factor thresholds rather than the wrapper checking
      noccupied_slots externally.

### Examples
- [X] [TODO-P2] `examples/adapt_demo.cpp` — insert, trigger false positive,
      adapt, verify correction.
- [X] [TODO-P2] `examples/Makefile` — builds examples against libqf.a.
- [ ] [TODO-P2] `examples/basic_insert_query.c` — insert keys, query, check
      false positives.
- [ ] [TODO-P2] `examples/bulk_load.c` — bulk insert + measure throughput.

### Benchmark driver
- [ ] [TODO-P2] `bench/cpp/bench_filters.cc` — Single executable,
      parameterized by CLI args:
      `--load-factor`, `--num-inserts`, `--num-queries`,
      `--num-adapt-keys`, `--adversarial`, `--adaptive`, `--db-type`.
- [ ] [TODO-P2] Add `cxxopts` or similar CLI parsing dependency.
- [ ] [TODO-P2] Output: throughput, false-positive rate, adapt latency,
      filter load factor, memory usage.
- [ ] [TODO-P2] `bench/cpp/workload_gen.cc` — Workload generation.

### Build system (initial)
- [ ] [TODO-P2] Add compile targets: `bench_filters`, `workload_gen`.
- [ ] [TODO-P2] Use `sponge/build/` for output artifacts.

---

## PHASE 3 — Delete Old Tests

- [X] [TODO-P3] Delete all `src/test_*.c`, `src/test_*.cc` files.
- [X] [TODO-P3] Delete all `test/test_*.c` files.
- [X] [TODO-P3] Delete `src/test_driver.c` and `include/test_driver.h`.
- [X] [TODO-P3] Remove all test targets from the Makefile.
- [ ] [TODO-P3] Update `CHANGELOG.md` if any new bugs are discovered during
      deletion (records what each deleted test covered).

---

## PHASE 4 — Complete the AQF and Remove gqf.c

### Goal: `aqf.c` replaces `gqf.c` entirely. No caller depends on gqf.c.

### AQF API gaps (features in gqf.h that aqf.h doesn't expose)

#### Lifecycle
- [ ] [TODO-P4] Port `aqf_init` (buffer-based init, like `qf_init`)
- [ ] [TODO-P4] Port `aqf_destroy` (return buffer, like `qf_destroy`)
- [ ] [TODO-P4] Port `aqf_use` (deserialize from buffer, like `qf_use`)
- [ ] [TODO-P4] Port `aqf_resize` / `aqf_resize_malloc`
- [ ] [TODO-P4] Port `aqf_set_auto_resize`

#### Core operations
- [ ] [TODO-P4] Port `aqf_remove` (remove a fingerprint) — needed for
      SplinterDB reverse-map adapt path
- [ ] [TODO-P4] Port `aqf_delete_key_value` (remove all instances)
- [ ] [TODO-P4] Port `aqf_set_count` (set counter)

#### Metadata accessors
- [ ] [TODO-P4] Add `aqf_get_bits_per_slot`
- [ ] [TODO-P4] Add `aqf_get_quotient_bits`
- [ ] [TODO-P4] Add `aqf_get_key_remainder_bits`
- [ ] [TODO-P4] Add `aqf_get_hash_mode`
- [ ] [TODO-P4] Add `aqf_get_seed`
- [ ] [TODO-P4] Add `aqf_get_num_distinct`
- [ ] [TODO-P4] Add `aqf_get_nblocks`
- [ ] [TODO-P4] Add `aqf_get_num_key_bits`
- [ ] [TODO-P4] Add `aqf_get_num_value_bits`

#### Iterators
- [ ] [TODO-P4] Port iterator API (`aqf_iterator_from_position`,
      `aqf_iterator_from_key_value`, `aqfi_get_key`, `aqfi_get_hash`,
      `aqfi_next`, `aqfi_end`)

#### Misc
- [ ] [TODO-P4] Port `aqf_reset` (clear filter)
- [ ] [TODO-P4] Port `aqf_copy`
- [ ] [TODO-P4] Port `aqf_merge` / `aqf_multi_merge`
- [ ] [TODO-P4] Port `aqf_inner_product` / `aqf_magnitude`
- [ ] [TODO-P4] Port `aqf_hash_cmp`

### Remove gqf.c dependency from callers

- [ ] [TODO-P4] `reverse_maps/splinterdb/splinter_util.c` — currently calls
      `qf_insert_using_ll_table`. Switch to `aqf_insert`.
- [ ] [TODO-P4] `examples/adapt_demo.cpp` — legacy path uses `extern` declaration
      for `qf_insert_using_ll_table`. Remove once gqf.c is gone.
- [ ] [TODO-P4] Remove `insert_aqf_copy` experiment from `gqf.c` (added during
      debugging, no longer needed).
- [ ] [TODO-P4] Remove `extern` declarations from `aqf.c` that reference gqf.c
      functions (`insert_and_extend`, `qf_get_count_using_ll_table`,
      `qf_adapt_using_ll_table`).

### Cleanup

- [ ] [TODO-P4] Delete `src/gqf.c`
- [ ] [TODO-P4] Delete `include/gqf.h`
- [ ] [TODO-P4] Delete `include/gqf_int.h`
- [ ] [TODO-P4] Remove `gqf.o` from the Makefile build targets
- [ ] [TODO-P4] Remove `gqf_ski_rental.cpp` from the repo (or keep as reference)
- [ ] [TODO-P4] Remove SEVEN_BIT_OFFSET blocks from `gqf_ski_rental.cpp`
      (only remaining file with them)
- [ ] [TODO-P4] Remove `insert_and_extend` from aqf.c's insert (AQF always
      uses count=1, the `count > 1` path is dead code)
- [ ] [TODO-P4] Encapsulate `qfblock` — replace direct struct access
      with accessor functions throughout aqf.c
- [ ] [TODO-P4] Remove `#include "ll_table.h"` from aqf_int.h (no longer needed
      once the SplinterDB reverse map is migrated)

### Resolve `src/gqf_backup.c`
- [ ] [TODO-P4] [DEAD] Run: `diff src/gqf.c src/gqf_backup.c`
- [ ] [TODO-P4] Determine if any logic in gqf_backup.c is missing from aqf.c
- [ ] [TODO-P4] Delete `src/gqf_backup.c`

### Clean vendored code
- [ ] [TODO-P4] [VENDORED] Audit sglib.h usage. Which macros are called?
      Likely candidates: `SGLIB_LIST_*`, `SGLIB_HASH_*` macros.
      Replace with `std::list`, `std::unordered_map`, etc.
- [ ] [TODO-P4] [DEAD] Remove `other_filters/` directory or move to
      `vendor/` with attribution.

---

## PHASE 5 — Clean Build System

### Makefile overhaul
- [ ] [TODO-P5] Current: 57+ explicit targets, repetitive rules.
      Target: 4–5 targets with pattern rules.
- [ ] [TODO-P5] Use `sponge/build/` directory structure:
    - `sponge/build/obj/` — `.o` files
    - `sponge/build/bin/` — executables
    - `sponge/build/dep/` — `.d` dependency files
- [ ] [TODO-P5] Add automatic dependency generation:
      `CFLAGS += -MMD -MP`
      `DEPS = $(OBJS:.o=.d)`
      `-include $(DEPS)`
- [ ] [TODO-P5] Reduce targets to:
    - `all` → `bench_filters` `libqf.a`
    - `bench_filters` — unified benchmark
    - `libqf.a` — static library
    - `workload_gen` — (if needed)
    - `clean`
- [ ] [TODO-P5] Match variable conventions from SkiRental/Deamortized:
      `D=1` (debug), `P=1` (profile), `NH=1` (no Haswell).
- [ ] [TODO-P5] Remove stale SplinterDB targets or make them optional
      with a `SPLINTER=1` flag.

### Formatting & conventions
- [ ] [TODO-P5] Consistent indentation throughout Makefile (tabs, not spaces).
- [ ] [TODO-P5] Consistent variable naming: `CC`, `CXX`, `CFLAGS`,
      `CXXFLAGS`, `LDFLAGS`, `OBJDIR`, `BINDIR`.

---

## Checklist Summary

| Phase | Items | Status |
|---|---|---|
| P1 — Document | File inventory, dependency map | ░░░░░ |
| P2 — C++ Layer | aqf_wrapper.hpp, examples/bench | ▒▒▒▒░ |
| P3 — Delete Tests | Legacy test files removed | ▒▒▒▒▒ |
| P4 — Complete AQF | Port remaining gqf.c features, remove gqf.c | ▒▒▒░░ |
| P5 — Build | Makefile cleanup, dependency tracking | ░░░░░ |