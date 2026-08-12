# Refactoring TODO — adaptiveqf/

This file lists every code-level work item for refactoring `adaptiveqf/`.
Use it as a checklist. Items are tagged by phase: `[TODO-P1]` through
`[TODO-P5]`. Dead-code items are tagged `[DEAD]`. Vendored items are
tagged `[VENDORED]`.

---

## PHASE 1 — Document & Map

### File Inventory

#### Core Implementation
- [ ] [TODO-P1] `src/gqf.c` (3081 lines) — Core QF implementation.
  Contains: insert, query, adapt, iterate, metadata, debug, locking,
  recording. ALL tangled in one file. See Phase 4 for modularization.
- [ ] [TODO-P1] [DEAD] `src/gqf_backup.c` (3352 lines) — Near-duplicate
  of gqf.c. Run `diff src/gqf.c src/gqf_backup.c` and document the
  actual differences. Leave in place until Phase 4.
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

## PHASE 2 — Add C++ Filter Layer

### Create `filters/`
- [ ] [TODO-P2] `filters/qf_filter.hpp` — Config/result structs.
      Model after `SkiRentalAdaptiveQF/variants/qf_filter.hpp`:
      `QFilterConfig`, `QFilterQueryResult`, `BenchmarkParams`.
- [ ] [TODO-P2] `filters/adaptiveqf.hpp` — Wraps C AQF API:
      `qf_insert_using_ll_table`, `qf_query_using_ll_table`,
      `qf_adapt_using_ll_table`. Should be ~35 lines.
      Model after `SkiRentalAdaptiveQF/variants/adaptiveqf.hpp`.
- [ ] [TODO-P2] Add a static factory method `AdaptiveQF::loadFromFile("dump.qf", cfg)`
      that calls `qf_use()` instead of `qf_malloc()` to recreate a filter from
      a saved file.
- [ ] [TODO-P2] `base_fingerprint_bits` is derived from `qbits + rbits` and passed
      to `reverseMap.init()`. It's an internal detail the fingerprint uses to
      encode rank. Find a way to handle this inside the reverse map so the
      filter doesn't need to pass it explicitly.
- [ ] [TODO-P2] `fullPoint` (max_load_factor * nslots) is computed in the filter
      wrapper but logically belongs inside gqf.c. The filter internals should
      track load factor thresholds rather than the wrapper checking
      noccupied_slots externally.

### Create `examples/`
- [ ] [TODO-P2] `examples/basic_insert_query.c` — insert keys, query, check
      false positives.
- [ ] [TODO-P2] `examples/adapt_demo.c` — insert, trigger a false positive,
      adapt, verify correction.
- [ ] [TODO-P2] `examples/bulk_load.c` — bulk insert + measure throughput.
- [ ] [TODO-P2] `examples/Makefile` — builds examples against the AQF library.

### Create benchmark driver
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
- [ ] [TODO-P2] Add `filters/` and `reverse_maps/*/` to include paths.
- [ ] [TODO-P2] Use `sponge/build/` for output artifacts.

---

## PHASE 3 — Delete Old Tests

- [ ] [TODO-P3] Delete all `src/test_*.c`, `src/test_*.cc` files.
- [ ] [TODO-P3] Delete all `test/test_*.c` files.
- [ ] [TODO-P3] Delete `src/test_driver.c` and `include/test_driver.h`.
- [ ] [TODO-P3] Remove all test targets from the Makefile.
- [ ] [TODO-P3] Update `CHANGELOG.md` if any new bugs are discovered during
      deletion (records what each deleted test covered).

---

## PHASE 4 — Decompose the C Core

### Modularization of `src/gqf.c`

Section boundaries in `gqf.c` (by line number):

| Lines | Module | Proposed file |
|---|---|---|
| 1–53 | Headers, macros, constants | (shared header) |
| 54–260 | Locking | `qf_lock.c` |
| 263–540 | Bit utilities (popcnt, select, rank) | `qf_bitutil.c` |
| 540–827 | Slot operations (get/set, run_end, shift) | `qf_slot.c` |
| 764–968 | Debug/dump | `qf_debug.c` |
| 968–987 | `insert_one_slot`, `get_slot_info` (first def) | `qf_slot.c` |
| 1021–1140 | Recording | `qf_record.c` |
| 1141–1418 | `insert`, `insert_using_ll_table` (internal) | `qf_insert.c` |
| 1419–1553 | Public query/insert-with-ll-table | `qf_query.c`, `qf_insert.c` |
| 1554–1693 | `insert_and_extend` | `qf_insert.c` |
| 1694–1936 | Init/use/destroy/malloc/free/copy/reset | `qf_metadata.c` |
| 1937–2090 | Resize, auto-resize, set_count, remove | `qf_resize.c`, `qf_insert.c` |
| 2091–2283 | Query, adapt (public) | `qf_query.c`, `qf_adapt.c` |
| 2285–2340 | Metadata accessors | `qf_metadata.c` |
| 2342–2556 | Iterators | `qf_iter.c` |
| 2557–2966 | Merge, bulk insert, multi-merge | `qf_merge.c` |
| 2967–3081 | Inner product, intersect, magnitude | `qf_misc.c` |

- [ ] [TODO-P4] Extract `qf_lock.c` (lines 54–260)
- [ ] [TODO-P4] Extract `qf_bitutil.c` (lines 263–540)
- [ ] [TODO-P4] Extract `qf_slot.c` (lines 540–827, plus 968–987)
- [ ] [TODO-P4] Extract `qf_debug.c` (lines 764–968)
- [ ] [TODO-P4] Extract `qf_record.c` (lines 1021–1140)
- [ ] [TODO-P4] Extract `qf_insert.c` (lines 1141–1418 plus 1554–2090)
- [ ] [TODO-P4] Extract `qf_query.c` (lines 1419–1553, 2091–2139)
- [ ] [TODO-P4] Extract `qf_adapt.c` (lines 2140–2283)
- [ ] [TODO-P4] Extract `qf_metadata.c` (lines 1694–1936, 2285–2340)
- [ ] [TODO-P4] Extract `qf_resize.c` (lines 1937–2090)
- [ ] [TODO-P4] Extract `qf_iter.c` (lines 2342–2556)
- [ ] [TODO-P4] Extract `qf_merge.c` (lines 2557–2966)
- [ ] [TODO-P4] Extract `qf_misc.c` (lines 2967–3081)
- [ ] [TODO-P4] Create `qf_internal.h` — shared internal macros/decls
      that the extracted modules all need (METADATA_WORD, BITMASK, etc.)

### Remove SEVEN_BIT_OFFSET compile-time flag
- [ ] [TODO-P4] Remove all `#ifndef SEVEN_BIT_OFFSET` / `#else` / `#endif`
      blocks from `gqf.c` (5+ locations: `block_offset`, `shift_runends`,
      `insert_one_slot`, `adapt`, `remove_replace_...`, `increment_block_counter`).
      Keep the full `uint8_t` offset path and delete the 7-bit variant.
- [ ] [TODO-P4] Remove `increment_block_counter()` — it exists only for the
      `SEVEN_BIT_OFFSET` path (sets bit 7 of the offset).
- [ ] [TODO-P4] Validate that no variant or test depends on the 7-bit behavior.
      See `CHANGELOG.md` for the full discussion.

### Resolve `src/gqf_backup.c`
- [ ] [TODO-P4] [DEAD] Run: `diff src/gqf.c src/gqf_backup.c`
- [ ] [TODO-P4] Determine which file is canonical. Options:
    - gqf_backup.c is strictly a subset → delete.
    - gqf_backup.c has diverging logic → merge differences or keep both.
- [ ] [TODO-P4] Document the result in a comment at the top of whichever
      file survives.

### Clean vendored code
- [ ] [TODO-P4] [VENDORED] Audit sglib.h usage. Which macros are called?
      Likely candidates: `SGLIB_LIST_*`, `SGLIB_HASH_*` macros.
      Replace with `std::list`, `std::unordered_map`, etc.
- [ ] [TODO-P4] [DEAD] Remove `other_filters/` directory or move to
      `vendor/` with attribution.

### API cleanup
- [ ] [TODO-P4] Remove duplicate function declarations from `gqf.h`
      (e.g., `qf_lock`/`qf_unlock` should probably not be public API).
- [ ] [TODO-P4] Add `extern "C"` guards to all headers (already done
      for gqf.h and gqf_int.h, but check others).
- [ ] [TODO-P4] Encapsulate `qfblock` — replace direct struct access
      with accessor functions throughout the codebase.
- [ ] [TODO-P4] Simplify the AQF interface: `bulkInsert(keys, numKeys)` should
      pass raw keys to the QF via `qf_insert_using_ll_table(key, ...)` without
      `QF_KEY_IS_HASH`. The QF handles hashing internally and returns the
      fingerprint, rank, and minirun info in the result struct. This removes
      the fingerprint pre-computation from `bulkInsert` and the `hash_mode`
      branching in `filters/adaptiveqf.hpp`.

### Shared core strategy
- [ ] [TODO-P4] Decide how `SkiRentalAdaptiveQF/` and
      `DeamortizedRehashing/` will share the C core once it's modularized:
      symlinks, shared `libqf/` directory, or sync script.

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
| P2 — C++ Layer | filters/, bench_filters.cc, examples/ | ░░░░░ |
| P3 — Delete Tests | Remove all legacy test files | ░░░░░ |
| P4 — Decompose C | gqf.c → modules, gqf_backup, vendored cleanup | ░░░░░ |
| P5 — Build | Makefile cleanup, dependency tracking | ░░░░░ |