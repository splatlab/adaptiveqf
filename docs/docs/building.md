<!-- BEGIN author: opencode -->
# Building

## Prerequisites

- gcc, g++ (C11 / C++17)
- OpenSSL (`libssl-dev`)
- Linux with `libaio-dev` and `libxxhash-dev`
- (Optional) SplinterDB — built automatically via submodule

## Build

```bash
cd adaptiveqf

# Build SplinterDB (first time only)
git submodule update --init external/splinterdb
cd external/splinterdb && make && cd ../..

# Build the AQF library
make

# Build examples
make -C examples
```

## Run the demos

```bash
# In-memory reverse map (no external deps)
./examples/adapt_demo

# SplinterDB-backed reverse map
./examples/adapt_demo_splinter
```

## Run the tests

```bash
# MockDB (no SplinterDB needed)
make -C tests test_e2e && ./tests/test_e2e

# SplinterDB
make -C tests test_e2e_splinter && ./tests/test_e2e_splinter --rm splinter
```

Both should report **0 false negatives** and identical false-positive counts.

<!-- END author: opencode -->