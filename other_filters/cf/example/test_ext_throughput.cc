#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <vector>
#include <openssl/rand.h>
#include <stxxl/unordered_map>
#include <stxxl/map>

using cuckoofilter::CuckooFilter;

#define DATA_NODE_BLOCK_SIZE (4096)
#define DATA_LEAF_BLOCK_SIZE (4096)

#define SUB_BLOCK_SIZE (8192)
#define SUB_BLOCKS_PER_BLOCK (256)

//! [hash]
struct HashFunctor {
        size_t operator () (uint64_t key) const {
               return (size_t)(key * 2654435761u);
        }
};
//! [hash]

//! [comparator]
struct CompareGreater {
        bool operator () (const uint64_t& a, const uint64_t& b) const {
                return a > b;
        }
        static uint64_t max_value() {
                return 0ULL;
        }
};
//! [comparator]

typedef stxxl::map<uint64_t, uint64_t, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> ordered_map_t;
typedef std::pair<uint64_t, uint64_t> pair_t;

typedef stxxl::unordered_map<uint64_t, uint64_t, HashFunctor, CompareGreater, SUB_BLOCK_SIZE, SUB_BLOCKS_PER_BLOCK> unordered_map_t;

#define USE_UNORDERED_MAP 0
#if USE_UNORDERED_MAP
#define BACKING_MAP_T unordered_map_t
#define BACKING_MAP_INSERT(X, Y, Z) X.insert(std::make_pair(Y, Z))
#else
#define BACKING_MAP_T ordered_map_t
#define BACKING_MAP_INSERT(X, Y, Z) X.insert(pair_t(Y, Z))
#endif

int main(int argc, char **argv) {
	uint64_t seed;
	size_t total_queries = 1000000;
	uint64_t total_churns = 100000;
	size_t nslots = 1 << 16;
	if (argc == 1) {
		printf("Main usage: ./test_ext_throughput [log of nslots] [num queries]\n");
		printf("            Using default values of 20 and 1000000\n");
	}
	else if (argc >= 3) {
		nslots = 1 << atoi(argv[1]);
		total_queries = strtoull(argv[2], NULL, 10);
	}
	if (argc >= 4) {
		seed = strtoull(argv[3], NULL, 10);
	}
	else {
		seed = time(NULL);
	}
	size_t total_inserts = nslots * 0.9;

	printf("running on seed %lu\n", seed);
	srand(seed);

	printf("generating inserts...\n");
	uint64_t *inserts = new uint64_t[total_inserts];
	RAND_bytes((unsigned char*)(inserts), total_inserts * sizeof(uint64_t));

	printf("generating queries...\n");
	uint64_t *queries = new uint64_t[total_queries];
	RAND_bytes((unsigned char*)(queries), total_queries * sizeof(uint64_t));

	printf("generating churns...\n");
	uint64_t *churns = new uint64_t[2 * total_churns];
	RAND_bytes((unsigned char*)(churns), 2 * total_churns * sizeof(uint64_t));

	// Create a cuckoo filter where each item is of type size_t and
	// use 12 bits for each item:
	//    CuckooFilter<size_t, 12> filter(total_items);
	// To enable semi-sorting, define the storage of cuckoo filter to be
	// PackedTable, accepting keys of size_t type and making 13 bits
	// for each key:
	//   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
	printf("initializing data structures...\n");
	CuckooFilter<size_t, 12> filter(total_inserts);
#if USE_UNORDERED_MAP
	unordered_map_t backing_map;
	//backing_map.max_buffer_size(SUB_BLOCK_SIZE);
#else
	ordered_map_t backing_map((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
#endif
	ordered_map_t database((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);

	printf("starting inserts...\n");

	// Insert items to this cuckoo filter
	double measure_interval = 0.05f;
	double current_interval = measure_interval;
	uint64_t measure_point = nslots * current_interval, last_point = 0;
	size_t num_inserted = 0;
	clock_t start_time = clock(), end_time, interval_time = start_time;

	for (size_t i = 0; i <= total_inserts; i++, num_inserted++) {
		if (filter.Add(inserts[i]) != cuckoofilter::Ok) {
			//filter.Add(inserts[i], backing_map);
			printf("insert failed at fill rate %f\n", (double)i/nslots);
			break;
		}
		BACKING_MAP_INSERT(backing_map, inserts[i], i);
		database.insert(pair_t(inserts[i], i));
		if (i >= measure_point) {
			printf("throughput for interval %f: \t%f\n", current_interval, 1000000. * (num_inserted - last_point) / (clock() - interval_time));
			current_interval += measure_interval;
			last_point = measure_point;
			measure_point = nslots * current_interval;
			interval_time = clock();
		}
	}
	end_time = clock();

	printf("made %lu inserts\n", num_inserted);
	printf("time per insert:      %f us\n", (double)(end_time - start_time) / num_inserted);
	printf("insert throughput:    %f ops/sec\n\n", 1000000. * num_inserted / (end_time - start_time));

	// Check non-existing items, a few false positives expected
	size_t fp_queries = 0;
	size_t p_queries = 0;


	std::cout << "performing queries... 0/" << total_queries << std::flush;
	measure_point = 1000;
	start_time = interval_time = clock();
	for (size_t i = last_point = 0; i < total_queries; i++) {
		uint64_t key = queries[i];
		if (filter.Contain(key)) {
			p_queries++;
			ordered_map_t::iterator item = database.find(key);
			if (item == database.end()) {
				fp_queries++;
				BACKING_MAP_T::iterator orig_key = backing_map.find(key);
			}
		}
		if (i >= measure_point) {
			std::cout << "\rperforming queries... " << i << "/" << total_queries << " - " << (double)(i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time) << " ops/sec" << std::flush;
			last_point = measure_point;
			measure_point += (i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time);
			interval_time = clock();
		}
	}
	end_time = clock();

	std::cout << std::endl;
	printf("made %lu queries\n", total_queries);
	printf("time per query:       %f us\n", (double)(end_time - start_time) / total_queries);
	printf("query throughput:     %f ops/sec\n", (double)total_queries * CLOCKS_PER_SEC / (end_time - start_time));

	printf("false positive rate:  %f%%\n\n", 100. * fp_queries / total_queries);


	std::cout << "performing churns... 0/" << total_churns << std::flush;
	measure_point = 1000;
	start_time = interval_time = clock();
	for (size_t i = 0; i < total_churns; i++) {
		int r = churns[2 * i] % total_inserts;
		filter.Delete(inserts[r]);
		database.erase(inserts[r]);
		backing_map.erase(inserts[r]);
		inserts[r] = churns[2 * i + 1];
		filter.Add(inserts[r]);
		database.insert(pair_t(inserts[r], i));
		BACKING_MAP_INSERT(backing_map, inserts[r], i);
		
		if (i >= measure_point) {
			std::cout << "\rperforming churns... " << i << "/" << total_churns << " - " << (double)(i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time) << " ops/sec" << std::flush;
			last_point = measure_point;
			measure_point += (i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time);
			interval_time = clock();
		}
	}
	end_time = clock();

	std::cout << std::endl;
	printf("made %lu churns\n", total_churns);
	printf("churn throughput:     %f ops/sec\n", (double)total_churns * CLOCKS_PER_SEC / (end_time - start_time));

	return 0;
}
