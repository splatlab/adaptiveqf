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

typedef stxxl::map<uint64_t, uint64_t, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> ordered_map_t;
typedef std::pair<uint64_t, uint64_t> pair_t;

int main(int argc, char **argv) {
	uint64_t seed;
	size_t total_queries = 1000000;
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

	// Create a cuckoo filter where each item is of type size_t and
	// use 12 bits for each item:
	//    CuckooFilter<size_t, 12> filter(total_items);
	// To enable semi-sorting, define the storage of cuckoo filter to be
	// PackedTable, accepting keys of size_t type and making 13 bits
	// for each key:
	//   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
	printf("initializing data structures...\n");
	CuckooFilter<size_t, 12> filter(total_inserts);
	unordered_map_t backing_map;
	//backing_map.max_buffer_size(SUB_BLOCK_SIZE);
	ordered_map_t database((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);

	printf("generating inserts...\n");
	uint64_t *inserts = new uint64_t[total_inserts];
	RAND_bytes((unsigned char*)(inserts), total_inserts * sizeof(uint64_t));

	printf("starting inserts...\n");

	// Insert items to this cuckoo filter
	double measure_interval = 0.05f;
	double current_interval = measure_interval;
	uint64_t measure_point = nslots * current_interval, last_point = 0;
	size_t num_inserted = 0;
	clock_t start_time = clock(), end_time, interval_time = start_time;

	for (size_t i = 0; i <= total_inserts; i++, num_inserted++) {
		if (filter.Add(inserts[i], backing_map) != cuckoofilter::Ok) {
			//filter.Add(inserts[i], backing_map);
			printf("insert failed at fill rate %f\n", (double)i/nslots);
			break;
		}
		database.insert(pair_t(inserts[i], i));
		if (i >= measure_point) {
			printf("throughput for interval %f: \t%f\n", current_interval, 1000000. * (num_inserted - last_point) / (clock() - interval_time));
			current_interval += measure_interval;
			last_point = measure_point;
			measure_point = nslots * current_interval;
			interval_time = clock();
			printf("inserts: %d\tkickouts: %d\tadapts: %d\n", filter.map_inserts, filter.map_kickouts, filter.map_adapts);
		}
	}
	end_time = clock();

	printf("made %lu inserts\n", num_inserted);
	printf("time per insert:      %f us\n", (double)(end_time - start_time) / num_inserted);
	printf("insert throughput:    %f ops/sec\n", 1000000. * num_inserted / (end_time - start_time));
	//return 0;

	/*uint64_t *query_set = (uint64_t*)calloc(total_items, sizeof(uint64_t));
	for (size_t i = 0; i < total_items; i++) {
		//query_set[i] = rand_zipfian(1.5f, 1lu << 30);
		query_set[i] = i + total_items;
	}*/

	// Check non-existing items, a few false positives expected
	size_t fp_queries = 0;
	size_t p_queries = 0;
	/*for (size_t i = total_items; i < 2 * total_items; i++) {
	  if (filter.Contain(i) == cuckoofilter::Ok) {
	  false_queries++;
	  }
	  total_queries++;
	  }*/

	printf("generating queries...\n");
	uint64_t *queries = new uint64_t[total_queries];
	/*for (size_t i = 0; i < total_queries; i++) {
		queries[i] = rand_zipfian(1.5f, 1ull << 30);
		//queries[i] = rand_uniform();
	}*/
	RAND_bytes((unsigned char*)queries, total_queries * sizeof(uint64_t));

	printf("performing queries...\n");
	start_time = clock();
	for (size_t i = 0; i < total_queries; i++) {
		uint64_t key = queries[i];
		//uint64_t key = query_set[i];
		uint64_t contain_ret = filter.ContainReturn(key);
		if (contain_ret) {
			p_queries++;
			ordered_map_t::iterator item = database.find(key);
			if (item == database.end()) {
				fp_queries++;
				unordered_map_t::iterator orig_key = backing_map.find(contain_ret);
				if (filter.Adapt(orig_key->second, backing_map) != cuckoofilter::Ok) printf("error: adapt failed to find previously queried item\n");
			}
		}
	}
	end_time = clock();
	printf("made %lu queries\n", total_queries);
	printf("time per query:       %f us\n", (double)(end_time - start_time) / total_queries);
	printf("query throughput:     %f ops/sec\n", 1000000. * total_queries / (end_time - start_time));

	printf("false positive rate:  %f%%\n", 100. * fp_queries / total_queries);

	return 0;
}
