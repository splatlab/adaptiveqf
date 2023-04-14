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
	size_t churning_inserts = 10000;
	size_t nslots = 1 << 16;
	if (argc == 1) {
		printf("Main usage: ./test_ext_throughput [log of nslots] [num churns]\n");
		printf("            Using default values of 16 and 10000\n");
	}
	else if (argc >= 3) {
		nslots = 1 << atoi(argv[1]);
		churning_inserts = strtoull(argv[2], NULL, 10);
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
#if USE_UNORDERED_MAP
	unordered_map_t backing_map;
	//backing_map.max_buffer_size(SUB_BLOCK_SIZE);
#else
	ordered_map_t backing_map((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
#endif
	ordered_map_t database((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);

	printf("generating inserts...\n");
	uint64_t *inserts = new uint64_t[total_inserts];
	RAND_bytes((unsigned char*)(inserts), total_inserts * sizeof(uint64_t));

	std::cout << "performing inserts... 0%" << std::flush;

	// Insert items to this cuckoo filter
	double measure_interval = 0.01f;
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
			std::cout << "\rperforming inserts... " << int(current_interval * 100) << '%' << " - " << 1000000. * (num_inserted - last_point) / (clock() - interval_time) << " ops/sec" << std::flush;
			current_interval += measure_interval;
			last_point = measure_point;
			measure_point = nslots * current_interval;
			interval_time = clock();
		}
	}
	end_time = clock();
	std::cout << std::endl;

	printf("made %lu inserts\n", num_inserted);
	printf("time per insert:      %f us\n", (double)(end_time - start_time) / num_inserted);
	printf("insert throughput:    %f ops/sec\n\n", (double)CLOCKS_PER_SEC * num_inserted / (end_time - start_time));

	std::cout << "performing churn... 0%" << std::flush;
	current_interval = measure_interval;
	measure_point = churning_inserts * current_interval;
	last_point = 0;

	start_time = interval_time = clock();
	for (size_t i = 0; i < churning_inserts; i++) {
		int r = rand() % total_inserts;
		filter.Delete(inserts[r], backing_map);
		database.erase(inserts[r]);
		RAND_bytes((unsigned char*)&(inserts[r]), sizeof(uint64_t));
		filter.Add(inserts[r], backing_map);
		database.insert(pair_t(inserts[r], r));

		if (churning_inserts >= measure_point) {
			std::cout << "\rperforming churn... " << int(current_interval * 100) << '%' << " - " << 1000000. * (i - last_point) / (clock() - interval_time) << " ops/sec" << std::flush;
			current_interval += measure_interval;
			last_point = measure_point;
			measure_point = churning_inserts * current_interval;
			interval_time = clock();
		}
	}
	end_time = clock();
	std::cout << std::endl;

	printf("churned %lu times\n", churning_inserts);
	printf("time per churn:     %f us\n", (double)(end_time - start_time) / churning_inserts);
	printf("churn throughput:   %f ops/sec\n", (double)CLOCKS_PER_SEC * churning_inserts / (end_time - start_time));

	return 0;
}
