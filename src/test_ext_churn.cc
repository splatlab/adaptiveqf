/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

extern "C" {
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <openssl/rand.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"
#include "include/hashutil.h"
}

#include <ctime>
#include <stxxl/unordered_map>
#include <stxxl/map>

int bp2() {
	return 0;
}

uint64_t rand_uniform(uint64_t max) {
	/*if (max <= RAND_MAX) return rand() % max;
	uint64_t a = rand();
	uint64_t b = rand();
	a |= (b << 32);
	return a % max;*/
	uint64_t a;
	RAND_bytes((unsigned char*)(&a), sizeof(uint64_t));
	if (max) a %= max;
	return a;
}

double rand_normal(double mean, double sd) {
	double a = (double)rand() / RAND_MAX;
	double b = (double)rand() / RAND_MAX;
	double c = sqrt(-2.0 * log(a)) * cos(2 * M_PI * b) * sd;
	return c + mean;
}

double rand_zipfian(double s, double max) {
	double p = (double)rand() / RAND_MAX;
	
	double pD = p * (12 * (pow(max, -s + 1) - 1) / (1 - s) + 6 + 6 * pow(max, -s) + s - s * pow(max, -s + 1));
	double x = max / 2;
	while (true) {
		double m = pow(x, -s - 2);
		double mx = m * x;
		double mxx = mx * x;
		double mxxx = mxx * x;
		
		double b = 12 * (mxxx - 1) / (1 - s) + 6 + 6 * mxx + s - (s * mx) - pD;
		double c = 12 * mxx - (6 * s * mx) + (m * s * (s + 1));
		double newx = x - b / c > 1 ? x - b / c : 1;
		if (abs(newx - x) <= 0.01) { // this is the tolerance for approximation
			return newx;
		}
		x = newx;
	}
}


#define USE_UNORDERED_MAP 0
#if USE_UNORDERED_MAP
#define BACKING_MAP_T unordered_map_t
#define BACKING_MAP_INSERT(X, Y, Z) X.insert(std::make_pair(Y, Z));
#else
#define BACKING_MAP_T ordered_map_t
#define BACKING_MAP_INSERT(X, Y, Z) X.insert(pair_t(Y, Z));
#endif

#define SUB_BLOCK_SIZE 8192
#define SUB_BLOCKS_PER_BLOCK 256

#define DATA_NODE_BLOCK_SIZE (4096)
#define DATA_LEAF_BLOCK_SIZE (4096)

//! [hash]
struct HashFunctor
{
	size_t operator () (int key) const {
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

typedef stxxl::unordered_map<uint64_t, uint64_t, HashFunctor, CompareGreater, SUB_BLOCK_SIZE, SUB_BLOCKS_PER_BLOCK> unordered_map_t;
typedef stxxl::map<uint64_t, uint64_t, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> ordered_map_t;
typedef std::pair<uint64_t, uint64_t> pair_t;

int map_inserts = 0;
int map_queries = 0;

int insert_key(QF *qf, BACKING_MAP_T& map, uint64_t key, int count) {
        uint64_t ret_index, ret_hash, ret_other_hash;
        int ret_hash_len;
        int ret = qf_insert_ret(qf, key, count, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
        if (ret == QF_NO_SPACE) {
                return 0;
        }
        else if (ret == 0) {
		BACKING_MAP_T::iterator item = map.find(ret_hash | (1ull << ret_hash_len));
		map_queries++;
                if (item == map.end()) {
                        printf("error:\tfilter claimed to have fingerprint %lu but hashtable could not find it\n", ret_hash);
                }
                else if (item->second == key) {
                        insert_and_extend(qf, ret_index, key, count, key, &ret_hash, &ret_other_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
                }
                else {
                        int ext_len = insert_and_extend(qf, ret_index, key, count, item->second, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
                        if (ext_len == QF_NO_SPACE) {
                                printf("filter is full after insert_and_extend\n");
                                return 0;
                        }

			map.erase(item);
			BACKING_MAP_INSERT(map, ret_other_hash | (1ull << ret_hash_len), item->second);
			BACKING_MAP_INSERT(map, ret_hash | (1ull << ret_hash_len), key);
			map_inserts++;
			map_inserts++;
                }
        }
        else if (ret == 1) {
		BACKING_MAP_INSERT(map, ret_hash | (1ull << ret_hash_len), key);
		map_inserts++;
        }
        else {
                printf("other error: errno %d\n", ret);
                return 0;
        }
        return 1;
}


int main(int argc, char **argv)
{
	if (argc < 4) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF\nthe number of remainder bits in the QF\nthe number of churns\n");
		fprintf(stderr, "./test_ext_churn 16 9 10000 0");
		exit(1);
	}
	if (argc >= 5) {
		srand(strtol(argv[4], NULL, 10));
		printf("running test on seed %ld\n", strtol(argv[5], NULL, 10));
	}
	else {
		time_t seed = time(NULL);
		printf("running test on seed %ld\n", seed);
		srand(seed);
	}
	printf("reading inputs...\n");
	uint64_t qbits = atoi(argv[1]);
	uint64_t rbits = atoi(argv[2]);
	int max_slots_per_item = (64 - qbits + rbits - 1) / rbits;
	int lbits = 0;
	while (max_slots_per_item > 0) {
		lbits++;
		max_slots_per_item >>= 1;
	}

	uint64_t nhashbits = qbits + rbits;
	uint64_t nslots = (1ULL << qbits);
	double load_factor = 0.9f;
	uint64_t num_inserts = nslots * load_factor;//strtoull(argv[3], NULL, 10);
	uint64_t num_churns = strtoull(argv[3], NULL, 10);


	printf("initializing hash table...\n");
#if USE_UNORDERED_MAP
	unordered_map_t backing_map;
	//backing_map.max_buffer_size(SUB_BLOCK_SIZE);
#else
	ordered_map_t backing_map((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
#endif	
	ordered_map_t database((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
	

	printf("initializing filter...\n");
	QF qf;
	if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
		fprintf(stderr, "Can't allocate CQF.\n");
		abort();
	}
	qf_set_auto_resize(&qf, false);


	printf("generating insert set of size %lu...\n", num_inserts);
	uint64_t i;
	uint64_t *insert_set = new uint64_t[num_inserts];
        RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));


	// PERFORM INSERTS
	printf("performing inserts...\n");
	uint64_t ret_hash;
	int ret_hash_len;

	uint64_t target_fill = nslots * load_factor;

	double measure_interval = 0.01f;
        double current_interval = measure_interval;
        uint64_t measure_point = nslots * current_interval, last_point = 0;

	std::cout << "performing inserts... 0" << '%' << std::flush;
	clock_t start_time = clock(), end_time, interval_time = start_time;
	for (i = 0; qf.metadata->noccupied_slots <= target_fill; i++) {
		if (!insert_key(&qf, backing_map, insert_set[i], 1)) break;
		database.insert(pair_t(insert_set[i], i));
		if (qf.metadata->noccupied_slots >= measure_point) {
			std::cout << "\rperforming inserts... " << int(current_interval * 100) << '%' << " - " << (double)(i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time) << " ops/sec" << std::flush;
                        current_interval += measure_interval;
                        last_point = i;
                        measure_point = nslots * current_interval;
                        interval_time = clock();
                }
	}
	end_time = clock();

	printf("time for inserts:      %f\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
	printf("avg insert throughput: %f ops/sec\n", (double)i * CLOCKS_PER_SEC / (end_time - start_time));

	current_interval = measure_interval;
	measure_point = num_churns * current_interval;
	last_point = 0;
	
	std::cout << "performing churns... 0" << '%' << std::flush;
	start_time = clock();
	for (i = 0; i < num_churns; i++) {
		int r = rand() % num_inserts;
		qf_remove(&qf, insert_set[r], &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK);
		database.erase(insert_set[r]);
		backing_map.erase(ret_hash | (1ull << ret_hash_len));
		RAND_bytes((unsigned char*)&(insert_set[r]), sizeof(uint64_t));
		insert_key(&qf, backing_map, insert_set[r], 1);
		database.insert(pair_t(insert_set[r], r));

		if (i >= measure_point) {
			std::cout << "\rperforming churns... " << int(current_interval * 100) << '%' << " - " << (double)(i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time) << " ops/sec" << std::flush;
                        current_interval += measure_interval;
                        last_point = i;
                        measure_point = num_churns * current_interval;
                        interval_time = clock();
                }

	}
	end_time = clock();
	std::cout << std::endl;

	printf("time for churns:   %f s\n", (double)(end_time - start_time));
	printf("churn throughput:  %f ops/sec\n", (double)num_churns * CLOCKS_PER_SEC / (end_time - start_time));
}

