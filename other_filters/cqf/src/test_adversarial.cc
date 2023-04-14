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
}

#include <time.h>
#include <stxxl/unordered_map>
#include <stxxl/map>

uint64_t MurmurHash64A(const void* key, int len, uint64_t seed) {
	const uint64_t m = 0xc6a4a7935bd1e995LLU;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t* data = (const uint64_t*)key;
	const uint64_t* end = (len >> 3) + data;

	while (data != end) {
		uint64_t k = *data++;

		k *= m; 
		k ^= k >> r; 
		k *= m; 

		h ^= k;
		h *= m; 
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch (len & 7) {
		case 7: h ^= (uint64_t)(data2[6]) << 48;
		case 6: h ^= (uint64_t)(data2[5]) << 40;
		case 5: h ^= (uint64_t)(data2[4]) << 32;
		case 4: h ^= (uint64_t)(data2[3]) << 24;
		case 3: h ^= (uint64_t)(data2[2]) << 16;
		case 2: h ^= (uint64_t)(data2[1]) << 8;
		case 1: h ^= (uint64_t)(data2[0]);
			h *= m;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
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

#define SUB_BLOCK_SIZE (8192)
#define SUB_BLOCKS_PER_BLOCK (256)

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

int main(int argc, char **argv)
{
#if 0
	uint64_t hash_seed_temp = rand();
	RAND_bytes((unsigned char*)(&hash_seed_temp), sizeof(uint64_t));
	for (int n = 0; n < 50; n++) {
		uint64_t z = rand_zipfian(1.5f, 100);
		printf("%lu\n", MurmurHash64A((unsigned char*)(&z), sizeof(uint64_t), hash_seed_temp));
	}
	return 0;
#endif

	QF qf;
	uint64_t qbits = 20;
	uint64_t rbits = 9;
	uint64_t num_queries = 10000000;

	if (argc < 4) {
		printf("provide qbits, rbits, num_queries\n");
		exit(0);
	}
	else {
		qbits = atoi(argv[1]);
		rbits = atoi(argv[2]);
		num_queries = strtoull(argv[3], NULL, 10);
	}

	uint64_t nhashbits = qbits + rbits;
	uint64_t nslots = (1ULL << qbits);
	uint64_t num_inserts = nslots * 0.96;
	uint64_t *insert_set;
	printf("qbits=%lu\trbits=%lu\n", qbits, rbits);

	/* Initialise the CQF */
	/*if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {*/
	/*fprintf(stderr, "Can't allocate CQF.\n");*/
	/*abort();*/
	/*}*/
	printf("initializing filter\n");
	if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
		fprintf(stderr, "Can't allocate CQF.\n");
		abort();
	}
	qf_set_auto_resize(&qf, false);

	printf("initializing maps\n");
#if USE_UNORDERED_MAP
	unordered_map_t backing_map;
        //backing_map.max_buffer_size(SUB_BLOCK_SIZE);
#else
	ordered_map_t backing_map((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
#endif
        ordered_map_t database((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);

	printf("generating insert set\n");
	/* Generate random values */
	//srand(time(NULL));
	srand(time(NULL));
	insert_set = new uint64_t[num_inserts];
	RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));

	printf("performing %lu inserts\n", num_inserts);
	double measure_interval = 0.05f;
        double current_interval = measure_interval;
        uint64_t measure_point = nslots * current_interval, last_point = 0;

	uint64_t i;
        clock_t start_time = clock(), end_time, interval_time = start_time;
        for (i = 0; i < num_inserts/*qf.metadata->noccupied_slots <= target_fill*/; i++) {
		int ret = qf_insert(&qf, insert_set[i], 0, 1, QF_NO_LOCK);
		if (ret < 0) {
			fprintf(stderr, "failed insertion for key: %lx %d.\n", insert_set[i], 50);
			if (ret == QF_NO_SPACE) {
				fprintf(stderr, "CQF is full.\n");
				break;
			}
			else if (ret == QF_COULDNT_LOCK)
				fprintf(stderr, "TRY_ONCE_LOCK failed.\n");
			else
				fprintf(stderr, "Does not recognise return value.\n");
			abort();
		}
#if USE_UNORDERED_MAP
		backing_map.insert(std::make_pair(i, insert_set[i]));
#else
		backing_map.insert(pair_t(i, insert_set[i]));
#endif
		map_inserts++;
                database.insert(pair_t(insert_set[i], i));
                if (i >= measure_point) {
                        printf("throughput for interval %f: \t%f\n", current_interval, (double)CLOCKS_PER_SEC * (i - last_point) / (clock() - interval_time));
                        printf("map inserts: %d\tmap_queries: %d\n", map_inserts, map_queries);
			printf("time since last measurement: %ld\n", clock() - interval_time);
                        current_interval += measure_interval;
                        last_point = measure_point;
                        measure_point = nslots * current_interval;
                        interval_time = clock();
                }
        }
        end_time = clock();

	if (0) {
		start_time = clock();
		uint64_t num_deletes = 0.05f * nslots;
		for (i = 0; i < num_deletes; i++) {
			ordered_map_t::iterator map_iter = database.find(insert_set[i]);
			database.erase(map_iter);
		}
		printf("delete throughput: %f\n", 1000000.f * num_deletes / (clock() - start_time));
		return 0;
	}

	printf("generating query set\n");
	uint64_t *query_set = new uint64_t[num_queries];
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));

	printf("performing %lu queries\n", num_queries);
	uint64_t *fp_set = new uint64_t[100];
	uint64_t fp_set_len = 0;

	uint64_t fp_attack_frequencies[10] = {2, 3, 4, 5, 10, 20, 50, 100, 1000, 1000000};
	uint64_t i_fp_attack_frequency = 0;

	start_time = clock();
	/* Lookup inserted keys and counts. */
	for (i_fp_attack_frequency = 0; i_fp_attack_frequency < sizeof(fp_attack_frequencies) / sizeof(uint64_t); i_fp_attack_frequency++) {
		uint64_t fp_attack_frequency = fp_attack_frequencies[i_fp_attack_frequency];
		uint64_t fp_count = 0, value;
		start_time = clock();
		for (i = 0; i < num_queries; i++) {
			if (i % fp_attack_frequency || !fp_set_len || !fp_attack_frequency) {
				if (qf_query(&qf, query_set[i], &value, QF_NO_LOCK)) {
#if USE_UNORDERED_MAP
					unordered_map_t::iterator orig_key = backing_map.find(query_set[i]);
#else
					ordered_map_t::iterator orig_key = backing_map.find(query_set[i]);
#endif
					if (orig_key == backing_map.end() || orig_key->first == 0 || orig_key->second != query_set[i]) {
						fp_count++;
						if (fp_set_len < 100) {
							fp_set[fp_set_len++] = query_set[i];
						}
					}
				}
			}
			else {
				uint64_t i_fp_query = query_set[i] % fp_set_len;
				uint64_t fp_query = fp_set[i_fp_query];
				if (qf_query(&qf, fp_query, &value, QF_NO_LOCK)) {
#if USE_UNORDERED_MAP
					unordered_map_t::iterator orig_key = backing_map.find(fp_query);
#else
					ordered_map_t::iterator orig_key = backing_map.find(fp_query);
#endif
					if (orig_key == backing_map.end()) {
						fp_count++;
					}
					else {
						fp_set[i_fp_query] = fp_set[fp_set_len--];
					}
				}
				else {
					fp_set[i_fp_query] = fp_set[fp_set_len--];
				}
			}
		}
		end_time = clock();
		printf("==query test with repeated false positive every %lu queries==\n", fp_attack_frequency);
		printf("query throughput:     %f ops/sec\n", (double)CLOCKS_PER_SEC * num_queries / (end_time - start_time));
		printf("false positive rate:  %f\n", (double)fp_count / num_queries);
	}
}

