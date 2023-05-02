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
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <openssl/rand.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"
}

#include <iostream>

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

double rand_zipfian(double s, double max, uint64_t rnd) {
        double p = (double)(rnd % RAND_MAX) / RAND_MAX;

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
#define ZIPFIAN_QUERIES 0

#if USE_UNORDERED_MAP
#define BACKING_MAP_T unordered_map_t
#else
#define BACKING_MAP_T ordered_map_t
#endif

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
	std::cout << "progress: 0" << std::flush;
	for (int i = 1; i <= 10; i++) {
		
		std::cout << "\rprogress: " << i << std::flush;
	}
	std::cout << std::endl;
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
	printf("qbits: %lu\trbits: %lu\n", qbits, rbits);

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
	printf("CLOCKS_PER_SEC: %ld\n", CLOCKS_PER_SEC);
        clock_t start_time = clock(), end_time, interval_time = start_time;
        for (i = 0; i < num_inserts/*qf.metadata->noccupied_slots <= target_fill*/; i++) {
		int ret = qf_insert(&qf, insert_set[i], 0, 1, QF_NO_LOCK);
		if (ret < 0) {
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
                        printf("throughput for interval %f: \t%f\n", current_interval, (double)(i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time));
                        printf("map inserts: %d\tmap_queries: %d\n", map_inserts, map_queries);
                        current_interval += measure_interval;
                        last_point = measure_point;
                        measure_point = nslots * current_interval;
                        interval_time = clock();
                }
        }
        end_time = clock();

	if (0) {
		time(&start_time);
		uint64_t num_deletes = 0.05f * nslots;
		for (i = 0; i < num_deletes; i++) {
			ordered_map_t::iterator map_iter = database.find(insert_set[i]);
			database.erase(map_iter);
		}
		printf("delete throughput: %f\n", (double)num_deletes / difftime(time(NULL), start_time));
		return 0;
	}

	printf("generating query set\n");
	uint64_t *query_set = new uint64_t[num_queries];
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));

#if ZIPFIAN_QUERIES
	std::cout << "turning query set zipfian... 0%" << std::flush;
	int zipfianfy_progress = 5;
	uint64_t zipfianfy_progress_point = 0.01f * zipfianfy_progress * num_queries;
	for (i = 0; i < num_queries; i++) {
		query_set[i] = (uint64_t)rand_zipfian(1.5f, 1000000ull, query_set[i]);
		if (i >= zipfianfy_progress_point) {
			std::cout << "\rturning query set zipfian... " << zipfianfy_progress << "%" << std::flush;
			zipfianfy_progress += 5;
			zipfianfy_progress_point = 0.01f * zipfianfy_progress * num_queries;
		}
	}
	std::cout << std::endl;

	int num_trials = 100;
	double *trial_throughputs = new double[num_trials + 1];
	double *trial_fprates = new double[num_trials + 1];
	for (int i_trial = 0; i_trial < num_trials; i_trial++) {
		uint64_t hash_seed = rand();
		for (i = 0; i < num_queries; i++) {
			query_set[i] = MurmurHash64A((unsigned char*)(&query_set[i]), sizeof(uint64_t), hash_seed);
		}

		uint64_t fp_count = 0, value;
		start_time = clock();
		for (i = 0; i < num_queries; i++) {
			if (qf_query(&qf, query_set[i], &value, QF_NO_LOCK)) {
				ordered_map_t::iterator orig_val = database.find(query_set[i]);
				if (orig_val == database.end()) {
					fp_count++;
				}
			}
		}
		end_time = clock();

		double trial_throughput = (double)num_queries * CLOCKS_PER_SEC / (end_time - start_time);
		double trial_fprate = (double)fp_count / num_queries;
		printf("trial %d throughput: %f\n", i_trial + 1, trial_throughput);
		printf("trial %d fp rate:    %f\n", i_trial + 1, trial_fprate);
		int ii;
		for (ii = 0; ii < i_trial; ii++) if (trial_throughputs[ii] > trial_throughput) break;
		for (int jj = i_trial; jj > ii; jj--) trial_throughputs[jj] = trial_throughputs[jj - 1];
		trial_throughputs[ii] = trial_throughput;
		for (ii = 0; ii < i_trial; ii++) if (trial_fprates[ii] > trial_fprate) break;
		for (int jj = i_trial; jj > ii; jj--) trial_fprates[jj] = trial_fprates[jj - 1];
		trial_fprates[ii] = trial_fprate;
	}

	double avg_throughput = 0;
	for (int ii = 0; ii < num_trials; ii++) avg_throughput += trial_throughputs[ii];
	avg_throughput /= num_trials;
	printf("avg throughput: %f\n", avg_throughput);
	printf("percentile 0:   %f\n", trial_throughputs[0]);
	printf("percentile 25:  %f\n", trial_throughputs[(int)(0.25f * num_trials)]);
	printf("percentile 50:  %f\n", trial_throughputs[(int)(0.50f * num_trials)]);
	printf("percentile 75:  %f\n", trial_throughputs[(int)(0.75f * num_trials)]);
	printf("percentile 100: %f\n", trial_throughputs[num_trials - 1]);

	double avg_fprate = 0;
	for (int ii = 0; ii < num_trials; ii++) avg_fprate += trial_fprates[ii];
	avg_fprate /= num_trials;
	printf("avg fp rate: %f\n", avg_fprate);
	printf("percentile 0:   %f\n", trial_fprates[0]);
	printf("percentile 25:  %f\n", trial_fprates[(int)(0.25f * num_trials)]);
	printf("percentile 50:  %f\n", trial_fprates[(int)(0.50f * num_trials)]);
	printf("percentile 75:  %f\n", trial_fprates[(int)(0.75f * num_trials)]);
	printf("percentile 100: %f\n", trial_fprates[num_trials - 1]);

	printf("throughputs:\t");
	for (int ii = 0; ii < num_trials; ii++) printf("%f,", trial_throughputs[ii]);
	printf("\nfprates:\t");
	for (int ii = 0; ii < num_trials; ii++) printf("%f,", trial_fprates[ii]);
	printf("\n");
#else
	printf("performing %lu queries\n", num_queries);
	start_time = clock();
	uint64_t accesses = 0, value;
	/* Lookup inserted keys and counts. */
	for (i = 0; i < num_queries; i++) {
		if (qf_query(&qf, query_set[i], &value, QF_NO_LOCK)) {
			BACKING_MAP_T::iterator orig_key = backing_map.find(query_set[i]);
                        if (orig_key == backing_map.end() || orig_key->first == 0 || orig_key->second != query_set[i]) {
				accesses++;
			}
		}
	}
	end_time = clock();

	printf("query throughput: %f ops/sec\n", (double)num_queries * CLOCKS_PER_SEC / (end_time - start_time));
	printf("false positive rate: %f\n", (double)accesses / num_queries);
#endif

	int num_churns = 10000;
	uint64_t *churn_set = new uint64_t[2 * num_churns];
	RAND_bytes((unsigned char*)churn_set, 2 * num_churns * sizeof(uint64_t));
	clock_t churn_start_time = clock(), churn_end_time;
	for (int ii = 0; ii < num_churns; ii++) {
		int r = churn_set[2 * ii] % num_inserts;
		qf_remove(&qf, insert_set[r], 0, 1, QF_NO_LOCK);
		backing_map.erase(r);
		database.erase(insert_set[r]);
		insert_set[r] = churn_set[2 * ii + 1];
		qf_insert(&qf, insert_set[r], 0, 1, QF_NO_LOCK);
		backing_map.insert(pair_t(i, insert_set[r]));
                database.insert(pair_t(insert_set[r], i));
	}
	churn_end_time = clock();
	printf("churn throughput: %f ops/sec\n", (double)num_churns * CLOCKS_PER_SEC / (churn_end_time - churn_start_time));
}

