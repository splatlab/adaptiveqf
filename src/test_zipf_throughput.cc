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

void csv_get(char* buffer, int col) {
	int i, j;
	for (i = 0; buffer[i] != '\0' && col > 0; i++) {
		if (buffer[i] == ',') col--;
	}
	for (j = 0; buffer[i + j] != '\0' && buffer[i + j] != ','; j++) {
		buffer[j] = buffer[i + j];
	}
	buffer[j] = '\0';
}


typedef struct test_struct {
	int a;
	int b;
} test_struct;

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

typedef stxxl::map<uint64_t, test_struct, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> map_t2;

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
		BACKING_MAP_T::iterator item = map.find(ret_hash | (1 << ret_hash_len));
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
			BACKING_MAP_INSERT(map, ret_other_hash | (1 << ret_hash_len), item->second);
			BACKING_MAP_INSERT(map, ret_hash | (1 << ret_hash_len), key);
			map_inserts++;
			map_inserts++;
                }
        }
        else if (ret == 1) {
		BACKING_MAP_INSERT(map, ret_hash | (1 << ret_hash_len), key);
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
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF\nthe number of remainder bits in the QF\nthe number of queries\n");
		// ./test 16 7 $((1 << 15)) 1000000 0
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
	uint64_t num_queries = strtoull(argv[3], NULL, 10);


	printf("initializing hash table...\n");
#if USE_UNORDERED_MAP
	unordered_map_t backing_map;
	//backing_map.max_buffer_size(SUB_BLOCK_SIZE);
#else
	ordered_map_t backing_map((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
#endif	
	ordered_map_t database((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);

	std::cout << "generating zipfian distribution... 0%" << std::flush;
	uint64_t *query_set = (uint64_t*)calloc(num_queries, sizeof(uint64_t));
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));
	for (uint64_t ii = 0; ii < num_queries; ii++) {
		query_set[ii] = (uint64_t)rand_zipfian(1.5f, 1000000ull, query_set[ii]);
		if (ii % 100000 == 0) std::cout << "\rgenerating zipfian distribution... " << ii * 100 / num_queries << '%' << std::flush;
	}
	std::cout << "\rgenerating zipfian distribution... 100%" << std::endl;

	int num_trials = 50;
	double *trial_throughputs = new double[num_trials];
	double *trial_fprates = new double[num_trials];
	/*for (int i_trial = 0; i_trial < num_trials; i_trial++) {
		printf("== TRIAL %d\n", i_trial + 1);*/
		printf("initializing filter...\n");
		QF qf;
		if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
			fprintf(stderr, "Can't allocate CQF.\n");
			abort();
		}
		qf_set_auto_resize(&qf, false);


		printf("generating insert set...\n");
		uint64_t i, j;
		uint64_t *insert_set = new uint64_t[num_inserts];
		RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));


		// PERFORM INSERTS
		std::cout << "performing inserts... 0" << '%' << " load" << std::flush;
		uint64_t progress_interval = nslots / 100;
		uint64_t ret_index, ret_hash;
		int ret_hash_len;

		uint64_t target_fill = nslots * load_factor;

		clock_t start_time = clock(), end_time;
		for (i = 0; qf.metadata->noccupied_slots <= target_fill; i++) {
			if (!insert_key(&qf, backing_map, insert_set[i], 1)) break;
			database.insert(pair_t(insert_set[i], i));
			if (i % progress_interval == 0) std::cout << "\rperforming inserts... " << i * 100 / nslots << '%' << " load" << std::flush;
		}
		end_time = clock();
		std::cout << "\rperforming inserts... " << i * 100 / nslots << '%' << " load" << std::endl;


		// PERFORM QUERIES
		printf("generating query set...\n");
		int still_have_space = 1;
		if (qf.metadata->noccupied_slots >= qf.metadata->nslots * 0.95) {
			still_have_space = 0;
			printf("filter is full; skipping query adaptations\n");
		}

	for (int i_trial = 0; i_trial < num_trials; i_trial++) {
		uint64_t hash_seed = rand();
                for (uint64_t ii = 0; ii < num_queries; ii++) query_set[ii] = MurmurHash64A((unsigned char*)(&query_set[ii]), sizeof(uint64_t), hash_seed);

		printf("performing queries...\n");
		uint64_t fp_count = 0;
		start_time = clock();
		for (i = 0; i < num_queries; i++) {
			j = query_set[i];

			if (qf_query(&qf, j, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
				ordered_map_t::iterator item = database.find(j);
				if (item == database.end()) {
					fp_count++;
					if (still_have_space) {
						BACKING_MAP_T::iterator orig_key = backing_map.find(ret_hash | (1 << ret_hash_len));
						ret_hash_len = qf_adapt(&qf, ret_index, orig_key->second, j, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
						if (ret_hash_len > 0) {
							backing_map.erase(orig_key);
							backing_map.insert(std::make_pair(ret_hash | (1 << ret_hash_len), orig_key->second));
						}
						else if (ret_hash_len == QF_NO_SPACE) {
							still_have_space = 0;
							printf("filter is full after %lu queries\n", i);
						}
					}
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

		/*qf_free(&qf);
		database.clear();
		backing_map.clear();
	}*/

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
}

