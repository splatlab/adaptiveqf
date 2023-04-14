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

#define DATA_NODE_BLOCK_SIZE (4096)
#define DATA_LEAF_BLOCK_SIZE (4096)

typedef stxxl::map<uint64_t, uint64_t, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> ordered_map_t;
typedef std::pair<uint64_t, uint64_t> pair_t;

int main(int argc, char **argv) {
	uint64_t seed;
	size_t num_queries = 1000000;
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
	uint64_t *query_set = new uint64_t[num_queries];
	RAND_bytes((unsigned char*)queries, num_queries * sizeof(uint64_t));
	for (int ii = 0; ii < num_queries; ii++) query_set[ii] = rand_zipfian(1.5f, 1000000ull, query_set[ii]);


	for (int i_trial = 0; i_trial < num_trials; i_trial++) {
                uint64_t hash_seed = rand();
                for (uint64_t ii = 0; ii < num_queries; ii++) query_set[ii] = MurmurHash64A((unsigned char*)(&query_set[ii]), sizeof(uint64_t), hash_seed);

                printf("performing queries...\n");
                uint64_t fp_count = 0;
                start_time = clock();
                for (i = 0; i < num_queries; i++) {
			uint64_t key = query_set[i];
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

	return 0;
}
