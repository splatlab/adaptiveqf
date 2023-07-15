#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <fcntl.h>

#include <iostream>
#include <vector>
#include <openssl/rand.h>

extern "C" {
#include <splinterdb/splinterdb.h>
#include <splinterdb/data.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/default_data_config.h>
}

using cuckoofilter::CuckooFilter;

uint64_t rand_uniform(uint64_t max) {
        if (max <= RAND_MAX) return rand() % max;
        uint64_t a = rand();
        uint64_t b = rand();
        a |= (b << 31);
        return a % max;
}

double rand_zipfian(double s, double max, uint64_t source, uint64_t rand_max) {
        double p = (double)source / rand_max;

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

/*uint64_t MurmurHash64A(const void* key, int len, uint64_t seed) {
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
}*/

void bp() {
	return;
}

int main(int argc, char **argv) {
	uint64_t seed;
	size_t num_queries;
	size_t nslots;
	if (argc < 3) {
		printf("Usage: ./test_ext_throughput [log of nslots] [num queries]\n");
		exit(1);
	}
	else if (argc >= 3) {
		nslots = 1ull << atoi(argv[1]);
		num_queries = strtoull(argv[2], NULL, 10);
	}
	if (argc >= 4) {
		seed = strtoull(argv[3], NULL, 10);
		printf("Warning: seeding may not necessarily work due to openssl's own generator\n");
	}
	else {
		seed = time(NULL);
	}
	size_t num_inserts = nslots * 0.9;

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
	CuckooFilter<size_t, 12> filter(num_inserts);

	uint64_t set_len = num_inserts * 1.5f;
	set_node *set = new set_node[set_len];

	printf("generating inserts...\n");
	uint64_t *inserts = new uint64_t[num_inserts];
	/*for (size_t i = 0; i < num_inserts; i++) {
		inserts[i] = rand_uniform(-1);
	}*/
	RAND_bytes((unsigned char*)(inserts), num_inserts * sizeof(uint64_t));

	printf("starting inserts...\n");

	// Insert items to this cuckoo filter
	double measure_interval = 0.01f;
	double current_interval = measure_interval;
	uint64_t measure_point = nslots * current_interval;

	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;

	size_t i;
	for (i = 0; i <= num_inserts; i++) {
		if (filter.AddUsingSet(inserts[i], set, set_len) != cuckoofilter::Ok) {
			printf("insert failed at fill rate %f\n", (double)i / nslots);
			break;
		}

		if (i >= measure_point) {
			fprintf(stderr, "\r%d%%", (int)(100 * current_interval));
			current_interval += measure_interval;
			measure_point = nslots * current_interval;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("\n");
	printf("made %lu inserts\n", i);
	printf("time for inserts:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughput:    %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
	printf("cpu time for inserts: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	printf("cpu insert through:   %f ops/sec\n", (double)i * CLOCKS_PER_SEC / (end_clock - start_clock));
	
	size_t fp_queries = 0;

	printf("generating queries...\n");
	uint64_t *queries = new uint64_t[num_queries];
	RAND_bytes((unsigned char*)queries, num_queries * sizeof(uint64_t));
	unsigned int murmur_seed = rand();
	for (i = 0; i < num_queries; i++) {
		queries[i] = (uint64_t)rand_zipfian(1.5f, 10000000ull, queries[i], -1ull);
		//queries[i] = queries[i] % (1ull << 24);
		queries[i] = MurmurHash64A(&queries[i], sizeof(queries[i]), murmur_seed);
	}

	printf("performing queries...\n");
	current_interval = measure_interval;
	measure_point = num_queries * current_interval;

	start_clock = clock();
	gettimeofday(&timecheck, NULL);
	start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	for (i = 0; i < num_queries; i++) {
		uint64_t location_data, orig_key;
		if (filter.ContainReturn(queries[i], &location_data) == cuckoofilter::Ok) {
			set_query(set, set_len, location_data, &orig_key);

			if (queries[i] != orig_key) {
				fp_queries++;

				if (filter.AdaptUsingSet(orig_key, set, set_len) != cuckoofilter::Ok) printf("error: adapt failed to find previously queried item\n");
			}
		}

		if (i >= measure_point) {
			fprintf(stderr, "\r%d%%", (int)(current_interval * 100));

			current_interval += measure_interval;
			measure_point = num_queries * current_interval;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("\n");
	printf("made %lu queries\n", num_queries);
	printf("time for queries:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
	printf("cpu time for queries: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	printf("cpu query throughput: %f ops/sec\n", (double)num_queries * CLOCKS_PER_SEC / (end_clock - start_clock));

	printf("false positive rate:  %f%%\n", 100. * fp_queries / num_queries);

	return 0;
}
