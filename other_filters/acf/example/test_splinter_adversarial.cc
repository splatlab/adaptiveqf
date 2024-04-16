#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <vector>
#include <fcntl.h>
#include <openssl/rand.h>
#include <sys/time.h>

extern "C" {
#include <splinterdb/splinterdb.h>
#include <splinterdb/data.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/default_data_config.h>
}

#define Kilo (1024UL)
#define Mega (1024UL * Kilo)
#define Giga (1024UL * Mega)

#define TEST_DB_NAME "db"

#define MAX_KEY_SIZE 16
#define MAX_VAL_SIZE 16

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

void bp() {
	return;
}

void pad_data(void *dest, const void *src, const size_t dest_len, const size_t src_len, const int flagged) {
	assert(dest_len >= src_len);
	if (flagged) memset(dest, 0xff, dest_len);
	else bzero(dest, dest_len);
	memcpy(dest, src, src_len);
}

slice padded_slice(const void *data, const size_t dest_len, const size_t src_len, void *buffer, const int flagged) {
	pad_data(buffer, data, dest_len, src_len, flagged);

	return slice_create(dest_len, buffer);
}

int db_insert(splinterdb *database, const void *key_data, const size_t key_len, unsigned char *key_buffer, const void *val_data, const size_t val_len, unsigned char *val_buffer, const int flagged) {
	pad_data(key_buffer, key_data, MAX_KEY_SIZE, key_len, flagged);
	pad_data(val_buffer, val_data, MAX_VAL_SIZE, val_len, 0);
	slice key = slice_create(MAX_KEY_SIZE, key_buffer);
	slice val = slice_create(MAX_VAL_SIZE, val_buffer);

	return splinterdb_insert(database, key, val);
}

int mycmp(const void *a, const void *b) {
        return memcmp(a, b, sizeof(uint64_t));
}

int main(int argc, char **argv) {
	uint64_t seed;
	size_t num_queries;
	size_t nslots;
	if (argc < 4) {
		printf("Usage: ./test_ext_throughput [log of nslots] [num queries] [adv freq]\n");
		exit(1);
	}
	else {
		nslots = 1ull << atoi(argv[1]);
		num_queries = strtoull(argv[2], NULL, 10);
	}
	seed = time(NULL);
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

	size_t cache_size_in_mb = 256;
	data_config data_cfg;
        default_data_config_init(MAX_KEY_SIZE, &data_cfg);
        splinterdb_config splinterdb_cfg = (splinterdb_config){
                .filename   = "db",
                .cache_size = cache_size_in_mb * Mega,
                .disk_size  = 20 * Giga,
                .data_cfg   = &data_cfg,
		.io_flags   = O_RDWR | O_CREAT | O_DIRECT
        };
        splinterdb *database;
        if (splinterdb_create(&splinterdb_cfg, &database) != 0) {
                printf("Error creating database\n");
                exit(0);
        }
        splinterdb_lookup_result db_result;
        splinterdb_lookup_result_init(database, &db_result, 0, NULL);

	printf("generating inserts...\n");
	uint64_t *inserts = new uint64_t[num_inserts];
	/*for (size_t i = 0; i < num_inserts; i++) {
		inserts[i] = rand_uniform(-1);
	}*/
	RAND_bytes((unsigned char*)(inserts), num_inserts * sizeof(uint64_t));

	printf("starting inserts...\n");

	unsigned char buffer[MAX_KEY_SIZE + MAX_VAL_SIZE];

	FILE *inserts_fp = fopen("stats_adv_inserts.csv", "w");
	fprintf(inserts_fp, "fill through\n");


	// Insert items to this cuckoo filter
	double measure_interval = 0.01f;
	double current_interval = measure_interval;
	uint64_t measure_point = nslots * current_interval, last_point = 0;

	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time, interval_time = start_time;

	size_t i;
	for (i = 0; i <= num_inserts; i++) {
		if (filter.AddUsingBackingMap(inserts[i], database, MAX_KEY_SIZE, MAX_VAL_SIZE, &db_result, buffer) != cuckoofilter::Ok) {
			printf("insert failed at fill rate %f\n", (double)i / nslots);
			break;
		}

		if (i >= measure_point) {
			gettimeofday(&timecheck, NULL);
			fprintf(inserts_fp, "%f %f\n", current_interval * 100, (double)(i - last_point) * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time));

			fprintf(stderr, "\r%d%%", (int)(100 * current_interval));
			current_interval += measure_interval;
			last_point = measure_point;
			measure_point = nslots * current_interval;

			gettimeofday(&timecheck, NULL);
			interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("\n");
	printf("made %lu inserts\n", i);
	printf("time for inserts:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughptu:    %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
	printf("cpu time for inserts: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	uint64_t *queries = new uint64_t[num_queries];
	
	uint64_t max_adv_set_len = 1000000;
	uint64_t *adv_set = new uint64_t[max_adv_set_len];

	for (int trial = 3; trial < argc; trial++) {
		uint64_t adv_freq = strtoull(argv[trial], NULL, 10);

		// Check non-existing items, a few false positives expected
		size_t fp_queries = 0;

		printf("generating queries...\n");
		RAND_bytes((unsigned char*)queries, num_queries * sizeof(uint64_t));

		printf("performing queries...\n");

		char filename[100];
		sprintf(filename, "%d-12-%lu-%lu-adversarial.csv", atoi(argv[1]), cache_size_in_mb, adv_freq);
		FILE *adv_fp = fopen(filename, "w");
		fprintf(adv_fp, "queries through fprate\n");

		current_interval = measure_interval;
		measure_point = num_queries * current_interval;
		uint64_t last_point = 0;

		uint64_t adv_index = 0;
		uint64_t adv_set_len = 0;

		uint64_t adv_attempted = 0, adv_successful = 0, adv_failed = 0, adv_found = 0;

		start_clock = clock();
		gettimeofday(&timecheck, NULL);
		uint64_t interval_time = start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		for (i = 0; i < num_queries; i++) {
			uint64_t location_data;
			if (i % adv_freq == 0 && adv_set_len > 0) {
				adv_attempted++;
				if (adv_index >= adv_set_len) adv_index = 0;
				uint64_t key = adv_set[adv_index];
				if (filter.ContainReturn(key, &location_data) == cuckoofilter::Ok) {
					slice query_slice = padded_slice(&location_data, MAX_KEY_SIZE, sizeof(location_data), buffer, 0);
					splinterdb_lookup(database, query_slice, &db_result);

					slice result_val;
					splinterdb_lookup_result_value(&db_result, &result_val);
					if (memcmp(&key, slice_data(result_val), sizeof(uint64_t)) != 0) {
						fp_queries++;
						adv_successful++;

						uint64_t orig_key;
						memcpy(&orig_key, slice_data(result_val), sizeof(orig_key));
						filter.Adapt(orig_key, database, MAX_KEY_SIZE, MAX_VAL_SIZE, &db_result, buffer);
					}
					else {
						adv_set[adv_index] = adv_set[--adv_set_len];
						adv_failed++;
					}
				}
				else {
					adv_set[adv_index] = adv_set[--adv_set_len];
					adv_failed++;
				}
				adv_index++;
			}
			else {
				if (filter.ContainReturn(queries[i], &location_data) == cuckoofilter::Ok) {
					slice query_slice = padded_slice(&location_data, MAX_KEY_SIZE, sizeof(location_data), buffer, 0);
					splinterdb_lookup(database, query_slice, &db_result);

					slice result_val;
					splinterdb_lookup_result_value(&db_result, &result_val);
					if (memcmp(&queries[i], slice_data(result_val), sizeof(uint64_t)) != 0) {
						fp_queries++;
						if (adv_set_len < max_adv_set_len) {
							adv_set[adv_set_len++] = queries[i];
							adv_found++;
						}

						uint64_t orig_key;
						memcpy(&orig_key, slice_data(result_val), sizeof(orig_key));
						filter.Adapt(orig_key, database, MAX_KEY_SIZE, MAX_VAL_SIZE, &db_result, buffer);
					}
				}
			}

			if (i >= measure_point) {
				gettimeofday(&timecheck, NULL);
				fprintf(adv_fp, "%lu %f %f\n", i, (double)(i - last_point) * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time), (double)fp_queries / i);
				fprintf(stderr, "\r%d%%", (int)(current_interval * 100));

				current_interval += measure_interval;
				measure_point = num_queries * current_interval;
				last_point = i;

				gettimeofday(&timecheck, NULL);
				interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
			}
		}
		gettimeofday(&timecheck, NULL);
		end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		end_clock = clock();

		fclose(adv_fp);

		printf("\n");
		printf("made %lu queries\n", num_queries);
		printf("time for queries:     %f sec\n", (double)(end_time - start_time) / 1000000);
		printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
		printf("cpu time for queries: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

		printf("adv attempted:        %lu\n", adv_attempted);
		printf("adv successful:       %lu\n", adv_successful);
		printf("adv failed:           %lu\n", adv_failed);
		printf("adv found:            %lu\n", adv_found);

		printf("false positives:      %lu\n", fp_queries);
		printf("false positive rate:  %f%%\n", 100. * fp_queries / num_queries);
	}

	return 0;
}
