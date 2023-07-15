#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <iostream>
#include <vector>
#include <openssl/rand.h>

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

int double_cmp(const void *a, const void *b) {
        double x = *((double*)a);
        double y = *((double*)b);
        if (x == y) return 0;
        if (x < y) return -1;
        return 1;
}

int main(int argc, char **argv) {
	uint64_t seed;
	size_t num_queries;
	size_t nslots;
	size_t num_trials;
	if (argc < 4) {
		printf("Usage: ./test_ext_throughput [log of nslots] [num queries per trial] [num trials]\n");
		exit(1);
	}
	else if (argc >= 4) {
		nslots = 1ull << atoi(argv[1]);
		num_queries = strtoull(argv[2], NULL, 10);
		num_trials = atoi(argv[3]);
	}
	if (argc >= 5) {
		seed = strtoull(argv[4], NULL, 10);
		printf("Warning: seeding may not necessarily work due to openssl's own generator\n");
	}
	else {
		seed = time(NULL);
	}
	size_t num_inserts = nslots * 0.9;

	printf("running on seed %lu\n", seed);
	srand(seed);

	unsigned int murmur_seed = rand();

	// Create a cuckoo filter where each item is of type size_t and
	// use 12 bits for each item:
	//    CuckooFilter<size_t, 12> filter(total_items);
	// To enable semi-sorting, define the storage of cuckoo filter to be
	// PackedTable, accepting keys of size_t type and making 13 bits
	// for each key:
	//   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
	printf("initializing data structures...\n");
	CuckooFilter<size_t, 12> filter(num_inserts);

	data_config data_cfg;
        default_data_config_init(MAX_KEY_SIZE, &data_cfg);
        splinterdb_config splinterdb_cfg = (splinterdb_config){
                .filename   = "db",
                .cache_size = 64 * Mega,
                .disk_size  = 20 * Giga,
                .data_cfg   = &data_cfg
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

	// Insert items to this cuckoo filter
	double measure_interval = 0.01f;
	double current_interval = measure_interval;
	uint64_t measure_point = nslots * current_interval;

	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;

	size_t i;
	for (i = 0; i <= num_inserts; i++) {
		if (filter.AddUsingBackingMap(inserts[i], database, MAX_KEY_SIZE, MAX_VAL_SIZE, &db_result, buffer) != cuckoofilter::Ok) {
			printf("insert failed at fill rate %f\n", (double)i / nslots);
			break;
		}

		//fprintf(stderr, "\r%lu/%lu", i, num_inserts);
		//if (i == 5652) bp();

		if (i >= measure_point) {
			fprintf(stderr, "\r%d%%", (int)(100 * current_interval));
			current_interval += measure_interval;
			measure_point = nslots * current_interval;
		}
	}

	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

	printf("\n");
	printf("made %lu inserts\n", i);
	printf("time for inserts:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughput:    %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
	
	// Check non-existing items, a few false positives expected
	size_t fp_queries = 0;

	printf("generating zipfian queries...\n");
	uint64_t *queries = new uint64_t[num_queries];
	RAND_bytes((unsigned char*)queries, num_queries * sizeof(uint64_t));
	for (i = 0; i < num_queries; i++) {
		queries[i] = (uint64_t)rand_zipfian(1.5f, 10000000ull, queries[i], -1ull);
	}

	printf("performing queries...\n");

        FILE *data_fp = fopen("stats_splinter_zipfian_trials.csv", "w");
        fprintf(data_fp, "trial avg_through min_through, med_through, max_through fprate\n");

        uint64_t num_trial_segments = 10;
        measure_interval = 1. / num_trial_segments;

        double *trials_avg = new double[num_trials];
        double *trials_min = new double[num_trials];
        double *trials_max = new double[num_trials];
        double *trials_med = new double[num_trials];
        double *trial_segments = new double[num_trial_segments];

	fprintf(stderr, "trials completed:\t0/%lu", num_trials);
	for (size_t trial = 0; trial < num_trials; trial++) {
		murmur_seed = rand();
		for (i = 0; i < num_queries; i++) {
			queries[i] = MurmurHash64A(&queries[i], sizeof(queries[i]), murmur_seed);
		}

                uint64_t trial_segment = 0, last_point = 0;
                measure_point = num_queries / num_trial_segments - 1;

                fp_queries = 0;

		gettimeofday(&timecheck, NULL);
		uint64_t interval_time = start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

		for (i = 0; i < num_queries; i++) {
			uint64_t location_data;
			if (filter.ContainReturn(queries[i], &location_data) == cuckoofilter::Ok) {
				slice query_slice = padded_slice(&queries[i], MAX_KEY_SIZE, sizeof(queries[i]), buffer, 0);
				splinterdb_lookup(database, query_slice, &db_result);

				if (!splinterdb_lookup_found(&db_result)) {
					fp_queries++;
					slice location_slice = padded_slice(&location_data, MAX_KEY_SIZE, sizeof(location_data), buffer, 0);
					splinterdb_lookup(database, location_slice, &db_result);

					slice result_val;
					splinterdb_lookup_result_value(&db_result, &result_val);
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val), sizeof(orig_key));

					if (filter.Adapt(orig_key, database, MAX_KEY_SIZE, MAX_VAL_SIZE, &db_result, buffer) != cuckoofilter::Ok) printf("error: adapt failed to find previously queried item\n");
				}
			}

                        if (i >= measure_point) {
                                gettimeofday(&timecheck, NULL);
                                trial_segments[trial_segment++] = (double)(i - last_point) * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time);

                                measure_point = ((trial_segment + 1) * num_queries / num_trial_segments) - 1;
                                last_point = i;
                                gettimeofday(&timecheck, NULL);
                                interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
                        }
		}

		gettimeofday(&timecheck, NULL);
		end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

                trials_avg[trial] = (double)(num_queries) * 1000000 / (end_time - start_time);

                qsort(trial_segments, num_trial_segments, sizeof(double), double_cmp);
                trials_min[trial] = trial_segments[0];
                trials_max[trial] = trial_segments[num_trial_segments - 1];
                trials_med[trial] = trial_segments[num_trial_segments / 2];

                fprintf(data_fp, "%lu %f %f %f %f %f\n", trial, trials_avg[trial], trials_min[trial], trials_med[trial], trials_max[trial], (double)fp_queries / num_queries);

		fprintf(stderr, "\rtrials completed:\t%lu/%lu", trial + 1, num_trials);
	}
	end_time = clock();

	printf("\n");
	printf("sorting results...\n");
	fclose(data_fp);

	printf("writing data to file...\n");

        qsort(trials_avg, num_trials, sizeof(double), double_cmp);
        FILE *quartiles_avg_fp = fopen("stats_splinter_zipfian_quartiles_avg.csv", "w");
        fprintf(quartiles_avg_fp, "filter whisker_bottom box_bottom median box_top whisker_top\n");
        fprintf(quartiles_avg_fp, "ACF %f %f %f %f %f\n", trials_avg[0], trials_avg[num_trials / 4], trials_avg[num_trials / 2], trials_avg[num_trials * 3 / 4], trials_avg[num_trials - 1]);
        fclose(quartiles_avg_fp);

        qsort(trials_min, num_trials, sizeof(double), double_cmp);
        FILE *quartiles_min_fp = fopen("stats_splinter_zipfian_quartiles_min.csv", "w");
        fprintf(quartiles_min_fp, "filter whisker_bottom box_bottom median box_top whisker_top\n");
        fprintf(quartiles_min_fp, "ACF %f %f %f %f %f\n", trials_min[0], trials_min[num_trials / 4], trials_min[num_trials / 2], trials_min[num_trials * 3 / 4], trials_min[num_trials - 1]);
        fclose(quartiles_min_fp);

        qsort(trials_med, num_trials, sizeof(double), double_cmp);
        FILE *quartiles_med_fp = fopen("stats_splinter_zipfian_quartiles_med.csv", "w");
        fprintf(quartiles_med_fp, "filter whisker_bottom box_bottom median box_top whisker_top\n");
        fprintf(quartiles_med_fp, "ACF %f %f %f %f %f\n", trials_med[0], trials_med[num_trials / 4], trials_med[num_trials / 2], trials_med[num_trials * 3 / 4], trials_med[num_trials - 1]);
        fclose(quartiles_med_fp);

        qsort(trials_max, num_trials, sizeof(double), double_cmp);
        FILE *quartiles_max_fp = fopen("stats_splinter_zipfian_quartiles_max.csv", "w");
        fprintf(quartiles_max_fp, "filter whisker_bottom box_bottom median box_top whisker_top\n");
        fprintf(quartiles_max_fp, "ACF %f %f %f %f %f\n", trials_max[0], trials_max[num_trials / 4], trials_max[num_trials / 2], trials_max[num_trials * 3 / 4], trials_max[num_trials - 1]);
        fclose(quartiles_max_fp);

	printf("done\n");

	return 0;
}
