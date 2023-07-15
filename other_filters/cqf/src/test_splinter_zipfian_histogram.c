/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

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

#include <splinterdb/splinterdb.h>
#include <splinterdb/data.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/default_data_config.h>


#define Kilo (1024UL)
#define Mega (1024UL * Kilo)
#define Giga (1024UL * Mega)

#define TEST_DB_NAME "db"

#define MAX_KEY_SIZE 16
#define MAX_VAL_SIZE 16

void bp2() {
	return;
}

uint64_t rand_uniform(uint64_t max) {
	if (max <= RAND_MAX) return rand() % max;
	uint64_t a = rand();
	uint64_t b = rand();
	a |= (b << 31);
	return a % max;
}

double rand_normal(double mean, double sd) {
	double a = (double)rand() / RAND_MAX;
	double b = (double)rand() / RAND_MAX;
	double c = sqrt(-2.0 * log(a)) * cos(2 * M_PI * b) * sd;
	return c + mean;
}

double rand_zipfian(double s, double max, uint64_t source) {
	double p = (double)source / (-1ULL);
	
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

uint64_t map_inserts = 0;
uint64_t map_queries = 0;

void pad_data(void *dest, const void *src, const size_t dest_len, const size_t src_len, const int flagged) {
	assert(dest_len >= src_len);
	if (flagged) memset(dest, 0xf, dest_len);
	else bzero(dest, dest_len);
	memcpy(dest, src, src_len);
}

slice padded_slice(const void *data, const size_t dest_len, const size_t src_len, void *buffer, const int flagged) {
	pad_data(buffer, data, dest_len, src_len, flagged);
	
	return slice_create(dest_len, buffer);
}

int db_insert(splinterdb *database, const void *key_data, const size_t key_len, const void *val_data, const size_t val_len, const int flagged) {
	char key_padded[MAX_KEY_SIZE];
	char val_padded[MAX_VAL_SIZE];
	pad_data(key_padded, key_data, MAX_KEY_SIZE, key_len, flagged);
	pad_data(val_padded, val_data, MAX_VAL_SIZE, val_len, 0);
	slice key = slice_create(MAX_KEY_SIZE, key_padded);
	slice val = slice_create(MAX_VAL_SIZE, val_padded);

	return splinterdb_insert(database, key, val);
}

int double_cmp(const void *a, const void *b) {
        double x = *((double*)a);
        double y = *((double*)b);
        if (x == y) return 0;
        if (x < y) return -1;
        return 1;
}

splinterdb_lookup_result db_result;
int main(int argc, char **argv)
{
	if (argc < 5) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries per trial [eg. 100000000]\nthe number of trials [eg. 100]\n");
		// ./test 16 7 $((1 << 15)) 1000000 disk_throughput_stats.txt 0
		exit(1);
	}
	if (argc >= 6) {
		srand(strtol(argv[5], NULL, 10));
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
	uint64_t num_trials = strtoull(argv[4], NULL, 10);

	unsigned int murmur_seed = rand();


	printf("initializing hash table...\n");

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
	splinterdb_lookup_result_init(database, &db_result, 0, NULL);

	printf("initializing filter...\n");
	QF qf;
	if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
		fprintf(stderr, "Can't allocate CQF.\n");
		abort();
	}
	qf_set_auto_resize(&qf, false);


	printf("generating insert set of size %lu...\n", num_inserts);
	uint64_t i;
	uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
        RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));


	// PERFORM INSERTS
	uint64_t target_fill = nslots * load_factor;
	printf("performing insertions... 0/%lu", target_fill);
	//std::cout << "performing insertions... 0%%" << std::flush;
	char buffer[MAX_KEY_SIZE + MAX_VAL_SIZE];

	double measure_interval = 0.01f;
        double current_interval = measure_interval;
        uint64_t measure_point = target_fill * current_interval;

	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	clock_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;
	for (i = 0; qf.metadata->noccupied_slots < target_fill; i++) {
		int ret = qf_insert(&qf, insert_set[i], 0, 1, QF_NO_LOCK);
		if (ret < 0) abort();

		db_insert(database, &insert_set[i], sizeof(insert_set[i]), &i, sizeof(i), 0);

		if (qf.metadata->noccupied_slots >= measure_point) {
			fprintf(stderr, "\rperforming insertions... %f%%          ", current_interval * 100);

                        current_interval += measure_interval;
                        measure_point = nslots * current_interval;
                }
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

	printf("\n");
	printf("time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
	printf("avg insert throughput: %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));

	// PERFORM QUERIES
	printf("\ngenerating query set of size %lu...\n", num_queries);

	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));
	for (i = 0; i < num_queries; i++) {
		query_set[i] = (uint64_t)rand_zipfian(1.5f, 10000000ull, query_set[i]);
	}

	printf("performing queries...\n");

        FILE *data_fp = fopen("stats_splinter_zipfian_trials.csv", "w");
        fprintf(data_fp, "trial avg_through min_through med_through max_through fprate\n");

        uint64_t num_trial_segments = 10;
        measure_interval = 1. / num_trial_segments;

        double *trials_avg = malloc(num_trials * sizeof(double));
        double *trials_min = malloc(num_trials * sizeof(double));
        double *trials_med = malloc(num_trials * sizeof(double));
        double *trials_max = malloc(num_trials * sizeof(double));
        double *trial_segments = malloc(num_trial_segments * sizeof(double));

        fprintf(stderr, "trials completed:\t0/%lu", num_trials);

	uint64_t fp_count = 0;
        for (size_t trial = 0; trial < num_trials; trial++) {
                murmur_seed = rand();
                for (i = 0; i < num_queries; i++) {
                        query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
                }

                uint64_t trial_segment = 0, last_point = 0;
                measure_point = num_queries / num_trial_segments - 1;

		fp_count = 0;

		gettimeofday(&timecheck, NULL);
		uint64_t interval_time = start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

		for (i = 0; i < num_queries; i++) {
			uint64_t value;
			if (qf_query(&qf, query_set[i], &value, QF_NO_LOCK)) {
				slice query = padded_slice(&query_set[i], MAX_KEY_SIZE, sizeof(query_set[i]), buffer, 0);
				splinterdb_lookup(database, query, &db_result);
				if (!splinterdb_lookup_found(&db_result)) {
					fp_count++;
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
                trials_med[trial] = trial_segments[num_trial_segments / 2];
                trials_max[trial] = trial_segments[num_trial_segments - 1];

                fprintf(data_fp, "%lu %f %f %f %f %f\n", trial, trials_avg[trial], trials_min[trial], trials_med[trial], trials_max[trial], (double)fp_count / num_queries);

                fprintf(stderr, "\rtrials completed:\t%lu/%lu", trial + 1, num_trials);
	}

        printf("\n");
        printf("sorting results...\n");
        fclose(data_fp);

        printf("writing data to file...\n");

        qsort(trials_avg, num_trials, sizeof(double), double_cmp);
        FILE *quartiles_avg_fp = fopen("stats_splinter_zipfian_quartiles_avg.csv", "w");
        fprintf(quartiles_avg_fp, "filter whisker_bottom box_bottom median box_top whisker_top\n");
        fprintf(quartiles_avg_fp, "QF %f %f %f %f %f\n", trials_avg[0], trials_avg[num_trials / 4], trials_avg[num_trials / 2], trials_avg[num_trials * 3 / 4], trials_avg[num_trials - 1]);
        fclose(quartiles_avg_fp);

        qsort(trials_min, num_trials, sizeof(double), double_cmp);
        FILE *quartiles_min_fp = fopen("stats_splinter_zipfian_quartiles_min.csv", "w");
        fprintf(quartiles_min_fp, "filter whisker_bottom box_bottom median box_top whisker_top\n");
        fprintf(quartiles_min_fp, "QF %f %f %f %f %f\n", trials_min[0], trials_min[num_trials / 4], trials_min[num_trials / 2], trials_min[num_trials * 3 / 4], trials_min[num_trials - 1]);
        fclose(quartiles_min_fp);

        qsort(trials_med, num_trials, sizeof(double), double_cmp);
        FILE *quartiles_med_fp = fopen("stats_splinter_zipfian_quartiles_med.csv", "w");
        fprintf(quartiles_med_fp, "filter whisker_bottom box_bottom median box_top whisker_top\n");
        fprintf(quartiles_med_fp, "QF %f %f %f %f %f\n", trials_med[0], trials_med[num_trials / 4], trials_med[num_trials / 2], trials_med[num_trials * 3 / 4], trials_med[num_trials - 1]);
        fclose(quartiles_med_fp);

        qsort(trials_max, num_trials, sizeof(double), double_cmp);
        FILE *quartiles_max_fp = fopen("stats_splinter_zipfian_quartiles_max.csv", "w");
        fprintf(quartiles_max_fp, "filter whisker_bottom box_bottom median box_top whisker_top\n");
        fprintf(quartiles_max_fp, "QF %f %f %f %f %f\n", trials_max[0], trials_max[num_trials / 4], trials_max[num_trials / 2], trials_max[num_trials * 3 / 4], trials_max[num_trials - 1]);
        fclose(quartiles_max_fp);

        printf("done\n");
}

