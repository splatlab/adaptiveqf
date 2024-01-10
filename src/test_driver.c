#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <openssl/rand.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"
#include "include/hashutil.h"
#include "include/rand_util.h"
#include "include/splinter_util.h"

#include <splinterdb/data.h>
#include <splinterdb/default_data_config.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/public_util.h>
#include <splinterdb/splinterdb.h>

#include "include/test_driver.h"


int run_throughput_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, int verbose, char *inserts_outfile, char *queries_outfile) {
	size_t num_slots = 1ull << qbits;
	size_t minirun_id_bitmask = (1ull << (qbits + rbits)) - 1;

	data_config data_cfg = qf_data_config_init();
	splinterdb_config splinterdb_cfg = qf_splinterdb_config_init("db", &data_cfg);
	splinterdb *db;
	if (splinterdb_create(&splinterdb_cfg, &db)) {
		return -1;
	}
	splinterdb_lookup_result db_result;
	splinterdb_lookup_result_init(db, &db_result, 0, NULL);

	QF qf;
	if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
		return -1;
	}


	double target_load = 0.9f;
	size_t max_inserts = num_slots * target_load;
	size_t num_inserts = insert_set_len > max_inserts ? max_inserts : insert_set_len;
	size_t i;

	size_t measure_freq = 100, curr_interval = 0;
	size_t measure_point = num_inserts * (curr_interval + 1) / measure_freq, prev_point = 0;


	FILE *inserts_file = inserts_outfile ? fopen(inserts_outfile, "w") : NULL;
	fprintf(inserts_file, "fill through\n");

	clock_t start_clock = clock(), end_clock;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t start_time = tv.tv_sec * 1000000 + tv.tv_usec, end_time, interval_time = start_time;

	for (i = 0; qf.metadata->noccupied_slots < num_inserts; i++) {
		if (!qf_splinter_insert(&qf, db, insert_set[i], 1)) break;

		if (qf.metadata->noccupied_slots >= measure_point) {
			gettimeofday(&tv, NULL);
			if (inserts_file) fprintf(inserts_file, "%.2f %f\n", (double)qf.metadata->noccupied_slots / num_slots * 100, (double)(i - prev_point) * 1000000 / (tv.tv_sec * 1000000 + tv.tv_usec - interval_time));
			if (verbose) fprintf(stderr, "\rPerforming insertions... %.2f%%", (double)curr_interval / measure_freq * 100);

			curr_interval++;
			prev_point = i;
			measure_point = num_inserts * (curr_interval + 1) / measure_freq;

			gettimeofday(&tv, NULL);
			interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
		}
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;
	end_clock = clock();

	if (inserts_file) fprintf(inserts_file, "%.2f %f\n", (double)qf.metadata->noccupied_slots / num_slots * 100, (double)(i - prev_point) * 1000000 / (end_time - interval_time));
	if (verbose) fprintf(stderr, "\rPerforming insertions... 100.00%%\n");

	fclose(inserts_file);

	if (verbose) {
		printf("Time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
		printf("Insert throughput:     %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
		printf("CPU time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	}


	size_t num_warmup_queries = 0; // should make this a function input
	for (i = 0; i < num_warmup_queries; i++) {
		// Do warmup queries here
	}

	size_t num_queries = query_set_len > num_warmup_queries ? query_set_len - num_warmup_queries : 0;
	uint64_t *warm_queries = &(query_set[num_warmup_queries]);

	curr_interval = 0;
	measure_point = num_queries * (curr_interval + 1) / measure_freq;
	prev_point = 0;

	FILE *queries_file = queries_outfile ? fopen(queries_outfile, "w") : NULL;
	fprintf(queries_file, "queries through fprate\n");

	int still_have_space = 1;
	size_t full_point = num_slots * 0.95f;
	qf_query_result query_result;
	char buffer[10 * MAX_VAL_SIZE];
	uint64_t fp_count = 0;

	start_clock = clock();
	gettimeofday(&tv, NULL);
	start_time = interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
	for (i = 0; i < num_queries; i++) {
		if (qf_query_using_ll_table(&qf, warm_queries[i], &query_result, QF_KEY_IS_HASH)) {
			uint64_t temp = query_result.hash & minirun_id_bitmask;
			slice query = padded_slice(&temp, MAX_KEY_SIZE, sizeof(temp), buffer, 0);
			splinterdb_lookup(db, query, &db_result);
			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);

			if (memcmp(&warm_queries[i], slice_data(result_val) + query_result.minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t)) != 0) {
				fp_count++;
				if (still_have_space) {
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val) + query_result.minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t));
					qf_adapt_using_ll_table(&qf, orig_key, warm_queries[i], query_result.minirun_rank, QF_KEY_IS_HASH);
					if (qf.metadata->noccupied_slots >= full_point) {
						still_have_space = 0;
						if (verbose) fprintf(stderr, "\rFilter is full after %lu queries\n", i);
					}
				}
			}
		}

		if (i >= measure_point) {
			gettimeofday(&tv, NULL);

			if (queries_file) fprintf(queries_file, "%lu %f %f\n", i, (double)(i - prev_point) * 1000000 / (tv.tv_sec * 1000000 + tv.tv_usec - interval_time), (double)fp_count / i);
			if (verbose) fprintf(stderr, "\rPerforming queries... %f%%", (double)curr_interval / measure_freq * 100);

			curr_interval++;
			prev_point = i;
			measure_point = num_queries * (curr_interval + 1) / measure_freq;

			gettimeofday(&tv, NULL);
			interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
		}
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;
	end_clock = clock();

	if (queries_file) fprintf(queries_file, "%lu %f %f\n", i, (double)(i - prev_point) * 1000000 / (end_time - interval_time), (double)fp_count / i);
	if (verbose) fprintf(stderr, "\rPerforming queries... 100.00%%\n");

	fclose(queries_file);

	if (verbose) {
		printf("Time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
		printf("Query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
		printf("CPU time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

		printf("False positives:      %lu\n", fp_count);
		printf("False positive rate:  %f%%\n", 100. * fp_count / num_queries);
	}
	
	splinterdb_close(&db);
	qf_free(&qf);
	return 0;
}
