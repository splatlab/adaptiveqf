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
#include <pthread.h>
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


void init_test_results(test_results_t *results) {
	results->exit_code = 0;
	results->insert_throughput = results->query_throughput = results->final_query_throughput = results->false_positive_rate = 0;
}

void warm_up_filter(const QF *qf, size_t num_query_set) {
	uint64_t hash;
	uint64_t *query_set = malloc(num_query_set * sizeof(uint64_t));
	for (int i = 0; i < num_query_set; i++) {
		qf_query_using_ll_table(qf, query_set[i], &hash, QF_KEY_IS_HASH);
	}
	free(query_set);
}

void bp() {
}

test_results_t run_throughput_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, int verbose, char *inserts_outfile, char *queries_outfile) {
	test_results_t results;
	init_test_results(&results);

	size_t num_slots = 1ull << qbits;
	size_t minirun_id_bitmask = (1ull << (qbits + rbits)) - 1;

	data_config data_cfg = qf_data_config_init();
	splinterdb_config splinterdb_cfg = qf_splinterdb_config_init("db", &data_cfg);
	remove(splinterdb_cfg.filename);
	splinterdb *db;
	if (splinterdb_create(&splinterdb_cfg, &db)) {
		results.exit_code = -1;
		return results;
	}
	splinterdb_lookup_result db_result;
	splinterdb_lookup_result_init(db, &db_result, 0, NULL);

	QF qf;
	if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
		results.exit_code = -1;
		return results;
	}


	double target_load = 0.9f;
	size_t max_inserts = num_slots * target_load;
	size_t num_inserts = insert_set_len > max_inserts ? max_inserts : insert_set_len;
	size_t i;

	size_t measure_freq = 100, curr_interval = 0;
	size_t measure_point = num_inserts * (curr_interval + 1) / measure_freq, prev_point = 0;


	FILE *inserts_file = inserts_outfile ? fopen(inserts_outfile, "w") : NULL;
	if (inserts_file) fprintf(inserts_file, "fill through\n");

	if (verbose) fprintf(stderr, "Performing insertions... 0.00%%");
	uint64_t num_updates = 0;
	clock_t start_clock = clock(), end_clock;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t start_time = tv.tv_sec * 1000000 + tv.tv_usec, end_time, interval_time = start_time;
	for (i = 0; qf.metadata->noccupied_slots < num_inserts; i++) {
		int ret = qf_splinter_insert(&qf, db, insert_set[i], 1);
		if (ret == 1) continue;
		if (ret == 0) break;
		num_updates++;

		//if (!qf_splinter_insert(&qf, db, insert_set[i], 1)) break;

		if (qf.metadata->noccupied_slots >= measure_point) {
			gettimeofday(&tv, NULL);
			if (inserts_file) fprintf(inserts_file, "%.2f %f\n", (double)qf.metadata->noccupied_slots / num_slots * 100, (double)(i - prev_point) * 1000000 / (tv.tv_sec * 1000000 + tv.tv_usec - interval_time));
			if (verbose) fprintf(stderr, "\rPerforming insertions... %.2f%%", (double)(curr_interval + 1) / measure_freq * 100);

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

	if (verbose) {
		printf("Number of inserts:     %lu\n", i);
		printf("Number of updates:     %lu\n", num_updates);
		printf("Time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
		printf("Insert throughput:     %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
		printf("CPU time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	}
	results.insert_throughput = (double)i * 1000000 / (end_time - start_time);


	curr_interval = 0;
	measure_point = query_set_len * (curr_interval + 1) / measure_freq;
	prev_point = 0;

	FILE *queries_file = queries_outfile ? fopen(queries_outfile, "w") : NULL;
	if (queries_file) fprintf(queries_file, "queries through fprate\n");

	int still_have_space = 1;
	size_t full_point = num_slots * 0.95f;
	char buffer[10 * MAX_VAL_SIZE];
	uint64_t fp_count = 0;
	uint64_t hash;
	int minirun_rank;

	if (verbose) fprintf(stderr, "Performing queries... 0.00%%");
	start_clock = clock();
	gettimeofday(&tv, NULL);
	start_time = interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
	for (i = 0; i < query_set_len; i++) {
		if ((minirun_rank = qf_query_using_ll_table(&qf, query_set[i], &hash, QF_KEY_IS_HASH)) >= 0) {
			hash = (hash & minirun_id_bitmask) << (64 - qf.metadata->quotient_remainder_bits);
			slice query = padded_slice(&hash, MAX_KEY_SIZE, sizeof(hash), buffer, 0);
			splinterdb_lookup(db, query, &db_result);

			/*slice query = padded_slice(&query_set[i], MAX_KEY_SIZE, sizeof(query_set[i]), buffer, 0);
			splinterdb_lookup(db, query, &db_result);
			if (!splinterdb_lookup_found(&db_result)) {
				fp_count++;
			}*/

			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);

			/*uint64_t orig_key = *((uint64_t*)(slice_data(result_val) + minirun_rank * MAX_VAL_SIZE));
			if (query_set[i] != orig_key) {
				fp_count++;
				if (still_have_space) {
					qf_adapt_using_ll_table(&qf, orig_key, query_set[i], minirun_rank, QF_KEY_IS_HASH);
					if (qf.metadata->noccupied_slots >= full_point) {
						still_have_space = 0;
						if (verbose) fprintf(stderr, "\rFilter is full after %lu queries\n", i);
					}
				}
			}*/

			if (memcmp(&query_set[i], slice_data(result_val) + minirun_rank * MAX_VAL_SIZE, sizeof(uint64_t)) != 0) {
				fp_count++;
				if (still_have_space) {
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val) + minirun_rank * MAX_VAL_SIZE, sizeof(uint64_t));
					qf_adapt_using_ll_table(&qf, orig_key, query_set[i], minirun_rank, QF_KEY_IS_HASH);
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
			if (verbose) fprintf(stderr, "\rPerforming queries... %.2f%%", (double)(curr_interval + 1) / measure_freq * 100);

			curr_interval++;
			prev_point = i;
			measure_point = query_set_len * (curr_interval + 1) / measure_freq;

			gettimeofday(&tv, NULL);
			interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
		}
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;
	end_clock = clock();

	if (queries_file) fprintf(queries_file, "%lu %f %f\n", i, (double)(i - prev_point) * 1000000 / (end_time - interval_time), (double)fp_count / i);
	if (verbose) fprintf(stderr, "\rPerforming queries... 100.00%%\n");

	if (verbose) {
		printf("Time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
		printf("Query throughput:     %f ops/sec\n", (double)query_set_len * 1000000 / (end_time - start_time));
		printf("CPU time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

		printf("False positives:      %lu\n", fp_count);
		printf("False positive rate:  %f%%\n", 100. * fp_count / query_set_len);
	}
	results.query_throughput = (double)i * 1000000 / (end_time - start_time);
	results.false_positive_rate = (double)fp_count / query_set_len;
	
	splinterdb_close(&db);
	qf_free(&qf);
	return results;
}

test_results_t run_adversarial_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, uint64_t cache_size, uint64_t adv_freq, uint64_t adv_set_max_size, int verbose, char *queries_outfile) {
	test_results_t results;
	init_test_results(&results);

	size_t num_slots = 1ull << qbits;
	size_t minirun_id_bitmask = (1ull << (qbits + rbits)) - 1;

	data_config data_cfg = qf_data_config_init();
	splinterdb_config splinterdb_cfg = qf_splinterdb_config_init("db", &data_cfg);
	splinterdb_cfg.cache_size = cache_size * Mega;
	remove(splinterdb_cfg.filename);
	splinterdb *db;
	if (splinterdb_create(&splinterdb_cfg, &db)) {
		results.exit_code = -1;
		return results;
	}
	splinterdb_lookup_result db_result;
	splinterdb_lookup_result_init(db, &db_result, 0, NULL);

	QF qf;
	if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
		results.exit_code = -1;
		return results;
	}


	double target_load = 0.9f;
	size_t max_inserts = num_slots * target_load;
	size_t num_inserts = insert_set_len > max_inserts ? max_inserts : insert_set_len;
	size_t i;

	size_t measure_freq = 100, curr_interval = 0;
	size_t measure_point = num_inserts * (curr_interval + 1) / measure_freq, prev_point = 0;

	
	if (verbose) fprintf(stderr, "Performing insertions... 0.00%%");
	clock_t start_clock = clock(), end_clock;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t start_time = tv.tv_sec * 1000000 + tv.tv_usec, end_time, interval_time = start_time;
	for (i = 0; qf.metadata->noccupied_slots < num_inserts; i++) {
		if (!qf_splinter_insert(&qf, db, insert_set[i], 1)) break;

		if (qf.metadata->noccupied_slots >= measure_point) {
			gettimeofday(&tv, NULL);
			if (verbose) fprintf(stderr, "\rPerforming insertions... %.2f%%", (double)(curr_interval + 1) / measure_freq * 100);

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

	if (verbose) fprintf(stderr, "\rPerforming insertions... 100.00%%\n");

	if (verbose) {
		printf("Number of inserts:     %lu\n", i);
		printf("Time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
		printf("Insert throughput:     %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
		printf("CPU time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	}
	results.insert_throughput = (double)i * 1000000 / (end_time - start_time);


	size_t adv_set_size = 0, curr_adv_query = 0;
	uint64_t *adv_queries = malloc(adv_set_max_size * sizeof(uint64_t));

	double final_throughput = 0;

	curr_interval = 0;
	measure_point = query_set_len * (curr_interval + 1) / measure_freq;
	prev_point = 0;

	FILE *queries_file = queries_outfile ? fopen(queries_outfile, "w") : NULL;
	fprintf(queries_file, "queries through fprate\n");

	int still_have_space = 1;
	size_t full_point = num_slots * 0.95f;
	char buffer[10 * MAX_VAL_SIZE];
	uint64_t hash;
	int minirun_rank;
	uint64_t fp_count = 0;
	uint64_t adv_succ = 0, adv_fail = 0, tp_count = 0;

	if (verbose) fprintf(stderr, "Performing queries... 0.00%%");
	start_clock = clock();
	gettimeofday(&tv, NULL);
	start_time = interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
	for (i = 0; i < query_set_len; i++) {
		if (i % adv_freq == 0 && adv_set_size > 0) {
			if (curr_adv_query >= adv_set_size) curr_adv_query = 0;
			if ((minirun_rank = qf_query_using_ll_table(&qf, adv_queries[curr_adv_query], &hash, QF_KEY_IS_HASH)) >= 0) {
				hash = (hash & minirun_id_bitmask) << (64 - qf.metadata->quotient_remainder_bits);
				slice query = padded_slice(&hash, MAX_KEY_SIZE, sizeof(hash), buffer, 0);
				splinterdb_lookup(db, query, &db_result);
				slice result_val;
				splinterdb_lookup_result_value(&db_result, &result_val);

				if (memcmp(&adv_queries[curr_adv_query], slice_data(result_val) + minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t)) != 0) {
					fp_count++;
					adv_succ++;
					if (still_have_space) {
						uint64_t orig_key;
						memcpy(&orig_key, slice_data(result_val) + minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t));
						qf_adapt_using_ll_table(&qf, orig_key, adv_queries[curr_adv_query], minirun_rank, QF_KEY_IS_HASH);
						if (qf.metadata->noccupied_slots >= full_point) {
							still_have_space = 0;
							if (verbose) fprintf(stderr, "\rFilter is full after %lu queries\n", i);
						}
					}
				}
				else tp_count++;

				curr_adv_query++;
			}
			else {
				adv_fail++;
				adv_queries[curr_adv_query] = adv_queries[--adv_set_size];
			}
		}
		else {
			if ((minirun_rank = qf_query_using_ll_table(&qf, query_set[i], &hash, QF_KEY_IS_HASH)) >= 0) {
				hash = (hash & minirun_id_bitmask) << (64 - qf.metadata->quotient_remainder_bits);
				slice query = padded_slice(&hash, MAX_KEY_SIZE, sizeof(hash), buffer, 0);
				splinterdb_lookup(db, query, &db_result);
				slice result_val;
				splinterdb_lookup_result_value(&db_result, &result_val);

				if (memcmp(&query_set[i], slice_data(result_val) + minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t)) != 0) {
					fp_count++;
					if (still_have_space) {
						uint64_t orig_key;
						memcpy(&orig_key, slice_data(result_val) + minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t));
						qf_adapt_using_ll_table(&qf, orig_key, query_set[i], minirun_rank, QF_KEY_IS_HASH);
						if (qf.metadata->noccupied_slots >= full_point) {
							still_have_space = 0;
							if (verbose) fprintf(stderr, "\rFilter is full after %lu queries\n", i);
						}
					}
				}
				else tp_count++;

				if (adv_set_size < adv_set_max_size) adv_queries[adv_set_size++] = query_set[i];
			}
		}

		if (i >= measure_point) {
			gettimeofday(&tv, NULL);

			final_throughput = (double)(i - prev_point) * 1000000 / (tv.tv_sec * 1000000 + tv.tv_usec - interval_time);
			if (queries_file) fprintf(queries_file, "%lu %f %f\n", i, final_throughput, (double)fp_count / i);
			if (verbose) fprintf(stderr, "\rPerforming queries... %.2f%%", (double)curr_interval / measure_freq * 100);

			curr_interval++;
			prev_point = i;
			measure_point = query_set_len * (curr_interval + 1) / measure_freq;

			gettimeofday(&tv, NULL);
			interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
		}
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;
	end_clock = clock();

	if (queries_file) fprintf(queries_file, "%lu %f %f\n", i, (double)(i - prev_point) * 1000000 / (end_time - interval_time), (double)fp_count / i);
	if (verbose) fprintf(stderr, "\rPerforming queries... 100.00%%\n");

	if (queries_file) fclose(queries_file);

	if (verbose) {
		printf("Time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
		printf("Query throughput:     %f ops/sec\n", (double)query_set_len * 1000000 / (end_time - start_time));
		printf("CPU time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

		printf("False positives:      %lu\n", fp_count);
		printf("False positive rate:  %f%%\n", 100. * fp_count / query_set_len);

		printf("Adversarial attacks:  %lu\n", adv_succ + adv_fail);
		printf("Successful attacks:   %lu\n", adv_succ);
		printf("Failed attacks:       %lu\n", adv_fail);
		printf("True positives:       %lu\n", tp_count);
	}
	results.query_throughput = (double)i * 1000000 / (end_time - start_time);
	results.false_positive_rate = (double)fp_count / query_set_len;
	results.final_query_throughput = final_throughput;
	
	free(adv_queries);
	splinterdb_close(&db);
	qf_free(&qf);
	return results;
}

test_results_t run_micro_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, int verbose) {
	test_results_t results;
	init_test_results(&results);

	size_t num_slots = 1ull << qbits;
	size_t minirun_id_bitmask = (1ull << (qbits + rbits)) - 1;

	QF qf;
	if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
		results.exit_code = -1;
		return results;
	}


	double target_load = 0.9f;
	size_t max_inserts = num_slots * target_load;
	size_t num_inserts = insert_set_len > max_inserts ? max_inserts : insert_set_len;

	qf_insert_result insert_result;
	size_t i;

	if (verbose) printf("Performing insertions...\n");
	clock_t start_clock = clock(), end_clock;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t start_time = tv.tv_sec * 1000000 + tv.tv_usec, end_time, interval_time = start_time;
	for (i = 0; qf.metadata->noccupied_slots < num_inserts; i++) {
		if (qf_insert_using_ll_table(&qf, insert_set[i], 1, &insert_result, QF_NO_LOCK | QF_KEY_IS_HASH) < 0) break;
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;
	end_clock = clock();

	if (verbose) {
		printf("Number of inserts:     %lu\n", i);
		printf("Time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
		printf("Insert throughput:     %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
		printf("CPU time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	}
	results.insert_throughput = (double)i * 1000000 / (end_time - start_time);


	//warm_up_filter(&qf, 10000000ull);

	uint64_t hash = 0;
	int minirun_rank = 0;
	uint64_t fp_count = 0;
	uint64_t dummy = 0;

	if (verbose) printf("Performing queries...\n");
	start_clock = clock();
	gettimeofday(&tv, NULL);
	start_time = interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
	for (i = 0; i < query_set_len; i++) {
		if ((minirun_rank = qf_query_using_ll_table(&qf, query_set[i], &hash, QF_KEY_IS_HASH)) >= 0) {
			fp_count++;
			qf_adapt_using_ll_table(&qf, hash & minirun_id_bitmask, query_set[i], minirun_rank, QF_KEY_IS_HASH);
		}
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;
	end_clock = clock();

	if (verbose) {
		printf("Time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
		printf("Query throughput:     %f ops/sec\n", (double)query_set_len * 1000000 / (end_time - start_time));
		printf("CPU time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

		printf("False positives:      %lu\n", fp_count);
		printf("False positive rate:  %f%%\n", 100. * fp_count / query_set_len);

		printf("Dummy: %lu\n", dummy);
	}
	results.query_throughput = (double)i * 1000000 / (end_time - start_time);
	results.false_positive_rate = (double)fp_count / query_set_len;
	
	qf_free(&qf);
	return results;
}

test_results_t run_split_throughput_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, int verbose, char *inserts_outfile, char *queries_outfile) {
	test_results_t results;
	init_test_results(&results);

	size_t num_slots = 1ull << qbits;
	size_t minirun_id_bitmask = (1ull << (qbits + rbits)) - 1;

	data_config db_data_cfg = qf_data_config_init();
	splinterdb_config splinterdb_cfg = qf_splinterdb_config_init("db", &db_data_cfg);
	remove(splinterdb_cfg.filename);
	splinterdb *db;
	if (splinterdb_create(&splinterdb_cfg, &db)) {
		results.exit_code = -1;
		return results;
	}
	splinterdb_lookup_result db_result;
	splinterdb_lookup_result_init(db, &db_result, 0, NULL);

	data_config bm_data_cfg = qf_data_config_init();
	splinterdb_config backing_cfg = qf_splinterdb_config_init("bm", &bm_data_cfg);
	remove(backing_cfg.filename);
	splinterdb *bm;
	if (splinterdb_create(&backing_cfg, &bm)) {
		results.exit_code = -1;
		return results;
	}
	splinterdb_lookup_result bm_result;
	splinterdb_lookup_result_init(bm, &bm_result, 0, NULL);

	QF qf;
	if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
		results.exit_code = -1;
		return results;
	}


	double target_load = 0.9f;
	size_t max_inserts = num_slots * target_load;
	size_t num_inserts = insert_set_len > max_inserts ? max_inserts : insert_set_len;
	size_t i;

	size_t measure_freq = 100, curr_interval = 0;
	size_t measure_point = num_inserts * (curr_interval + 1) / measure_freq, prev_point = 0;


	FILE *inserts_file = inserts_outfile ? fopen(inserts_outfile, "w") : NULL;
	if (inserts_file) fprintf(inserts_file, "fill through\n");

	if (verbose) fprintf(stderr, "Performing insertions... 0.00%%");
	uint64_t num_updates = 0;
	clock_t start_clock = clock(), end_clock;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t start_time = tv.tv_sec * 1000000 + tv.tv_usec, end_time, interval_time = start_time;
	for (i = 0; qf.metadata->noccupied_slots < num_inserts; i++) {
		int ret = qf_splinter_insert_split(&qf, db, bm, insert_set[i], i);
		//int ret = qf_splinter_insert(&qf, bm, insert_set[i], 1);
		if (ret == 1) continue;
		if (ret == 0) break;
		num_updates++;

		if (qf.metadata->noccupied_slots >= measure_point) {
			gettimeofday(&tv, NULL);
			if (inserts_file) fprintf(inserts_file, "%.2f %f\n", (double)qf.metadata->noccupied_slots / num_slots * 100, (double)(i - prev_point) * 1000000 / (tv.tv_sec * 1000000 + tv.tv_usec - interval_time));
			if (verbose) fprintf(stderr, "\rPerforming insertions... %.2f%%", (double)(curr_interval + 1) / measure_freq * 100);

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

	if (verbose) {
		printf("Number of inserts:     %lu\n", i);
		printf("Number of updates:     %lu\n", num_updates);
		printf("Time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
		printf("Insert throughput:     %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
		printf("CPU time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	}
	results.insert_throughput = (double)i * 1000000 / (end_time - start_time);


	curr_interval = 0;
	measure_point = query_set_len * (curr_interval + 1) / measure_freq;
	prev_point = 0;

	FILE *queries_file = queries_outfile ? fopen(queries_outfile, "w") : NULL;
	if (queries_file) fprintf(queries_file, "queries through fprate\n");

	int still_have_space = 1;
	size_t full_point = num_slots * 0.95f;
	char buffer[10 * MAX_VAL_SIZE];
	uint64_t fp_count = 0;
	uint64_t hash;
	int minirun_rank;

	if (verbose) fprintf(stderr, "Performing queries... 0.00%%");
	start_clock = clock();
	gettimeofday(&tv, NULL);
	start_time = interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
	for (i = 0; i < query_set_len; i++) {
		if ((minirun_rank = qf_query_using_ll_table(&qf, query_set[i], &hash, QF_KEY_IS_HASH)) >= 0) {
			slice db_query = padded_slice(&query_set[i], MAX_KEY_SIZE, sizeof(query_set[i]), buffer, 0);
			splinterdb_lookup(db, db_query, &db_result);
			if (!splinterdb_lookup_found(&db_result)) {
			//if (true) {
				fp_count++;

				hash = (hash & minirun_id_bitmask) << (64 - qf.metadata->quotient_remainder_bits);
				slice bm_query = padded_slice(&hash, MAX_KEY_SIZE, sizeof(hash), buffer, 0);
				splinterdb_lookup(bm, bm_query, &bm_result);

				slice result_val;
				splinterdb_lookup_result_value(&bm_result, &result_val);

				if (still_have_space) {
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val) + minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t));
					assert((orig_key & minirun_id_bitmask) == hash);
					qf_adapt_using_ll_table(&qf, orig_key, query_set[i], minirun_rank, QF_KEY_IS_HASH);
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
			if (verbose) fprintf(stderr, "\rPerforming queries... %.2f%%", (double)(curr_interval + 1) / measure_freq * 100);

			curr_interval++;
			prev_point = i;
			measure_point = query_set_len * (curr_interval + 1) / measure_freq;

			gettimeofday(&tv, NULL);
			interval_time = tv.tv_sec * 1000000 + tv.tv_usec;
		}
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;
	end_clock = clock();

	if (queries_file) fprintf(queries_file, "%lu %f %f\n", i, (double)(i - prev_point) * 1000000 / (end_time - interval_time), (double)fp_count / i);
	if (verbose) fprintf(stderr, "\rPerforming queries... 100.00%%\n");

	if (verbose) {
		printf("Time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
		printf("Query throughput:     %f ops/sec\n", (double)query_set_len * 1000000 / (end_time - start_time));
		printf("CPU time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

		printf("False positives:      %lu\n", fp_count);
		printf("False positive rate:  %f%%\n", 100. * fp_count / query_set_len);
	}
	results.query_throughput = (double)i * 1000000 / (end_time - start_time);
	results.false_positive_rate = (double)fp_count / query_set_len;
	
	splinterdb_close(&db);
	//splinterdb_close(&bm);
	qf_free(&qf);
	return results;
}


typedef struct {
	QF *qf;
	splinterdb *db;
	uint64_t *inserts;
	size_t num_inserts;
	uint64_t num_inserted;
} thread_insert_args;

void *thread_splinter_insert_keys(void *args) {
	thread_insert_args *targs = (thread_insert_args*)args;
	splinterdb_register_thread(targs->db);

	for (int i = 0; i < targs->num_inserts; i++) {
		uint64_t hash_bucket_index = (targs->inserts[i] & ((1ULL << (targs->qf->metadata->quotient_bits + targs->qf->metadata->bits_per_slot)) - 1)) >> targs->qf->metadata->bits_per_slot;

		while (!qf_lock(targs->qf, hash_bucket_index, false, QF_KEY_IS_HASH));
		bool stop = (qf_splinter_insert(targs->qf, targs->db, targs->inserts[i], 1) == 0);
		//qf_insert_result result;
		//bool stop = (qf_insert_using_ll_table(targs->qf, targs->inserts[i], 1, &result, QF_NO_LOCK | QF_KEY_IS_HASH) != 0);
		qf_unlock(targs->qf, hash_bucket_index, false);
		
		if (stop) break;
		targs->num_inserted++;
	}

	splinterdb_deregister_thread(targs->db);
	pthread_exit(NULL);
}
/*
typedef struct {
	QF *qf;
	splinterdb *db;
	uint64_t *queries;
	size_t num_queries;
	uint64_t num_queried;
} thread_query_args;

void *thread_query_keys(void *args) {
	thread_query_keys *targs = (thread_query_args*)args;
	splinterdb_register_thread(targs->db);

	for (i = 0; i < query_set_len; i++) {
		if ((minirun_rank = qf_query_using_ll_table(&qf, query_set[i], &hash, QF_KEY_IS_HASH)) >= 0) {
			hash = (hash & minirun_id_bitmask) << (64 - qf.metadata->quotient_remainder_bits);
			slice query = padded_slice(&hash, MAX_KEY_SIZE, sizeof(hash), buffer, 0);
			splinterdb_lookup(db, query, &db_result);

			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);

			if (memcmp(&query_set[i], slice_data(result_val) + minirun_rank * MAX_VAL_SIZE, sizeof(uint64_t)) != 0) {
				fp_count++;
				if (still_have_space) {
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val) + minirun_rank * MAX_VAL_SIZE, sizeof(uint64_t));
					qf_adapt_using_ll_table(&qf, orig_key, query_set[i], minirun_rank, QF_KEY_IS_HASH);
					if (qf.metadata->noccupied_slots >= full_point) {
						still_have_space = 0;
						fprintf(stderr, "\rFilter is full after %lu queries\n", i);
					}
				}
			}
		}
	}

	for (int i = 0; i < targs->num_inserts; i++) {
		uint64_t hash_bucket_index = (targs->inserts[i] & ((1ULL << (targs->qf->metadata->quotient_bits + targs->qf->metadata->bits_per_slot)) - 1)) >> targs->qf->metadata->bits_per_slot;

		while (!qf_lock(targs->qf, hash_bucket_index, false, QF_KEY_IS_HASH));
		bool stop = (qf_splinter_insert(targs->qf, targs->db, targs->inserts[i], 1) == 0);
		//qf_insert_result result;
		//bool stop = (qf_insert_using_ll_table(targs->qf, targs->inserts[i], 1, &result, QF_NO_LOCK | QF_KEY_IS_HASH) != 0);
		qf_unlock(targs->qf, hash_bucket_index, false);
		
		if (stop) break;
		targs->num_inserted++;
	}

	splinterdb_deregister_thread(targs->db);
	pthread_exit(NULL);
	
}*/

test_results_t run_parallel_splinter_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, size_t num_threads) {
	test_results_t results;
	init_test_results(&results);

	size_t num_slots = 1ull << qbits;
	size_t minirun_id_bitmask = (1ull << (qbits + rbits)) - 1;

	data_config data_cfg = qf_data_config_init();
	splinterdb_config splinterdb_cfg = qf_splinterdb_config_init("db", &data_cfg);
	remove(splinterdb_cfg.filename);
	splinterdb *db;
	if (splinterdb_create(&splinterdb_cfg, &db)) {
		results.exit_code = -1;
		return results;
	}
	splinterdb_lookup_result db_result;
	splinterdb_lookup_result_init(db, &db_result, 0, NULL);

	QF qf;
	if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
		results.exit_code = -1;
		return results;
	}


	double target_load = 0.9f;
	size_t max_inserts = num_slots * target_load;
	size_t num_inserts = insert_set_len > max_inserts ? max_inserts : insert_set_len;
	size_t i;

	size_t measure_freq = 100, curr_interval = 0;
	size_t measure_point = num_inserts * (curr_interval + 1) / measure_freq, prev_point = 0;

	pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
	thread_insert_args *args = malloc(num_threads * sizeof(thread_insert_args));

	printf("Performing insertions... \n");
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t start_time = tv.tv_sec * 1000000 + tv.tv_usec, end_time;
	for (i = 0; i < num_threads; i++) {
		args[i].qf = &qf;
		args[i].db = db;
		args[i].inserts = &insert_set[num_inserts * i / num_threads];
		args[i].num_inserts = (num_inserts * (i + 1) / num_threads) - (num_inserts * i / num_threads);
		args[i].num_inserted = 0;
		if (pthread_create(&threads[i], NULL, thread_splinter_insert_keys, (void*)&args[i])) {
			fprintf(stderr, "Error creating thread %lu\n", i);
			results.exit_code = 1;
			return results;
		}
	}

	size_t num_inserted = 0;
	for (i = 0; i < num_threads; i++) {
		pthread_join(threads[i], NULL);
		printf("Thread %lu performed %lu successful insertions\n", i, args[i].num_inserted);
		num_inserted += args[i].num_inserted;
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;

	free(threads);
	free(args);

	printf("Number of inserts:     %lu\n", num_inserted);
	printf("Time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
	printf("Insert throughput:     %f ops/sec\n", (double)num_inserted * 1000000 / (end_time - start_time));
	results.insert_throughput = (double)num_inserted * 1000000 / (end_time - start_time);

	return results;


	int still_have_space = 1;
	size_t full_point = num_slots * 0.95f;
	char buffer[10 * MAX_VAL_SIZE];
	uint64_t fp_count = 0;
	uint64_t hash;
	int minirun_rank;

	printf("Performing queries... \n");
	gettimeofday(&tv, NULL);
	start_time = tv.tv_sec * 1000000 + tv.tv_usec;
	for (i = 0; i < query_set_len; i++) {
		if ((minirun_rank = qf_query_using_ll_table(&qf, query_set[i], &hash, QF_KEY_IS_HASH)) >= 0) {
			hash = (hash & minirun_id_bitmask) << (64 - qf.metadata->quotient_remainder_bits);
			slice query = padded_slice(&hash, MAX_KEY_SIZE, sizeof(hash), buffer, 0);
			splinterdb_lookup(db, query, &db_result);

			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);

			if (memcmp(&query_set[i], slice_data(result_val) + minirun_rank * MAX_VAL_SIZE, sizeof(uint64_t)) != 0) {
				fp_count++;
				if (still_have_space) {
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val) + minirun_rank * MAX_VAL_SIZE, sizeof(uint64_t));
					qf_adapt_using_ll_table(&qf, orig_key, query_set[i], minirun_rank, QF_KEY_IS_HASH);
					if (qf.metadata->noccupied_slots >= full_point) {
						still_have_space = 0;
						fprintf(stderr, "\rFilter is full after %lu queries\n", i);
					}
				}
			}
		}
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;

	printf("Time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
	printf("Query throughput:     %f ops/sec\n", (double)query_set_len * 1000000 / (end_time - start_time));

	printf("False positives:      %lu\n", fp_count);
	printf("False positive rate:  %f%%\n", 100. * fp_count / query_set_len);
	results.query_throughput = (double)i * 1000000 / (end_time - start_time);
	results.false_positive_rate = (double)fp_count / query_set_len;
	
	splinterdb_close(&db);
	qf_free(&qf);
	return results;
}

void *thread_insert_keys(void *args) {
	thread_insert_args *targs = (thread_insert_args*)args;

	qf_insert_result insert_result;
	int i;
	for (i = 0; i < targs->num_inserts; i++) {
		//uint64_t hash_bucket_index = (targs->inserts[i] & ((1ULL << (targs->qf->metadata->quotient_bits + targs->qf->metadata->bits_per_slot)) - 1)) >> targs->qf->metadata->bits_per_slot;

		//while (!qf_lock(targs->qf, hash_bucket_index, false, QF_KEY_IS_HASH));
		if (qf_insert_using_ll_table(targs->qf, targs->inserts[i], 1, &insert_result, QF_WAIT_FOR_LOCK | QF_KEY_IS_HASH) < 0) break;
		//qf_unlock(targs->qf, hash_bucket_index, false);
		
		//targs->num_inserted++;
	}
	targs->num_inserted = i;

	pthread_exit(NULL);
}

test_results_t run_parallel_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, size_t num_threads) {
	test_results_t results;
	init_test_results(&results);

	size_t num_slots = 1ull << qbits;
	size_t minirun_id_bitmask = (1ull << (qbits + rbits)) - 1;

	QF qf;
	if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
		results.exit_code = -1;
		return results;
	}


	double target_load = 0.9f;
	size_t max_inserts = num_slots * target_load;
	size_t num_inserts = insert_set_len > max_inserts ? max_inserts : insert_set_len;
	size_t i;

	size_t measure_freq = 100, curr_interval = 0;
	size_t measure_point = num_inserts * (curr_interval + 1) / measure_freq, prev_point = 0;

	pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
	thread_insert_args *args = malloc(num_threads * sizeof(thread_insert_args));

	printf("Performing insertions... \n");
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t start_time = tv.tv_sec * 1000000 + tv.tv_usec, end_time;
	for (i = 0; i < num_threads; i++) {
		args[i].qf = &qf;
		args[i].inserts = &insert_set[num_inserts * i / num_threads];
		args[i].num_inserts = (num_inserts * (i + 1) / num_threads) - (num_inserts * i / num_threads);
		args[i].num_inserted = 0;
		if (pthread_create(&threads[i], NULL, thread_insert_keys, (void*)&args[i])) {
			fprintf(stderr, "Error creating thread %lu\n", i);
			results.exit_code = 1;
			return results;
		}
	}

	size_t num_inserted = 0;
	for (i = 0; i < num_threads; i++) {
		pthread_join(threads[i], NULL);
		printf("Thread %lu performed %lu successful insertions\n", i, args[i].num_inserted);
		num_inserted += args[i].num_inserted;
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;

	free(threads);
	free(args);

	printf("Number of inserts:     %lu\n", num_inserted);
	printf("Time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
	printf("Insert throughput:     %f ops/sec\n", (double)num_inserted * 1000000 / (end_time - start_time));
	results.insert_throughput = (double)num_inserted * 1000000 / (end_time - start_time);

	return results;


	int still_have_space = 1;
	size_t full_point = num_slots * 0.95f;
	char buffer[10 * MAX_VAL_SIZE];
	uint64_t fp_count = 0;
	uint64_t hash;
	int minirun_rank;

	printf("Performing queries... \n");
	gettimeofday(&tv, NULL);
	start_time = tv.tv_sec * 1000000 + tv.tv_usec;
	for (i = 0; i < query_set_len; i++) {
		/*if ((minirun_rank = qf_query_using_ll_table(&qf, query_set[i], &hash, QF_KEY_IS_HASH)) >= 0) {
			hash = (hash & minirun_id_bitmask) << (64 - qf.metadata->quotient_remainder_bits);
			slice query = padded_slice(&hash, MAX_KEY_SIZE, sizeof(hash), buffer, 0);
			splinterdb_lookup(db, query, &db_result);

			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);

			if (memcmp(&query_set[i], slice_data(result_val) + minirun_rank * MAX_VAL_SIZE, sizeof(uint64_t)) != 0) {
				fp_count++;
				if (still_have_space) {
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val) + minirun_rank * MAX_VAL_SIZE, sizeof(uint64_t));
					qf_adapt_using_ll_table(&qf, orig_key, query_set[i], minirun_rank, QF_KEY_IS_HASH);
					if (qf.metadata->noccupied_slots >= full_point) {
						still_have_space = 0;
						fprintf(stderr, "\rFilter is full after %lu queries\n", i);
					}
				}
			}
		}*/
	}
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec * 1000000 + tv.tv_usec;

	printf("Time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
	printf("Query throughput:     %f ops/sec\n", (double)query_set_len * 1000000 / (end_time - start_time));

	printf("False positives:      %lu\n", fp_count);
	printf("False positive rate:  %f%%\n", 100. * fp_count / query_set_len);
	results.query_throughput = (double)i * 1000000 / (end_time - start_time);
	results.false_positive_rate = (double)fp_count / query_set_len;
	
	qf_free(&qf);
	return results;
}

