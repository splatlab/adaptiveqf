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

int bp_count = 0;
void bp() {
	printf("%d\n", bp_count);
	bp_count++;
}

void csv_get_col(char* buffer, int col) {
	int i, j;
	for (i = 0; buffer[i] != '\0' && col > 0; i++) {
		if (buffer[i] == ',') col--;
	}
	for (j = 0; buffer[i + j] != '\0' && buffer[i + j] != ','; j++) {
		buffer[j] = buffer[i + j];
	}
	buffer[j] = '\0';
}

uint64_t hash_str(char *str) {
	uint64_t hash = 5381;
	int c;
	while ((c = *str++)) {
		hash = ((hash << 5) + hash) + c;
	}
	return hash;
}

int main(int argc, char **argv)
{
	if (argc < 5) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries [eg. 100000000]\nthe number of incremental queries [eg. 1000000]\n");
		// eg ./test_splinter_throughput 26 9 100000000 1000000
		// results will be output to stats_splinter_inserts.csv, stats_splinter_inc_queries.csv, stats_splinter_queries.csv, stats_splinter_fprates.csv
		exit(1);
	}
	if (argc >= 6) { // optional seed
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
	uint64_t num_inc_queries = strtoull(argv[4], NULL, 10);


	printf("initializing database...\n");
	data_config data_cfg = qf_data_config_init();
	splinterdb_config splinterdb_cfg = qf_splinterdb_config_init("db", &data_cfg);
	splinterdb *database;
	if (splinterdb_create(&splinterdb_cfg, &database)) {
		printf("Error creating database\n");
		exit(0);
	}
	splinterdb_lookup_result db_result;
	splinterdb_lookup_result_init(database, &db_result, 0, NULL);

	printf("initializing filter...\n");
	QF qf;
	if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
		printf("Can't allocate CQF.\n");
		exit(0);
	}
	qf_set_auto_resize(&qf, false);


	printf("generating insert set of size %lu...\n", num_inserts);
	uint64_t i;
	uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
	//RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));
	for (uint64_t j = 0; j < num_inserts; j++) {
		insert_set[j] = rand_uniform(-1ULL);
	}

	// PERFORM INSERTS
	uint64_t target_fill = nslots * load_factor;
	uint64_t ret_index, ret_hash;
	int ret_hash_len;
	char buffer[MAX_KEY_SIZE + MAX_VAL_SIZE];

	double measure_interval = 0.01f;
	double current_interval = measure_interval;
	uint64_t measure_point = target_fill * current_interval, last_point = 0;

	FILE *inserts_fp = fopen("stats_splinter_inserts.csv", "w");
	fprintf(inserts_fp, "fill through\n");
	FILE *inc_queries_fp = fopen("stats_splinter_inc_queries.csv", "w");
	fprintf(inc_queries_fp, "fill through\n");

	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time, interval_time = start_time;
	for (i = 0; qf.metadata->noccupied_slots < target_fill; i++) {
		if (!qf_splinter_insert(&qf, database, insert_set[i], 1)) break;

		if (qf.metadata->noccupied_slots >= measure_point) {
			gettimeofday(&timecheck, NULL);
			fprintf(inserts_fp, "%f %f\n", (double)qf.metadata->noccupied_slots / nslots * 100, (double)(i - last_point) * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time));

			uint64_t *inc_query_set = malloc(num_inc_queries * sizeof(uint64_t));
			RAND_bytes((unsigned char*)inc_query_set, num_inc_queries * sizeof(uint64_t));

			gettimeofday(&timecheck, NULL);
			interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
			for (int j = 0; j < num_inc_queries; j++) {
				if (qf_query(&qf, inc_query_set[j], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
					uint64_t temp = ret_hash | (1ull << ret_hash_len);
					slice query = padded_slice(&temp, MAX_KEY_SIZE, sizeof(temp), buffer, 0);
					splinterdb_lookup(database, query, &db_result);
				}
			}
			gettimeofday(&timecheck, NULL);
			fprintf(inc_queries_fp, "%f %f\n", (double)qf.metadata->noccupied_slots / nslots * 100, (double)(num_inc_queries) * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time));
			free(inc_query_set);

			fprintf(stderr, "\rperforming insertions... %f%%          ", current_interval * 100);

			current_interval += measure_interval;
			last_point = i;
			measure_point = nslots * current_interval;

			gettimeofday(&timecheck, NULL);
			interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	fclose(inserts_fp);
	fclose(inc_queries_fp);

	printf("\n");

	printf("time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughput:     %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
	printf("cpu time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	// PERFORM QUERIES
	printf("\ngenerating query set of size %lu...\n", num_queries);
	int still_have_space = 1;
	uint64_t full_point = 0.95 * nslots;
	if (qf.metadata->noccupied_slots >= qf.metadata->nslots * 0.95) {
		still_have_space = 0;
		printf("filter is full; skipping query adaptations\n");
	}

	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	for (int i = 0; i < num_queries; i++) {
		query_set[i] = rand_uniform(-1ull);
	}
	//RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));
	//unsigned int murmur_seed = rand();
	/*for (i = 0; i < num_queries; i++) { // making the distrubution uniform from a limited universe
		query_set[i] = query_set[i] % (1ull << 24);
		query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
	}*/
	/*char dataset_buffer[256];
	//FILE *shalla = fopen("data/shalla.txt", "r");
	FILE *caida = fopen("data/20140619-140100.csv", "r");
	fgets(dataset_buffer, sizeof(dataset_buffer), caida);
	for (int q = 0; q < num_queries; q++) {
		//fgets(buffer, sizeof(buffer), shalla);
		//query_set[q] = hash_str(buffer);
		fgets(dataset_buffer, sizeof(dataset_buffer), caida);
		csv_get_col(dataset_buffer, 3);
		query_set[q] = hash_str(dataset_buffer);
	}
	//fclose(shalla);
	fclose(caida);*/

	printf("performing queries... 0%%");
	uint64_t warmup_queries = 0;//49999999ull;

	current_interval = measure_interval;
	measure_point = num_queries * current_interval;
	last_point = warmup_queries;

	FILE *queries_fp = fopen("stats_splinter_queries.csv", "w");
	fprintf(queries_fp, "queries through\n");
	FILE *fprates_fp = fopen("stats_splinter_fprates.csv", "w");
	fprintf(fprates_fp, "queries fprate\n");

	qf_query_result query_result;
	uint64_t fp_count = 0;
	//uint64_t faulty_fps = 0;
	start_clock = clock();
	gettimeofday(&timecheck, NULL);
	start_time = interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	for (i = 0; i < num_queries; i++) {
		if (qf_query_using_ll_table(&qf, query_set[i], &query_result, QF_KEY_IS_HASH) >= 0) {
			uint64_t temp = query_result.hash & ((1ull << (qbits + rbits)) - 1);
			slice query = padded_slice(&temp, MAX_KEY_SIZE, sizeof(temp), buffer, 0);
			splinterdb_lookup(database, query, &db_result);
			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);

			if (i >= warmup_queries && memcmp(&query_set[i], slice_data(result_val) + query_result.minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t)) != 0) {
				fp_count++;
				if (still_have_space) {
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val) + query_result.minirun_rank * MAX_KEY_SIZE, sizeof(uint64_t));
					qf_adapt_using_ll_table(&qf, orig_key, query_set[i], query_result.minirun_rank, QF_KEY_IS_HASH | QF_NO_LOCK);
					if (qf.metadata->noccupied_slots >= full_point) {
						still_have_space = 0;
						fprintf(stderr, "\rfilter is full after %lu queries\n", i);
					}
				}
			}
		}

		if (i >= measure_point) {
			gettimeofday(&timecheck, NULL);

			if (i >= warmup_queries) {
				fprintf(queries_fp, "%lu %f\n", i - warmup_queries, (double)(i - last_point) * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time));
				fprintf(fprates_fp, "%lu %f\n", i - warmup_queries, (double)fp_count / i);
			}
			fprintf(stderr, "\rperforming queries... %f%%           ", current_interval * 100);

			current_interval += measure_interval;
			last_point = i;
			measure_point = num_queries * current_interval;

			gettimeofday(&timecheck, NULL);
			interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	fclose(queries_fp);
	fclose(fprates_fp);

	printf("\n");

	printf("time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
	printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
	printf("cpu time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	printf("false positives:      %lu\n", fp_count);
	printf("false positive rate:  %f%%\n", 100. * fp_count / num_queries);
	//printf("faulty false pos:     %lu\n", faulty_fps);
}

