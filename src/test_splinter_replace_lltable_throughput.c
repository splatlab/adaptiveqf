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
#include "include/ll_table.h"

#include <splinterdb/data.h>
#include <splinterdb/default_data_config.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/public_util.h>
#include <splinterdb/splinterdb.h>


#define Kilo (1024UL)
#define Mega (1024UL * Kilo)
#define Giga (1024UL * Mega)

#define TEST_DB_NAME "db"

#define MAX_KEY_SIZE 16
#define MAX_VAL_SIZE 16

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

int merge_tuples(const data_config *cfg, slice key, message old_message, merge_accumulator *new_message) {
	assert(merge_accumulator_length(new_message) == MAX_VAL_SIZE);
	if (!merge_accumulator_resize(new_message, message_length(old_message) + MAX_VAL_SIZE)) return -1;
	memcpy(merge_accumulator_data(new_message) + MAX_VAL_SIZE, message_data(old_message), message_length(old_message));
	return 0;
}

int merge_tuples_final(const data_config *cfg, slice key, merge_accumulator *oldest_message) {
	merge_accumulator_set_class(oldest_message, MESSAGE_TYPE_INSERT);
	return 0;
}

void pad_data(void *dest, const void *src, const size_t dest_len, const size_t src_len, const int flagged) {
	if (dest_len >= src_len) {
		if (flagged) memset(dest, 0xff, dest_len);
		else bzero(dest, dest_len);
	}
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
	slice key_slice = slice_create(MAX_KEY_SIZE, key_padded);
	slice val_slice = slice_create(MAX_VAL_SIZE, val_padded);

	return splinterdb_insert(database, key_slice, val_slice);
}

splinterdb_lookup_result db_result;
int insert_key(QF *qf, splinterdb *db, uint64_t key, int count) {
	uint64_t ret_hash;
	int ret = qf_insert_using_ll_table(qf, key, count, &ret_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
	if (ret < 0) {
		return 0;
	}
	else {
		uint64_t minirun_id = ret_hash | ((1ull << (qf->metadata->quotient_bits + qf->metadata->bits_per_slot)) - 1);
		if (ret == 0) {
			db_insert(db, &minirun_id, sizeof(minirun_id), &key, sizeof(key), 0);
		}
		else {
			char buffer[MAX_VAL_SIZE * 8]; // just a safe guess for buffer size
			slice minirun_id_slice = padded_slice(&minirun_id, MAX_KEY_SIZE, sizeof(minirun_id), buffer, 0);
			splinterdb_lookup(db, minirun_id_slice, &db_result);
			assert(splinterdb_lookup_found(&db_result));
			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);
			assert(slice_length(result_val) + MAX_VAL_SIZE <= MAX_VAL_SIZE * 8);

			bzero(buffer, MAX_VAL_SIZE * 8);
			memcpy(buffer, &minirun_id, sizeof(minirun_id));
			memcpy(buffer + MAX_VAL_SIZE, slice_data(result_val), slice_length(result_val));

			splinterdb_delete(db, minirun_id_slice);
			slice new_linked_list_slice = slice_create(slice_length(result_val) + MAX_VAL_SIZE, buffer);
			splinterdb_insert(db, minirun_id_slice, new_linked_list_slice);
		}
		return 1;
	}
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

	unsigned int murmur_seed = rand();


	printf("initializing hash table...\n");
	const char *db_path = "db";//"/ssd1/richard/aqf/db";
	remove(db_path);
	data_config data_cfg;
	default_data_config_init(MAX_KEY_SIZE, &data_cfg);
	data_cfg.merge_tuples = merge_tuples;
	data_cfg.merge_tuples_final = merge_tuples_final;
	splinterdb_config splinterdb_cfg = (splinterdb_config){
		.filename   = db_path,
		.cache_size = 64 * Mega,
		.disk_size  = 20 * Giga,
		.data_cfg   = &data_cfg,
		.io_flags   = O_RDWR | O_CREAT | O_DIRECT
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

	/*
	// UNIT TESTING
	uint64_t k = 1, v = 2;
	db_insert(database, &k, sizeof(k), &v, sizeof(v), 0);

	char unit_buffer[128];
	slice unit_query = padded_slice(&k, MAX_KEY_SIZE, sizeof(k), unit_buffer, 0);
	splinterdb_lookup(database, unit_query, &db_result);
	slice unit_result;
	splinterdb_lookup_result_value(&db_result, &unit_result);
	assert(slice_length(unit_result) == MAX_VAL_SIZE);
	memcpy(&v, slice_data(unit_result), sizeof(v));
	assert(v == 2);

	v = 3;
	db_insert(database, &k, sizeof(k), &v, sizeof(v), 0);
	
	splinterdb_lookup(database, unit_query, &db_result);
	splinterdb_lookup_result_value(&db_result, &unit_result);
	assert(slice_length(unit_result) == MAX_VAL_SIZE * 2);
	memcpy(&v, slice_data(unit_result), MAX_KEY_SIZE);
	assert(v == 3);
	memcpy(&v, slice_data(unit_result) + MAX_KEY_SIZE, MAX_KEY_SIZE);
	assert(v == 2);

	exit(0);*/

	printf("generating insert set of size %lu...\n", num_inserts);
	uint64_t i;
	uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
	RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));

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
		if (!insert_key(&qf, database, insert_set[i], 1)) break;

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
	/*for (i = 0; i < num_queries; i++) { // making the distrubution uniform from a limited universe
		query_set[i] = query_set[i] % (1ull << 24);
		query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
	}*/

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
		if (qf_query_using_ll_table(&qf, query_set[i], &query_result, QF_KEY_IS_HASH)) {
			uint64_t temp = query_result.hash | ((1ull << (qbits + rbits)) - 1);
			slice query = padded_slice(&temp, MAX_KEY_SIZE, sizeof(temp), buffer, 0);
			splinterdb_lookup(database, query, &db_result);
			/*if (!splinterdb_lookup_found(&db_result)) {
				fp_count++;
				faulty_fps++;
				continue;
			}*/
			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);

			if (i >= warmup_queries && memcmp(&query_set[i], slice_data(result_val) + query_result.intralist_rank * MAX_KEY_SIZE, sizeof(uint64_t)) != 0) {
				fp_count++;
				if (still_have_space) {
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val) + query_result.intralist_rank * MAX_KEY_SIZE, sizeof(uint64_t));
					qf_query_using_ll_table(&qf, query_set[i], &query_result, QF_KEY_IS_HASH);
					//qf_adapt_using_ll_table(&qf, ret_index, orig_key, query_set[i], &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
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

