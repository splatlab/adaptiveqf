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

double rand_zipfian(double s, double max) {
	double p = (double)rand() / RAND_MAX;
	
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

uint64_t hash_str(char *str) {
	uint64_t hash = 5381;
	int c;
	while ((c = *str++)) {
		hash = ((hash << 5) + hash) + c;
	}
	return hash;
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

splinterdb_lookup_result db_result;
int insert_key(QF *qf, splinterdb *db, uint64_t key, int count, void *buffer) {
        uint64_t ret_index, ret_hash, ret_other_hash;
        int ret_hash_len;
        int ret = qf_insert_ret(qf, key, count, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
        if (ret == QF_NO_SPACE) {
                return 0;
        }
        else if (ret == 0) {
		uint64_t fingerprint_data = ret_hash | (1ull << ret_hash_len);
		slice fingerprint = padded_slice(&fingerprint_data, MAX_KEY_SIZE, sizeof(fingerprint_data), buffer, 1);
		splinterdb_lookup(db, fingerprint, &db_result);
		map_queries++;
                if (!splinterdb_lookup_found(&db_result)) {
                        printf("error:\tfilter claimed to have fingerprint %lu but hashtable could not find it\n", ret_hash);
                }
		else {
			slice result_val;
			splinterdb_lookup_result_value(&db_result, &result_val);
			if (memcmp(slice_data(result_val), &key, sizeof(uint64_t)) == 0) {
				insert_and_extend(qf, ret_index, key, count, key, &ret_hash, &ret_other_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
			}
			else {
				uint64_t orig_key;
				memcpy(&orig_key, slice_data(result_val), sizeof(uint64_t));
				int ext_len = insert_and_extend(qf, ret_index, key, count, orig_key, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
				if (ext_len == QF_NO_SPACE) {
					printf("filter is full after insert_and_extend\n");
					return 0;
				}

				splinterdb_delete(db, result_val);
				uint64_t temp = ret_other_hash | (1ull << ext_len) | (1ull << (sizeof(uint64_t) - 1));
				db_insert(db, &temp, sizeof(temp), slice_data(result_val), MAX_VAL_SIZE, 1);
				temp = ret_hash | (1ull << ext_len);
				db_insert(db, &temp, sizeof(temp), &key, sizeof(key), 1);
				map_inserts++;
				map_inserts++;
			}
                }
		return 1;
        }
        else if (ret == 1) {
		uint64_t temp = ret_hash | (1ull << ret_hash_len);
		db_insert(db, &temp, sizeof(temp), &key, sizeof(key), 1);
		map_inserts++;
		return 1;
        }
        printf("other error: errno %d\n", ret);
        return 0;
}


int main(int argc, char **argv)
{
	if (argc < 5) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries [eg. 100000000]\nthe number of incremental queries [eg. 1000000]\n");
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
	uint64_t num_inc_queries = strtoull(argv[4], NULL, 10);


	printf("initializing hash table...\n");
	/*splinterdb *database;
	db_init(&database, "db", 64 * Mega, 128 * Mega);
	splinterdb_lookup_result_init(database, &db_result, 0, NULL);*/
	/*splinterdb *backing_map;
	if (db_init(&backing_map, "bm", 64 * Mega, 4 * Giga) != 0) {
		printf("error initializing database\n");
		exit(0);
	}
	splinterdb_lookup_result_init(backing_map, &bm_result, 0, NULL);*/

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
        //RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));
	for (i = 0; i < num_inserts; i++) {
		insert_set[i] = rand_uniform(-1);
	}


	// PERFORM INSERTS
	uint64_t target_fill = nslots * load_factor;
	printf("performing insertions... 0/%lu", target_fill);
	//std::cout << "performing insertions... 0%%" << std::flush;
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

	clock_t start_time = clock(), end_time, interval_time = start_time;
	for (i = 0; qf.metadata->noccupied_slots < target_fill; i++) {
		//if (insert_set[i] == 2205438568271438432ULL) bp2();
		if (!insert_key(&qf, database, insert_set[i], 1, buffer)) break;
		//assert(qf_query(&qf, insert_set[i], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH) != 0);
		//char query_buffer[MAX_KEY_SIZE];
		//uint64_t query_data = ret_hash | (1ULL << ret_hash_len);
		//slice query = padded_slice(&query_data, MAX_KEY_SIZE, sizeof(query_data), query_buffer);
		//splinterdb_lookup(backing_map, query, &bm_result);
		//assert(splinterdb_lookup_found(&bm_result));

		db_insert(database, &insert_set[i], sizeof(insert_set[i]), &i, sizeof(i), 0);
		//fprintf(stderr, "\rperforming insertions... %lu/%lu          ", qf.metadata->noccupied_slots, target_fill);

		if (qf.metadata->noccupied_slots >= measure_point) {
			fprintf(stderr, "\rperforming insertions... %f%%          ", current_interval * 100);
			fprintf(inserts_fp, "%f %f\n", (double)qf.metadata->noccupied_slots / nslots, (double)(i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time));

			uint64_t *inc_query_set = (uint64_t*)calloc(num_inc_queries, sizeof(uint64_t));
			RAND_bytes((unsigned char*)inc_query_set, num_inc_queries * sizeof(uint64_t));
			interval_time = clock();
			
			for (int j = 0; j < num_inc_queries; j++) {
				if (qf_query(&qf, inc_query_set[j], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
					slice query = padded_slice(&inc_query_set[j], MAX_KEY_SIZE, sizeof(inc_query_set[j]), buffer, 0);
					splinterdb_lookup(database, query, &db_result);
				}
			}
			fprintf(inc_queries_fp, "%f %f\n", (double)qf.metadata->noccupied_slots / nslots, (double)(num_inc_queries) * CLOCKS_PER_SEC / (clock() - interval_time));
			free(inc_query_set);

                        current_interval += measure_interval;
                        last_point = i;
                        measure_point = nslots * current_interval;
			interval_time = clock();
                }
	}
	end_time = clock();
	printf("\rperforming insertions... %lu/%lu          \n", target_fill, target_fill);
	//std::cout << "\rperforming insertions... 100%%         " << std::endl;
	//fprintf(outfp, "%f\t%f\n", (double)target_fill / nslots, (double)(i - last_point) * CLOCKS_PER_SEC / (end_time - interval_time));
	fclose(inserts_fp);
	fclose(inc_queries_fp);

	printf("time for inserts:      %f\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
	printf("avg insert throughput: %f ops/sec\n", (double)i * CLOCKS_PER_SEC / (end_time - start_time));

	//backing_map.print_statistics(std::cout);
	//backing_map.print_load_statistics(std::cout);

	// PERFORM QUERIES
	printf("\ngenerating query set of size %lu...\n", num_queries);
	int still_have_space = 1;
	uint64_t full_point = 0.95 * nslots;
	if (qf.metadata->noccupied_slots >= qf.metadata->nslots * 0.95) {
		still_have_space = 0;
		printf("filter is full; skipping query adaptations\n");
	}

	uint64_t *query_set = (uint64_t*)calloc(num_queries, sizeof(uint64_t));
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));

	printf("performing queries... 0%%");
	current_interval = measure_interval;
	measure_point = num_queries * current_interval;
	last_point = 0;

	FILE *queries_fp = fopen("stats_splinter_queries.csv", "w");
	fprintf(queries_fp, "queries through\n");

	uint64_t fp_count = 0;
	start_time = clock();
	for (i = 0; i < num_queries; i++) {
		if (qf_query(&qf, query_set[i], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
			slice query = padded_slice(&query_set[i], MAX_KEY_SIZE, sizeof(query_set[i]), buffer, 0);
			splinterdb_lookup(database, query, &db_result);
			if (!splinterdb_lookup_found(&db_result)) {
				fp_count++;
				if (still_have_space) {
					uint64_t fingerprint_data = ret_hash | (1ull << ret_hash_len);
					slice fingerprint = padded_slice(&fingerprint_data, MAX_KEY_SIZE, sizeof(fingerprint_data), buffer, 1);
					splinterdb_lookup(database, fingerprint, &db_result);
					slice result_val;
					splinterdb_lookup_result_value(&db_result, &result_val);
					uint64_t orig_key;
					memcpy(&orig_key, slice_data(result_val), sizeof(uint64_t));
					ret_hash_len = qf_adapt(&qf, ret_index, orig_key, query_set[i], &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
					if (ret_hash_len > 0) {
						splinterdb_delete(database, result_val);
						uint64_t temp = ret_hash | (1ull << ret_hash_len);
						db_insert(database, &temp, sizeof(temp), &orig_key, sizeof(orig_key), 1);
					}
					else if (ret_hash_len == QF_NO_SPACE) {
						still_have_space = 0;
						fprintf(stderr, "\rfilter is full after %lu queries\n", i);
						continue;
					}
					if (qf.metadata->noccupied_slots >= full_point) {
						still_have_space = 0;
						fprintf(stderr, "\rfilter is full after %lu queries\n", i);
					}
				}
			}
		}
		if (i >= measure_point) {
			fprintf(stderr, "\rperforming queries... %f%%           ", current_interval * 100);
			fprintf(queries_fp, "%lu %f\n", i, (double)(i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time));
                        current_interval += measure_interval;
                        last_point = i;
                        measure_point = num_queries * current_interval;
                        interval_time = clock();
                }

	}
	end_time = clock();
	fclose(queries_fp);
	printf("\rperforming queries... 100%%           \n");

	printf("time for queries:  %f s\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
	printf("query throughput:  %f ops/sec\n", (double)num_queries * CLOCKS_PER_SEC / (end_time - start_time));

	printf("false positive rate: %f%%\n", 100. * fp_count / num_queries);
}

