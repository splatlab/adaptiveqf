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
	if (flagged) memset(dest, 0xff, dest_len);
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
		slice fingerprint = padded_slice(&fingerprint_data, MAX_KEY_SIZE, sizeof(fingerprint_data), buffer, 0);
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
				uint64_t temp = ret_other_hash | (1ull << ext_len);
				db_insert(db, &temp, sizeof(temp), slice_data(result_val), MAX_VAL_SIZE, 0);
				temp = ret_hash | (1ull << ext_len);
				db_insert(db, &temp, sizeof(temp), &key, sizeof(key), 0);
				map_inserts++;
				map_inserts++;
			}
                }
		return 1;
        }
        else if (ret == 1) {
		uint64_t temp = ret_hash | (1ull << ret_hash_len);
		db_insert(db, &temp, sizeof(temp), &key, sizeof(key), 0);
		map_inserts++;
		return 1;
        }
        printf("other error: errno %d\n", ret);
        return 0;
}


int main(int argc, char **argv)
{
	if (argc < 5) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries [eg. 100000000]\nthe adversarial frequency [eg. 10]\n");
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
	double load_factor = 0.85f;
	uint64_t num_inserts = nslots * load_factor;//strtoull(argv[3], NULL, 10);
	uint64_t num_queries = strtoull(argv[3], NULL, 10);
	uint64_t adv_freq = strtoull(argv[4], NULL, 10);

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


	printf("generating insert set of size %lu...\n", num_inserts);
	uint64_t i;
	uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
        RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));
	/*for (i = 0; i < num_inserts; i++) {
		insert_set[i] = rand_uniform(-1);
	}*/


	// PERFORM INSERTS
	uint64_t target_fill = nslots * load_factor;
	uint64_t ret_index, ret_hash;
	int ret_hash_len;
	char buffer[MAX_KEY_SIZE + MAX_VAL_SIZE];

	double measure_interval = 0.01f;
        double current_interval = measure_interval;
        uint64_t measure_point = target_fill * current_interval;

	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;
	clock_t start_clock = clock(), end_clock;
	for (i = 0; qf.metadata->noccupied_slots < target_fill; i++) {
		if (!insert_key(&qf, database, insert_set[i], 1, buffer)) break;

		if (qf.metadata->noccupied_slots >= measure_point) {
			fprintf(stderr, "\rperforming insertions... %f%%          ", current_interval * 100);

			current_interval += measure_interval;
			measure_point = nslots * current_interval;
		}
	}
	end_clock = clock();
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

	printf("\n");

	printf("time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
	printf("avg insert throughput: %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
	printf("cpu time for inserts:  %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	// PERFORM QUERIES
	printf("\ngenerating query set of size %lu...\n", num_queries);
	int still_have_space = 1;
	uint64_t full_point = 0.95 * nslots;
	if (qf.metadata->noccupied_slots >= qf.metadata->nslots * 0.95) {
		still_have_space = 0;
		printf("filter is full; skipping query adaptations\n");
	}

	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));

	current_interval = measure_interval;
	measure_point = num_queries * current_interval;
	uint64_t last_point = 0;

	//double adversarial_clock = 0;
	uint64_t adv_index = 0;
	int max_adv_set_len = 1000000;
	uint64_t *adv_set = malloc(max_adv_set_len * sizeof(uint64_t));
	int adv_set_len = 0;

	FILE *adv_fp = fopen("stats_splinter_adversarial.csv", "w");
	fprintf(adv_fp, "queries through fprate slots\n");

	uint64_t fp_count = 0;
	uint64_t adv_attempted = 0, adv_successful = 0, adv_failed = 0, adv_found = 0;

	gettimeofday(&timecheck, NULL);
	uint64_t interval_time = start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	start_clock = clock();
	for (i = 0; i < num_queries; i++) {
		if (i % adv_freq == 0 && adv_set_len > 0) {
			adv_attempted++;
			if (adv_index >= adv_set_len) adv_index = 0;
			uint64_t key = adv_set[adv_index];
			if (qf_query(&qf, key, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
				uint64_t temp = ret_hash & BITMASK(ret_hash_len);
				slice query = padded_slice(&temp, MAX_KEY_SIZE, sizeof(temp), buffer, 0);
				splinterdb_lookup(database, query, &db_result);
				slice result_val;
				splinterdb_lookup_result_value(&db_result, &result_val);

				if (memcmp(&key, slice_data(result_val), sizeof(uint64_t)) != 0) {
					fp_count++;
					adv_successful++;
					if (still_have_space) {
						uint64_t orig_key;
						memcpy(&orig_key, slice_data(result_val), sizeof(uint64_t));
						ret_hash_len = qf_adapt(&qf, ret_index, orig_key, key, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
						if (ret_hash_len > 0) {
							splinterdb_delete(database, result_val);
							uint64_t temp = ret_hash & BITMASK(ret_hash_len);
							db_insert(database, &temp, sizeof(temp), &orig_key, sizeof(orig_key), 0);
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
			if (qf_query(&qf, query_set[i], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
				uint64_t temp = ret_hash & BITMASK(ret_hash_len);
				slice query = padded_slice(&temp, MAX_KEY_SIZE, sizeof(temp), buffer, 0);
				splinterdb_lookup(database, query, &db_result);
				slice result_val;
				splinterdb_lookup_result_value(&db_result, &result_val);

				if (memcmp(&query_set[i], slice_data(result_val), sizeof(uint64_t)) != 0) {
					fp_count++;
					if (adv_set_len < max_adv_set_len) {
						adv_set[adv_set_len++] = query_set[i];
						adv_found++;
					}
					if (still_have_space) {
						uint64_t orig_key;
						memcpy(&orig_key, slice_data(result_val), sizeof(uint64_t));
						ret_hash_len = qf_adapt(&qf, ret_index, orig_key, query_set[i], &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
						if (ret_hash_len > 0) {
							splinterdb_delete(database, result_val);
							uint64_t temp = ret_hash | (1ull << ret_hash_len);
							db_insert(database, &temp, sizeof(temp), &orig_key, sizeof(orig_key), 0);
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
		}

		if (i >= measure_point) {
			gettimeofday(&timecheck, NULL);

			fprintf(adv_fp, "%lu %f %f %lu\n", i, (double)(i - last_point) * 1000000 / (timecheck.tv_sec * 1000000 + timecheck.tv_usec - interval_time), (double)fp_count / i, qf.metadata->noccupied_slots);

			fprintf(stderr, "\rperforming queries... %f%%           ", current_interval * 100);
                        current_interval += measure_interval;
                        measure_point = num_queries * current_interval;
			last_point = i;

			gettimeofday(&timecheck, NULL);
			interval_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
                }

	}
	end_clock = clock();
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;

	printf("\n");

	printf("time for queries:      %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("query throughput:      %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
	printf("cpu time for queries:  %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	printf("adv found:             %lu\n", adv_found);
	printf("adv attempted:         %lu\n", adv_attempted);
	printf("adv successful:        %lu\n", adv_successful);
	printf("adv failed:            %lu\n", adv_failed);

	printf("false positives:     %lu\n", fp_count);
	printf("false positive rate: %f%%\n", 100. * fp_count / num_queries);
}

