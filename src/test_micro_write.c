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

FILE *insert_trace_fp;
FILE *query_trace_fp;

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

#define HASH_SET_SEED 26571997
struct _set_node {
        struct _set_node *next;
        uint64_t key;
	uint64_t value;
} typedef set_node;

int set_insert(set_node *set, int set_len, uint64_t key, uint64_t value) {
        if (!key) key++;
        uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->key) {
                ptr->key = key;
		ptr->value = value;
        }
	else {
                while (ptr->next) {
                        if (ptr->key == key) return 0;
                        ptr = ptr->next;
                }
                if (ptr->key == key) return 0;
                set_node *node = malloc(sizeof(set_node));
                ptr->next = node;

                node->next = NULL;
                node->key = key;
		node->value = value;
        }
        return 1;
}

int set_query(set_node *set, int set_len, uint64_t key, uint64_t *value) {
        if (!key) key++;
        uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->key) {
                return 0;
        }
	else {
                while (ptr->next) {
                        if (ptr->key == key) {
				*value = ptr->value;
                                return 1;
                        }
                        ptr = ptr->next;
                }
                if (ptr->key == key){
			*value = ptr->value;
                        return 1;
                }
		else return 0;
        }
}

int set_delete(set_node *set, int set_len, uint64_t key) {
	if (!key) key++;
	uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
	set_node *ptr = &set[hash % set_len];
	if (!ptr->key) {
		return 0;
	}
	else if (ptr->key == key) {
		if (ptr->next) {
			set_node *temp = ptr->next;
			ptr->key = ptr->next->key;
			ptr->value = ptr->next->value;
			ptr->next = ptr->next->next;
			free(temp);
		}
		else {
			ptr->key = 0;
		}
		return 1;
	}
	else if (!ptr->next) {
		return 0;
	}
	else {	
		do {
			if (ptr->next->key == key) {
				set_node *temp = ptr->next;
				ptr->next = ptr->next->next;
				free(temp);
				return 1;
			}
			ptr = ptr->next;
		} while (ptr->next);
		return 0;
	}
}

int insert_key(QF *qf, set_node *set, uint64_t set_len, uint64_t key, int count) {
        uint64_t ret_index, ret_hash, ret_other_hash;
        int ret_hash_len;
        int ret = qf_insert_ret(qf, key, count, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
	fprintf(insert_trace_fp, "i %lu 0 0\n", key);
        if (ret == QF_NO_SPACE) {
                return 0;
        }
        else if (ret == 0) {
		uint64_t fingerprint = ret_hash | (1ull << ret_hash_len), orig_key;
		ret = set_query(set, set_len, fingerprint, &orig_key);
                if (!ret) {
                        printf("error:\tfilter claimed to have fingerprint %lu but hashtable could not find it\n", ret_hash);
                }
		else {
			if (key == orig_key) {
				insert_and_extend(qf, ret_index, key, count, key, &ret_hash, &ret_other_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
				fprintf(insert_trace_fp, "e %lu %lu %lu\n", ret_index, key, key);
			}
			else {
				int ext_len = insert_and_extend(qf, ret_index, key, count, orig_key, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
				fprintf(insert_trace_fp, "e %lu %lu %lu\n", ret_index, key, orig_key);
				if (ext_len == QF_NO_SPACE) {
					printf("filter is full after insert_and_extend\n");
					return 0;
				}

				set_delete(set, set_len, fingerprint);
				set_insert(set, set_len, ret_other_hash | (1ull << ext_len), orig_key);
				set_insert(set, set_len, ret_hash | (1ull << ext_len), key);
			}
                }
		return 1;
        }
        else if (ret == 1) {
		set_insert(set, set_len, ret_hash | (1ull << ret_hash_len), key);
		return 1;
        }
        printf("other error: errno %d\n", ret);
        return 0;
}


int main(int argc, char **argv)
{
	uint64_t seed = 0;
	if (argc < 4) {
                fprintf(stderr, "Please specify \nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries [eg. 100000000]\n");
                // ./test 16 7 $((1 << 15)) 1000000 disk_throughput_stats.txt 0
                exit(1);
        }
        if (argc >= 5) {
                srand(strtoull(argv[4], NULL, 10));
                printf("running test on seed %llu\n", strtoull(argv[4], NULL, 10));
        }
        else {
                seed = time(NULL);
                printf("running test on seed %lu\n", seed);
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

	insert_trace_fp = fopen("sim_trace_insert.txt", "w");
	fprintf(insert_trace_fp, "%lu %lu %lu %lu\n", seed, qbits, rbits, num_inserts);
	query_trace_fp = fopen("sim_trace_query.txt", "w");
	fprintf(query_trace_fp, "%lu\n", num_queries);

	printf("initializing hash table...\n");
	uint64_t set_len = num_inserts * 1.5f;
	set_node *set = calloc(set_len, sizeof(set_node));

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

	double measure_interval = 0.01f;
        double current_interval = measure_interval;
        uint64_t measure_point = target_fill * current_interval;

	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;
	for (i = 0; qf.metadata->noccupied_slots < target_fill; i++) {
		if (!insert_key(&qf, set, set_len, insert_set[i], 1)) break;

		if (qf.metadata->noccupied_slots >= measure_point) {
			fprintf(stderr, "\rperforming insertions... %f%%          ", current_interval * 100);

                        current_interval += measure_interval;
                        measure_point = nslots * current_interval;
                }
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("\n");

	printf("time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughput:     %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
	printf("cpu time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	// PERFORM QUERIES
	printf("\ngenerating query set of size %lu...\n", num_queries);
	int still_have_space = 1;

	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));
	/*for (i = 0; i < num_queries; i++) {
		query_set[i] = rand_uniform(-1);
	}*/
	//unsigned int murmur_seed = rand();
	/*for (i = 0; i < num_queries; i++) { // making the distrubution uniform from a limited universe
		query_set[i] = rand_zipfian(1.5f, 10000000ull, query_set[i]);
		//query_set[i] = query_set[i] % (1ull << 24);
		//query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
	}*/

	current_interval = measure_interval;
	measure_point = num_queries * current_interval;

	uint64_t full_point = nslots * 0.95f;

	uint64_t fp_count = 0;
	start_clock = clock();
	gettimeofday(&timecheck, NULL);
	start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	for (i = 0; i < num_queries; i++) {
		fprintf(query_trace_fp, "q %lu 0 0\n", query_set[i]);
		if (qf_query(&qf, query_set[i], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
			uint64_t temp = ret_hash | (1ull << ret_hash_len), orig_key = 0;
			set_query(set, set_len, temp, &orig_key);

			if (query_set[i] != orig_key) {
				fp_count++;
				if (still_have_space) {
					fprintf(query_trace_fp, "a %lu %lu %lu\n", ret_index, orig_key, query_set[i]);
					ret_hash_len = qf_adapt(&qf, ret_index, orig_key, query_set[i], &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
					if (ret_hash_len > 0) {
						int ret = set_delete(set, set_len, temp);
						if (ret == 0) {
							printf("%lu\n", i);
							abort();
						}
						set_insert(set, set_len, ret_hash | (1ull << ret_hash_len), orig_key);
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

                        current_interval += measure_interval;
                        measure_point = num_queries * current_interval;
                }

	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("\rperforming queries... 100%%           \n");

	printf("time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
	printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
	printf("cpu time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	printf("false positives:      %lu\n", fp_count);
	printf("false positive rate:  %f%%\n", 100. * fp_count / num_queries);
}

