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
#include <fcntl.h>
#include <pthread.h>
#include <openssl/rand.h>

//#define INC_PC

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"
#include "include/hashutil.h"

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

int cmp(const void *a, const void *b) {
	return *(uint64_t*)a - *(uint64_t*)b;
}

typedef struct {
	QF *qf;
	uint64_t *inserts;
	size_t num_inserts;
	uint64_t num_inserted;
} thread_insert_args;

void *thread_insert_keys(void *args) {
	thread_insert_args *targs = (thread_insert_args*)args;

	for (int i = 0; i < targs->num_inserts; i++) {
		//uint64_t hash_bucket_index = (targs->inserts[i] % targs->qf->metadata->range) >> targs->qf->metadata->bits_per_slot;

		//while (!qf_lock(targs->qf, hash_bucket_index, false, QF_KEY_IS_HASH));
		bool stop = (qf_insert(targs->qf, targs->inserts[i], 0, 1, QF_WAIT_FOR_LOCK | QF_KEY_IS_HASH) == QF_NO_SPACE);
		//qf_unlock(targs->qf, hash_bucket_index, false);
		
		if (stop) break;
		targs->num_inserted++;
	}

	pthread_exit(NULL);
}


int main(int argc, char **argv)
{
	uint64_t seed = 0;
	if (argc < 5) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries [eg. 100000000]\nthe number of threads [eg. 16]\n");
		// ./test 16 7 $((1 << 15)) 1000000 disk_throughput_stats.txt 0
		exit(1);
	}
	if (argc >= 6) {
		seed = strtoull(argv[5], NULL, 10);
		printf("running test on seed %lu\n", seed);
		srand(seed);
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
	uint64_t num_threads = strtoull(argv[4], NULL, 10);

	/*printf("initializing hash table...\n");
	uint64_t set_len = num_inserts * 1.5f;
	set_node *set = calloc(set_len, sizeof(set_node));*/

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

	/*printf("sorting insert set...\n");
	qsort(insert_set, num_inserts, sizeof(uint64_t), cmp);*/

	printf("generating query set of size %lu...\n", num_queries);
	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));
	//unsigned int murmur_seed = rand();
	for (i = 0; i < num_queries; i++) { // making the distrubution uniform from a limited universe
		//query_set[i] = rand_zipfian(1.5f, 10000000ull, query_set[i]);
		//query_set[i] = query_set[i] % (1ull << 24);
		//query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
	}

	// PERFORM INSERTS
	printf("performing inserts...\n");

	pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
	thread_insert_args *args = malloc(num_threads * sizeof(thread_insert_args));

	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	clock_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;

	for (i = 0; i < num_threads; i++) {
		args[i].qf = &qf;
		args[i].inserts = &insert_set[num_inserts * i / num_threads];
		args[i].num_inserts = (num_inserts * (i + 1) / num_threads) - (num_inserts * i / num_threads);
		args[i].num_inserted = 0;
		if (pthread_create(&threads[i], NULL, thread_insert_keys, (void*)&args[i])) {
			fprintf(stderr, "Error creating thread %lu\n", i);
			return 1;
		}
	}

	size_t num_inserted = 0;
	for (i = 0; i < num_threads; i++) {
		pthread_join(threads[i], NULL);
		printf("Thread %lu performed %lu successful insertions\n", i, args[i].num_inserted);
		num_inserted += args[i].num_inserted;
	}
	
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("num inserts:           %lu\n", num_inserted);
	printf("time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
	printf("avg insert throughput: %f ops/sec\n", (double)num_inserted * 1000000 / (end_time - start_time));
	printf("cpu time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	// PERFORM QUERIES

	uint64_t fp_count = 0, value;
	printf("performing queries...\n");

	start_clock = clock();
	gettimeofday(&timecheck, NULL);
	start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	for (i = 0; i < num_queries; i++) {
		if (qf_query(&qf, query_set[i], &value, QF_KEY_IS_HASH)) {
			fp_count++;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
	printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
	printf("cpu time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	printf("false positives:      %lu\n", fp_count);
	printf("false positive rate:  %f%%\n", 100. * fp_count / num_queries);
}

