
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <openssl/rand.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"

#include "sglib.h"

int bp2() {
	return 0;
}

int clear_log() {
	FILE *fp = fopen("data/frames.txt", "w");
	if (fp == NULL) return 0;
	fclose(fp);
	return 1;
}

uint64_t rand_uniform(uint64_t max) {
  if (max <= RAND_MAX) return rand() % max;
  uint64_t a = rand();
  uint64_t b = rand();
  a |= (b << 32);
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

#define HASH_TABLE_SIZE  88000000

struct _ilist {
	uint64_t val; // the value of the actual item
	uint64_t rem; // the remainder representing the item in the filter including extensions and the quotient
	uint64_t len; // the number of bits used by the remainder and quotient
	struct _ilist *next;
} typedef ilist;
typedef struct sglib_hashed_ilist_iterator ilist_iter;
ilist *htab[HASH_TABLE_SIZE];

//#define ILIST_COMPARATOR(e1, e2)    ((e1->rem + (1 << e1->len)) - (e2->rem + (1 << e2->len)))
#define ILIST_COMPARATOR(e1, e2)    (e1->len != e2->len ? e2->len - e1->len : e2->rem - e1->rem)

unsigned int ilist_hash_function(ilist *e) {
	return e->rem + (1 << e->len);
}

SGLIB_DEFINE_LIST_PROTOTYPES(ilist, ILIST_COMPARATOR, next)
SGLIB_DEFINE_LIST_FUNCTIONS(ilist, ILIST_COMPARATOR, next)
SGLIB_DEFINE_HASHED_CONTAINER_PROTOTYPES(ilist, HASH_TABLE_SIZE, ilist_hash_function)
SGLIB_DEFINE_HASHED_CONTAINER_FUNCTIONS(ilist, HASH_TABLE_SIZE, ilist_hash_function)

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


int main(int argc, char **argv)
{
	if (argc < 6) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF\nthe number of remainder bits in the QF\nthe number of inserts\nthe number of queries\nthe number of trials\n");
		// ./test 16 7 $((1 << 15)) 1000000 1 0
		exit(1);
	}
	if (argc >= 7) {
		srand(strtol(argv[6], NULL, 10));
		printf("running test on seed %ld\n", strtol(argv[6], NULL, 10));
	}
	else {
		time_t seed = time(NULL);
		printf("running test on seed %ld\n", seed);
		srand(seed);
	}
	
	clear_log();
	
	uint64_t qbits = atoi(argv[1]);
	uint64_t rbits = atoi(argv[2]);
	uint64_t nhashbits = qbits + rbits;
	uint64_t nslots = (1ULL << qbits);
	uint64_t universe = -1;
	uint64_t num_inserts = strtoull(argv[3], NULL, 10);
	uint64_t num_queries = strtoull(argv[4], NULL, 10);
	uint64_t expected_fill = 1lu << (qbits + rbits) < num_inserts ? 1lu << (qbits + rbits) : num_inserts;
	if (universe < expected_fill) {
		printf("warning: universe may be too small to fill filter to completion\n");
	}
	else if (universe < expected_fill * 10) {
		printf("warning: filling the filter may take longer than necessary because of low universe size\n");
	}

	char buffer[1024];
	//FILE *shalla = fopen("data/shalla.txt", "r");
	//FILE *caida = fopen("data/20140619-140100.csv", "r");
	//fgets(buffer, sizeof(buffer), caida);
	
	double avgInsTime = 0, avgInsPer = 0, avgQryTime = 0, avgQryPer = 0, avgFP = 0, avgFill = 0, minFP = INFINITY, maxFP = 0;
	double avgInsSlots = 0, avgQrySlots = 0;
	int numtrials = atoi(argv[5]);
	int trials = 0;
	for (; trials < numtrials; trials++) {
		QF qf;

		if (!qf_initfile(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0,
										 "mycqf.file")) {
			fprintf(stderr, "Can't allocate CQF.\n");
			abort();
		}

		qf_set_auto_resize(&qf, false);
		
		uint64_t count_fp = 0, count_p = 0;
		
		//uint64_t *nodes = malloc(sizeof(node) * num_inserts);
		//uint64_t *tree = nodes[0]; // simple reverse map for testing - index equals hash
		//uint64_t* values = malloc(sizeof(uint64_t) * num_inserts);
		uint64_t ret_index, ret_hash, ret_other_hash;
		int ret_hash_len;

		/*qf_insert_ret(&qf, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		assert(qf_query(&qf, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH));
		qf_remove(&qf, 1, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		assert(!qf_query(&qf, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH));
		qf_insert_ret(&qf, 3, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		qf_insert_ret(&qf, 2, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		qf_insert_ret(&qf, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		qf_remove(&qf, 2, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		assert(!qf_query(&qf, 2, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH));
		qf_remove(&qf, 3, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		assert(!qf_query(&qf, 3, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH));
		qf_insert_ret(&qf, 2, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		qf_insert_ret(&qf, 1 + (1 << rbits), 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		assert(!qf_insert_ret(&qf, 1 + (1 << (rbits + qbits)), 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH));
		insert_and_extend(&qf, ret_index, 1 + (1 << (rbits + qbits)), 1, 1, &ret_hash, &ret_other_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
		qf_remove(&qf, 1, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
		snapshot(&qf);
		abort();*/
		
		sglib_hashed_ilist_init(htab);
		ilist ii, *nn, *ll;
		ilist_iter it;
		
		if (universe < (double)qf.metadata->range * 0.8) {
			printf("warning: universe may be too small to fill filter to completion\n");
		}
		
		clock_t start_time, end_time;
		start_time = clock();

		// PERFORM INSERTS
		printf("starting %lu inserts...\n", num_inserts);
		uint64_t i, j, k = 0;
		for (i = 0; i < num_inserts; i++) {
			j = rand_uniform(universe); // pick a random number
			//j = rand_zipfian(1.01, 100000000);
			//fgets(buffer, sizeof(buffer), shalla);
			//j = hash_str(buffer);
			int ret = qf_insert_ret(&qf, j, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH); // attempt to insert
			if (ret == QF_NO_SPACE) { // if filter is full, stop
				printf("filter is full after %lu inserts\n", i);
				break;
			}
			else if (ret == 0) { // if a matching fingerprint is found, search hash table for original item
				ii.rem = ret_hash;
				ii.len = ret_hash_len;
				ilist *item_in_table = sglib_hashed_ilist_find_member(htab, &ii);
				if (item_in_table != NULL) {
					int ext_len = insert_and_extend(&qf, ret_index, j, 1, item_in_table->val, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
					if (ext_len == QF_NO_SPACE) {
						printf("filter is full after %lu inserts\n", i);
						break;
					}
					else if (ext_len == 0) continue;

					sglib_hashed_ilist_delete(htab, item_in_table);
					item_in_table->rem = ret_other_hash;
					item_in_table->len = ext_len;
					sglib_hashed_ilist_add(htab, item_in_table);

					nn = malloc(sizeof(ilist));
					nn->val = j;
					nn->rem = ret_hash;
					nn->len = ext_len;
					sglib_hashed_ilist_add(htab, nn);
					k++;
				}
			}
			else if (ret == 1) {
				nn = malloc(sizeof(ilist));
				nn->val = j;
				nn->rem = ret_hash;
				nn->len = ret_hash_len;
				sglib_hashed_ilist_add(htab, nn);
			}
			else {
				printf("other error: errno %d\n", ret);
				break;
			}
		}

		avgFill += (double)i / nslots;
		printf("nslots = %lu\n", qf_get_nslots(&qf));
		printf("made %lu inserts\n", i);
		printf("number of slots used after inserts: %lu\n", qf_get_num_occupied_slots(&qf));
		printf("extended %lu times\n", k);
		avgInsTime += (clock() - start_time);
		avgInsPer += (double)(clock() - start_time) / i;
		avgInsSlots += qf_get_num_occupied_slots(&qf);
		
		// PERFORM QUERIES
		printf("starting %lu queries...\n", num_queries);
		int still_have_space = 1;
		if (qf_get_num_occupied_slots(&qf) >= qf.metadata->nslots * 0.95) {
			still_have_space = 0;
			printf("filter is full; skipping query adaptations\n");
		}

		uint64_t *query_set = calloc(num_queries, sizeof(uint64_t));
		for (i = 0; i < num_queries; i++) {
			query_set[i] = rand_uniform(-1);
			//query_set[i] = rand_zipfian(1.5f, 1lu << 30);
		}

		start_time = clock();

		for (i = 0; i < num_queries; i++) {
			//j = rand_uniform(universe);
			//j = rand_zipfian(1.5f, (1lu << 30));
			j = query_set[i];
			//fgets(buffer, sizeof(buffer), shalla);
			//j = hash_str(buffer);
			/*fgets(buffer, sizeof(buffer), caida);
			csv_get(buffer, 3);
			j = hash_str(buffer);*/

			if (qf_query(&qf, j, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
				count_p++;
				ii.rem = ret_hash;
				ii.len = ret_hash_len;
				ilist *potential_item = sglib_hashed_ilist_find_member(htab, &ii);
				if (potential_item == NULL) bp2();
				else if (potential_item->val != j) {
					count_fp++;
					if (still_have_space) {
						if (qf_get_num_occupied_slots(&qf) >= qf.metadata->nslots * 0.95) {
							still_have_space = 0;
							printf("out of space to adapt after %lu queries\n", i);
						}
						else {
							ret_hash_len = qf_adapt(&qf, ret_index, potential_item->val, j, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
							if (ret_hash_len > 0) {
								sglib_hashed_ilist_delete(htab, potential_item);
								potential_item->rem = ret_hash;
								potential_item->len = ret_hash_len;
								sglib_hashed_ilist_add(htab, potential_item);
							}
						}
					}
				}
			}
			if (i % 100000 == 0 && i != 0) {
				printf("churning...\n");
				uint64_t current_fill = qf_get_num_occupied_slots(&qf);
				for (int r = 0; r < 0.01 * current_fill;) {
					j = rand_uniform(-1);
					if (qf_query(&qf, j, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
						ii.rem = ret_hash;
						ii.len = ret_hash_len;
						ilist **potential_items = NULL;
						sglib_hashed_ilist_delete_if_member(htab, &ii, potential_items);
						qf_remove(&qf, j, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK);
					}
				}
			}
		}
		
		end_time = clock();
		printf("completed in time %ld us\n", end_time - start_time);
		printf("performed %lu queries uniformly from a universe of %lu\n", num_queries, universe);
		printf("number of slots used after queries: %lu\n", qf_get_num_occupied_slots(&qf));
		printf("false positive rate: %f\n", (double)count_fp / num_queries);
		
		for(ll=sglib_hashed_ilist_it_init(&it,htab); ll!=NULL; ll=sglib_hashed_ilist_it_next(&it)) {
			free(ll);
		}
		avgQryTime += (end_time - start_time);
		avgQryPer += (double)(end_time - start_time) / num_queries;
		double fp_rate = (double)count_fp / num_queries;
		avgFP += fp_rate;
		if (fp_rate > maxFP) maxFP = fp_rate;
		if (fp_rate < minFP) minFP = fp_rate;
		avgQrySlots += qf_get_num_occupied_slots(&qf);
	}
	printf("\nperformed %d trials\n", numtrials);
	printf("min false positive rate: %f\n", minFP);
	printf("avg false positive rate: %f\n", avgFP / numtrials);
	printf("max false positive rate: %f\n", maxFP);
	printf("avg fill rate: %f\n", avgFill / numtrials);
	printf("avg slots used after inserts: %f\n", avgInsSlots / numtrials);
	printf("avg slots used after queries: %f\n", avgQrySlots / numtrials);
	printf("avg total insert time: %f\n", avgInsTime / numtrials);
	printf("avg insert time per item: %f\n", avgInsPer / numtrials);
	printf("avg total query time: %f\n", avgQryTime / numtrials);
	printf("avg query time per item: %f\n", avgQryPer / numtrials);
}

