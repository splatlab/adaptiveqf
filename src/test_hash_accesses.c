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

#define HASH_TABLE_SIZE  87000000

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

int insert_key(QF *qf, ilist **htab, uint64_t key, int count) {
        uint64_t ret_index, ret_hash, ret_other_hash;
        int ret_hash_len;
        int ret = qf_insert_ret(qf, key, count, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH); // attempt to insert
        if (ret == QF_NO_SPACE) {
                return 0;
        }
        else if (ret == 0) { // if a matching fingerprint is found, search hash table for original item
                ilist ii, *nn;
                ii.rem = ret_hash;
                ii.len = ret_hash_len;
                ilist *item_in_table = sglib_hashed_ilist_find_member(htab, &ii);
                if (item_in_table == NULL) {
                        printf("error:\tfilter claimed to have fingerprint %lu but hashtable could not find it\n", ii.rem);
                }
                else if (item_in_table->val == key) {
                        insert_and_extend(qf, ret_index, key, count, key, &ret_hash, &ret_other_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
                }
                else if (1){
                        int ext_len = insert_and_extend(qf, ret_index, key, count, item_in_table->val, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
                        if (ext_len == QF_NO_SPACE) {
                                printf("filter is full after insert_and_extend\n");
                                return 0;
                        }
                        sglib_hashed_ilist_delete(htab, item_in_table);
                        item_in_table->rem = ret_other_hash;
                        item_in_table->len = ext_len;
                        sglib_hashed_ilist_add(htab, item_in_table);

                        nn = malloc(sizeof(ilist));
                        nn->val = key;
                        nn->rem = ret_hash;
                        nn->len = ext_len;
                        sglib_hashed_ilist_add(htab, nn);
                }
        }
        else if (ret == 1) {
                ilist *nn = malloc(sizeof(ilist));
                nn->val = key;
                nn->rem = ret_hash;
                nn->len = ret_hash_len;
                sglib_hashed_ilist_add(htab, nn);
        }
        else {
                printf("other error: errno %d\n", ret);
                return 0;
        }
        return 1;
}

//ilist hash_space[1 << 20];
//int hash_space_used = 0;

int main(int argc, char **argv)
{
	if (argc < 5) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF\nthe number of remainder bits in the QF\nthe load factor\nthe number of queries\n");
		// ./test 16 7 $((1 << 15)) 1000000 1 0
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
	
	uint64_t qbits = atoi(argv[1]);
	uint64_t rbits = atoi(argv[2]);
	uint64_t nhashbits = qbits + rbits;
	uint64_t nslots = (1ULL << qbits);
	double load_factor = atof(argv[3]);
	uint64_t num_inserts = nslots * load_factor;//strtoull(argv[3], NULL, 10);
	uint64_t num_queries = strtoull(argv[4], NULL, 10);

	QF qf;
	if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
		fprintf(stderr, "Can't allocate CQF.\n");
		abort();
	}

	qf_set_auto_resize(&qf, false);

	//start_recording();

	//uint64_t *nodes = malloc(sizeof(node) * num_inserts);
	//uint64_t *tree = nodes[0]; // simple reverse map for testing - index equals hash
	//uint64_t* values = malloc(sizeof(uint64_t) * num_inserts);
	uint64_t ret_index, ret_hash;
	int ret_hash_len;

	sglib_hashed_ilist_init(htab);
	ilist ii;

	uint64_t i, j;

	// PERFORM INSERTS
	printf("starting %lu inserts...\n", num_inserts);

	uint64_t target_fill = nslots * load_factor;
	for (i = 0; qf.metadata->noccupied_slots < target_fill; i++) {
		if (!insert_key(&qf, htab, rand_uniform(-1), 1)) break;
	}

	//snapshot(&qf);
	//stop_recording();

	// PERFORM QUERIES
	printf("starting %lu queries...\n", num_queries);
	int still_have_space = 1;
	if (qf.metadata->noccupied_slots >= qf.metadata->nslots * 0.95) {
		still_have_space = 0;
		printf("filter is full; skipping query adaptations\n");
	}

	FILE *fp = fopen("target/hash_accesses.txt", "w");
	fprintf(fp, "queries\taccesses\n");
	fclose(fp);
	fp = fopen("target/hash_accesses.txt", "a");

	uint64_t step_size = 100000;
	uint64_t hash_accesses = 0;
	for (i = 0; i < num_queries; i++) {
		j = rand_zipfian(1.5f, 1ull << 30);

		if (qf_query(&qf, j, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
			ii.rem = ret_hash;
			ii.len = ret_hash_len;
			ilist *potential_item = sglib_hashed_ilist_find_member(htab, &ii);
			hash_accesses++;
			if (potential_item == NULL) bp2();
			else if (potential_item->val != j) {
				if (still_have_space) {
					ret_hash_len = qf_adapt(&qf, ret_index, potential_item->val, j, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
					if (ret_hash_len > 0) {
						sglib_hashed_ilist_delete(htab, potential_item);
						potential_item->rem = ret_hash;
						potential_item->len = ret_hash_len;
						sglib_hashed_ilist_add(htab, potential_item);
					}
					else if (ret_hash_len == QF_NO_SPACE) {
						still_have_space = 0;
						printf("filter is full after %lu queries\n", i);
					}
				}
			}
		}

		if (i % step_size == 0) {
			fprintf(fp, "%lu\t%lu\n", i, hash_accesses);
		}
	}

	fclose(fp);

	printf("done\n");
}

