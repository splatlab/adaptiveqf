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

void print_bin_hash(uint64_t num, int qbits, int rbits) {
	for (int i = 63; i >= qbits + rbits; i--) {
		printf("%c", num & (1 << i) ? '1' : '0');
	}
	printf("|");
	for (int i = qbits + rbits - 1; i >= rbits; i--) {
		printf("%c", num & (1 << i) ? '1' : '0');
	}
	printf("|");
	for (int i = rbits - 1; i >= 0; i--) {
		printf("%c", num & (1 << i) ? '1' : '0');
	}
	printf("\n");
}

uint64_t rand_uniform(uint64_t max) {
	if (max <= RAND_MAX) return rand() % max;
	uint64_t a = rand();
	a *= ((uint64_t)RAND_MAX) + 1;
	a += rand();
	if (max == -1) return a;
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

#define HASH_TABLE_SIZE  1300000

struct _ilist {
	uint64_t val; // the value of the actual item
	uint64_t rem; // the remainder representing the item in the filter including extensions and the quotient
	int len; // the number of bits used by the remainder and quotient
	struct _ilist *next;
} typedef ilist;
typedef struct sglib_hashed_ilist_iterator ilist_iter;
ilist *htab[HASH_TABLE_SIZE];

//#define ILIST_COMPARATOR(e1, e2)    ((e1->rem + (1 << e1->len)) - (e2->rem + (1 << e2->len)))
#define ILIST_COMPARATOR(e1, e2)    (e1->len == e2->len ? (e1->rem == e2->rem ? 0 : (e1->rem > e2->rem ? 1 : -1)) : (e1->len > e2->len ? 1 : -1))

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

void bp2() {}

int insert_key(QF *qf, ilist **htab, uint64_t key, int count, uint8_t flags) {
	uint64_t ret_index, ret_hash, ret_other_hash;
	int ret_hash_len;
	int ret = qf_insert_ret(qf, key, count, &ret_index, &ret_hash, &ret_hash_len, flags); // attempt to insert
	if (ret == QF_NO_SPACE) {
		return 0;
	}
	else if (ret == 0) { // if a matching fingerprint is found, search hash table for original item
		ilist ii, *nn;
		ii.rem = ret_hash;
		ii.len = ret_hash_len;
		ilist *item_in_table = sglib_hashed_ilist_find_member(htab, &ii);
		if (item_in_table == NULL) {
			//printf("error:\tfilter claimed to have fingerprint %lu but hashtable could not find it\n", ii.rem);
		}
		else if (item_in_table->val == key) {
			insert_and_extend(qf, ret_index, key, count, key, &ret_hash, &ret_other_hash, flags);
		}
		else {
			int ext_len = insert_and_extend(qf, ret_index, key, count, item_in_table->val, &ret_hash, &ret_other_hash, flags);
			if (ext_len == QF_NO_SPACE) {
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

QF qf;
int cmp(const void *elem1, const void *elem2) {
	uint64_t hash1 = *((uint64_t*)elem1);
	uint64_t hash2 = *((uint64_t*)elem2);
	return qf_hash_cmp(&qf, hash1, hash2);
}

int main(int argc, char **argv)
{
	if (argc < 5) {
                fprintf(stderr, "Please specify \nthe log of the number of slots in the QF\nthe number of remainder bits in the QF\nthe blacklist size\nthe whitelist size\n");
                // ./test 16 7 $((1 << 15)) 1000000 1 0
                exit(1);
        }
        if (argc > 5) {
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
	uint64_t blacksize;
	uint64_t whitesize = atoi(argv[4]);
	uint64_t max_entries = 0.94 * (1 << qbits);
	if (whitesize == -1) {
		blacksize = max_entries * (atof(argv[3]));
		whitesize = max_entries - blacksize;
	}
	else {
		blacksize = atoi(argv[3]);
	}

	qf_malloc(&qf, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0);
	qf_set_auto_resize(&qf, false);

	sglib_hashed_ilist_init(htab);
	ilist ii, *nn;

	//start_recording();
	uint64_t ret_index, ret_hash;
	int ret_hash_len;

	printf("making blacklist inserts...\n");
	uint64_t i, j, k;
	for (i = 0; i < blacksize; i++) {
		//printf("%lu\n", i);
		/*if (i == 130435) {
			printf("here\n");
		}*/
		uint64_t key = rand_uniform(-1);
		if (!insert_key(&qf, htab, key, 1, QF_NO_LOCK)) {
			printf("out of space to insert after %lu insertions\n", i);
			abort();
		}
		assert(qf_query(&qf, key, &ret_index, &ii.rem, &ii.len, QF_NO_LOCK));
		nn = sglib_hashed_ilist_find_member(htab, &ii);
		assert(nn != NULL && nn->val == key);
	}
	printf("inserted %lu items\n", i);

	//stop_recording();

	int exts = 0;
	printf("making whitelist extensions...\n");
	for (j = 0; j < whitesize; j++) {
		uint64_t key = rand_uniform(-1);
		if (qf_query(&qf, key, &ret_index, &ii.rem, &ii.len, QF_NO_LOCK)) {
			nn = sglib_hashed_ilist_find_member(htab, &ii);
			if (nn == NULL) printf("hash table could not find item %lu\n", key);
			else if (nn->val != key) {
				if (qf_get_num_occupied_slots(&qf) >= qf.metadata->nslots * 0.95) {
					printf("out of space to adapt after %lu queries\n", j);
					abort();
				}
				else {
					ret_hash_len = qf_adapt(&qf, ret_index, nn->val, key, &ret_hash, QF_NO_LOCK);
					if (nn->len > 0) {
						exts++;
						sglib_hashed_ilist_delete(htab, nn);
						nn->rem = ret_hash;
						nn->len = ret_hash_len;
						sglib_hashed_ilist_add(htab, nn);
					}
				}
				assert(!qf_query(&qf, key, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK));
				assert(qf_query(&qf, nn->val, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK));
			}
		}
	}
	printf("extended from %lu items\n", j);
	printf("made %d extensions\n", exts);

	/*printf("verifying results...\n");
	uint64_t ret_index, ret_hash;
	int ret_hash_len;
	for (i = 0; i < nkeys; i++) {
		assert(qf_query(&qf, keys[i], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	}*/

	printf("testing false positive rate...\n");
	int fps = 0;
	uint64_t nqueries = 10000000;
	for (k = 0; k < nqueries; k++) {
		uint64_t key = rand_uniform(-1);
		if (qf_query(&qf, key, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK)) {
			ii.rem = ret_hash;
			ii.len = ret_hash_len;
			nn = sglib_hashed_ilist_find_member(htab, &ii);
			if (nn == NULL) printf("hash table could not find item %lu\n", key);
			else if (nn->val != key) fps++;
		}
	}
	printf("false positive rate: %f\n", (double)fps / nqueries);
	printf("space usage (bytes): %lu\n", (uint64_t)(qf.metadata->noccupied_slots * (3.125f + rbits) / 8));
}

