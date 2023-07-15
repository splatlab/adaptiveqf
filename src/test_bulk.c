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

#define HASH_TABLE_SIZE  11000000

struct _ilist {
	uint64_t val; // the value of the actual item
	uint64_t rem; // the remainder representing the item in the filter including extensions and the quotient
	uint64_t len; // the number of bits used by the remainder and quotient
	struct _ilist *next;
} typedef ilist;
typedef struct sglib_hashed_ilist_iterator ilist_iter;
ilist *htaba[HASH_TABLE_SIZE];
ilist *htabb[HASH_TABLE_SIZE];

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

void bp2() {}

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
			//printf("error:\tfilter claimed to have fingerprint %lu but hashtable could not find it\n", ii.rem);
		}
		else if (item_in_table->val == key) {
			insert_and_extend(qf, ret_index, key, count, key, &ret_hash, &ret_other_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
		}
		else {
			int ext_len = insert_and_extend(qf, ret_index, key, count, item_in_table->val, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
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
	if (argc < 3) {
                fprintf(stderr, "Please specify \nthe log of the number of slots in the QF\nthe number of remainder bits in the QF\n");
                // ./test 16 7 $((1 << 15)) 1000000 1 0
                exit(1);
        }
        if (argc >= 4) {
                srand(strtol(argv[3], NULL, 10));
                printf("running test on seed %ld\n", strtol(argv[3], NULL, 10));
        }
        else {
                time_t seed = time(NULL);
                printf("running test on seed %ld\n", seed);
                srand(seed);
        }

        uint64_t qbits = atoi(argv[1]);
        uint64_t rbits = atoi(argv[2]);

	qf_malloc(&qf, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0);
	qf_set_auto_resize(&qf, false);

	printf("generating inserts...\n");
	int nkeys = 0.9 * (1 << qbits);
	uint64_t *keys = malloc(nkeys * sizeof(uint64_t));
	/*for (int i = 0; i < nkeys; i++) {
		keys[i] = rand();
	}*/
	RAND_bytes((unsigned char*)keys, nkeys * sizeof(uint64_t));

	clock_t start_time, end_time;
	printf("sorting inserts...\n");
	start_time = clock();
	qsort(keys, nkeys, sizeof(uint64_t), cmp);

	end_time = clock();
	uint64_t sort_time = end_time - start_time;
	printf("time for sorting: %f sec\n", (double)(sort_time) / CLOCKS_PER_SEC);
	printf("time per item:    %f usec\n", (double)(sort_time) / nkeys);

	printf("bulk inserting...\n");
	start_time = clock();
	qf_bulk_insert(&qf, keys, nkeys);

	end_time = clock();
	uint64_t bulk_time = end_time - start_time;
	printf("time for bulk:    %f sec\n", (double)(bulk_time) / CLOCKS_PER_SEC);
	printf("time per item:    %f usec\n", (double)(bulk_time) / nkeys);

	printf("total time per item: %f usec\n", (double)(sort_time + bulk_time) / nkeys);
}

