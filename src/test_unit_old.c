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
#include "include/hashutil.h"

#include "sglib.h"


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

#define HASH_TABLE_SIZE  85000

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


void test_insertions(int qbits, int rbits) {
	QF qf;
	assert(qf_malloc(&qf, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));

	uint64_t ret_index, ret_hash;
	int ret_hash_len;
	assert(qf_insert_ret(&qf, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_ret(&qf, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) == 0);

	assert(qf_free(&qf));

	printf("insertion test successful\n");
}

void test_deletions(int qbits, int rbits) {
	QF qf;
	assert(qf_malloc(&qf, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));

	uint64_t ret_index, ret_hash, ret_other_hash;
	int ret_hash_len;
	assert(qf_insert_ret(&qf, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_ret(&qf, 2, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_ret(&qf, 1 + (1 << rbits), 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_ret(&qf, 1 + (2 << rbits), 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_ret(&qf, 2 + (2 << rbits), 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);

	assert(qf_query(&qf, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1 + (1 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1 + (2 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2 + (2 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(!qf_query(&qf, 3, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(!qf_query(&qf, 1 + (3 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));

	// delete item from run
	qf_remove(&qf, 2, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK);

	assert(!qf_query(&qf, 2, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2 + (2 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));

	// delete only item in run
	qf_remove(&qf, 1 + (1 << rbits), &ret_hash, &ret_hash_len, QF_KEY_IS_HASH);

	assert(!qf_query(&qf, 1 + (1 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));

	// delete runend
	qf_remove(&qf, 1 + (2 << rbits), &ret_hash, &ret_hash_len, QF_KEY_IS_HASH);
	
	assert(!qf_query(&qf, 1 + (2 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2 + (2 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));

	assert(!qf_insert_ret(&qf, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK));
	assert(insert_and_extend(&qf, ret_index, 1 + (1 << (rbits + qbits)), 2, 1, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK));
	assert(qf_query(&qf, 1 + (1 << (rbits + qbits)), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH) == 2);

	qf_remove(&qf, 1, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH);

	assert(!qf_query(&qf, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1 + (1 << (rbits + qbits)), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));

	qf_remove(&qf, 1 + (1 << (rbits + qbits)), &ret_hash, &ret_hash_len, QF_KEY_IS_HASH);

	assert(!qf_query(&qf, 1 + (1 << (rbits + qbits)), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2 + (2 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));

	snapshot(&qf);

	assert(qf_free(&qf));
	printf("deletion test successful\n");
}

void test_merge() {
	int qbits = 6, rbits = 6;
	QF qfa, qfb, qfc;
	assert(qf_malloc(&qfa, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));
	assert(qf_malloc(&qfb, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));
	assert(qf_malloc(&qfc, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));

	uint64_t ret_index, ret_hash, ret_other_hash;
	int ret_hash_len;
	
	assert(qf_insert_ret(&qfa, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_ret(&qfb, 2, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_ret(&qfb, 1 + (1 << rbits), 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_ret(&qfa, 1 + (1 << (rbits + qbits)), 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) == 0);
	assert(insert_and_extend(&qfa, ret_index, 1 + (1 << (rbits + qbits)), 1, 1, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK));
	assert(qf_insert_ret(&qfa, 2, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	
	qf_merge(&qfa, &qfb, &qfc);

	assert(qf_query(&qfc, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qfc, 2, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH) == 2);
	assert(qf_query(&qfc, 1 + (1 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qfc, 1 + (1 << (rbits + qbits)), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(!qf_query(&qfc, 1 + (2 << (rbits + qbits)), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	
	snapshot(&qfc);

	assert(qf_free(&qfa) && qf_free(&qfb) && qf_free(&qfc));
	printf("merge test successful\n");
}

void test_merge_2() {
	int qbits = 6, rbits = 7;
	QF qfa, qfb, qfc;
	assert(qf_malloc(&qfa, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));
	assert(qf_malloc(&qfb, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));
	assert(qf_malloc(&qfc, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));

	uint64_t ret_index, ret_hash;
	int ret_hash_len;

	double load = 0.4f;
	int fill = load * (1 << qbits);

	uint64_t *keys_a = calloc(fill, sizeof(uint64_t));
	uint64_t *keys_b = calloc(fill, sizeof(uint64_t));
	for (int i = 0; i < fill; i++) {
		uint64_t key = rand_uniform(-1);
		if (!qf_insert_ret(&qfa, key, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) i--;
		else keys_a[i] = key;
	}
	for (int i = 0; i < fill; i++) {
		uint64_t key = rand_uniform(-1);
		if (!qf_insert_ret(&qfb, key, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) i--;
		else keys_b[i] = key;
	}

	qf_merge(&qfa, &qfb, &qfc);
	
	for (int i = 0; i < fill; i++) {
		assert(qf_query(&qfc, keys_a[i], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
		assert(qf_query(&qfc, keys_b[i], &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	}
	for (int i = 0; i < 100000; i++) {
		uint64_t key = rand_uniform(-1);
		if (!qf_query(&qfc, key, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
			assert(!qf_query(&qfa, key, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH) && !qf_query(&qfb, key, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
		}
	}

	assert(qf_free(&qfa) && qf_free(&qfb) && qf_free(&qfc));
	printf("merge test 2 successful\n");
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

/*uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
        const uint64_t m = 0xc6a4a7935bd1e995;
        const int r = 47;

        uint64_t h = seed ^ (len * m);

        const uint64_t * data = (const uint64_t *)key;
        const uint64_t * end = data + (len/8);

        while(data != end)
        {
                uint64_t k = *data++;

                k *= m;
                k ^= k >> r;
                k *= m;

                h ^= k;
                h *= m;
        }

        const unsigned char * data2 = (const unsigned char*)data;

        switch(len & 7)
        {
                case 7: h ^= (uint64_t)data2[6] << 48;
                case 6: h ^= (uint64_t)data2[5] << 40;
                case 5: h ^= (uint64_t)data2[4] << 32;
                case 4: h ^= (uint64_t)data2[3] << 24;
                case 3: h ^= (uint64_t)data2[2] << 16;
                case 2: h ^= (uint64_t)data2[1] << 8;
                case 1: h ^= (uint64_t)data2[0];
                                                h *= m;
        };

        h ^= h >> r;
        h *= m;
        h ^= h >> r;

        return h;
}*/

struct _ht_node {
	uint64_t key;
	uint64_t tag;
	struct _ht_node* next;
} typedef ht_node;

ht_node** ht_init(uint64_t size) {
	return calloc(size, sizeof(ht_node*));
}

#define HASH_SET_SEED 26571997
void ht_insert(ht_node** ht, int len, uint64_t tag, uint64_t key) {
	uint64_t hash = MurmurHash64A((void*)&tag, sizeof(uint64_t), HASH_SET_SEED) % len;
	ht_node *node = malloc(sizeof(node));
	node->key = key;
	node->tag = tag;
	node->next = ht[hash];
	ht[hash] = node;
}

void ht_free(ht_node** ht, int len) {
	int i;
	for (i = 0; i < len; i++) {
		if (ht[i] == NULL) continue;
		ht_node *ptr = ht[i];
		ht_node *next;
		do {
			next = ptr->next;
			free(ptr);
			ptr = next;
		} while (ptr != NULL);
	}
	free(ht);
}

int main(int argc, char **argv)
{
	int ht_len = 1500000;
	ht_node **ht = ht_init(ht_len);

	int i;
	for (i = 0; i < 1000000; i++) {
		uint64_t r = rand();
		ht_insert(ht, ht_len, r & 0b11111111, r);
	}

	ht_free(ht, ht_len);

	/*uint64_t seed = time(NULL);
	if (argc > 1) {
		seed = strtoull(argv[1], NULL, 10);
	}
	printf("running tests on seed %lu\n", seed);
	srand(seed);

	if (argc > 3) {
		uint64_t key = strtoull(argv[2], NULL, 10);
		uint64_t mask = (1ULL << atoi(argv[3])) - 1;

		uint64_t hash = hash_64(key, mask);
		printf("%lu -> %lu\n", key, hash);
		printf("%lu <- %lu\n", hash_64i(hash, mask), hash);
	}*/

	//test_insertions(16, 7);
	//test_deletions(16, 7);
	//test_merge(16, 7);
	//test_merge_2(16, 7);
}

