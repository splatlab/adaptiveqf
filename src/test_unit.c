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
#include "include/rand_util.h"
#include "include/splinter_util.h"


void test_insertions(int qbits, int rbits) {
	QF qf;
	assert(qf_malloc(&qf, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));

	uint64_t ret_hash;
	assert(qf_insert_using_ll_table(&qf, 1, 1, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qf, 1, 1, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qf, 2, 2, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qf, 1 + (1 << rbits), 1, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);

	assert(qf_free(&qf));

	printf("insertion test successful\n");
}

void test_deletions(int qbits, int rbits) {
	QF qf;
	assert(qf_malloc(&qf, 1 << qbits, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0));

	uint64_t ret_hash;
	assert(qf_insert_using_ll_table(&qf, 1, 1, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qf, 2, 1, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qf, 1 + (1 << rbits), 1, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qf, 1 + (2 << rbits), 1, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qf, 2 + (2 << rbits), 1, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);

	assert(qf_query(&qf, 1, &ret_hash, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2, &ret_hash, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1 + (1 << rbits), &ret_hash, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1 + (2 << rbits), &ret_hash, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2 + (2 << rbits), &ret_hash, QF_KEY_IS_HASH));
	assert(!qf_query(&qf, 3, &ret_hash, QF_KEY_IS_HASH));
	assert(!qf_query(&qf, 1 + (3 << rbits), &ret_hash, QF_KEY_IS_HASH));

	// delete item from run
	qf_remove(&qf, 2, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK);

	assert(!qf_query(&qf, 2, &ret_hash, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1, &ret_hash, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2 + (2 << rbits), &ret_hash, QF_KEY_IS_HASH));

	// delete only item in run
	qf_remove(&qf, 1 + (1 << rbits), &ret_hash, &ret_hash_len, QF_KEY_IS_HASH);

	assert(!qf_query(&qf, 1 + (1 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));

	// delete runend
	qf_remove(&qf, 1 + (2 << rbits), &ret_hash, &ret_hash_len, QF_KEY_IS_HASH);
	
	assert(!qf_query(&qf, 1 + (2 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));
	assert(qf_query(&qf, 2 + (2 << rbits), &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH));

	assert(!qf_insert_using_ll_table(&qf, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK));
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
	
	assert(qf_insert_using_ll_table(&qfa, 1, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qfb, 2, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qfb, 1 + (1 << rbits), 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	assert(qf_insert_using_ll_table(&qfa, 1 + (1 << (rbits + qbits)), 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) == 0);
	assert(insert_and_extend(&qfa, ret_index, 1 + (1 << (rbits + qbits)), 1, 1, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK));
	assert(qf_insert_using_ll_table(&qfa, 2, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH | QF_NO_LOCK) > 0);
	
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
		if (!qf_insert_using_ll_table(&qfa, key, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) i--;
		else keys_a[i] = key;
	}
	for (int i = 0; i < fill; i++) {
		uint64_t key = rand_uniform(-1);
		if (!qf_insert_using_ll_table(&qfb, key, 1, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) i--;
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

	test_insertions(20, 7);
	test_deletions(20, 7);
	test_merge(20, 7);
	test_merge_2(20, 7);
}

