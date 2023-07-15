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
#include <assert.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"

#include "include/ll_table.h"


/*int insert_key(QF *qf, ll_table *htab, uint64_t key, int count) {
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
}*/

//ilist hash_space[1 << 20];
//int hash_space_used = 0;

int main(int argc, char **argv) {
	uint64_t seed;
	if (argc > 1) {
		seed = strtoull(argv[1], NULL, 10);
	}
	else {
		seed = time(NULL);
	}
	printf("running test on seed %lu\n", seed);
	srand(seed);

	uint64_t pbits = 20;
	uint64_t starting_qbits = 8;

	QF qf;
	if (!qf_malloc(&qf, 1ULL << starting_qbits, pbits, 0, QF_HASH_INVERTIBLE, 0)) {
		fprintf(stderr, "Can't allocate CQF.\n");
		abort();
	}
	qf_set_auto_resize(&qf, true);

	uint64_t ret_index, ret_hash;
	int ret_hash_len;

	ll_table table;
	ll_table_init(&table, 1ULL << starting_qbits);

	uint64_t keya = rand(), keyb = rand(), keyc = rand();

	qf_insert_using_ll_table(&qf, &table, keya, 1, QF_NO_LOCK);
	qf_insert_using_ll_table(&qf, &table, keyb, 1, QF_NO_LOCK);
	qf_insert_using_ll_table(&qf, &table, keyc, 1, QF_NO_LOCK);

	assert(qf_query(&qf, keya, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK) > 0);
	assert(qf_query(&qf, keyb, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK) > 0);
	assert(qf_query(&qf, keyc, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK) > 0);

	qf_free(&qf);

	printf("test completed\n");
}

