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

int snapshot(QF *qf) {
	FILE *fp = fopen("data/snapshot.txt", "w");
	if (fp == NULL) return 0;
	char buffer1[128];
	char buffer2[256];
	bzero(buffer1, 128);
	bzero(buffer2, 256);
	int i, j;
	for (i = 0; i * 64 < qf->metadata->xnslots; i++) {
		uint64_t occupied = get_block(qf, i)->occupieds[0];
		for (j = 0; j < 64; j++) {
			buffer1[63 - j] = '0' + occupied % 2;
			occupied >>= 1;
		}
		sprintf(buffer2, "%d\t%s\n", i, buffer1);
		//printf("%s", buffer2);
		fputs(buffer2, fp);
		uint64_t runend = get_block(qf, i)->runends[0];
		for (j = 0; j < 64; j++) {
			buffer1[63 - j] = '0' + runend % 2;
			runend >>= 1;
		}
		sprintf(buffer2, "\t%s\n", buffer1);
		//printf("%s", buffer2);
		fputs(buffer2, fp);
		uint64_t extension = get_block(qf, i)->extensions[0];
		for (j = 0; j < 64; j++) {
			buffer1[63 - j] = '0' + extension % 2;
			extension >>= 1;
		}
		sprintf(buffer2, "\t%s\n", buffer1);
		fputs(buffer2, fp);
	}
	fclose(fp);
	return 1;
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

struct _keyValuePair {
	uint64_t key;
	uint64_t val;
	struct _keyValuePair *left;
	struct _keyValuePair *right;
} typedef keyValuePair;

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

int find(uint64_t* array, int len, uint64_t item) {
	uint64_t i;
	for (i = 0; i < len; i++)
		if (array[i] == item) return 1;
	return 0;
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

// returns the number of low order bits on which hash1 and hash2 match
int hashCmp(uint64_t hash1, uint64_t hash2) {
	//printf("hashCmp: %lu, %lu\n", hash1, hash2);
	int i;
	for (i = 0; i < 64; i++) {
		if ((hash1 & 1) != (hash2 & 1)) break;
		hash1 >>= 1;
		hash2 >>= 1;
	}
	return i;
}

keyValuePair *getItem(keyValuePair *root, uint64_t hash) {
	if (root == NULL) return root;
	//printf("getItem: %lu, %lu\n", root->key, hash);
	int cmp = hashCmp(root->key, hash);
	if (cmp == 64) return root;
	keyValuePair *ret = NULL;
	if (((hash >> cmp) & 1) > ((root->key >> cmp) & 1)) {
		if (root->right != NULL) ret = getItem(root->right, hash);
		else return root;
	}
	else {
		if (root->left != NULL) ret = getItem(root->left, hash);
		else return root;
	}
	//printf("%lu, %lu\n", ret->key, hash);
	if (hashCmp(ret->key, hash) > cmp) return ret;
	else return root;
}

keyValuePair *insertItem(keyValuePair* root, keyValuePair *item) {
	if (root == NULL || item == NULL) return item;
	//printf("insItem: %lu, %lu\n", root->key, item->key);
	int cmp = hashCmp(root->key, item->key);
	if (cmp / 128 > 0) {
		if (root->right == NULL) root->right = item;
		else insertItem(root->right, item);
	}
	else {
		if (root->left == NULL) root->left = item;
		else insertItem(root->left, item);
	}
	return root;
}

int matchpart(uint64_t item, uint64_t hash, uint64_t hash_len, uint64_t qbits, uint64_t rbits) {
	if ((item & ((2 << rbits) - 1)) != (hash & ((2 << rbits) - 1))) return 0;
	hash_len -= rbits;
	hash >>= rbits;
	item >>= qbits + rbits;
	while (hash_len > 0) {
		if ((item & ((2 << rbits) - 1)) != (hash & ((2 << rbits) - 1))) return 0;
		hash_len -= rbits;
		hash >>= rbits;
		item >>= rbits;
	}
	return 1;
}

int findpart(uint64_t* array, int len, uint64_t hash, uint64_t hash_len, uint64_t qbits, uint64_t rbits) {
	int i;
	for (i = 0; i < len; i++)
		if (matchpart(array[i], hash, hash_len, qbits, rbits)) return 1;
	return -1;
}

void printbin(uint64_t val) {
	int i;
	for (i = 63; i >= 0; i--) {
		printf("%lu", (val >> i) % 2);
	}
	printf("\n");
}

ilist *find_ilist(uint64_t hash) {
	return NULL;
}

int main(int argc, char **argv)
{
	if (argc < 7) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF\nthe number of remainder bits in the QF\nthe universe size\nthe number of inserts\nthe number of queries\nthe number of trials\n");
		// eg. ./test 8 7 100000000 20000000 1000000 20
		// ./test 16 7 100000000 20000000 1000000 1 0
    // ./test 16 7 -1 $((1 << 15)) 1000000 1 0
		exit(1);
	}
	if (argc >= 8) {
		srand(strtol(argv[7], NULL, 10));
		printf("running test on seed %ld\n", strtol(argv[7], NULL, 10));
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
	uint64_t universe = strtoull(argv[3], NULL, 10);
	uint64_t num_inserts = strtoull(argv[4], NULL, 10);
	uint64_t num_queries = strtoull(argv[5], NULL, 10);
	uint64_t expected_fill = 1lu << (qbits + rbits) < num_inserts ? 1lu << (qbits + rbits) : num_inserts;
	if (universe < expected_fill) {
		printf("warning: universe may be too small to fill filter to completion\n");
	}
	else if (universe < expected_fill * 10) {
		printf("warning: filling the filter may take longer than necessary because of low universe size\n");
	}

	FILE *shalla = fopen("data/shalla.txt", "r");
	char buffer[1024];
	FILE *caida = fopen("data/20140619-140100.csv", "r");
	fgets(buffer, sizeof(buffer), caida);
	
	double avgInsTime = 0, avgInsPer = 0, avgQryTime = 0, avgQryPer = 0, avgFP = 0, avgFill = 0, maxFP = 0;
	double avgInsSlots = 0, avgQrySlots = 0;
	int numtrials = atoi(argv[6]);
	int trials = 0;
	for (; trials < numtrials; trials++) {
		QF qf;

		/* Initialise the CQF */
		/*if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {*/
		/*fprintf(stderr, "Can't allocate CQF.\n");*/
		/*abort();*/
		/*}*/
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
		uint64_t *ret_index = malloc(sizeof(uint64_t));
		uint64_t *ret_hash = malloc(sizeof(uint64_t));
		uint64_t *ret_other_hash = malloc(sizeof(uint64_t));
		int *ret_hash_len = malloc(sizeof(int));
		
		/*printf("insert returned %lu\n", qf_insert(&qf, 1, 0, 1, QF_NO_LOCK | QF_KEY_IS_HASH));
		printf("insert returned %lu\n", qf_insert(&qf, 2, 0, 1, QF_NO_LOCK | QF_KEY_IS_HASH));
		
		printf("query returned %lu\n", qf_query(&qf, 1, 0, ret_index, ret_hash, QF_KEY_IS_HASH));
		printf("query returned %lu\n", qf_query(&qf, 2, 0, ret_index, ret_hash, QF_KEY_IS_HASH));
		printf("query returned %lu\n", qf_query(&qf, 3, 0, ret_index, ret_hash, QF_KEY_IS_HASH));
		
		printf("adapt returned %d\n", qf_adapt(&qf, 1, 1, (1ULL << 32) | 1ULL, QF_KEY_IS_HASH));
		printf("adapt returned %d\n", qf_adapt(&qf, 0, 2, (1ULL << 32) | 2ULL, QF_KEY_IS_HASH));
		
		printf("query returned %lu\n", qf_query(&qf, 1, 0, ret_index, ret_hash, QF_KEY_IS_HASH));
		printf("query returned %lu\n", qf_query(&qf, (1ULL << 32) | 1ULL, 0, ret_index, ret_hash, QF_KEY_IS_HASH));
		printf("query returned %lu\n", qf_query(&qf, 2, 0, ret_index, ret_hash, QF_KEY_IS_HASH));
		printf("query returned %lu\n", qf_query(&qf, (1ULL << 32) | 2ULL, 0, ret_index, ret_hash, QF_KEY_IS_HASH));*/
		
		/*qf_insert_ret(&qf, 1, 0, 1, QF_KEY_IS_HASH, ret_index, ret_hash, ret_hash_len);
		printf("collided query expected, got %d\n", qf_insert_ret(&qf, 1 + qf.metadata->range, 0, 1, QF_KEY_IS_HASH | QF_NO_LOCK, ret_index, ret_hash, ret_hash_len));
		printf("extended length: %d\n", insert_and_extend(&qf, *ret_index, 1 + qf.metadata->range, 0, 1, 1, 0, QF_KEY_IS_HASH | QF_NO_LOCK));*/
		
		// testing extending on insert
		/*assert(qf_insert_ret(&qf, 1, 0, 1, QF_KEY_IS_HASH, ret_index, ret_hash, ret_hash_len)); // insert 1
		assert(qf_query(&qf, 1, 0, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(qf_query(&qf, 1 + (2 << (qbits + rbits)), 0, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_insert_ret(&qf, 1 + (1 << (qbits + rbits)), 0, 1, QF_KEY_IS_HASH, ret_index, ret_hash, ret_hash_len)); // try insert 11
		assert(insert_and_extend(&qf, *ret_index, 1 + (1 << (qbits + rbits)), 0, 1, 1, 0, ret_hash, ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0); // insert 11
		assert(qf_query(&qf, 1, 0, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(qf_query(&qf, 1 + (1 << (qbits + rbits)), 0, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(qf_query(&qf, 1 + (1lu << (qbits + 2*rbits)), 0, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_insert_ret(&qf, 1 + (1lu << (qbits + 2*rbits)), 0, 1, QF_KEY_IS_HASH, ret_index, ret_hash, ret_hash_len)); // try insert 101
		assert(insert_and_extend(&qf, *ret_index, 1 + (1lu << (qbits + 2*rbits)), 0, 1, 1, 0, ret_hash, ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK) > 0); // insert 101
		assert(!qf_query(&qf, 1 + (2 << (qbits + rbits)), 0, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(qf_insert_ret(&qf, 1 + (2 << (qbits + rbits)), 0, 1, QF_KEY_IS_HASH, ret_index, ret_hash, ret_hash_len)); // insert 21
		assert(qf_query(&qf, 1 + (2 << (qbits + rbits)), 0, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		abort();*/
		
		// testing extending on false positive query
		/*qf_insert_ret(&qf, 1, 1, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH);
		assert(qf_query(&qf, 1, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_query(&qf, 2, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_query(&qf, 1 + (1 << rbits), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(qf_query(&qf, 1 + (1 << (qbits + rbits)), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		
		assert(qf_adapt(&qf, *ret_index, 1, 1 + (1 << (qbits + rbits)), ret_hash, QF_KEY_IS_HASH) > 0);
		assert(qf_query(&qf, 1, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_query(&qf, 2, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_query(&qf, 1 + (1 << (qbits + rbits)), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_query(&qf, 1 + (2 << (qbits + rbits)), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(qf_query(&qf, 1 + (1 << (qbits + 2*rbits)), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		
		assert((*ret_hash_len = qf_adapt(&qf, *ret_index, 1, 1 + (1 << (qbits + 2*rbits)), ret_hash, QF_KEY_IS_HASH)) > 0);
		assert(qf_query(&qf, 1, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_query(&qf, 1 + (1 << (qbits + rbits)), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(!qf_query(&qf, 1 + (1 << (qbits + 2*rbits)), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		assert(qf_query(&qf, 1 + (1lu << (qbits + 3*rbits)), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
		abort();*/

    // testing insertion of duplicate items
    /*qf_insert_ret(&qf, 1, 2, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH);
    assert(qf_query(&qf, 1, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH) == 2);
    assert(!qf_insert_ret(&qf, 1, 3, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
    assert(insert_and_extend(&qf, *ret_index, 1, 3, 1, ret_hash, ret_other_hash, QF_KEY_IS_HASH) == 0);
    assert(qf_query(&qf, 1, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH) == 5);
    
    assert(!qf_insert_ret(&qf, 1 + (1 << (qbits + rbits)), 4, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH));
    assert(insert_and_extend(&qf, *ret_index, 1 + (1 << (qbits + rbits)), 4, 1, ret_hash, ret_other_hash, QF_KEY_IS_HASH));
    assert(qf_query(&qf, 1, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH) == 5);
    assert(qf_query(&qf, 1 + (1 << (qbits + rbits)), ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH) == 4);
    abort();*/
		
		if (0) {
			printf("%d\n", QF_COULDNT_LOCK);
			uint64_t q = 0, p = 0, m;
			for (m = 0; m < num_queries; m++) {
				q = (uint64_t)(rand() % universe);
				if (q % (1 << rbits) < (1 << qbits) && q % (1 << rbits) == 0) {
					p++;
				}
			}
			printf("false positive rate: %f\n", (double) p / num_queries);
			abort();
		}
		
		sglib_hashed_ilist_init(htab);
		ilist ii, *nn, *ll;
		ilist_iter it;
		
		if (universe < (double)qf.metadata->range * 0.8) {
			printf("warning: universe may be too small to fill filter to completion\n");
		}
		
		clock_t start_time, end_time;
		start_time = clock();
		
		//FILE *fp = fopen("timing.txt", "w");
		// PERFORM INSERTS
		printf("starting %lu inserts...\n", num_inserts);
		uint64_t i, j, k = 0, l = 0;
		for (i = 0; i < num_inserts;) {
			//uint64_t num_occupied_slots = qf_get_num_occupied_slots(&qf);
			/*if (qf_get_num_occupied_slots(&qf) >= 100000) {
				fclose(fp);
				bp2();
			}
			if (clock() - start_time > 500000) {
				start_time = clock();
				printf("%lu\t%lu\n", qf_get_num_occupied_slots(&qf), i);
			}*/
			/*if (qf_get_num_occupied_slots(&qf) == 97264) {
				bp2();
			}*/
			j = rand_uniform(universe); // pick a random number
			//j = rand_zipfian(1.01, 100000000);
			//fgets(buffer, sizeof(buffer), shalla);
			//j = hash_str(buffer);
			if (j == -1804289383) {
				bp2();
			}
			int ret = qf_insert_ret(&qf, j, 1, ret_index, ret_hash, ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH); // attempt to insert
      //if (*ret_hash == 3311657) bp2();
			if (ret == QF_NO_SPACE) { // if filter is full, stop
				printf("filter is full after qf_insert_ret\n");
				break;
			}
			else if (ret == 0) { // if a matching fingerprint is found, search hash table for original item
				if (0) {
					continue;
				}
				if (1) {
					//printf("wanted to extend, skipping for testing\n");
					k++;
					l += 3;
					i++;
					continue;
				}
				ii.rem = *ret_hash;
				ii.len = *ret_hash_len;
				ilist *item_in_table = sglib_hashed_ilist_find_member(htab, &ii);
				if (item_in_table == NULL) {
          printf("error:\tfilter claimed to have fingerprint %lu but hashtable could not find it\n", ii.rem);
          bp2();
        }
				else if (item_in_table->val == j) {
					//printf("log:\trandom insert %lu would have been duplicate, skipping\n", ii.rem);
				}
				else if (1){
					uint64_t new_rem, upd_rem;
					int ext_len = insert_and_extend(&qf, *ret_index, j, 1, item_in_table->val, &new_rem, &upd_rem, QF_KEY_IS_HASH | QF_NO_LOCK);
					if (ext_len == QF_NO_SPACE) {
						//printf("filter is full after insert_and_extend at %lu\n", num_occupied_slots);
						printf("filter is full after insert_and_extend\n");
						break;
					}
          sglib_hashed_ilist_delete(htab, item_in_table);
          item_in_table->rem = upd_rem;
          item_in_table->len = ext_len;
          sglib_hashed_ilist_add(htab, item_in_table);

					nn = malloc(sizeof(ilist));
					nn->val = j;
					nn->rem = new_rem;
					nn->len = ext_len;
					sglib_hashed_ilist_add(htab, nn);
					//printf("extended to length %d\n", insert_and_extend(&qf, *ret_index, j, 0, 1, findpart(vals, i - 1, *ret_hash, *ret_hash_len, qbits, rbits), 0, QF_KEY_IS_HASH | QF_NO_LOCK));
					k++;
					i++;
				}
			}
			else if (ret == 1) {
				nn = malloc(sizeof(ilist));
				nn->val = j;
				nn->rem = *ret_hash;
				nn->len = *ret_hash_len;
				sglib_hashed_ilist_add(htab, nn);
				/*val_mem[val_cnt].key = val_mem[val_cnt].val = j & BITMASK(rbits);
				val_mem[val_cnt].left = val_mem[val_cnt].right = NULL;
				values = insertItem(values, &(val_mem[val_cnt]));
				val_cnt++;*/
				//l++;
				i++;
			}
			/*else if (ret == -5) {
				continue;
			}*/
			else {
				printf("other error: errno %d\n", ret);
				break;
			}
			/*if (get_block(&qf, 23)->offset > 0) {
				bp2();
				//printf(" ");
			}*/
		}
		avgFill += (double)i / nslots;
		printf("nslots = %lu\n", qf_get_nslots(&qf));
		printf("made %lu inserts\n", i);
		printf("number of slots used after inserts: %lu\n", qf_get_num_occupied_slots(&qf));
		printf("extended %lu times\n", k);
		avgInsTime += (clock() - start_time);
    avgInsPer += (double)(clock() - start_time) / i;
		start_time = clock();
		avgInsSlots += qf_get_num_occupied_slots(&qf);
		
		// PERFORM QUERIES
		printf("starting %lu queries...\n", num_queries);
		int still_have_space = 1;
		if (qf_get_num_occupied_slots(&qf) >= qf.metadata->nslots * 0.95) {
			still_have_space = 0;
			printf("filter is full; skipping query adaptations\n");
		}
		for (i = 0; i < num_queries; i++) {
			//printf("%lu\n", i);
			
			j = rand_uniform(universe);
			//if (j == 6907711552940998528) bp2();
			//j = rand_zipfian(1.5f, (1lu << 30));
			//fgets(buffer, sizeof(buffer), shalla);
			//j = hash_str(buffer);
			/*fgets(buffer, sizeof(buffer), caida);
			csv_get(buffer, 3);
			j = hash_str(buffer);*/

      //if (j == 582125609) bp2();
			if (qf_query(&qf, j, ret_index, ret_hash, ret_hash_len, QF_KEY_IS_HASH)) {
				if (*ret_index / 64 == 489) frame(&qf, *ret_index / 64, *ret_index % 64, i);
				count_p++;
				ii.rem = *ret_hash;
				ii.len = *ret_hash_len;
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
							*ret_hash_len = qf_adapt(&qf, *ret_index, potential_item->val, j, ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
							if (*ret_hash_len > 0) {
								sglib_hashed_ilist_delete(htab, potential_item);
								potential_item->rem = *ret_hash;
								potential_item->len = *ret_hash_len;
								sglib_hashed_ilist_add(htab, potential_item);
							}
						}
					}
				}
			}
		}
		
		end_time = clock();
		printf("completed in time %ld us\n", end_time - start_time);
		printf("performed %lu queries uniformly from a universe of %lu\n", num_queries, universe);
		printf("number of slots used after queries: %lu\n", qf_get_num_occupied_slots(&qf));
		printf("false positive rate: %f\n", (double)count_fp / num_queries);
		
		free(ret_index);
		free(ret_hash);
		for(ll=sglib_hashed_ilist_it_init(&it,htab); ll!=NULL; ll=sglib_hashed_ilist_it_next(&it)) {
			free(ll);
		}
		avgQryTime += (end_time - start_time);
    avgQryPer += (double)(end_time - start_time) / num_queries;
		double fp_rate = (double)count_fp / num_queries;
		avgFP += fp_rate;
		if (fp_rate > maxFP) maxFP = fp_rate;
		avgQrySlots += qf_get_num_occupied_slots(&qf);
	}
	printf("\nperformed %d trials\n", numtrials);
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

