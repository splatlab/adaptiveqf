/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

extern "C" {
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <openssl/rand.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"
#include "include/hashutil.h"
}

#include <ctime>
#include <stxxl/unordered_map>
#include <stxxl/map>

int bp2() {
	return 0;
}

uint64_t rand_uniform(uint64_t max) {
	/*if (max <= RAND_MAX) return rand() % max;
	uint64_t a = rand();
	uint64_t b = rand();
	a |= (b << 32);
	return a % max;*/
	uint64_t a;
	RAND_bytes((unsigned char*)(&a), sizeof(uint64_t));
	if (max) a %= max;
	return a;
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

#define HASH_SET_SEED 26571997

struct _set_node {
        struct _set_node *next;
        uint64_t key;
	uint64_t rem;
	int len;
} typedef set_node;

int set_insert(set_node *set, int set_len, set_node *new_node) {
        uint64_t hash = MurmurHash64A((void*)(&(new_node->key)), sizeof(new_node->key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->key) {
                ptr->key = new_node->key;
        } else {
                while (ptr->next) {
                        if (ptr->key == new_node->key) {
                                return 0;
                        }
                        ptr = ptr->next;
                }
                if (ptr->key == new_node->key) {
                        return 0;
                }
                ptr->next = new_node;
        }
        return 1;
}

int set_query(set_node *set, int set_len, uint64_t key) {
        uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->key) {
                return 0;
        } else {
                while (ptr->next){
                        if (ptr->key == key) {
                                return 1;
                        }
                        ptr = ptr->next;
                }
                if (ptr->key == key){
                        return 1;
                } else {
                        return 0;
                }
        }
}

set_node *set_get(set_node *set, int set_len, set_node *node) {
        uint64_t hash = MurmurHash64A((void*)(&(node->key)), sizeof(node->key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->key) {
                return NULL;
        } else {
                while (ptr->next){
                        if (ptr->key == node->key) {
                                return ptr;
                        }
                        ptr = ptr->next;
                }
                if (ptr->key == node->key){
                        return ptr;
                } else {
                        return NULL;
                }
        }
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


typedef struct test_struct {
	int a;
	int b;
} test_struct;

#define USE_UNORDERED_MAP 1
#if USE_UNORDERED_MAP
#define BACKING_MAP_T unordered_map_t
#define BACKING_MAP_INSERT(X, Y, Z) X.insert(std::make_pair(Y, Z));
#else
#define BACKING_MAP_T ordered_map_t
#define BACKING_MAP_INSERT(X, Y, Z) X.insert(pair_t(Y, Z));
#endif

#define SUB_BLOCK_SIZE 256
#define SUB_BLOCKS_PER_BLOCK 256

#define DATA_NODE_BLOCK_SIZE (4096)
#define DATA_LEAF_BLOCK_SIZE (4096)

//! [hash]
struct HashFunctor
{
	size_t operator () (int key) const {
		return (size_t)(key * 2654435761u);
	}
};
//! [hash]

//! [comparator]
struct CompareGreater {
	bool operator () (const uint64_t& a, const uint64_t& b) const {
		return a > b;
	}
	static uint64_t max_value() {
		return 0ULL;
	}
};
//! [comparator]

typedef stxxl::map<uint64_t, test_struct, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> map_t2;

typedef stxxl::unordered_map<uint64_t, uint64_t, HashFunctor, CompareGreater, SUB_BLOCK_SIZE, SUB_BLOCKS_PER_BLOCK> unordered_map_t;
typedef stxxl::map<uint64_t, uint64_t, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> ordered_map_t;
typedef std::pair<uint64_t, uint64_t> pair_t;

int map_inserts = 0;
int map_queries = 0;

int insert_key(QF *qf, BACKING_MAP_T& map, uint64_t key, int count) {
        uint64_t ret_index, ret_hash, ret_other_hash;
        int ret_hash_len;
        int ret = qf_insert_ret(qf, key, count, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
        if (ret == QF_NO_SPACE) {
                return 0;
        }
        else if (ret == 0) {
		BACKING_MAP_T::iterator item = map.find(ret_hash | (1 << ret_hash_len));
		map_queries++;
                if (item == map.end()) {
                        printf("error:\tfilter claimed to have fingerprint %lu but hashtable could not find it\n", ret_hash);
                }
                else if (item->second == key) {
                        insert_and_extend(qf, ret_index, key, count, key, &ret_hash, &ret_other_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
                }
                else {
                        int ext_len = insert_and_extend(qf, ret_index, key, count, item->second, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
                        if (ext_len == QF_NO_SPACE) {
                                printf("filter is full after insert_and_extend\n");
                                return 0;
                        }

			map.erase(item);
			BACKING_MAP_INSERT(map, ret_other_hash | (1 << ret_hash_len), item->second);
			BACKING_MAP_INSERT(map, ret_hash | (1 << ret_hash_len), key);
			map_inserts++;
			map_inserts++;
                }
        }
        else if (ret == 1) {
		BACKING_MAP_INSERT(map, ret_hash | (1 << ret_hash_len), key);
		map_inserts++;
        }
        else {
                printf("other error: errno %d\n", ret);
                return 0;
        }
        return 1;
}


int main(int argc, char **argv)
{
	if (argc < 4) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF\nthe number of remainder bits in the QF\nthe number of queries\n");
		// ./test 16 7 $((1 << 15)) 1000000 0
		exit(1);
	}
	if (argc >= 5) {
		srand(strtol(argv[4], NULL, 10));
		printf("running test on seed %ld\n", strtol(argv[5], NULL, 10));
	}
	else {
		time_t seed = time(NULL);
		printf("running test on seed %ld\n", seed);
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
	double load_factor = 0.95f;
	uint64_t num_inserts = nslots * load_factor;//strtoull(argv[3], NULL, 10);
	uint64_t num_queries = strtoull(argv[3], NULL, 10);


	printf("initializing hash table...\n");
#if USE_UNORDERED_MAP
	unordered_map_t backing_map;
	//backing_map.max_buffer_size(SUB_BLOCK_SIZE);
	unordered_map_t::allocator_type allocator = backing_map.get_allocator();
#else
	ordered_map_t backing_map((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
#endif	
	ordered_map_t database((ordered_map_t::node_block_type::raw_size) * 3, (ordered_map_t::leaf_block_type::raw_size) * 3);
	

	printf("initializing filter...\n");
	QF qf;
	if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
		fprintf(stderr, "Can't allocate CQF.\n");
		abort();
	}
	qf_set_auto_resize(&qf, false);


	printf("generating insert set of size %lu...\n", num_inserts);
	uint64_t i, j;
	uint64_t *insert_set = new uint64_t[num_inserts];
        RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));
	uint64_t *query_set = (uint64_t*)calloc(num_queries, sizeof(uint64_t));


	// PERFORM INSERTS
	printf("performing inserts...\n");
	uint64_t ret_index, ret_hash;
	int ret_hash_len;

	uint64_t target_fill = nslots * load_factor;

	double measure_interval = 0.05f;
        double current_interval = measure_interval;
        uint64_t measure_point = nslots * current_interval, last_point = 0;
	int still_have_space = 1;

	printf("CLOCKS_PER_SEC: %ld\n", CLOCKS_PER_SEC);
	clock_t start_time = clock(), end_time, interval_time = start_time;
	for (i = 0; current_interval != 1./*qf.metadata->noccupied_slots <= target_fill*/; i++) {
		if (!insert_key(&qf, backing_map, insert_set[i], 1)) break;
		//database.insert(pair_t(insert_set[i], i));
		if (qf.metadata->noccupied_slots >= measure_point) {
                        printf("throughput for interval %f: \t%f\n", current_interval, (double)(i - last_point) * CLOCKS_PER_SEC / (clock() - interval_time));
			printf("map inserts: %d\tmap_queries: %d\n", map_inserts, map_queries);

			RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));

			uint64_t k;
			uint64_t fp_count = 0;
			clock_t query_start_time = clock();
			for (k = 0; k < num_queries; k++) {
				j = query_set[k];

				if (qf_query(&qf, j, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
					ordered_map_t::iterator item = database.find(j);
					if (item == database.end()) {
						fp_count++;
						if (0 && still_have_space) {
							BACKING_MAP_T::iterator orig_key = backing_map.find(ret_hash | (1 << ret_hash_len));
							ret_hash_len = qf_adapt(&qf, ret_index, orig_key->second, j, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
							if (ret_hash_len > 0) {
								backing_map.erase(orig_key);
								backing_map.insert(std::make_pair(ret_hash | (1 << ret_hash_len), orig_key->second));
							}
							else if (ret_hash_len == QF_NO_SPACE) {
								still_have_space = 0;
								printf("filter is full after %lu queries\n", i);
							}
						}
					}
				}
			}
			printf("query throughput at %f fill: \t%f\n", current_interval, (double)(num_queries) * CLOCKS_PER_SEC / (clock() - query_start_time));
			printf("false positive rate: \t%f\n\n", (double)(fp_count) / num_queries);

                        current_interval += measure_interval;
                        last_point = i;
                        measure_point = nslots * current_interval;
			interval_time = clock();
                }
	}
	end_time = clock();

	printf("time for inserts:      %f\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
	printf("avg insert throughput: %f ops/sec\n", (double)i * CLOCKS_PER_SEC / (end_time - start_time));
}

