#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <vector>

#include <openssl/rand.h>

using cuckoofilter::CuckooFilter;

#define HASH_SET_SEED 26571997

uint64_t MurmurHash64A (const void *key, int len, unsigned int seed) {
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
}

uint64_t rand_uniform() {
	uint64_t a = rand();
	a *= (uint64_t)RAND_MAX + 1;
	a += rand();
	return a;
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

struct _set_node {
        struct _set_node *next;
        uint64_t value;
} typedef set_node;

int set_insert(set_node *set, int set_len, uint64_t key) {
        if (!key) key++;
        uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->value) {
                ptr->value = key;
        } else {
                while (ptr->next) {
                        if (ptr->value == key) {
                                return 0;
                        }
                        ptr = ptr->next;
                }
                if (ptr->value == key) {
                        return 0;
                }
                set_node *node = new set_node;
                ptr->next = node;

                node->next = NULL;
                node->value = key;
        }
        return 1;
}

//look up word in the set
//returns sources if it is found; returns 0 otherwise
int set_query(set_node *set, int set_len, uint64_t key) {
        if (!key) key++;
        uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->value) {
                return 0;
        } else {
                while (ptr->next){
                        if (ptr->value == key) {
                                return 1;
                        }
                        ptr = ptr->next;
                }
                if (ptr->value == key){
                        return 1;
                } else {
                        return 0;
                }
        }
}

int main(int argc, char **argv) {
  size_t total_items = 1 << 24;
  size_t total_inserts = total_items * 0.9;

  // Create a cuckoo filter where each item is of type size_t and
  // use 12 bits for each item:
  //    CuckooFilter<size_t, 12> filter(total_items);
  // To enable semi-sorting, define the storage of cuckoo filter to be
  // PackedTable, accepting keys of size_t type and making 13 bits
  // for each key:
  //   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
  CuckooFilter<size_t, 12> filter(total_inserts);
  size_t set_len = total_inserts * 1.3;
  set_node *set = new set_node[set_len];

	uint64_t *inserts = new uint64_t[total_inserts];
	for (size_t i = 0; i < total_inserts; i++) {
		inserts[i] = rand_uniform();
	}

  // Insert items to this cuckoo filter
  clock_t start_time = clock(), end_time;
  size_t num_inserted = 0;
  for (size_t i = 0; i < total_items; i++, num_inserted++) {
    if (filter.Add(rand_uniform()) != cuckoofilter::Ok) {
      break;
    }
    //set_insert(set, set_len, inserts[i]);
  }
  end_time = clock();
  printf("time per insert: %f\n", (double)(end_time - start_time) / num_inserted);

  uint64_t total_queries = 100000000ull;
  uint64_t *query_set = new uint64_t[total_queries];
  RAND_bytes((unsigned char*)query_set, total_items * sizeof(uint64_t));

  printf("starting queries...\n");

  // Check non-existing items, a few false positives expected
  start_time = clock();
  size_t false_queries = 0;
  for (size_t i = 0; i < total_queries; i++) {
    if (filter.Contain(query_set[i]) == cuckoofilter::Ok) {
      false_queries++;
    }
    total_queries++;
  }
  end_time = clock();
  printf("query throughput: %f\n", (double)(total_queries) * CLOCKS_PER_SEC / (end_time - start_time));

  // Output the measured false positive rate
  std::cout << "false positive rate is "
            << 100.0 * false_queries / total_queries << "%\n";

  return 0;
}
