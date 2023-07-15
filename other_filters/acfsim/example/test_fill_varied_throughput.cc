#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <vector>
#include <openssl/rand.h>

using cuckoofilter::CuckooFilter;

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
	a += ((uint64_t)RAND_MAX + 1) * rand();
	return a;
}

int zipfian_seed = 0;
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
			uint64_t result = newx;
                        return MurmurHash64A((void*)(&result), sizeof(result), zipfian_seed);
                }
                x = newx;
        }
}


#define HASH_SET_SEED 26571997

//linked list node for hashing with chaining
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
	//srand(time(NULL));
	uint64_t seed = 1669602171;
	seed = time(NULL);
	seed = 0;
	printf("running on seed %lu\n", seed);
	srand(seed);
	zipfian_seed = rand();

	size_t nslots = 1 << 24;
	size_t total_inserts = nslots * 0.9;

	printf("initializing filter...\n");
	// Create a cuckoo filter where each item is of type size_t and
	// use 12 bits for each item:
	//    CuckooFilter<size_t, 12> filter(total_items);
	// To enable semi-sorting, define the storage of cuckoo filter to be
	// PackedTable, accepting keys of size_t type and making 13 bits
	// for each key:
	//   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
	CuckooFilter<size_t, 12> filter(total_inserts);

	printf("initializing hash table...\n");
	int set_len = total_inserts * 1.3;
	set_node *set = new set_node[set_len];
	for (int i = 0; i < set_len; i++) {
		set[i].value = 0;
		set[i].next = NULL;
	}

	printf("generating inserts...\n");
	uint64_t *inserts = new uint64_t[total_inserts];
	/*for (size_t i = 0; i < total_inserts; i++) {
		inserts[i] = rand_uniform();
	}*/
	RAND_bytes((unsigned char*)inserts, total_inserts * sizeof(uint64_t));

	printf("starting inserts...\n");
	// Insert items to this cuckoo filter
	clock_t start_time = clock();
	size_t num_inserted = 0;

	double measurement_interval = 0.05f;
	double current_measurement_point = measurement_interval;
	uint64_t next_measurement = current_measurement_point * nslots;
	start_time = clock();
	for (size_t i = 0; i < total_inserts; i++, num_inserted++) {
		//uint64_t key = rand_uniform();
		if (filter.Add(inserts[i]) != cuckoofilter::Ok) {
			filter.Add(inserts[i - 1]);
			break;
		}
		set_insert(set, set_len, inserts[i]);
		if (i > next_measurement) {
			printf("%f\t%f\n", current_measurement_point, 1000000.0f * (measurement_interval * nslots) / (clock() - start_time));
			current_measurement_point += measurement_interval;
			next_measurement = current_measurement_point * nslots;
			start_time = clock();
		}
	}
	printf("%f\t%f\n", current_measurement_point, 1000000.0f * (measurement_interval * nslots) / (clock() - start_time));
	printf("made %lu inserts\n", num_inserted);

	return 0;
}
