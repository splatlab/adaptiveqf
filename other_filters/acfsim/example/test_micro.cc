#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <vector>

using cuckoofilter::CuckooFilter;

uint64_t rand_uniform() {
	uint64_t a = rand();
	a += ((uint64_t)RAND_MAX + 1) * rand();
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

uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
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
}

#define HASH_SET_SEED 26571997

//linked list node for hashing with chaining
struct _set_node {
	struct _set_node *next;
	uint64_t value;
} typedef set_node;

int set_insert(set_node *set, int set_len, uint64_t key) {
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
			return 0;
		} else {
			return 1;
		}
	}
}

int main(int argc, char **argv) {
	size_t nslots = 1 << 24;
	size_t total_inserts = nslots * 0.9;
	size_t total_queries = 10000000;

	if (argc > 1) {
		FILE *fp = fopen(argv[1], "r");
		fscanf(fp, "%lu\n%lu\n", &nslots, &total_queries);
		total_inserts = nslots * 0.9;

		CuckooFilter<size_t, 12> filter(total_inserts);
		
		uint64_t *insert_set = new uint64_t[total_inserts];
		uint64_t *query_set = new uint64_t[total_queries];
		char *oper_set = new char[total_queries];

		for (size_t i = 0; i < total_inserts; i++) {
			fscanf(fp, "%lu\n", &insert_set[i]);
		}
		for (size_t i = 0; i < total_queries; i++) {
			fscanf(fp, "%c\t%lu\n", &oper_set[i], &query_set[i]);
		}

		uint64_t num_inserts = 0;
		clock_t start_time = clock(), end_time;
		for (size_t i = 0; i < total_inserts; i++, num_inserts++) {
			if (filter.Add(insert_set[i]) == cuckoofilter::Ok);
			/*if (num_inserts % 1000000 == 0) {
				printf("%lu\n", num_inserts);
			}*/
		}
		end_time = clock();
		printf("total insert time: %ld\n", end_time - start_time);
		printf("total num inserts: %lu\n", num_inserts);
		printf("insert throughput: %f\n\n", 1000000. * num_inserts / (end_time - start_time));

		uint64_t queries_made = 0;
		start_time = clock();
		for (size_t i = 0; i < total_queries; i++) {
			filter.Contain(query_set[i]);
			if (1 || oper_set[i] == '0') {
				if (filter.Contain(query_set[i]) != cuckoofilter::Ok);
			}
			else if (oper_set[i] == '1') {
				//if (filter.Contain(query_set[i]) != cuckoofilter::Ok) printf("error: unexpected query result\n");
				//filter.Adapt(filter.Adapt(query_set[i]));
			}
			queries_made++;
		}
		end_time = clock();
		printf("total query time:  %ld\n", end_time - start_time);
		printf("total num queries: %lu\n", total_queries);
		printf("query throughput:  %f\n", 1000000. * total_queries / (end_time - start_time));
		assert(queries_made == total_queries);

		fclose(fp);
		return 0;
	}

	// Create a cuckoo filter where each item is of type size_t and
	// use 12 bits for each item:
	//    CuckooFilter<size_t, 12> filter(total_items);
	// To enable semi-sorting, define the storage of cuckoo filter to be
	// PackedTable, accepting keys of size_t type and making 13 bits
	// for each key:
	//   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
	CuckooFilter<size_t, 12> filter(total_inserts);
	int set_len = total_inserts * 1.3;
	set_node *set = new set_node[set_len];

	FILE *fp = fopen("sim_recording.txt", "w");
	fprintf(fp, "%lu\n%lu\n", nslots, total_queries);

	// Insert items to this cuckoo filter
	size_t num_inserted = 0;
	for (size_t i = 0; i < total_inserts; i++, num_inserted++) {
		uint64_t key = rand_uniform();
		if (filter.Add(key) != cuckoofilter::Ok) {
			printf("fatal error: inserts terminated early\n");
			abort();
		}
		set_insert(set, set_len, key);
		fprintf(fp, "%lu\n", key);
	}

	/*uint64_t *query_set = (uint64_t*)calloc(total_items, sizeof(uint64_t));
	for (size_t i = 0; i < total_items; i++) {
		//query_set[i] = rand_zipfian(1.5f, 1lu << 30);
		query_set[i] = i + total_items;
	}*/

	// Check non-existing items, a few false positives expected
	size_t fp_queries = 0;
	/*for (size_t i = total_items; i < 2 * total_items; i++) {
	  if (filter.Contain(i) == cuckoofilter::Ok) {
	  false_queries++;
	  }
	  total_queries++;
	  }*/
	for (size_t i = 0; i < total_queries; i++) {
		//uint64_t key = rand_uniform();
		uint64_t key = rand_zipfian(1.5f, 1ull << 30);
		//uint64_t key = query_set[i];
		if (filter.Contain(key) == cuckoofilter::Ok) {
			if (!set_query(set, set_len, key)) {
				fp_queries++;
				filter.Adapt(key);
				fprintf(fp, "1\t%lu\n", key);
				continue;
			}
		}
		fprintf(fp, "0\t%lu\n", key);
	}
	fclose(fp);

	printf("done\n");

	return 0;
}
