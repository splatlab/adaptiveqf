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

uint64_t rand_uniform(uint64_t max) {
	if (max <= RAND_MAX) return rand() % max;
	uint64_t a = rand();
	a *= (uint64_t)RAND_MAX + 1;
	a += rand();
	if (max != -1) return a % max;
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

uint64_t hash_str(char *str) {
	uint64_t hash = 5381;
	int c;
	while ((c = *str++)) {
		hash = ((hash << 5) + hash) + c;
	}
	return hash;
}

void csv_get_col(char* buffer, int col) {
        int i, j;
        for (i = 0; buffer[i] != '\0' && col > 0; i++) {
                if (buffer[i] == ',') col--;
        }
        for (j = 0; buffer[i + j] != '\0' && buffer[i + j] != ','; j++) {
                buffer[j] = buffer[i + j];
        }
        buffer[j] = '\0';
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


uint64_t query_set_size = 2900000;
uint64_t step_size = 20000;
uint64_t step_sample = 1000000;

int main(int argc, char **argv)
{
	if (argc < 6) {
		fprintf(stderr, "Please specify\nnumber of quotient bits\nnumber of remainder bits\nload factor\nnumber of queries\nnumber of trials\n");
		exit(1);
	}
	if (argc >= 7) {
		srand(strtol(argv[6], NULL, 10));
		printf("running test on seed %ld\n", strtol(argv[6], NULL, 10));
	}
	else {
		time_t seed = time(NULL);
		printf("running test on seed %ld\n", seed);
		srand(seed);
	}
	
	uint64_t qbits = atoi(argv[1]);
	uint64_t rbits = atoi(argv[2]);
	uint64_t nhashbits = qbits + rbits;
	uint64_t nslots = (1ULL << qbits);
	double load_factor = atof(argv[3]);
	uint64_t num_inserts = nslots * load_factor;
	uint64_t num_queries = strtoull(argv[4], NULL, 10);
	

	char buffer[256];
	//FILE *shalla = fopen("data/shalla.txt", "r");
	//FILE *caida = fopen("data/20140619-140100.csv", "r");
	//fgets(buffer, sizeof(buffer), caida);
	
	FILE *outfp = fopen("data/progress.csv", "w");
	sprintf(buffer, "queries,false positive rate\n");
	fputs(buffer, outfp);
	fclose(outfp);
	outfp = fopen("data/progress.csv", "a");
	
	uint64_t set_len = num_inserts * 1.3;
	set_node *set = calloc(set_len, sizeof(set_node));

	uint64_t *query_set = calloc(query_set_size, sizeof(uint64_t));
	for (int q = 0; q < query_set_size; q++) {
		/*fgets(buffer, sizeof(buffer), shalla);
		query_set[q] = hash_str(buffer);*/
		//fgets(buffer, sizeof(buffer), caida);
		//csv_get_col(buffer, 3);
		//query_set[q] = hash_str(buffer);
	}

	double *fp_rate_set = calloc(num_queries / step_size, sizeof(double));

	double avgInsTime = 0, avgInsPer = 0, avgQryTime = 0, avgQryPer = 0, avgFP = 0, avgFill = 0, maxFP = 0;
	double avgInsSlots = 0, avgQrySlots = 0;
	int num_trials = atoi(argv[5]);
	int trials;
	for (trials = 0; trials < num_trials; trials++) {
		QF qf;

		if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0,
										 "mycqf.file")) {
			fprintf(stderr, "Can't allocate CQF.\n");
			abort();
		}

		qf_set_auto_resize(&qf, false);
		
		uint64_t count_fp = 0;
		
		uint64_t ret_index, ret_hash, ret_other_hash;
		int ret_hash_len;
		
		sglib_hashed_ilist_init(htab);
		ilist ii, *nn, *ll;
		ilist_iter it;
		
		clock_t start_time = clock(), end_time;
		
		// PERFORM INSERTS
		printf("starting inserts...\n");
		uint64_t i, j, k = 0;
		for (i = 0; i < num_inserts; i++) {
			j = rand_uniform(INT_MAX);
			//fgets(buffer, sizeof(buffer), caida);
			//csv_get_col(buffer, 3);
			//j = hash_str(buffer);

			int ret = qf_insert_ret(&qf, j, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH); // attempt to insert
			if (ret == QF_NO_SPACE) {
				printf("filter is full after %lu inserts\n", i);
				break;
			}
			else if (ret == 0) { // if a matching fingerprint is found, search hash table for original item
				ii.rem = ret_hash;
				ii.len = ret_hash_len;
				ilist *item_in_table = sglib_hashed_ilist_find_member(htab, &ii);
				if (item_in_table != NULL) {
					int ext_len = insert_and_extend(&qf, ret_index, j, 1, item_in_table->val, &ret_hash, &ret_other_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
					if (ext_len == QF_NO_SPACE) {
						printf("filter is full after %lu inserts\n", i);
						break;
					}
					else if (ext_len == 0) {
						continue;
					}

					sglib_hashed_ilist_delete(htab, item_in_table);
					item_in_table->rem = ret_other_hash;
					item_in_table->len = ext_len;
					sglib_hashed_ilist_add(htab, item_in_table);

					nn = malloc(sizeof(ilist));
					nn->val = j;
					nn->rem = ret_hash;
					nn->len = ext_len;
					sglib_hashed_ilist_add(htab, nn);
					k++;
				}
			}
			else if (ret == 1) {
				nn = malloc(sizeof(ilist));
				nn->val = j;
				nn->rem = ret_hash;
				nn->len = ret_hash_len;
				sglib_hashed_ilist_add(htab, nn);
			}
			else {
				printf("other error: errno %d\n", ret);
				break;
			}
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

		double fpr;
		count_fp = 0;
		for (i = 0; i < num_queries; i++) {
			j = rand_uniform(INT_MAX);
			//j = rand_zipfian(1.5, 1lu << 30);
			//j = query_set[rand() % query_set_size];

			if (qf_query(&qf, j, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH)) {
				ii.rem = ret_hash;
				ii.len = ret_hash_len;
				ilist *potential_item = sglib_hashed_ilist_find_member(htab, &ii);
				if (potential_item != NULL && potential_item->val != j) {
					count_fp++;
					if (still_have_space) {
						if (qf_get_num_occupied_slots(&qf) >= qf.metadata->nslots * 0.95) {
							still_have_space = 0;
							printf("out of space to adapt after %lu queries\n", i);
						}
						else {
							ret_hash_len = qf_adapt(&qf, ret_index, potential_item->val, j, &ret_hash, QF_KEY_IS_HASH | QF_NO_LOCK);
							if (ret_hash_len > 0) {
								sglib_hashed_ilist_delete(htab, potential_item);
								potential_item->rem = ret_hash;
								potential_item->len = ret_hash_len;
								sglib_hashed_ilist_add(htab, potential_item);
							}
						}
					}
				}
			}
			if (i % step_size == 0) {
				fpr = qf_get_num_occupied_slots(&qf);
				fpr = (i == 0 ? 0 : fpr / num_inserts);
				fp_rate_set[i / step_size] += fpr;
				//printf("%f\n", fpr);
			}
		}

		
		end_time = clock();
		printf("completed in time %ld us\n", end_time - start_time);
		printf("performed %lu queries\n", num_queries);
		printf("number of slots used after queries: %lu\n", qf_get_num_occupied_slots(&qf));
		printf("false positive rate: %f\n", (double)count_fp / num_queries);
		
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

	for (uint64_t q = 0; q < num_queries / step_size; q++) {
		printf("%f\n", fp_rate_set[q]);
		//fp_rate_set[q] /= num_trials;
		sprintf(buffer, "%f\n", fp_rate_set[q]);
		fputs(buffer, outfp);
	}
	//fclose(shalla);
	//fclose(caida);
	fclose(outfp);

	printf("\nperformed %d trials\n", num_trials);
	printf("avg false positive rate: %f\n", avgFP / num_trials);
	printf("max false positive rate: %f\n", maxFP);
	printf("avg fill rate: %f\n", avgFill / num_trials);
	printf("avg slots used after inserts: %f\n", avgInsSlots / num_trials);
	printf("avg slots used after queries: %f\n", avgQrySlots / num_trials);
	printf("avg total insert time: %f\n", avgInsTime / num_trials);
	printf("avg insert time per item: %f\n", avgInsPer / num_trials);
	printf("avg total query time: %f\n", avgQryTime / num_trials);
	printf("avg query time per item: %f\n", avgQryPer / num_trials);
}

