#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <openssl/rand.h>

#include "include/hashutil.h"
#include "include/rand_util.h"
#include "include/splinter_util.h"
#include "include/test_driver.h"


int main(int argc, char **argv)
{
	if (argc < 5) {
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries [eg. 100000000]\nthe number of threads [eg. 16]\n");
		exit(1);
	}
	if (argc > 5) { // optional seed
		srand(strtol(argv[5], NULL, 10));
		printf("Running test on seed %ld\n", strtol(argv[4], NULL, 10));
	}
	else {
		time_t seed = time(NULL);
		printf("Running test on seed %ld\n", seed);
		srand(seed);
	}
	size_t qbits = atoi(argv[1]);
	size_t rbits = atoi(argv[2]);

	size_t num_inserts = (1ull << qbits) * 0.9f;//strtoull(argv[3], NULL, 10);
	size_t num_queries = strtoull(argv[3], NULL, 10);

	size_t num_threads = strtoull(argv[4], NULL, 10);

	test_results_t ret;

	
	uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
	RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));
	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));

	ret = run_parallel_test(qbits, rbits, insert_set, num_inserts, query_set, num_queries, num_threads);
	if (ret.exit_code) {
		printf("Test failed\n");
	}

	/*for (int i = 0; i < num_inserts; i++) {
		query_set[i] = rand_zipfian(1.5f, 10000000ull, query_set[i]);
		query_set[i] = MurmurHash64A(&query_set[i], sizeof(query_set[i]), murmur_seed);
	}
	
	ret = run_throughput_test(qbits, rbits, insert_set, num_inserts, query_set, num_queries, 1, "zipf_i.csv", "zipf_q.csv");
	if (ret.exit_code) {
		printf("Test failed\n");
	}*/

	free(insert_set);
	free(query_set);

	return 0;
}
