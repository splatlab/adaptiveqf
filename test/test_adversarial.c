#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
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
		fprintf(stderr, "Please specify \nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries [eg. 100000000]\nthe adversarial frequency [eg. 100]\n");
		exit(1);
	}
	if (argc > 5) { // optional seed
		srand(strtol(argv[5], NULL, 10));
		printf("Running test on seed %ld\n", strtol(argv[5], NULL, 10));
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
	size_t adv_freq = strtoull(argv[4], NULL, 10);
	size_t cache_size = 512;

	test_results_t ret;
	char buffer[100];

	uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
	RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));
	uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
	RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));

	//ret = run_adversarial_test(qbits, rbits, insert_set, num_inserts, query_set, num_queries, 1000, 10000, 1, "adv_q.csv");
	//if (ret.exit_code) printf("Test failed\n");

	printf("Running adversarial test with attacks every %lu queries\n", adv_freq);

	mkdir("logs", 0777);
	sprintf(buffer, "logs/adv_q_%lu_%lu_%lu.csv", qbits, cache_size, adv_freq);

	ret = run_adversarial_test(qbits, rbits, insert_set, num_inserts, query_set, num_queries, cache_size, adv_freq, 10000, 1, buffer);
	if (ret.exit_code) {
		printf("Test failed with exit code %d\n", ret.exit_code);
	}

	return 0;
}
