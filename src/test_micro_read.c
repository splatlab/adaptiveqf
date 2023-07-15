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

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"
#include "include/hashutil.h"

struct operation {
	char op;
	uint64_t a;
	uint64_t b;
	uint64_t c;
} typedef operation;


int main(int argc, char **argv)
{
	uint64_t i, qbits, rbits, seed, num_inserts, num_queries, num_lines;
	char buffer[1000];
	operation *ops;

	printf("reading insert trace...\n");
	FILE *insert_trace_fp = fopen("sim_trace_insert.txt", "r");
	fgets(buffer, 1000, insert_trace_fp);
	num_lines = 0;
	while (fgets(buffer, 1000, insert_trace_fp) != NULL) num_lines++;
	ops = malloc(num_lines * sizeof(operation));

	rewind(insert_trace_fp);
	fscanf(insert_trace_fp, "%lu %lu %lu %lu\n", &seed, &qbits, &rbits, &num_inserts);
	srand(seed);
	uint64_t nhashbits = qbits + rbits;
	uint64_t nslots = (1ull << qbits);
	for (i = 0; i < num_lines; i++) {
		fscanf(insert_trace_fp, "%c %lu %lu %lu\n", &ops[i].op, &ops[i].a, &ops[i].b, &ops[i].c);
	}
	fclose(insert_trace_fp);

        printf("initializing filter...\n");
        QF qf;
        if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
                fprintf(stderr, "Can't allocate CQF.\n");
                abort();
        }
        qf_set_auto_resize(&qf, false);

	printf("performing inserts...\n");
	uint64_t ret_index, ret_hash, ret_other_hash;
	int ret_hash_len;

	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;
	for (i = 0; i < num_lines; i++) {
		switch (ops[i].op) {
			case 'i':
				qf_insert_ret(&qf, ops[i].a, 1, &ret_index, &ret_hash, &ret_hash_len, QF_NO_LOCK | QF_KEY_IS_HASH);
				break;
			case 'e':
				insert_and_extend(&qf, ops[i].a, ops[i].b, 1, ops[i].c, &ret_hash, &ret_other_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
				break;
			default:
				break;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("time for inserts:      %f\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughput:     %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
	printf("cpu time for inserts:  %f\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);

	free(ops);

	// PERFORM QUERIES
	printf("reading query trace...\n");
	FILE *query_trace_fp = fopen("sim_trace_query.txt", "r");
	fgets(buffer, 1000, query_trace_fp);
	num_lines = 0;
	while (fgets(buffer, 1000, query_trace_fp) != NULL) num_lines++;
	ops = malloc(num_lines * sizeof(operation));

	rewind(query_trace_fp);
	fscanf(query_trace_fp, "%lu\n", &num_queries);
	for (i = 0; i < num_lines; i++) {
		fscanf(insert_trace_fp, "%c %lu %lu %lu\n", &ops[i].op, &ops[i].a, &ops[i].b, &ops[i].c);
	}
	fclose(query_trace_fp);

	start_clock = clock();
	gettimeofday(&timecheck, NULL);
	start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	for (i = 0; i < num_lines; i++) {
		switch(ops[i].op) {
			case 'q':
				qf_query(&qf, ops[i].a, &ret_index, &ret_hash, &ret_hash_len, QF_KEY_IS_HASH);
				break;
			case 'a':
				qf_adapt(&qf, ops[i].a, ops[i].b, ops[i].c, &ret_hash, QF_NO_LOCK | QF_KEY_IS_HASH);
				break;
			default:
				break;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("time for queries:     %f s\n", (double)(end_time - start_time) / 1000000);
	printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
	printf("cpu time for queries: %f s\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
}

