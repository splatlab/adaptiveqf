#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <fcntl.h>

#include <iostream>
#include <vector>
#include <openssl/rand.h>

using cuckoofilter::CuckooFilter;

void bp() {
	return;
}

int main(int argc, char **argv) {
	uint64_t seed;
	uint64_t i;
	size_t num_queries;
	size_t num_inserts;
	uint64_t qbits;
	uint64_t num_lines;

	char buffer[1000];
        operation *ops;

        printf("reading insert trace...\n");
        FILE *insert_trace_fp = fopen("sim_trace_insert.txt", "r");
        fgets(buffer, 1000, insert_trace_fp);
        num_lines = 0;
        while (fgets(buffer, 1000, insert_trace_fp) != NULL) num_lines++;
        ops = new operation[num_lines];

        rewind(insert_trace_fp);
        fscanf(insert_trace_fp, "%lu %lu %lu\n", &seed, &qbits, &num_inserts);
	uint64_t nslots = (1ull << qbits);
        srand(seed);
        for (i = 0; i < num_lines; i++) {
                fscanf(insert_trace_fp, "%c %lu %d\n", &ops[i].op, &ops[i].a, &ops[i].b);
        }
        fclose(insert_trace_fp);

	printf("initializing data structures...\n");
	CuckooFilter<size_t, 12> filter(nslots);

	printf("starting inserts...\n");
	clock_t start_clock = clock(), end_clock;
	struct timeval timecheck;
	gettimeofday(&timecheck, NULL);
	uint64_t start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec, end_time;

	for (i = 0; i <= num_lines; i++) {
		switch (ops[i].op) {
			case 'i':
				filter.AddFromRecording(ops, &i);
				break;
			default:
				break;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	printf("made %lu inserts\n", i);
	printf("time for inserts:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("insert throughput:    %f ops/sec\n", (double)i * 1000000 / (end_time - start_time));
	printf("cpu time for inserts: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	printf("cpu insert through:   %f ops/sec\n", (double)i * CLOCKS_PER_SEC / (end_clock - start_clock));

	free(ops);
	
	size_t fp_queries = 0;
	printf("reading query trace...\n");
        FILE *query_trace_fp = fopen("sim_trace_query.txt", "r");
        fgets(buffer, 1000, query_trace_fp);
        num_lines = 0;
        while (fgets(buffer, 1000, query_trace_fp) != NULL) num_lines++;
        ops = new operation[num_lines];

        rewind(query_trace_fp);
        fscanf(query_trace_fp, "%lu\n", &num_queries);
        for (i = 0; i < num_lines; i++) {
                fscanf(insert_trace_fp, "%c %lu %d\n", &ops[i].op, &ops[i].a, &ops[i].b);
        }
        fclose(query_trace_fp);

	printf("performing queries...\n");
	uint64_t location_data;
	start_clock = clock();
	gettimeofday(&timecheck, NULL);
	start_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	for (i = 0; i < num_lines; i++) {
		switch(ops[i].op) {
			case 'q':
				if (filter.ContainReturn(ops[i].a, &location_data) == cuckoofilter::Ok) {
					fp_queries++;
				}
				break;
			case 'a':
				if (filter.AdaptFromRecording(ops, &i) != cuckoofilter::Ok) printf("error: adapt failed to find previously queried item\n");
				break;
			default:
				break;
		}
	}
	gettimeofday(&timecheck, NULL);
	end_time = timecheck.tv_sec * 1000000 + timecheck.tv_usec;
	end_clock = clock();

	free(ops);

	printf("made %lu queries\n", num_queries);
	printf("time for queries:     %f sec\n", (double)(end_time - start_time) / 1000000);
	printf("query throughput:     %f ops/sec\n", (double)num_queries * 1000000 / (end_time - start_time));
	printf("cpu time for queries: %f sec\n", (double)(end_clock - start_clock) / CLOCKS_PER_SEC);
	printf("cpu query throughput: %f ops/sec\n", (double)num_queries * CLOCKS_PER_SEC / (end_clock - start_clock));

	printf("false positive rate:  %f%%\n", 100. * fp_queries / num_queries);

	return 0;
}
