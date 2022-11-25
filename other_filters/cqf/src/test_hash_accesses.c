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

uint64_t rand_uniform() {
	uint64_t a;
	RAND_bytes((void*)(&a), sizeof(uint64_t));
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

int main(int argc, char **argv)
{
	QF qf;
	uint64_t qbits = 23;
	uint64_t rbits = 8;
	uint64_t nhashbits = qbits + rbits;
	uint64_t nslots = (1ULL << qbits);
	uint64_t nvals = nslots * 0.9;
	uint64_t *vals;

	/* Initialise the CQF */
	/*if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {*/
	/*fprintf(stderr, "Can't allocate CQF.\n");*/
	/*abort();*/
	/*}*/
	if (!qf_malloc(&qf, nslots, nhashbits, 0, QF_HASH_INVERTIBLE, 0)) {
		fprintf(stderr, "Can't allocate CQF.\n");
		abort();
	}

	qf_set_auto_resize(&qf, true);

	/* Generate random values */
	//srand(time(NULL));
	srand(time(NULL));
	vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
	RAND_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	for (uint64_t i = 0; i < nvals; i++) {
		//vals[i] = rand_uniform();
		vals[i] = (1 * vals[i]) % qf.metadata->range;
		/*vals[i] = rand() % qf.metadata->range;*/
		/*fprintf(stdout, "%lx\n", vals[i]);*/
	}
	srand(time(NULL));

	clock_t start_time = clock(), end_time;
	/* Insert keys in the CQF */
	for (uint64_t i = 0; i < nvals; i++) {
		//int ret = qf_insert(&qf, vals[i], 0, key_count, QF_NO_LOCK);
		int ret = qf_insert(&qf, vals[i], 0, 1, QF_NO_LOCK);

		if (ret < 0) {
			fprintf(stderr, "failed insertion for key: %lx %d.\n", vals[i], 50);
			if (ret == QF_NO_SPACE)
				fprintf(stderr, "CQF is full.\n");
			else if (ret == QF_COULDNT_LOCK)
				fprintf(stderr, "TRY_ONCE_LOCK failed.\n");
			else
				fprintf(stderr, "Does not recognise return value.\n");
			abort();
		}
	}
	end_time = clock();
	printf("time per insert: %f\n", (double)(end_time - start_time) / nvals);

	/*uint64_t *query_set = calloc(nvals, sizeof(uint64_t));
	for (int i = 0; i < nvals; i++) {
		//uint64_t j = rand();
		//j <<= 32;
		//j |= rand();
		query_set[i] = rand_zipfian(1.5f, 1lu << 30);
	}*/

	FILE *fp = fopen("hash_accesses.txt", "w");
	fprintf(fp, "queries\taccesses\n");
	fclose(fp);
	fp = fopen("hash_accesses.txt", "a");

	start_time = clock();
	uint64_t j, accesses = 0, value;
	/* Lookup inserted keys and counts. */
	for (uint64_t i = 0; i < 10000000; i++) {
		//uint64_t count = qf_count_key_value(&qf, vals[i], 0, 0);
		/*j = rand();
		j <<= 32;
		j |= rand();*/
		//j = query_set[i];
		j = rand_zipfian(1.5f, 1ull << 30);
		if (qf_query(&qf, j, &value, QF_NO_LOCK)) accesses++;

		/*uint64_t count = qf_count_key_value(&qf, j, 0, 0);
		if (count > 0) {
			accesses++;
		}*/
		/*if (count < key_count) {
			fprintf(stderr, "failed lookup after insertion for %lx %ld.\n", vals[i],
							count);
			abort();
		}*/
		if (i % 100000 == 0) {
			fprintf(fp, "%lu\t%lu\n", i, accesses);
		}
	}
	end_time = clock();

	fclose(fp);
	printf("time per query: %f\n", (double)(end_time - start_time) / 10000000);
	printf("false positive rate: %f\n", (double)accesses / 10000000);

#if 0
	for (uint64_t i = 0; i < nvals; i++) {
		uint64_t count = qf_count_key_value(&qf, vals[i], 0, 0);
		if (count < key_count) {
			fprintf(stderr, "failed lookup during deletion for %lx %ld.\n", vals[i],
							count);
			abort();
		}
		if (count > 0) {
			/*fprintf(stdout, "deleting: %lx\n", vals[i]);*/
			qf_delete_key_value(&qf, vals[i], 0, QF_NO_LOCK);
			/*qf_dump(&qf);*/
			uint64_t cnt = qf_count_key_value(&qf, vals[i], 0, 0);
			if (cnt > 0) {
				fprintf(stderr, "failed lookup after deletion for %lx %ld.\n", vals[i],
								cnt);
				abort();
			}
		}
	}
#endif


}

