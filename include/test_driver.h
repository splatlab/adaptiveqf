#ifndef _TEST_DRIVER_H_
#define _TEST_DRIVER_H_

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
#include "include/rand_util.h"
#include "include/splinter_util.h"

#include <splinterdb/data.h>
#include <splinterdb/default_data_config.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/public_util.h>
#include <splinterdb/splinterdb.h>


struct _test_results_t {
	int exit_code;
	double insert_throughput;
	double query_throughput;
	double final_query_throughput;
	double false_positive_rate;
} typedef test_results_t;

void warm_up_filter(const QF *qf, uint64_t num_warmup_queries);

test_results_t run_throughput_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, int verbose, char *inserts_outfile, char *queries_outfile);

test_results_t run_split_throughput_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, int verbose, char *inserts_outfile, char *queries_outfile);

test_results_t run_adversarial_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, uint64_t cache_size, uint64_t adv_freq, uint64_t adv_set_max_size, int verbose, char *queries_outfile);

test_results_t run_micro_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, int verbose);

test_results_t run_parallel_splinter_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, size_t num_threads);

test_results_t run_parallel_test(size_t qbits, size_t rbits, uint64_t *insert_set, size_t insert_set_len, uint64_t *query_set, size_t query_set_len, size_t num_threads);

#endif
