#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <openssl/rand.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/gqf_file.h"
#include "include/hashutil.h"

#include <splinterdb/splinterdb.h>
#include <splinterdb/data.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/default_data_config.h>


#define TEST_DB_NAME "db"

#define Kilo (1024UL)
#define Mega (1024UL * Kilo)
#define Giga (1024UL * Mega)


#define MAX_KEY_SIZE 8
#define MAX_VAL_SIZE 8


int main(int argc, char **argv) {
	data_config data_cfg;
	default_data_config_init(MAX_KEY_SIZE, &data_cfg);
	splinterdb_config splinterdb_cfg = (splinterdb_config){
		.filename   = TEST_DB_NAME,
                .cache_size = 64 * Mega,
                .disk_size  = 127 * Mega,
                .data_cfg   = &data_cfg
	};

	splinterdb *database;
	int rc = splinterdb_create(&splinterdb_cfg, &database);
	assert(rc == 0);

	char key_val_data[MAX_KEY_SIZE + MAX_VAL_SIZE];
	RAND_bytes((void*)key_val_data, MAX_KEY_SIZE + MAX_VAL_SIZE);

	slice key = slice_create(MAX_KEY_SIZE, key_val_data);
	slice val = slice_create(MAX_VAL_SIZE, key_val_data + MAX_KEY_SIZE);

	rc = splinterdb_insert(database, key, val);
	assert(rc == 0);

	splinterdb_lookup_result result;
	splinterdb_lookup_result_init(database, &result, 0, NULL);

	rc = splinterdb_lookup(database, key, &result);
	assert(rc == 0);
	assert(splinterdb_lookup_found(&result));
	slice result_val;
	splinterdb_lookup_result_value(&result, &result_val);
	assert(slice_length(result_val) == MAX_VAL_SIZE);
	assert(strncmp(slice_data(result_val), key_val_data + MAX_KEY_SIZE, MAX_VAL_SIZE) == 0);

	rc = splinterdb_delete(database, key);
	assert(rc == 0);
	
	rc = splinterdb_lookup(database, key, &result);
	assert(rc == 0);
	assert(!splinterdb_lookup_found(&result));

	splinterdb_lookup_result_deinit(&result);
	splinterdb_close(&database);

	printf("Test completed\n");
}


