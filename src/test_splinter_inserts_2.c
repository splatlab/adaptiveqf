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


#define MAX_KEY_SIZE 16
#define MAX_VAL_SIZE 16


void pad_data(void *dest, const void *src, const size_t dest_len, const size_t src_len) {
        assert(dest_len >= src_len);
        memcpy(dest, src, src_len);
        bzero(dest + src_len, dest_len - src_len);
}

slice padded_slice(const void *data, const size_t dest_len, const size_t src_len, void *buffer) {
        pad_data(buffer, data, dest_len, src_len);

        return slice_create(dest_len, data);
}

int db_init(splinterdb **database, const char *db_name, uint64_t cache_size, uint64_t disk_size) {
        data_config data_cfg;
        default_data_config_init(MAX_KEY_SIZE, &data_cfg);
        splinterdb_config splinterdb_cfg = (splinterdb_config){
                //.filename   = db_name,
                .filename   = "db",
                .cache_size = cache_size,
                .disk_size  = disk_size,
                .data_cfg   = &data_cfg
        };

        return splinterdb_create(&splinterdb_cfg, database);
}

int db_insert(splinterdb *database, const void *key_data, const size_t key_len, const void *val_data, const size_t val_len) {
        char key_padded[MAX_KEY_SIZE];
        char val_padded[MAX_VAL_SIZE];
        pad_data(key_padded, key_data, MAX_KEY_SIZE, key_len);
        pad_data(val_padded, val_data, MAX_VAL_SIZE, val_len);
        slice key = slice_create(MAX_KEY_SIZE, key_padded);
        slice val = slice_create(MAX_VAL_SIZE, val_padded);

        return splinterdb_insert(database, key, val);
}


splinterdb_lookup_result bm_result;
int main(int argc, char **argv) {
	uint64_t num_inserts = 0;
	//uint64_t size_estimate = 0;
	if (argc < 2) {
		printf("Provide number of insertions\n");
		exit(0);
	}
	else {
		num_inserts = strtoull(argv[1], NULL, 10);
	}

	data_config data_cfg;
        default_data_config_init(MAX_KEY_SIZE, &data_cfg);
        splinterdb_config splinterdb_cfg = (splinterdb_config){
                .filename   = TEST_DB_NAME,
                .cache_size = 64 * Mega,
                .disk_size  = 8 * Giga,
                .data_cfg   = &data_cfg
        };

        splinterdb *backing_map;
        if (splinterdb_create(&splinterdb_cfg, &backing_map) != 0) {
                printf("Error creating database\n");
                exit(0);
        }

	/*splinterdb *backing_map;
        if (db_init(&backing_map, "db", 64 * Mega, 8 * Giga) != 0) {
                printf("error initializing database\n");
                exit(0);
        }
        //splinterdb_lookup_result_init(backing_map, &bm_result, 0, NULL);*/

	char *key_val_data = malloc((MAX_KEY_SIZE + MAX_VAL_SIZE) * num_inserts);
	RAND_bytes((void*)key_val_data, (MAX_KEY_SIZE + MAX_VAL_SIZE) * num_inserts);

	clock_t start_time = clock();
	for (uint64_t i = 0; i < num_inserts; i++) {
		uint64_t j = i * (MAX_KEY_SIZE + MAX_VAL_SIZE);
		slice key = slice_create(MAX_KEY_SIZE, key_val_data + j);
		slice val = slice_create(MAX_VAL_SIZE, key_val_data + j + MAX_KEY_SIZE);

		if (splinterdb_insert(backing_map, key, val) != 0) {
			printf("\nFailed insertion %lu", i);
			break;
		}
		fprintf(stderr, "\r%lu/%lu", i + 1, num_inserts);
	}
	clock_t end_time = clock();
	fprintf(stderr, "\n");

	printf("Test completed in %f sec\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
	printf("Average throughput: %f ops/sec\n", (double)num_inserts * CLOCKS_PER_SEC / (end_time - start_time));

	free(key_val_data);
}


