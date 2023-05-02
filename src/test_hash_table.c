#include <stdio.h>
#include <assert.h>

#include "hash_table.c"

uint64_t hash_str(char *str) {
        uint64_t hash = 5381;
        int c;
        while ((c = *str++)) {
                hash = ((hash << 5) + hash) + c;
        }
        return hash;
}

void csv_get(char* buffer, int col) {
        int i, j, k;
        for (i = 0; buffer[i] != '\0' && col > 0; i++) {
                if (buffer[i] == ',') col--;
        }
	int found_one_comma = 0;
        for (j = k = 0; buffer[i + k] != '\0' && (buffer[i + k] != ',' || !found_one_comma); k++) {
		if (buffer[i + k] == '"') continue;
		if (buffer[i + k] == ',') found_one_comma = 1;
                buffer[j++] = buffer[i + k];
        }
        buffer[j] = '\0';
}

int main() {
	int num_inserts = 12000000;
	int hash_table_size = num_inserts * 1.3;

	ht_node *hash_table = calloc(hash_table_size, sizeof(ht_node));

	uint64_t key = rand();
	hash_table_insert(hash_table, hash_table_size, key);

	assert(hash_table_query(hash_table, hash_table_size, key));
	assert(!hash_table_query(hash_table, hash_table_size, key + 1));

	FILE *fp = fopen("data/20140619-125911.csv", "r");
	if (fp == NULL) {
		printf("unable to open data file\n");
		abort();
	}

	char buffer[2048];
	fgets(buffer, sizeof(buffer) / sizeof(char), fp);

	uint64_t num_entries = 0;
	uint64_t num_duplicates = 0;
	while (fgets(buffer, sizeof(buffer) / sizeof(char), fp) != NULL) {
		num_entries++;
		csv_get(buffer, 2);
		uint64_t key = hash_str(buffer);
		if (hash_table_query(hash_table, hash_table_size, key)) {
			num_duplicates++;
			continue;
		}
		//printf("%s\n", buffer);
		hash_table_insert(hash_table, hash_table_size, key);
	}

	printf("total entries: %lu\n", num_entries);
	printf("duplicates: %lu\n", num_duplicates);
	printf("unique entries: %lu\n", num_entries - num_duplicates);
	printf("unique proportion: %f\n", (double)(num_entries - num_duplicates) / num_entries);

	return 0;
}
