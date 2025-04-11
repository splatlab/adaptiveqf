#include <stdio.h>

#include "include/splinter_util.h"

void test_assert(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "Assertion failed: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

// This test shows the existence of a temporary false negative.
int main(int argc, char **argv) {
    size_t qbits = 20;
    size_t rbits = 7;

    size_t num_slots = 1ull << qbits;

    data_config data_cfg = qf_data_config_init();
    splinterdb_config splinterdb_cfg = qf_splinterdb_config_init("db", &data_cfg);
    remove(splinterdb_cfg.filename);
    splinterdb *db;
    if (splinterdb_create(&splinterdb_cfg, &db)) {
        fprintf(stderr, "Failed to create database\n");
        return -1;
    }
    splinterdb_lookup_result db_result;
    splinterdb_lookup_result_init(db, &db_result, 0, NULL);

    QF qf;
    if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
        fprintf(stderr, "Failed to allocate QF\n");
        return -1;
    }

    test_assert(qf_splinter_insert(&qf, db, 1, 1) == 1, "insert 1"); // Normal insert
    test_assert(qf_splinter_insert(&qf, db, 1 + (1ull << (qbits + rbits)), 1) == 2, "insert 2"); // Insert into existing minirun

    test_assert(qf_splinter_query_and_adapt(&qf, db, 1 + (1ull << (qbits + rbits))) > 0, "query 1"); // Query finds key 2
    test_assert(qf_splinter_query_and_adapt(&qf, db, 1) < 0, "query 2"); // Query misses key 1 (false negative)
    test_assert(qf_splinter_query_and_adapt(&qf, db, 1) > 0, "query 3"); // Query finds key 1 (because of adapt)

    splinterdb_close(&db);
    qf_free(&qf);
    printf("Test successful\n");
    return 0;
}