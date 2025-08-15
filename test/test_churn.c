#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/ll_table.h"

void try_assert(int condition, const char *message) {
    if (!condition) {
        fprintf(stderr, "Assertion failed: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

void unit_test_ll_table() {
    printf("Running unit tests for ll_table...\n");

    ll_table inv_map;
    ll_table_init(&inv_map, 1024);
    
    // Inserts
    ll_table_insert(&inv_map, 1, 0, 100);
    ll_table_insert(&inv_map, 1, 1, 200);
    ll_table_insert(&inv_map, 2, 0, 300);
    ll_table_insert(&inv_map, 2, 0, 400); // Insert into front, should shift 300 to the front

    // Queries
    uint64_t *result = ll_table_query(&inv_map, 1, 0);
    try_assert(result != NULL && *result == 100, "Failed query for (1, 0)");

    result = ll_table_query(&inv_map, 1, 1);
    try_assert(result != NULL && *result == 200, "Failed query for (1, 1)");

    result = ll_table_query(&inv_map, 2, 0);
    try_assert(result != NULL && *result == 400, "Failed query for (2, 0)");

    result = ll_table_query(&inv_map, 2, 1);
    try_assert(result != NULL && *result == 300, "Failed query for (2, 1)");

    // Deletes
    ll_table_delete(&inv_map, 1, 0);
    result = ll_table_query(&inv_map, 1, 0);
    try_assert(result != NULL && *result == 200, "Failed to delete (1, 0)");
    result = ll_table_query(&inv_map, 1, 1);
    try_assert(result == NULL, "Failed to delete (1, 0)");

    ll_table_delete(&inv_map, 1, 0);
    result = ll_table_query(&inv_map, 1, 0);
    try_assert(result == NULL, "Failed to delete (1, 0)");

    // Clean up
    ll_free(&inv_map);
    
    printf("Test passed!\n");
}

void unit_test_qf_deletes() {
    printf("Running unit tests for QF deletes...\n");

    // Initialize a quotient filter
    uint64_t qbits = 16;
    uint64_t rbits = 7;
    uint64_t num_slots = 1ULL << qbits;
    QF qf;
    if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
        fprintf(stderr, "Failed to allocate quotient filter\n");
        return;
    }

    // Initialize the linked list table for inverse mapping
    ll_table inv_map;
    ll_table_init(&inv_map, num_slots);

    uint64_t insert_set[] = {
        (1ULL << rbits) + 1ULL, // Run 1
        (1ULL << rbits) + 1ULL + (1ULL << (qbits + rbits)), // Same minirun
        (1ULL << rbits) + 1ULL + (2ULL << (qbits + rbits)), // Same minirun
        (1ULL << rbits) + 2ULL, // Same run
        (1ULL << rbits) + 3ULL, // Same run
        (1ULL << rbits) + 3ULL + (1ULL << (qbits + rbits)), // Same minirun
        (2ULL << rbits) + 1ULL, // Run 2
        (2ULL << rbits) + 2ULL, // Same run
        (2ULL << rbits) + 2ULL + (1ULL << (qbits + rbits)), // Same minirun
        (4ULL << rbits) + 1ULL, // Run 4
        (5ULL << rbits) + 1ULL, // Run 5
        (5ULL << rbits) + 1ULL + (1ULL << (qbits + rbits)), // Same minirun
        (5ULL << rbits) + 2ULL, // Same run
        (5ULL << rbits) + 3ULL, // Same run
        (5ULL << rbits) + 4ULL, // Same run
    };
    // The canonical slots should be at
    // 0120450000000000
    // The runs should be at
    // 0111111222455555
    // The miniruns should be at
    // 0111233122111234
    for (size_t i = 0; i < sizeof(insert_set) / sizeof(insert_set[0]); i++) {
        qf_insert_result insert_result;
        int ret = qf_insert_using_ll_table(&qf, insert_set[i], 1, &insert_result, QF_NO_LOCK | QF_KEY_IS_HASH);
        try_assert(ret >= 0, "Failed to insert into quotient filter");
        try_assert(insert_result.minirun_id == (insert_set[i] & ((1ULL << (qbits + rbits)) - 1)), "Minirun ID mismatch");
        ll_table_insert(&inv_map, insert_result.minirun_id, 0, insert_set[i]);
    }

    uint64_t ret_hash;
    int minirun_rank = qf_query_using_ll_table(&qf, (1ULL << rbits) + 1ULL, &ret_hash, QF_KEY_IS_HASH);
    uint64_t *ll_result = ll_table_query(&inv_map, ret_hash & ((1ULL << (qbits + rbits)) - 1), minirun_rank);
    try_assert(ll_result != NULL && *ll_result == (1ULL << rbits) + 1ULL + (2ULL << (qbits + rbits)), "Failed query for minirun (1, 0)");
    int freed_slots = qf_remove_using_ll_table(&qf, *ll_result, minirun_rank, QF_KEY_IS_HASH);
    try_assert(freed_slots == 1, "Failed to delete first (1, 0) from QF");
    ll_table_delete(&inv_map, ret_hash & ((1ULL << (qbits + rbits)) - 1), minirun_rank);

    minirun_rank = qf_query_using_ll_table(&qf, (1ULL << rbits) + 1ULL, &ret_hash, QF_KEY_IS_HASH);
    ll_result = ll_table_query(&inv_map, ret_hash & ((1ULL << (qbits + rbits)) - 1), minirun_rank);
    try_assert(ll_result != NULL && *ll_result == (1ULL << rbits) + 1ULL + (1ULL << (qbits + rbits)), "Failed query for minirun (1, 0)");
    freed_slots = qf_remove_using_ll_table(&qf, *ll_result, minirun_rank, QF_KEY_IS_HASH);
    try_assert(freed_slots == 1, "Failed to delete second (1, 0) from QF");
    ll_table_delete(&inv_map, ret_hash & ((1ULL << (qbits + rbits)) - 1), minirun_rank);

    minirun_rank = qf_query_using_ll_table(&qf, (1ULL << rbits) + 1ULL, &ret_hash, QF_KEY_IS_HASH);
    ll_result = ll_table_query(&inv_map, ret_hash & ((1ULL << (qbits + rbits)) - 1), minirun_rank);
    try_assert(ll_result != NULL && *ll_result == (1ULL << rbits) + 1ULL, "Failed query for minirun (1, 0)");
    freed_slots = qf_remove_using_ll_table(&qf, *ll_result, minirun_rank, QF_KEY_IS_HASH);
    try_assert(freed_slots == 1, "Failed to delete third (1, 0) from QF");

    minirun_rank = qf_query_using_ll_table(&qf, (1ULL << rbits) + 1ULL, &ret_hash, QF_KEY_IS_HASH);
    try_assert(minirun_rank < 0, "Some of previous 3 deletions failed");

    qf_free(&qf);
    ll_free(&inv_map);

    printf("Test passed!\n");
}

int main(int argc, char **argv) {
    unit_test_ll_table();
    unit_test_qf_deletes();
    return 0;
}