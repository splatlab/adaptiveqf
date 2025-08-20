#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>

#include <openssl/rand.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/ll_table.h"
#include "include/hashutil.h"
#include "include/rand_util.h"

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

void breakpoint() {

}

#define ZIPFIAN_COEFF 1.5f
#define ZIPFIAN_MAX 10000000ULL

int main(int argc, char **argv) {
    // Unit tests for correctness
    // unit_test_ll_table();
    // unit_test_qf_deletes();

    // Read command line arguments
    int seeded;
    char mode = 'z';
    if (argc < 4) {
        fprintf(stderr, "Please specify\nthe log of the number of slots in the QF [eg. 20]\nthe number of remainder bits in the QF [eg. 9]\nthe number of queries [eg. 1000000]\n(optional) mode: z (default) - zipfian; u - uniform\n(optional) seed\n");
        return 1;
    }
    if (argc > 4) {
        mode = argv[4][0];
        if (mode != 'u' && mode != 'z') {
            fprintf(stderr, "Invalid mode specified. Defaulting to zipfian.\n");
            mode = 'z';
        }
    }
    if (argc > 5) {
        srand(strtol(argv[5], NULL, 10));
        printf("Running test on seed %ld\n", strtol(argv[5], NULL, 10));
        seeded = 1;
    }
    else {
        srand(time(NULL));
    }

    size_t qbits = atoi(argv[1]);
    size_t rbits = atoi(argv[2]);
    size_t num_slots = 1ull << qbits;
    size_t num_inserts = num_slots * 0.9f;
    size_t max_capacity = num_slots * 0.95f;
    size_t num_queries = strtoull(argv[3], NULL, 10);

    // Generate input data
    uint64_t *insert_set = malloc(num_inserts * sizeof(uint64_t));
    if (seeded) for (size_t i = 0; i < num_inserts; i++) insert_set[i] = rand_uniform(-1ull);
    else RAND_bytes((unsigned char*)insert_set, num_inserts * sizeof(uint64_t));

    uint64_t *query_set = malloc(num_queries * sizeof(uint64_t));
    if (seeded) for (size_t i = 0; i < num_queries; i++) query_set[i] = rand_uniform(-1ull);
    else RAND_bytes((unsigned char*)query_set, num_queries * sizeof(uint64_t));
    if (mode == 'z') {
        int hash_seed = rand();
        for (size_t i = 0; i < num_queries; i++) {
            query_set[i] = rand_zipfian(ZIPFIAN_COEFF, ZIPFIAN_MAX, query_set[i]);
            query_set[i] = MurmurHash64A(&query_set[i], sizeof(uint64_t), hash_seed);
        }
    }

    // Initialize data structures
    QF qf;
    if (!qf_malloc(&qf, num_slots, qbits + rbits, 0, QF_HASH_INVERTIBLE, 0)) {
        fprintf(stderr, "Failed to allocate memory for QF\n");
        return 1;
    }

    ll_table inv_map;
    ll_table_init(&inv_map, num_slots);

    // Perform inserts
    size_t i;
    for (i = 0; i < num_inserts; i++) {
        // printf("Inserting key %lu (%zu/%zu)\n", insert_set[i], i + 1, num_inserts);
        if (insert_set[i] == 16374950167423305250ULL) {
            printf("Found key\n");
            breakpoint();
        }
        qf_insert_result insert_result;
        int ret = qf_insert_using_ll_table(&qf, insert_set[i], 1, &insert_result, QF_KEY_IS_HASH | QF_NO_LOCK);
        if (ret < 0) break;
        ll_table_insert(&inv_map, insert_result.minirun_id, 0, insert_set[i]);

        // for (uint64_t j = 0; j < num_slots; j++) {
        //     ll_list *list_ptr = inv_map.buckets[j];
        //     while (list_ptr != NULL) {
        //         ll_node *node_ptr = list_ptr->head;
        //         while (node_ptr != NULL) {
        //             try_assert((node_ptr->key & ((1ULL << (qbits + rbits)) - 1)) == list_ptr->family, "Inverse map key does not match minirun ID");
        //             try_assert((MurmurHash64A(&list_ptr->family, sizeof(uint64_t), inv_map.seed) % num_slots) == j, "Inverse map key does not hash to correct bucket");
        //             node_ptr = node_ptr->next;
        //         }
        //         list_ptr = list_ptr->next;
        //     }
        // }
    }
    num_inserts = i;
    printf("Number of inserts performed: %zu\n", num_inserts);
    // snapshot(&qf);
    // for (i = 0; i < num_slots; i++) {
    //     fprintf(stderr, "\rTesting correctness of inverse map for slot %zu... ", i);
    //     ll_list *list_ptr = inv_map.buckets[i];
    //     while (list_ptr != NULL) {
    //         ll_node *node_ptr = list_ptr->head;
    //         while (node_ptr != NULL) {
    //             try_assert((node_ptr->key & ((1ULL << (qbits + rbits)) - 1)) == list_ptr->family, "Inverse map key does not match minirun ID");
    //             try_assert((MurmurHash64A(&list_ptr->family, sizeof(uint64_t), inv_map.seed) % num_slots) == i, "Inverse map key does not hash to correct bucket");
    //             node_ptr = node_ptr->next;
    //         }
    //         list_ptr = list_ptr->next;
    //     }
    // }
    fprintf(stderr, "\n");

    int inserted_test_key = 0;

    // Perform queries
    size_t num_churns = 10;
    size_t curr_churn = 1;
    size_t next_churn = num_queries * curr_churn / num_churns;
    double churn_factor = 0.1;
    size_t updates_per_churn = num_inserts * churn_factor;
    size_t false_positives = 0;
    for (i = 0; i < num_queries; i++) {
        // printf("Performing query %zu\n", i);
        // if (i == 64) {
        //     breakpoint();
        // }
        size_t ret_hash;
        int minirun_rank = qf_query_using_ll_table(&qf, query_set[i], &ret_hash, QF_KEY_IS_HASH);
        if (minirun_rank >= 0) {
            uint64_t *ll_result = ll_table_query(&inv_map, ret_hash & ((1ULL << (qbits + rbits)) - 1), minirun_rank);
            try_assert(ll_result != NULL, "Mismatch detected between filter and inverse map");
            if (*ll_result != query_set[i]) {
                false_positives++;
                if (qf.metadata->noccupied_slots < max_capacity) qf_adapt_using_ll_table(&qf, *ll_result, query_set[i], minirun_rank, QF_KEY_IS_HASH | QF_NO_LOCK);
            }
        }

        // Performing a churn
        if (i == next_churn) {
            printf("Churn %zu/%zu at query %zu: performing %zu deletions and inserts\n", curr_churn, num_churns, i, updates_per_churn);
            for (size_t j = 0; j < updates_per_churn; j++) {
                // if (inv_map.buckets[79] != NULL && inv_map.buckets[79]->head == NULL) {
                //     breakpoint();
                // }

                fprintf(stderr, "\rChurn delete %zu/%zu", j + 1, updates_per_churn);

                // Find a random key to delete
                size_t delete_index = rand() % num_slots;
                while (inv_map.buckets[delete_index] == NULL) {
                    delete_index = (delete_index + 1) % num_slots;
                }
                uint64_t delete_key = inv_map.buckets[delete_index]->head->key;
                
                if (i == 7 && j == 6) breakpoint();
                if (i == 7 && j == 49) {
                    printf("Key: %lu\n", delete_key);
                    breakpoint();
                }

                minirun_rank = qf_query_using_ll_table(&qf, delete_key, &ret_hash, QF_KEY_IS_HASH);
                try_assert(minirun_rank >= 0, "Item in inverse map not found in QF");
                // try_assert((ret_hash & ((1ULL << (qbits + rbits)) - 1)) == (delete_key & ((1ULL << (qbits + rbits)) - 1)), "Hashes don't match"); // Only should hold if QF_KEY_IS_HASH
                uint64_t *ll_result = ll_table_query(&inv_map, ret_hash & ((1ULL << (qbits + rbits)) - 1), minirun_rank);
                try_assert(ll_result != NULL, "Item in QF not found in inverse map");
                int freed_slots = qf_remove_using_ll_table(&qf, *ll_result, minirun_rank, QF_KEY_IS_HASH);
                try_assert(freed_slots > 0, "Failed to delete key from QF");
                
                ll_table_delete(&inv_map, ret_hash & ((1ULL << (qbits + rbits)) - 1), minirun_rank);

                try_assert(!inserted_test_key || (qf_query_using_ll_table(&qf, 16374950167423305250ULL, &ret_hash, QF_KEY_IS_HASH) >= 0), "Test key lost");

                // for (uint64_t k = 0; k < num_slots; k++) {
                //     fprintf(stderr, "\rTesting correctness of inverse map for slot %zu... ", k);
                //     ll_list *list_ptr = inv_map.buckets[k];
                //     while (list_ptr != NULL) {
                //         ll_node *node_ptr = list_ptr->head;
                //         while (node_ptr != NULL) {
                //             try_assert((node_ptr->key & ((1ULL << (qbits + rbits)) - 1)) == list_ptr->family, "Inverse map key does not match minirun ID");
                //             try_assert((MurmurHash64A(&list_ptr->family, sizeof(uint64_t), inv_map.seed) % num_slots) == k, "Inverse map key does not hash to correct bucket");
                //             node_ptr = node_ptr->next;
                //         }
                //         list_ptr = list_ptr->next;
                //     }
                // }
            }
            // fprintf(stderr, "\n");

            if (seeded) for (size_t j = 0; j < updates_per_churn; j++) insert_set[j] = rand_uniform(-1ULL);
            else RAND_bytes((unsigned char*)insert_set, updates_per_churn * sizeof(uint64_t));
            for (size_t j = 0; j < updates_per_churn; j++) {
                if (insert_set[j] == 16374950167423305250ULL) {
                    printf("Found key\n");
                    inserted_test_key = 1;
                    breakpoint();
                }

                qf_insert_result insert_result;
                int ret = qf_insert_using_ll_table(&qf, insert_set[j], 1, &insert_result, QF_KEY_IS_HASH | QF_NO_LOCK);
                if (ret < 0) {
                    fprintf(stderr, "\nFailed to insert key %zu into QF\n", insert_set[j]);
                    break;
                }
                ll_table_insert(&inv_map, insert_result.minirun_id, 0, insert_set[j]);

                try_assert(!inserted_test_key || (qf_query_using_ll_table(&qf, 16374950167423305250ULL, &ret_hash, QF_KEY_IS_HASH) >= 0), "Test key lost");

                // for (uint64_t k = 0; k < num_slots; k++) {
                //     fprintf(stderr, "\rTesting correctness of inverse map for slot %zu after insert %zu... ", k, j);
                //     ll_list *list_ptr = inv_map.buckets[k];
                //     while (list_ptr != NULL) {
                //         ll_node *node_ptr = list_ptr->head;
                //         while (node_ptr != NULL) {
                //             try_assert((node_ptr->key & ((1ULL << (qbits + rbits)) - 1)) == list_ptr->family, "Inverse map key does not match minirun ID");
                //             try_assert((MurmurHash64A(&list_ptr->family, sizeof(uint64_t), inv_map.seed) % num_slots) == k, "Inverse map key does not hash to correct bucket");
                //             node_ptr = node_ptr->next;
                //         }
                //         list_ptr = list_ptr->next;
                //     }
                // }
            }

            curr_churn++;
            next_churn = num_queries * curr_churn / num_churns;
        }
    }
    printf("False positive rate: %f%%\n", (false_positives / (double)num_queries) * 100);

    return 0;
}