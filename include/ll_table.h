#ifndef _LL_TABLE_H_
#define _LL_TABLE_H_

#include <sys/types.h>
#include <stdlib.h>
#include <stdint.h>

// a node in a linked list
// intended to store a key that was inserted into an AQF
struct ll_node {
	uint64_t key;
	struct ll_node *next;
} typedef ll_node;

// a linked list of keys whose fingerprints in the AQF share the same quotient and remainder
// family is the combination of quotient and remainder shared by all items in this list
// head is the start of the actual linked list
// next is a pointer to the next linked list (these linked lists are also stored in a linked list in a hash table)
struct ll_list {
	uint64_t family;
	ll_node *head;
	struct ll_list *next;
} typedef ll_list;

// a hash table mapping from signatures to linked lists of items sharing a family
//
struct ll_table {
	ll_list **buckets;
	unsigned int seed;
	uint64_t size;
	uint64_t num_keys;
	uint64_t num_families;
} typedef ll_table;

void ll_table_init(ll_table *table, uint64_t size);

void ll_free(ll_table *table);

void ll_table_insert(ll_table *table, uint64_t family, uint64_t rank, uint64_t key);

uint64_t *ll_table_query(ll_table *table, uint64_t family, uint64_t rank);

#endif
