#include "ll_table.h"
#include "hashutil.h"
#include <assert.h>

void ll_table_init(ll_table *table, uint64_t size) {
        table->size = size;
        table->num_keys = 0;
        table->num_families = 0;
        table->buckets = calloc(size, sizeof(ll_list*));
	table->seed = rand();
}

void ll_free(ll_table *table) {
	for (int i = 0; i < table->size; i++) {
		ll_list *list_ptr = table->buckets[i], *list_next = NULL;
		while (list_ptr != NULL) {
			list_next = list_ptr->next;
			ll_node *ptr = list_ptr->head, *next = NULL;
			while (ptr != NULL) {
				next = ptr->next;
				free(ptr);
				ptr = next;
			}
			free(list_ptr);
			list_ptr = list_next;
		}
		table->buckets[i] = NULL;
	}
	free(table->buckets);
	table->buckets = NULL;
}

void ll_table_insert(ll_table *table, uint64_t family, uint64_t rank, uint64_t key) {
	uint64_t index = MurmurHash64A((void*)(&family), sizeof(uint64_t), table->seed) % table->size;

	ll_list *family_list = NULL;
	if (table->buckets[index] == NULL) {
		family_list = table->buckets[index] = malloc(sizeof(ll_list));
		family_list->family = family;
		family_list->head = NULL;
		family_list->next = NULL;
	}
	else if (table->buckets[index]->family >= family) {
		if (table->buckets[index]->family == family) family_list = table->buckets[index];
		else {
			family_list = malloc(sizeof(ll_list));
			family_list->family = family;
			family_list->next = table->buckets[index];
			table->buckets[index] = family_list;
			table->num_families++;
		}
	}
	else {
		ll_list *ptr = table->buckets[index];
		while (ptr->next != NULL) {
			if (ptr->next->family >= family) {
				if (ptr->next->family == family) family_list = ptr->next;
				else {
					family_list = malloc(sizeof(ll_list));
					family_list->family = family;
					family_list->next = ptr->next;
					ptr->next = family_list;
					table->num_families++;
				}
				break;
			}
			ptr = ptr->next;
		}
		if (family_list == NULL) {
			family_list = malloc(sizeof(ll_list));
			family_list->family = family;
			ptr->next = family_list;
			table->num_families++;
		}
	}

	ll_node *new_node = malloc(sizeof(ll_node));
	new_node->key = key;
	if (rank == 0) {
		new_node->next = family_list->head;
		family_list->head = new_node;
	}
	else {
		ll_node *ptr = family_list->head;
		for (int i = 1; i < rank; i++) ptr = ptr->next;
		new_node->next = ptr->next;
		ptr->next = new_node;
	}
	assert(family_list->family == family);

	table->num_keys++;
}

uint64_t *ll_table_query(ll_table *table, uint64_t family, uint64_t rank) {
	uint64_t index = MurmurHash64A((void*)(&family), sizeof(uint64_t), table->seed) % table->size;
	
	ll_list *list_ptr = table->buckets[index];
	while (list_ptr != NULL) {
		if (list_ptr->family >= family) {
			if (list_ptr->family == family) {
				ll_node *ptr = list_ptr->head;
				while (rank > 0) {
					ptr = ptr->next;
					rank--;
				}
				return &(ptr->key);
			}
			else return NULL;
		}
	}
	return NULL;
}


