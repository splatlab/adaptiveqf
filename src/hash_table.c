#include "../include/hash_table.h"


uint64_t _murmur_64(const void *key, int len, uint64_t seed) {
        const uint64_t m = 0xc6a4a7935bd1e995;
        const int r = 47;

        uint64_t h = seed ^ (len * m);

        const uint64_t * data = (const uint64_t *)key;
        const uint64_t * end = data + (len/8);

        while(data != end)
        {
                uint64_t k = *data++;

                k *= m;
                k ^= k >> r;
                k *= m;

                h ^= k;
                h *= m;
        }

        const unsigned char * data2 = (const unsigned char*)data;

        switch(len & 7)
        {
                case 7: h ^= (uint64_t)data2[6] << 48;
                case 6: h ^= (uint64_t)data2[5] << 40;
                case 5: h ^= (uint64_t)data2[4] << 32;
                case 4: h ^= (uint64_t)data2[3] << 24;
                case 3: h ^= (uint64_t)data2[2] << 16;
                case 2: h ^= (uint64_t)data2[1] << 8;
                case 1: h ^= (uint64_t)data2[0];
                                                h *= m;
        };

        h ^= h >> r;
        h *= m;
        h ^= h >> r;

        return h;
}


#define HASH_SEED 26571997

int hash_table_insert(ht_node *ht, int ht_len, uint64_t key) {
        if (!key) key++;
        uint64_t hash = _murmur_64((void*)(&key), sizeof(uint64_t), HASH_SEED);
        ht_node *ptr = &ht[hash % ht_len];
        if (!ptr->value) {
                ptr->value = key;
        } else {
                while (ptr->next) {
                        if (ptr->value == key) {
                                return 0;
                        }
                        ptr = ptr->next;
                }
                if (ptr->value == key) {
                        return 0;
                }
                ht_node *node = malloc(sizeof(ht_node));
                ptr->next = node;

                node->next = NULL;
                node->value = key;
        }
        return 1;
}

int hash_table_query(ht_node *ht, int ht_len, uint64_t key) {
        if (!key) key++;
        uint64_t hash = _murmur_64((void*)(&key), sizeof(uint64_t), HASH_SEED);
        ht_node *ptr = &ht[hash % ht_len];
        if (!ptr->value) {
                return 0;
        } else {
                while (ptr->next){
                        if (ptr->value == key) {
                                return 1;
                        }
                        ptr = ptr->next;
                }
                if (ptr->value == key){
                        return 1;
                } else {
                        return 0;
                }
        }
}
