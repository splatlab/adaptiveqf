#ifndef CUCKOO_FILTER_CUCKOO_FILTER_H_
#define CUCKOO_FILTER_CUCKOO_FILTER_H_

#include <assert.h>
#include <algorithm>

#include "debug.h"
#include "hashutil.h"
#include "packedtable.h"
#include "printutil.h"
#include "singletable.h"
#include "mirroredtable.h"

extern "C" {
#include <splinterdb/splinterdb.h>
#include <splinterdb/data.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/default_data_config.h>
}

uint64_t MurmurHash64A(const void* key, int len, uint64_t seed) {
        const uint64_t m = 0xc6a4a7935bd1e995LLU;
        const int r = 47;

        uint64_t h = seed ^ (len * m);

        const uint64_t* data = (const uint64_t*)key;
        const uint64_t* end = (len >> 3) + data;

        while (data != end) {
                uint64_t k = *data++;

                k *= m;
                k ^= k >> r;
                k *= m;

                h ^= k;
                h *= m;
        }

        const unsigned char * data2 = (const unsigned char*)data;

        switch (len & 7) {
                case 7: h ^= (uint64_t)(data2[6]) << 48;
                case 6: h ^= (uint64_t)(data2[5]) << 40;
                case 5: h ^= (uint64_t)(data2[4]) << 32;
                case 4: h ^= (uint64_t)(data2[3]) << 24;
                case 3: h ^= (uint64_t)(data2[2]) << 16;
                case 2: h ^= (uint64_t)(data2[1]) << 8;
                case 1: h ^= (uint64_t)(data2[0]);
                        h *= m;
        };

        h ^= h >> r;
        h *= m;
        h ^= h >> r;

        return h;
}

#define HASH_SET_SEED 26571997
struct _set_node {
        struct _set_node *next;
        uint64_t key;
        uint64_t value;
} typedef set_node;

int set_insert(set_node *set, int set_len, uint64_t key, uint64_t value) {
        if (!key) key++;
        uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->key) {
                ptr->key = key;
                ptr->value = value;
        }
        else {
                while (ptr->next) {
                        if (ptr->key == key) return 0;
                        ptr = ptr->next;
                }
                if (ptr->key == key) return 0;
                set_node *node = new set_node;
                ptr->next = node;

                node->next = NULL;
                node->key = key;
                node->value = value;
        }
        return 1;
}

int set_query(set_node *set, int set_len, uint64_t key, uint64_t *value) {
        if (!key) key++;
        uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->key) {
                return 0;
        }
        else {
                while (ptr->next) {
                        if (ptr->key == key) {
                                *value = ptr->value;
                                return 1;
                        }
                        ptr = ptr->next;
                }
                if (ptr->key == key){
                        *value = ptr->value;
                        return 1;
                }
                else return 0;
        }
}

int set_delete(set_node *set, int set_len, uint64_t key) {
        if (!key) key++;
        uint64_t hash = MurmurHash64A((void*)(&key), sizeof(key), HASH_SET_SEED);
        set_node *ptr = &set[hash % set_len];
        if (!ptr->key) {
                return 0;
        }
        else if (ptr->key == key) {
                if (ptr->next) {
                        set_node *temp = ptr->next;
                        ptr->key = ptr->next->key;
                        ptr->value = ptr->next->value;
                        ptr->next = ptr->next->next;
                        free(temp);
                }
                else {
                        ptr->key = 0;
                }
                return 1;
        }
        else if (!ptr->next) {
                return 0;
        }
        else {
                do {
                        if (ptr->next->key == key) {
                                set_node *temp = ptr->next;
                                ptr->next = ptr->next->next;
                                free(temp);
                                return 1;
                        }
                        ptr = ptr->next;
                } while (ptr->next);
                return 0;
        }
}

struct operation {
        char op;
        uint64_t a;
        int b;
} typedef operation;

namespace cuckoofilter {
// status returned by a cuckoo filter operation
	enum Status {
		Ok = 0,
		NotFound = 1,
		NotEnoughSpace = 2,
		NotSupported = 3,
	};

	// maximum number of cuckoo kicks before claiming failure
	const size_t kMaxCuckooCount = 500;

	// A cuckoo filter class exposes a Bloomier filter interface,
	// providing methods of Add, Delete, Contain. It takes three
	// template parameters:
	//   ItemType:  the type of item you want to insert
	//   bits_per_item: how many bits each item is hashed into
	//   TableType: the storage of table, SingleTable by default, and
	// PackedTable to enable semi-sorting
	template <typename ItemType, size_t bits_per_item,
		 template <size_t> class TableType = MirroredTable,
		 typename HashFamily = TwoIndependentMultiplyShift>
			 class CuckooFilter {
				 // Storage of items
				 TableType<bits_per_item> *table_;
				 TableType<64> *remote;

				 // Number of items stored
				 size_t num_items_;

				 typedef struct {
					 size_t index;
					 uint32_t tag;
					 bool used;
				 } VictimCache;

				 VictimCache victim_;

				 HashFamily hasher_;
				 HashFamily alt_hasher_;
				 char *hash_sels;
				 HashFamily *hash_fcns;

				 inline size_t GenerateIndexHash(const ItemType item) const {
					 size_t index = hasher_(item);
					 return index % table_->NumBuckets();
				 }
				 
				 inline size_t GenerateAltIndexHash(const ItemType item, const uint64_t first_index) const {
					 size_t index = alt_hasher_(item) % table_->NumBuckets();
					 if (index == first_index) index = (index ? index - 1 : table_->NumBuckets() - 1);
					 return index;
				 }

				 inline uint64_t GenerateTagHash(const ItemType item, const size_t index) const {
					 size_t hash_sel;
					 if (index & 1) hash_sel = ((size_t)hash_sels[index / 2]) >> 4;
					 else hash_sel = ((size_t)hash_sels[index / 2]) & ((1 << 4) - 1);
					 uint64_t tag = (hash_fcns[hash_sel])(item);
					 tag &= ((1ULL << bits_per_item) - 1);
					 return tag + (tag == 0);
				 }

				 inline size_t AltIndexFromKey(const size_t index, const uint64_t key) const {
					 size_t alt_index = alt_hasher_(key) % table_->NumBuckets();
					 if (alt_index == index) {
						 alt_index = hasher_(key) % table_->NumBuckets();
						 if (alt_index == index) return (alt_index ? alt_index - 1 : table_->NumBuckets() - 1);
					 }
					 return alt_index;
				 }

				 inline void IncrHashSel(const size_t index) const {
					 size_t n0 = ((size_t)hash_sels[index / 2]) & ((1 << 4) - 1);
					 size_t n1 = ((size_t)hash_sels[index / 2]) >> 4;
					 if (index & 1) n1 = (n1 + 1) & ((1 << 4) - 1);
					 else n0 = (n0 + 1) & ((1 << 4) - 1);
					 hash_sels[index / 2] = (char)(n0 | (n1 << 4));
				 }

				 Status AddImpl(const uint64_t key);

				 // load factor is the fraction of occupancy
				 double LoadFactor() const { return 1.0 * Size() / table_->SizeInTags(); }

				 double BitsPerItem() const { return 8.0 * table_->SizeInBytes() / Size(); }

				 public: int map_inserts, map_kickouts, map_adapts;
				 public:
				 explicit CuckooFilter(const size_t max_num_keys) : num_items_(0), victim_(), hasher_() {
					 map_inserts = map_kickouts = map_adapts = 0;
					 size_t assoc = 4;
					 size_t num_buckets = upperpower2(std::max<uint64_t>(1, max_num_keys / assoc));
					 double frac = (double)max_num_keys / num_buckets / assoc;
					 if (frac > 0.96) {
						 num_buckets <<= 1;
					 }
					 victim_.used = false;
					 table_ = new TableType<bits_per_item>(num_buckets);
					 hash_sels = new char[num_buckets / 2];
					 hash_fcns = new HashFamily[1 << 4];
				 }

				 ~CuckooFilter() { delete table_; }

				 // Add an item to the filter.
				 Status Add(const ItemType item);

				 Status AddAndRecord(const uint64_t item, set_node *set, uint64_t set_len, FILE *fp);

				 Status AddFromRecording(const operation *ops, uint64_t *i);

				 // Report if the item is inserted, with false positive rate.
				 Status Contain(const ItemType item) const;

				 Status ContainReturn(const ItemType item, uint64_t *ret_location) const;


				 Status Adapt(const ItemType item, splinterdb *backing_map, const size_t bm_max_key_size, const size_t bm_max_val_size, splinterdb_lookup_result *bm_result, unsigned char* buffer); // TODO: rename this to AdaptWithSplinterDB and add back in the old Adapt so older experiments still work (do same with Delete)

				 Status AdaptAndRecord(const ItemType item, set_node *set, const uint64_t set_len, FILE *fp);

				 Status AdaptFromRecording(const operation *ops, uint64_t *i);

				 // Delete an key from the filter
				 Status Delete(const ItemType item, splinterdb *backing_map, const size_t bm_max_key_size, const size_t bm_max_val_size, splinterdb_lookup_result *bm_result, unsigned char* buffer);

				 /* methods for providing stats  */
				 // summary infomation
				 std::string Info() const;

				 // number of current inserted items;
				 size_t Size() const { return num_items_; }

				 // size of the filter in bytes.
				 size_t SizeInBytes() const { return table_->SizeInBytes(); }
			 };

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Add(const ItemType item) {
			if (victim_.used) {
				//victim_.used = false;
				return NotEnoughSpace;
			}

			return AddImpl(item);
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::AddImpl(const uint64_t key) {
			size_t curindex = GenerateIndexHash(key);
			uint64_t curkey = key;
			uint64_t curtag = GenerateTagHash(curkey, curindex);

			for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
				bool kickout = count > 0;
				int insert_index = table_->InsertTagToBucket(curindex, curtag, curkey, kickout);
				if (insert_index > 0) {
					//BACKING_MAP_INSERT(backing_map, curindex + (insert_index * table_->NumBuckets()), curkey);
					map_inserts++;
					num_items_++;
					return Ok;
				}
				insert_index = -insert_index; // TODO: FIX THIS - THIS IS A REMNANT FROM USING STXXL TO REPLACE THE MIRRORED_TABLE'S BACKING STORE
				//                                     CHANGE THE InsertTagToBucket FUNCTION BACK TO RETURNING THE OLD KEY?
				//BACKING_MAP_T::iterator item = backing_map.find(curindex + (insert_index * table_->NumBuckets()));

				map_kickouts++;
				//uint64_t nextkey = item->second;
				//item->second = curkey;
				//curkey = nextkey;
				curindex = AltIndexFromKey(curindex, curkey);
				curtag = GenerateTagHash(curkey, curindex);
			}

			victim_.index = curindex;
			victim_.tag = curtag;
			victim_.used = true;
			return Ok;
		}

	void bp_acf() {}

	void pad_data(void *dest, const void *src, const size_t dest_len, const size_t src_len, const int flagged) {
		assert(dest_len >= src_len);
		if (flagged) memset(dest, 0xff, dest_len);
		else bzero(dest, dest_len);
		memcpy(dest, src, src_len);
	}

	slice padded_slice(const void *data, const size_t dest_len, const size_t src_len, void *buffer, const int flagged) {
		pad_data(buffer, data, dest_len, src_len, flagged);

		return slice_create(dest_len, buffer);
	}

	int db_insert(splinterdb *database, const void *key_data, const size_t key_len, const size_t max_key_len, unsigned char *key_buffer, const void *val_data, const size_t val_len, const size_t max_val_len, unsigned char *val_buffer, const int flagged) {
		pad_data(key_buffer, key_data, max_key_len, key_len, flagged);
		pad_data(val_buffer, val_data, max_val_len, val_len, 0);
		slice key = slice_create(max_key_len, key_buffer);
		slice val = slice_create(max_val_len, val_buffer);

		return splinterdb_insert(database, key, val);
	}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::AddAndRecord(const uint64_t key, set_node *set, uint64_t set_len, FILE *fp) {
			if (victim_.used) {
				return NotEnoughSpace;
			}

			size_t curindex = GenerateIndexHash(key);
			uint64_t curkey = key;
			uint64_t curtag = GenerateTagHash(curkey, curindex);

			for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
				bool kickout = count < kMaxCuckooCount; // TODO: either remove this line or figure out if still needed
				int kickout_index = rand() % 4;
				int insert_index = table_->InsertTagToBucketWithSpecificKickout(curindex, curtag, curkey, kickout, kickout_index);
				fprintf(fp, "i %lu %d\n", curkey, kickout_index);
				if (insert_index == 0) break;
				if (insert_index > 0) {
					set_insert(set, set_len, curindex + ((insert_index - 1) * table_->NumBuckets()), curkey);

					num_items_++;
					return Ok;
				}
				insert_index = -insert_index;

				uint64_t location = curindex + ((insert_index - 1) * table_->NumBuckets()), nextkey;
				set_query(set, set_len, location, &nextkey);

				set_delete(set, set_len, location);
				set_insert(set, set_len, location, curkey);

				curkey = nextkey;
				curindex = AltIndexFromKey(curindex, curkey);
				curtag = GenerateTagHash(curkey, curindex);
			}

			victim_.index = curindex;
			victim_.tag = curtag;
			victim_.used = true;
			return Ok;
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::AddFromRecording(const operation *ops, uint64_t *i) {
			if (victim_.used) {
				return NotEnoughSpace;
			}

			uint64_t key = ops[*i].a;

			size_t curindex = GenerateIndexHash(key);
			uint64_t curkey = key;
			uint64_t curtag = GenerateTagHash(curkey, curindex);

			for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
				bool kickout = count < kMaxCuckooCount; // TODO: either remove this line or figure out if still needed
				int insert_index = table_->InsertTagToBucketWithSpecificKickout(curindex, curtag, curkey, kickout, ops[*i].b);
				if (insert_index == 0) break;
				if (insert_index > 0) {
					num_items_++;
					return Ok;
				}
				insert_index = -insert_index;

				*i = *i + 1;
				curkey = ops[*i].a;

				curindex = AltIndexFromKey(curindex, curkey);
				curtag = GenerateTagHash(curkey, curindex);
			}

			victim_.index = curindex;
			victim_.tag = curtag;
			victim_.used = true;
			*i = *i - 1;
			return Ok;
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Adapt(const ItemType key, splinterdb *backing_map, const size_t bm_max_key_size, const size_t bm_max_val_size, splinterdb_lookup_result *bm_result, unsigned char* buffer) {
			uint64_t i1, i2, index;
			uint64_t tag1, tag2;

			i1 = GenerateIndexHash(key);
			tag1 = GenerateTagHash(key, i1);
			i2 = GenerateAltIndexHash(key, i1);
			tag2 = GenerateTagHash(key, i2);
			if (victim_.used && ((tag1 == victim_.tag && i1 == victim_.index) || (tag2 == victim_.tag && i2 == victim_.index))) return Ok;

			if (table_->FindTagInBucket(i1, tag1)) index = i1;
			else {
				assert(table_->FindTagInBucket(i2, tag2));
				index = i2;
			}
			
			IncrHashSel(index);
			
			for (size_t j = 0; j < 4; j++) {
				uint64_t location_data = index + (j * table_->NumBuckets());
				slice location = padded_slice(&location_data, bm_max_key_size, sizeof(location_data), buffer, 0);
				splinterdb_lookup(backing_map, location, bm_result);

				if (!splinterdb_lookup_found(bm_result)) continue;

				map_adapts++;
				slice result_val;
				splinterdb_lookup_result_value(bm_result, &result_val);

				uint64_t updated_tag;
				memcpy(&updated_tag, slice_data(result_val), sizeof(uint64_t));
				table_->UpdateTag(index, j, GenerateTagHash(updated_tag, index));
			}

			/*uint64_t keys[4];
			table_->GetKeys(index, keys);
			for (size_t j = 0; j < 4; j++) {
				if (keys[j] == 0) break;
				table_->UpdateTag(index, j, GenerateTagHash(keys[j], index));
			}*/

			return Ok;
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::AdaptAndRecord(const ItemType key, set_node *set, uint64_t set_len, FILE *fp) {
			uint64_t i1, i2, index;
			uint64_t tag1, tag2;

			i1 = GenerateIndexHash(key);
			tag1 = GenerateTagHash(key, i1);
			i2 = GenerateAltIndexHash(key, i1);
			tag2 = GenerateTagHash(key, i2);
			if (victim_.used && ((tag1 == victim_.tag && i1 == victim_.index) || (tag2 == victim_.tag && i2 == victim_.index))) return Ok;

			if (table_->FindTagInBucket(i1, tag1)) index = i1;
			else {
				assert(table_->FindTagInBucket(i2, tag2));
				index = i2;
			}
			
			IncrHashSel(index);
			
			uint64_t i = 0;
			uint64_t keys[4];
			uint64_t indexes[4];
			for (size_t j = 0; j < 4; j++) {
				uint64_t orig_key;
				if (!set_query(set, set_len, index + (j * table_->NumBuckets()), &orig_key)) continue;
				keys[i] = orig_key;
				indexes[i] = j;
				i++;
				table_->UpdateTag(index, j, GenerateTagHash(orig_key, index));
			}
			fprintf(fp, "a %lu %lu\n", key, i);
			for (size_t j = 0; j < i; j++) {
				fprintf(fp, "c %lu %lu\n", keys[j], indexes[j]);
			}

			return Ok;
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::AdaptFromRecording(const operation *ops, uint64_t *i) {
			uint64_t i1, i2, index;
			uint64_t tag1, tag2;

			uint64_t key = ops[*i].a;

			i1 = GenerateIndexHash(key);
			tag1 = GenerateTagHash(key, i1);
			i2 = GenerateAltIndexHash(key, i1);
			tag2 = GenerateTagHash(key, i2);
			if (victim_.used && ((tag1 == victim_.tag && i1 == victim_.index) || (tag2 == victim_.tag && i2 == victim_.index))) return Ok;

			if (table_->FindTagInBucket(i1, tag1)) index = i1;
			else {
				assert(table_->FindTagInBucket(i2, tag2));
				index = i2;
			}
			
			IncrHashSel(index);
			
			uint64_t num_adapts = ops[*i].b;
			for (size_t j = 0; j < num_adapts; j++) {
				*i = *i + 1;
				uint64_t orig_key = ops[*i].a;
				table_->UpdateTag(index, ops[*i].b, GenerateTagHash(orig_key, index));
			}

			return Ok;
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Contain(const ItemType key) const {
			bool found = false;
			uint64_t i1, i2;
			uint64_t tag1, tag2;

			i1 = GenerateIndexHash(key);
			tag1 = GenerateTagHash(key, i1);
			i2 = GenerateAltIndexHash(key, i1);
			tag2 = GenerateTagHash(key, i2);

			assert(i1 == AltIndexFromKey(i2, key));

			found = victim_.used && ((tag1 == victim_.tag && i1 == victim_.index) ||
				(tag2 == victim_.tag && i2 == victim_.index));

			if (found || table_->FindTagInBucket(i1, tag1) || table_->FindTagInBucket(i2, tag2)) {
				return Ok;
			} else {
				return NotFound;
			}
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::ContainReturn(const ItemType key, uint64_t *ret_location) const {
			bool found = false;
			uint64_t i1, i2;
			uint64_t tag1, tag2;

			i1 = GenerateIndexHash(key);
			tag1 = GenerateTagHash(key, i1);
			i2 = GenerateAltIndexHash(key, i1);
			tag2 = GenerateTagHash(key, i2);

			assert(i1 == AltIndexFromKey(i2, key));

			found = victim_.used && ((tag1 == victim_.tag && i1 == victim_.index) ||
				(tag2 == victim_.tag && i2 == victim_.index));
			if (found) return NotFound; // The victim will not be able to adapt
			int insert_index;
			if ((insert_index = table_->FindTagLocationInBucket(i1, tag1)) >= 0) {
				*ret_location = i1 + (insert_index * table_->NumBuckets());
				return Ok;
			}
			if ((insert_index = table_->FindTagLocationInBucket(i2, tag2)) >= 0) {
				*ret_location = i2 + (insert_index * table_->NumBuckets());
				return Ok;
			}
			return NotFound;
		}


	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Delete(const ItemType key, splinterdb *backing_map, const size_t bm_max_key_size, const size_t bm_max_val_size, splinterdb_lookup_result *bm_result, unsigned char* buffer) {
			size_t i1, i2;
			uint64_t tag1, tag2;

			i1 = GenerateIndexHash(key);
			tag1 = GenerateTagHash(key, i1);
			i2 = GenerateAltIndexHash(key, i1);
			tag2 = GenerateTagHash(key, i2);

			int j;

			if ((j = table_->DeleteTagFromBucket(i1, tag1)) != 0) {
				num_items_--;
				uint64_t location_data = i1 + (j * table_->NumBuckets());
				slice location = padded_slice(&location_data, bm_max_key_size, sizeof(location_data), buffer, 0);
				splinterdb_delete(backing_map, location);
				goto TryEliminateVictim;
			} else if ((j = table_->DeleteTagFromBucket(i2, tag2)) != 0) {
				num_items_--;
				uint64_t location_data = i2 + (j * table_->NumBuckets());
				slice location = padded_slice(&location_data, bm_max_key_size, sizeof(location_data), buffer, 0);
				splinterdb_delete(backing_map, location);
				goto TryEliminateVictim;
			} else if (victim_.used && ((tag1 == victim_.tag && i1 == victim_.index) ||
					(tag2 == victim_.tag && i2 == victim_.index))) {
				// num_items_--;
				victim_.used = false;
				return Ok;
			} else {
				return NotFound;
			}
TryEliminateVictim:
			/*if (victim_.used) {
				victim_.used = false;
				AddImpl(key, backing_map);
			}*/
			return Ok;
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		std::string CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Info() const {
			std::stringstream ss;
			ss << "CuckooFilter Status:\n"
				<< "\t\t" << table_->Info() << "\n"
				<< "\t\tKeys stored: " << Size() << "\n"
				<< "\t\tLoad factor: " << LoadFactor() << "\n"
				<< "\t\tHashtable size: " << (table_->SizeInBytes() >> 10) << " KB\n";
			if (Size() > 0) {
				ss << "\t\tbit/key:   " << BitsPerItem() << "\n";
			} else {
				ss << "\t\tbit/key:   N/A\n";
			}
			return ss.str();
		}
}  // namespace cuckoofilter
#endif  // CUCKOO_FILTER_CUCKOO_FILTER_H_
