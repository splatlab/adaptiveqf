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

				 public:
				 explicit CuckooFilter(const size_t max_num_keys) : num_items_(0), victim_(), hasher_() {
					 size_t assoc = 4;
					 size_t num_buckets = upperpower2(std::max<uint64_t>(1, max_num_keys / assoc));
					 double frac = (double)max_num_keys / num_buckets / assoc;
					 if (frac > 0.96) {
						 num_buckets <<= 1;
					 }
					 victim_.used = false;
					 table_ = new TableType<bits_per_item>(num_buckets);
					 hash_sels = new char[num_buckets / 2];
					 hash_fcns = new HashFamily[num_buckets];
				 }

				 ~CuckooFilter() { delete table_; }

				 // Add an item to the filter.
				 Status Add(const ItemType item);

				 // Report if the item is inserted, with false positive rate.
				 Status Contain(const ItemType item) const;


				 Status Adapt(const ItemType item) const;

				 // Delete an key from the filter
				 Status Delete(const ItemType item);

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
				if (table_->InsertTagToBucket(curindex, curtag, curkey, kickout)) {
					num_items_++;
					return Ok;
				}
				curindex = AltIndexFromKey(curindex, key);
				curtag = GenerateTagHash(curkey, curindex);
			}

			victim_.index = curindex;
			victim_.tag = curtag;
			victim_.used = true;
			return Ok;
		}

	template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Adapt(const ItemType key) const {
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
			uint64_t keys[4];
			table_->GetKeys(index, keys);
			for (size_t j = 0; j < 4; j++) {
				if (keys[j] == 0) break;
				table_->UpdateTag(index, j, GenerateTagHash(keys[j], index));
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
		Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Delete(const ItemType key) {
			/*size_t i1, i2;
			uint64_t tag;

			GenerateIndexTagHash(key, &i1, &tag);
			i2 = AltIndex(i1, tag);

			if (table_->DeleteTagFromBucket(i1, tag)) {
				num_items_--;
				goto TryEliminateVictim;
			} else if (table_->DeleteTagFromBucket(i2, tag)) {
				num_items_--;
				goto TryEliminateVictim;
			} else if (victim_.used && tag == victim_.tag &&
					(i1 == victim_.index || i2 == victim_.index)) {
				// num_items_--;
				victim_.used = false;
				return Ok;
			} else {
				return NotFound;
			}
TryEliminateVictim:
			if (victim_.used) {
				victim_.used = false;
				size_t i = victim_.index;
				uint32_t tag = victim_.tag;
				AddImpl(i, tag);
			}
			return Ok;*/
			return NotFound;
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
