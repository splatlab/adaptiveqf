#include "include/splinter_util.h"
#include <assert.h>

int merge_tuples(const data_config *cfg, slice key, message old_message, merge_accumulator *new_message) {
	int new_message_len = merge_accumulator_length(new_message);
	if (!merge_accumulator_resize(new_message, message_length(old_message) + new_message_len)) return -1;
	memcpy(merge_accumulator_data(new_message) + new_message_len, message_data(old_message), message_length(old_message));
	/*if (merge_accumulator_length(new_message) > 2000) {
		printf("%lu\n", merge_accumulator_length(new_message));
	}*/
	return 0;
}

int merge_tuples_final(const data_config *cfg, slice key, merge_accumulator *oldest_message) {
	merge_accumulator_set_class(oldest_message, MESSAGE_TYPE_INSERT);
	return 0;
}

data_config qf_data_config_init() {
	data_config data_cfg;
	default_data_config_init(MAX_KEY_SIZE, &data_cfg);

	data_cfg.merge_tuples = merge_tuples;
	data_cfg.merge_tuples_final = merge_tuples_final;

	return data_cfg;
}

splinterdb_config qf_splinterdb_config_init(char *db_path, data_config *data_cfg) {
	splinterdb_config splinterdb_cfg = (splinterdb_config){
		.filename   = db_path,
		.cache_size = 64 * Mega,
		.disk_size  = 20 * Giga,
		.data_cfg   = data_cfg,
		.io_flags   = O_RDWR | O_CREAT | O_DIRECT
	};
	
	return splinterdb_cfg;
}

void pad_data(void *dest, const void *src, const size_t dest_len, const size_t src_len, const int flagged) {
	if (dest_len < src_len) return;
	if (flagged) memset(dest, 0xff, dest_len);
	else bzero(dest, dest_len);
	memcpy(dest, src, src_len);
}

slice padded_slice(const void *data, const size_t dest_len, const size_t src_len, void *buffer, const int flagged) {
	pad_data(buffer, data, dest_len, src_len, flagged);
	
	return slice_create(dest_len, buffer);
}

int db_insert(splinterdb *database, const void *key_data, const size_t key_len, const void *val_data, const size_t val_len, const int update, const int flagged) {
	char key_padded[MAX_KEY_SIZE];
	char val_padded[MAX_VAL_SIZE];
	pad_data(key_padded, key_data, MAX_KEY_SIZE, key_len, flagged);
	pad_data(val_padded, val_data, MAX_VAL_SIZE, val_len, 0);
	slice key_slice = slice_create(MAX_KEY_SIZE, key_padded);
	slice val_slice = slice_create(MAX_VAL_SIZE, val_padded);

	if (update) return splinterdb_update(database, key_slice, val_slice);
	else return splinterdb_insert(database, key_slice, val_slice);
}

int qf_splinter_insert(QF *qf, splinterdb *db, uint64_t key, int count) {
	qf_insert_result result;
	int ret = qf_insert_using_ll_table(qf, key, count, &result, QF_NO_LOCK | QF_KEY_IS_HASH);
	if (ret < 0) {
		return 0;
	}
	else {
		//uint64_t id_hashed = MurmurHash64A(&result.minirun_id, sizeof(result.minirun_id), 26571997);
		//if (!db_insert(db, &id_hashed, sizeof(id_hashed), &key, sizeof(key), result.minirun_existed, 0)) return result.minirun_existed + 1;

		result.minirun_id <<= 64 - qf->metadata->quotient_remainder_bits;
		if (!db_insert(db, &result.minirun_id, sizeof(result.minirun_id), &key, sizeof(key), result.minirun_existed, 0)) return result.minirun_existed + 1;

		//if (!db_insert(db, &result.minirun_id, sizeof(result.minirun_id), &key, sizeof(key), result.minirun_existed, 0)) return result.minirun_existed + 1;
		return 0;
	}
}

int qf_splinter_insert_split(QF *qf, splinterdb *db, splinterdb *bm, uint64_t key, uint64_t val) {
	qf_insert_result result;
	int ret = qf_insert_using_ll_table(qf, key, 1, &result, QF_NO_LOCK | QF_KEY_IS_HASH);
	if (ret < 0) {
		return 0;
	}
	else {
		if (db_insert(db, &key, sizeof(key), &val, sizeof(val), 0, 0)) return 0;
		assert(result.minirun_id == key % (1 << (qf->metadata->quotient_bits + qf->metadata->bits_per_slot)));
		assert(sizeof(result.minirun_id) <= MAX_KEY_SIZE);
		if (db_insert(bm, &result.minirun_id, sizeof(result.minirun_id), &key, sizeof(key), result.minirun_existed, 0)) return 0;
		return result.minirun_existed + 1;
	}
}
