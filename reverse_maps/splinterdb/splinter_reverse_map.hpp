#ifndef SPLINTER_REVERSE_MAP_HPP
#define SPLINTER_REVERSE_MAP_HPP

#include <cstdint>
#include <cstring>
#include <cstdio>
#include <string>
#include <thread>
#include <chrono>
#include <fcntl.h>

extern "C" {
#include <splinterdb/data.h>
#include <splinterdb/default_data_config.h>
#include <splinterdb/public_platform.h>
#include <splinterdb/public_util.h>
#include <splinterdb/splinterdb.h>
}

#include "reverse_map_config.hpp"

#define SPLINTER_MAX_KEY_SIZE 16
#define SPLINTER_MAX_VAL_SIZE 16
#define SPLINTER_Mega (1024UL * 1024UL)
#define SPLINTER_Giga (1024UL * SPLINTER_Mega)

class SplinterDB {
public:
  SplinterDB() : db(nullptr), quotient_remainder_bits(0), query_delay_us(0) {}

  int init(ReverseMapConfig config, int base_fingerprint_bits) {
    quotient_remainder_bits = base_fingerprint_bits;

    data_cfg = qf_data_config_init();
    splinterdb_cfg = qf_splinterdb_config_init(
        const_cast<char *>(config.storage_path.c_str()), &data_cfg);
    splinterdb_cfg.cache_size = config.cache_size_bytes;
    splinterdb_cfg.disk_size  = 20 * SPLINTER_Giga;

    if (config.overwrite) {
      remove(splinterdb_cfg.filename);
      if (splinterdb_create(&splinterdb_cfg, &db)) {
        fprintf(stderr, "SplinterDB: failed to create database at %s\n",
                config.storage_path.c_str());
        return -1;
      }
    } else {
      if (splinterdb_open(&splinterdb_cfg, &db)) {
        fprintf(stderr, "SplinterDB: failed to open database at %s\n",
                config.storage_path.c_str());
        return -1;
      }
    }

    splinterdb_lookup_result_init(db, &db_result, SPLINTERDB_LOOKUP_VALUE, 0, NULL);
    return 0;
  }

  int insert(uint64_t fingerprint, uint64_t rank, uint64_t key) {
    uint64_t encoded = encodeKey(fingerprint, rank);
    int isUpdate = (rank > 0);
    return db_insert(db, &encoded, sizeof(encoded), &key, sizeof(key),
                     isUpdate, 0) ? 0 : -1;
  }

  void bulkInsert(uint64_t *keys, uint64_t *fingerprints, uint64_t *ranks,
                  uint64_t numKeys) {
    for (uint64_t i = 0; i < numKeys; i++) {
      insert(fingerprints[i], ranks[i], keys[i]);
    }
  }

  int getFingerprint(uint64_t fingerprint, int rank, uint64_t *value) {
    if (query_delay_us > 0) {
      std::this_thread::sleep_for(
          std::chrono::microseconds(query_delay_us));
    }

    uint64_t encoded = encodeKey(fingerprint, rank);
    char buffer[2 * SPLINTER_MAX_KEY_SIZE];
    slice query = padded_slice(&encoded, SPLINTER_MAX_KEY_SIZE,
                               sizeof(encoded), buffer, 0);
    splinterdb_lookup(db, query, &db_result);
    if (splinterdb_lookup_found(&db_result)) {
      slice result_val;
      splinterdb_lookup_result_value(&db_result, &result_val);
      memcpy(value, slice_data(result_val), sizeof(uint64_t));
      return 0;
    }
    return -1;
  }

  int close() {
    if (db) {
      splinterdb_close(&db);
      db = nullptr;
    }
    return 0;
  }

private:
  uint64_t encodeKey(uint64_t fingerprint, int rank) {
    uint64_t mask = (1ULL << quotient_remainder_bits) - 1;
    uint64_t key = (fingerprint & mask) << (64 - quotient_remainder_bits);
    key += rank;
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    key = __builtin_bswap64(key);
#endif
    return key;
  }

  static int merge_tuples(const data_config *cfg, slice key,
                           message old_message, merge_accumulator *new_message) {
    int new_message_len = merge_accumulator_length(new_message);
    if (!merge_accumulator_resize(
            new_message, message_length(old_message) + new_message_len))
      return -1;
    memcpy(merge_accumulator_data(new_message) + new_message_len,
           message_data(old_message), message_length(old_message));
    return 0;
  }

  static int merge_tuples_final(const data_config *cfg, slice key,
                                 merge_accumulator *oldest_message) {
    merge_accumulator_set_class(oldest_message, MESSAGE_TYPE_INSERT);
    return 0;
  }

  data_config qf_data_config_init() {
    data_config data_cfg;
    default_data_config_init(&data_cfg);
    data_cfg.merge_tuples = merge_tuples;
    data_cfg.merge_tuples_final = merge_tuples_final;
    return data_cfg;
  }

  splinterdb_config qf_splinterdb_config_init(char *db_path,
                                               data_config *data_cfg) {
    splinterdb_config cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.filename   = db_path;
    cfg.cache_size = 64 * SPLINTER_Mega;
    cfg.disk_size  = 20 * SPLINTER_Giga;
    cfg.data_cfg   = data_cfg;
    cfg.io_flags   = O_RDWR | O_CREAT | O_DIRECT;
    return cfg;
  }

  static void pad_data(void *dest, const void *src, size_t dest_len,
                        size_t src_len, int flagged) {
    if (dest_len < src_len) return;
    if (flagged)
      memset(dest, 0xff, dest_len);
    else
      memset(dest, 0, dest_len);
    memcpy(dest, src, src_len);
  }

  static slice padded_slice(const void *data, size_t dest_len,
                             size_t src_len, void *buffer, int flagged) {
    pad_data(buffer, data, dest_len, src_len, flagged);
    return slice_create(dest_len, buffer);
  }

  static int db_insert(splinterdb *database, const void *key_data,
                        size_t key_len, const void *val_data, size_t val_len,
                        int update, int flagged) {
    char key_padded[SPLINTER_MAX_KEY_SIZE];
    char val_padded[SPLINTER_MAX_VAL_SIZE];
    pad_data(key_padded, key_data, SPLINTER_MAX_KEY_SIZE, key_len, flagged);
    pad_data(val_padded, val_data, SPLINTER_MAX_VAL_SIZE, val_len, 0);
    slice key_slice = slice_create(SPLINTER_MAX_KEY_SIZE, key_padded);
    slice val_slice = slice_create(SPLINTER_MAX_VAL_SIZE, val_padded);

    if (update)
      return splinterdb_update(database, key_slice, val_slice, NULL);
    return splinterdb_insert(database, key_slice, val_slice, NULL);
  }

  data_config data_cfg;
  splinterdb_config splinterdb_cfg;
  splinterdb *db;
  splinterdb_lookup_result db_result;
  int quotient_remainder_bits;
  uint64_t query_delay_us;
};

#endif