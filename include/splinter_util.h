#ifndef _SPLINTER_UTIL_H_
#define _SPLINTER_UTIL_H_

#include <string.h>
#include <fcntl.h>

#include "include/gqf.h"
#include "include/gqf_int.h"
#include "include/hashutil.h"

#include <splinterdb/data.h> // for data config
#include <splinterdb/default_data_config.h> // for default data config
#include <splinterdb/public_platform.h>
#include <splinterdb/public_util.h>
#include <splinterdb/splinterdb.h> // for insert/update/delete

#define Kilo (1024UL)
#define Mega (1024UL * Kilo)
#define Giga (1024UL * Mega)

#define MAX_KEY_SIZE 16
#define MAX_VAL_SIZE 16

int merge_tuples(const data_config *cfg, slice key, message old_message, merge_accumulator *new_message);
int merge_tuples_final(const data_config *cfg, slice key, merge_accumulator *oldest_message);

data_config qf_data_config_init();
splinterdb_config qf_splinterdb_config_init(char *db_path, data_config *data_cfg);

void pad_data(void *dest, const void *src, const size_t dest_len, const size_t src_len, const int flagged);
slice padded_slice(const void *data, const size_t dest_len, const size_t src_len, void *buffer, const int flagged);
int db_insert(splinterdb *database, const void *key_data, const size_t key_len, const void *val_data, const size_t val_len, const int update, const int flagged);

int qf_splinter_insert(QF *qf, splinterdb *db, uint64_t key, int count);
int qf_splinter_insert_split(QF *qf, splinterdb *db, splinterdb *bm, uint64_t key, uint64_t val);

#endif
