/**
 * @file bamdb_lmdb.h
 * @brief Functions for creating a bamdb index
 *
 */
#ifndef BAMDB_LMDB_H
#define BAMDB_LMDB_H

#include <stdbool.h>

#include <lmdb.h>

#include "bam_api.h"
#include "bamdb.h"

/* HTSLib */
#include "bgzf.h"
#include "hts.h"
#include "sam.h"

char *get_default_dbname(const char *filename);

int get_lmdb_env(MDB_env **env, const char *full_db_path, bool read_only);

int commit_lmdb_transaction(MDB_txn *txn);

/** @brief Return matching bam offsets from an LMDB based index
 *
 * @param[out] output Offset list to populate with the results
 * @param[in] db_path Top-level directory of the index database
 * @param[in] index_name Name of the field to search in
 * @param[in] key Specific index value to search for
 * @return 0 on success or a non-zero error value on failure
 */
int get_offsets_lmdb(offset_list_t *offset_list, const char *db_path,
                     const char *index_name, const char *key);

/**
 * Get a list of the available indices in an existing lmdb database
 */
bamdb_indices_t *get_available_indices(const char *db_path);

bool is_index_present(const char *db_path, const char *index_name);

#endif
