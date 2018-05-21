#ifndef BAMLMDB_H
#define BAMLMDB_H

#include <lmdb.h>
#include <stdbool.h>

#include "bam_api.h"
#include "htslib/sam.h"

typedef struct indices {
  bool includes_qname;
  size_t num_key_indices;  // Does not include qname index
  char **key_indices;
} indices_t;

int convert_to_lmdb(samFile *input_file, char *db_path,
                    indices_t *target_indices);

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
indices_t *get_available_indices(const char *db_path);

bool is_index_present(const char *db_path, const char *index_name);

/** @brief Find all rows matching a key in an indexed bam file
 *
 * This function will allocate space for the resulting rows; it is up to the
 * caller to free the results. Callers should NOT pass a preallocated or
 * existing row set to this function.
 *
 * @param[out] output Location to store the resulting records
 * @param[in] input_file_name Path of the bam file to query
 * @param[in] db_path Top-level directory of the index database
 * @param[in] index_name Name of the field to search in
 * @param[in] key Specific index value to search for
 * @return 0 on success or a non-zero error value on failure
 */
int get_bam_rows(bam_row_set_t **output, const char *input_file_name,
                 const char *db_path, const char *index_name, const char *key);

#endif
