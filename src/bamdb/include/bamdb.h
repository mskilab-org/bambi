#ifndef BAMDB_H
#define BAMDB_H

#include <stdbool.h>
#include <stddef.h>

#include "bamdb_status.h"
#include "bam_api.h"

#define MAX_FILENAME 1024

typedef struct bamdb_indices {
  bool includes_qname;
  size_t num_key_indices;  // Does not include qname index
  char **key_indices;
} bamdb_indices_t;


#ifdef BUILD_BAMDB_WRITER
/** @brief Create an index for a given bam file
 *
 * Currently all bamdb indices are LMDB based, though this function may change
 * in the future to support multiple index types. This function will spawn
 * a deserialize thread and an additional thread per desired index column.
 *
 * @param[in] input_file The path of the bam file to index
 * @param[in] db_path Optional path of the generated index; a default path
 * based on the input filename will be used if this is NULL
 * @param[in] target_indices Struct containing the desired fields to be indexed
 * @return 0 on success or a non-zero error value on failure
 */
int generate_index_file(char *input_file_name, char *output_file_name,
                        bamdb_indices_t *target_indices);
#endif

void print_bamdb_rows(const char *input_file_name, const char *db_path,
                      const char *index_name, const char *key);

int write_row_subset(char *input_file_name, offset_list_t *offset_list,
                     char *out_filename);

/** @brief Find all rows matching a key in an indexed bam file
 *
 * This function will allocate space for the resulting rows; it is up to the
 * caller to free the results. In the event of an error we will return an
 * empty row set object. Callers should NOT pass a preallocated or
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
