#ifndef BAMLMDB_H
#define BAMLMDB_H

#include <lmdb.h>
#include <stdbool.h>

#include "bam_api.h"
#include "htslib/sam.h"


typedef struct indices {
	bool includes_qname;
	size_t num_key_indices; // Does not include qname index
	char **key_indices;
} indices_t;

int convert_to_lmdb(samFile *input_file, char *db_path, indices_t *target_indices);

/* Return number of offsets on success, -1 on failure */
int get_offsets(offset_list_t *offset_list, const char *db_path, const char *index_name, const char *key);

/**
 * Get a list of the available indices in an existing lmdb database
 */
indices_t *get_available_indices(const char *db_path);

bool is_index_present(const char *db_path, const char *index_name);

bam_row_set_t *get_bam_rows(const char *input_file_name, const char *db_path, const char *key_type, const char *key_value);

#endif
