#ifndef BAMDB_H
#define BAMDB_H

#include <stdbool.h>
#include <stddef.h>

#define MAX_FILENAME 1024

typedef struct indices {
  bool includes_qname;
  size_t num_key_indices;  // Does not include qname index
  char **key_indices;
} indices_t;

enum bamdb_convert_to {
  BAMDB_CONVERT_TO_TEXT,
  BAMDB_CONVERT_TO_SQLITE,
  BAMDB_CONVERT_TO_LMDB
};

typedef struct bamdb_args {
  enum bamdb_convert_to convert_to;
  char input_file_name[MAX_FILENAME];
  char *index_file_name;
  char *output_file_name;
  char *bx;
} bam_args_t;

int generate_index_file(char *input_file_name, char *output_file_name,
                        indices_t *target_indices);

void print_bx_rows(char **input_file_name, char **db_path, char **bx);

int write_row_subset(char *input_file_name, offset_list_t *offset_list,
                     char *out_filename);

#endif
