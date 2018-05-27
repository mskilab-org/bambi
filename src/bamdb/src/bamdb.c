#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* HTSlib */
#include "sam.h"

#include "bam_api.h"
#include "bamdb.h"
#include "bamdb_index_writer.h"
#include "bamdb_lmdb.h"

// Return number of characters an unsigned int takes when represented in base 10
#define get_int_chars(i) ((i == 0) ? 1 : floor(log10(i)) + 1)

int write_row_subset(char *input_file_name, offset_list_t *offset_list,
                     char *out_filename) {
  int rc = 0;
  int64_t src = 0;
  bam1_t *bam_row;
  offset_node_t *offset_node;
  bam_hdr_t *new_header = NULL;
  bam_hdr_t *header = NULL;
  BGZF *fp = bgzf_open(out_filename, "w");
  samFile *input_file = 0;

  if ((input_file = sam_open(input_file_name, "r")) == 0) {
    fprintf(stderr, "Unable to open file %s\n", input_file_name);
    return 1;
  }

  header = sam_hdr_read(input_file);
  if (header == NULL) {
    fprintf(stderr, "Unable to read the header from %s\n", input_file->fn);
    rc = 1;
    return 1;
  }

  new_header = bam_hdr_dup(header);
  rc = bam_hdr_write(fp, new_header);
  if (rc != 0) {
    fprintf(stderr, "Unable to write header for %s\n", out_filename);
    return rc;
  }

  bam_row = bam_init1();
  if (offset_list != NULL) {
    offset_node = offset_list->head;
    while (offset_node != NULL) {
      src = bgzf_seek(input_file->fp.bgzf, offset_node->offset, SEEK_SET);
      if (src != 0) {
        fprintf(stderr, "Error seeking to file offset\n");
        return 1;
      }

      rc = sam_read1(input_file, header, bam_row);
      rc = bam_write1(fp, bam_row);
      if (rc < 0) {
        fprintf(stderr, "Error writing row to file %s\n", out_filename);
        return rc;
      }

      offset_node = offset_node->next;
    }
  }

  bgzf_close(fp);
  return rc;
}

#ifdef BUILD_BAMDB_WRITER
int generate_index_file(char *input_file_name, char *output_file_name,
                        bamdb_indices_t *target_indices) {
  samFile *input_file = 0;

  if ((input_file = sam_open(input_file_name, "r")) == 0) {
    fprintf(stderr, "Unable to open file %s\n", input_file_name);
    return 1;
  }

  return generate_lmdb_index(input_file, output_file_name, target_indices);
}
#endif

void print_bamdb_rows(const char *input_file_name, const char *db_path,
                      const char *index_name, const char *key) {
  int rc = 0;
  bam_row_set_t *row_set = NULL;
  rc = get_bam_rows(&row_set, input_file_name, db_path, index_name, key);

  if (row_set != NULL) {
    for (size_t j = 0; j < row_set->num_entries; ++j) {
      print_sequence_row(row_set->rows[j]);
    }
  }
}
