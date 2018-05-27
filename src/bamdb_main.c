#ifdef BUILD_BAMDB_WRITER
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "bam_api.h"
#include "bamdb.h"
#include "bamdb_lmdb.h"

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

int main(int argc, char *argv[]) {
  int rc = 0;
  int c;
  bam_args_t bam_args;
  int max_rows = 0;

  bam_args.index_file_name = NULL;
  bam_args.bx = NULL;
  bam_args.output_file_name = NULL;
  bam_args.convert_to = BAMDB_CONVERT_TO_TEXT;
  while ((c = getopt(argc, argv, "t:f:n:i:b:o:")) != -1) {
    switch (c) {
      case 't':
        if (strcmp(optarg, "lmdb") == 0) {
          bam_args.convert_to = BAMDB_CONVERT_TO_LMDB;
        } else if (strcmp(optarg, "text") == 0) {
          bam_args.convert_to = BAMDB_CONVERT_TO_TEXT;
        } else {
          fprintf(stderr, "Invalid output format %s\n", optarg);
          return 1;
        }
        break;
      case 'f':
        strcpy(bam_args.input_file_name, optarg);
        break;
      case 'n':
        max_rows = atoi(optarg);
        break;
      case 'i':
        bam_args.index_file_name = strdup(optarg);
        break;
      case 'b':
        bam_args.bx = strdup(optarg);
        break;
      case 'o':
        bam_args.output_file_name = strdup(optarg);
        break;
      default:
        fprintf(stderr, "Unknown argument\n");
        return 1;
    }
  }

  if (bam_args.convert_to == BAMDB_CONVERT_TO_LMDB) {
    bamdb_indices_t target_indices = {.includes_qname = true,
                                      .num_key_indices = 1,
                                      .key_indices = malloc(sizeof(char *))};

    target_indices.key_indices[0] = calloc(1, 3);
    /* Get key name from first non optional argument */
    if (optind < argc) {
      strncpy(target_indices.key_indices[0], argv[optind], 2);
    }

    return generate_index_file(bam_args.input_file_name,
                               bam_args.output_file_name, &target_indices);
  }

  /* Get filename from first non optional argument */
  if (optind < argc) {
    strcpy(bam_args.input_file_name, argv[optind]);
  }

  if (bam_args.bx != NULL && bam_args.index_file_name != NULL) {
    if (bam_args.output_file_name != NULL) {
      /* Write resulting rows to file */
      offset_list_t *offset_list = calloc(1, sizeof(offset_list_t));

      rc = get_offsets_lmdb(offset_list, bam_args.index_file_name, "BX",
                            bam_args.bx);
      rc = write_row_subset(bam_args.input_file_name, offset_list,
                            bam_args.output_file_name);
      free(offset_list);
    } else {
      /* Print rows in tab delim format */
      bam_row_set_t *row_set = NULL;
      rc = get_bam_rows(&row_set, bam_args.input_file_name,
                        bam_args.index_file_name, "BX", bam_args.bx);

      if (row_set != NULL) {
        for (size_t j = 0; j < row_set->num_entries; ++j) {
          print_sequence_row(row_set->rows[j]);
        }
        free_bamdb_row_set(row_set);
      }
    }
  }
}
#endif
