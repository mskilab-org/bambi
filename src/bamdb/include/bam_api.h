/* Some additional functions to simplify pulling rows from BAM data */
#ifndef BAMAPI_H
#define BAMAPI_H

/* HTSlib */
#include "bgzf.h"
#include "sam.h"

const char *bam_get_rname(const bam1_t *row, const bam_hdr_t *header);
const char *bam_get_rnext(const bam1_t *row, const bam_hdr_t *header);

/**
 * Return pointer to beginning of formatted string in work_buffer,
 * advance work_buffer to the element after the formatted string
 */
char *bam_cigar_str(const bam1_t *row, char *work_buffer);
char *bam_seq_str(const bam1_t *row, char *work_buffer);
char *bam_qual_str(const bam1_t *row, char *work_buffer);

/**
 * Return the value associated with an optional string key on a BAM row.
 * This will advance the work buffer to the element after the value.
 * Returns a default value of '*' if they key is not found.
 */
char *bam_str_key(const bam1_t *row, const char *key, char *work_buffer);

typedef struct aux_elm_key {
  char key[2];
  char type;
  char subtype;
} aux_elm_key_t;

typedef struct aux_elm {
  aux_elm_key_t key;
  uint32_t val_size;
  void *val;
  struct aux_elm *next;
} aux_elm_t;

typedef struct aux_list {
  size_t n_tags;
  aux_elm_t *head;
  aux_elm_t *tail;
} aux_list_t;

typedef struct bam_sequence_row {
  char *qname;
  int flag;
  char *rname;
  int pos;
  int mapq;
  char *cigar;
  char *rnext;
  int pnext;
  int tlen;
  char *seq;
  char *qual;
  aux_list_t aux_list;
} bam_sequence_row_t;

typedef struct bam_aux_header {
  aux_elm_key_t key;
  struct bam_aux_header *next;
} bam_aux_header_t;

typedef struct bam_aux_header_list {
  bam_aux_header_t *head;
  bam_aux_header_t *tail;
} bam_aux_header_list_t;

typedef struct bam_row_set {
  size_t num_entries;
  bam_sequence_row_t **rows;
  bam_aux_header_list_t aux_tags;
} bam_row_set_t;

typedef struct offset_node {
  int64_t offset;
  struct offset_node *next;
} offset_node_t;

typedef struct offset_list {
  size_t num_entries;
  struct offset_node *head;
  struct offset_node *tail;
} offset_list_t;

int deserialize_bam_row(bam_sequence_row_t **out,
                        bam_aux_header_list_t *tag_list, const bam1_t *row,
                        const bam_hdr_t *header);
int get_bam_row(bam_sequence_row_t **out, bam_aux_header_list_t *tag_list,
                const int64_t offset, samFile *input_file, bam_hdr_t *header);
void print_sequence_row(bam_sequence_row_t *row);
void destroy_bam_sequence_row(bam_sequence_row_t *row);

/** @brief Destroy a bamdb row set including all sub objects
 *
 * @param[in] row_set The row set to be destroyed
 */
void free_bamdb_row_set(bam_row_set_t *row_set);

/* Return 0 on success */
int write_row_set_to_file(bam_row_set_t *row_set, bam_hdr_t *header,
                          char *out_filename);

#define bam_row_size(b) (sizeof(bam1_t) + b->l_data * sizeof(uint8_t))
#define bam_row_max_size(b) (sizeof(bam1_t) + b->m_data * sizeof(uint8_t))

/* I really hope we don't have sequences longer than this */
#define WORK_BUFFER_SIZE 65536

#endif
