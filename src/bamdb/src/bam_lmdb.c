#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lmdb.h>

#include "bam_api.h"
#include "bam_lmdb.h"
#include "bamdb_status.h"

#define LMDB_POSTFIX "_lmdb"
#define WORK_BUFFER_SIZE 65536
/* Making this huge because so we don't have to resize while running*/
#define LMDB_INIT_MAPSIZE 100000000000
#define MAX_PATH_CHARS 2048

char *get_default_dbname(const char *filename) {
  char *db_name;
  char *dot = strrchr(filename, '.');

  /* This needs to be freed by the caller */
  db_name = malloc(dot - filename + sizeof(LMDB_POSTFIX) + 1);
  strncpy(db_name, filename, dot - filename);
  strcpy(db_name + (dot - filename), LMDB_POSTFIX);

  return db_name;
}

int commit_lmdb_transaction(MDB_txn *txn) {
  int rc;

  rc = mdb_txn_commit(txn);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error commiting transaction: %s\n", mdb_strerror(rc));
    return BAMDB_DB_ERROR;
  }

  return BAMDB_SUCCESS;
}

int get_lmdb_env(MDB_env **env, const char *full_db_path, bool read_only) {
  int rc;
  int flags;
  int ret = BAMDB_SUCCESS;

  rc = mdb_env_create(env);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error creating env: %s\n", mdb_strerror(rc));
    return BAMDB_DB_ERROR;
  }

  rc = mdb_env_set_maxdbs(*env, 1);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error setting maxdbs: %s\n", mdb_strerror(rc));
    return BAMDB_DB_ERROR;
  }

  /* We don't ever read while writing */
  flags = MDB_NOLOCK;
  if (read_only) {
    flags = flags | MDB_RDONLY;
  }

  rc = mdb_env_open(*env, full_db_path, flags, 0664);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error opening env: %s\n", mdb_strerror(rc));
    return BAMDB_DB_ERROR;
  }

  rc = mdb_env_set_mapsize(*env, LMDB_INIT_MAPSIZE);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error setting map size: %s\n", mdb_strerror(rc));
    return BAMDB_DB_ERROR;
  }

  return ret;
}

static int open_ro_handle(MDB_env *output, MDB_dbi *dbi, MDB_cursor **cur,
                          MDB_txn *txn, const char *key_name,
                          const char *db_path) {
  char target_path[MAX_PATH_CHARS];
  int rc = 0;

  snprintf(target_path, MAX_PATH_CHARS, "%s/%s", db_path, key_name);
  rc = get_lmdb_env(&output, target_path, true);
  if (rc != BAMDB_SUCCESS) {
    return BAMDB_DB_ERROR;
  }

  rc = mdb_txn_begin(output, NULL, MDB_RDONLY, &txn);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error beginning LMDB transaction: %s\n", mdb_strerror(rc));
    return BAMDB_DB_ERROR;
  }

  rc = mdb_dbi_open(txn, NULL, MDB_DUPSORT | MDB_DUPFIXED, dbi);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error opening LMDB database handle: %s\n",
            mdb_strerror(rc));
    return BAMDB_DB_ERROR;
  }

  rc = mdb_cursor_open(txn, *dbi, cur);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error creating LMDB cursor: %s\n", mdb_strerror(rc));
    return BAMDB_DB_ERROR;
  }

  return BAMDB_SUCCESS;
}

int get_offsets_lmdb(offset_list_t *offset_list, const char *db_path,
                     const char *index_name, const char *key) {
  MDB_env *env = NULL;
  MDB_dbi dbi;
  MDB_cursor *cur = NULL;
  MDB_txn *txn = NULL;
  MDB_val db_key, data;
  int rc;

  db_key.mv_size = strlen(key);
  db_key.mv_data = strdup(key);

  rc = open_ro_handle(env, &dbi, &cur, txn, index_name, db_path);
  if (rc != MDB_SUCCESS) {
    return BAMDB_DB_ERROR;
  }

  rc = mdb_cursor_get(cur, &db_key, &data, MDB_SET);

  if (rc == MDB_NOTFOUND) {
    /* No matching rows for the given index query */
    return BAMDB_SUCCESS;
  } else if (rc != MDB_SUCCESS) {
    return BAMDB_DB_ERROR;
  }
  if ((rc = mdb_cursor_get(cur, &db_key, &data, MDB_FIRST_DUP)) == 0) {
    offset_node_t *new_node = calloc(1, sizeof(offset_node_t));
    new_node->offset = *(int64_t *)data.mv_data;

    offset_list->head = new_node;
    offset_list->tail = new_node;
    offset_list->num_entries++;

    while ((rc = mdb_cursor_get(cur, &db_key, &data, MDB_NEXT_DUP)) == 0) {
      offset_node_t *new_node = calloc(1, sizeof(offset_node_t));
      new_node->offset = *(int64_t *)data.mv_data;

      offset_list->tail->next = new_node;
      offset_list->tail = new_node;
      offset_list->num_entries++;
    }
  }

  mdb_cursor_close(cur);
  mdb_txn_abort(txn);
  return BAMDB_SUCCESS;
}

int get_bam_rows(bam_row_set_t **output, const char *input_file_name,
                 const char *db_path, const char *index_name, const char *key) {
  samFile *input_file = 0;
  bam_hdr_t *header = NULL;
  offset_node_t *offset_node;
  offset_list_t *offsets;
  int i = 0;
  int rc = 0;
  int ret = BAMDB_SUCCESS;

  /* Always create an object for the caller */
  *output = calloc(1, sizeof(bam_row_set_t));

  if ((input_file = sam_open(input_file_name, "r")) == 0) {
    return BAMDB_SEQUENCE_FILE_ERROR;
  }

  header = sam_hdr_read(input_file);
  if (header == NULL) {
    return BAMDB_DB_ERROR;
  }

  offsets = calloc(1, sizeof(offset_list_t));
  rc = get_offsets_lmdb(offsets, db_path, index_name, key);
  if (rc != 0) {
    free(offsets);
    return rc;
  }

  (*output)->num_entries = offsets->num_entries;
  (*output)->rows = malloc((*output)->num_entries * sizeof(bam_row_set_t));

  offset_node = offsets->head;
  while (offset_node != NULL) {
    /* TODO: make sure we don't overrun row_set */
    ret = get_bam_row(&(*output)->rows[i], &(*output)->aux_tags,
                      offset_node->offset, input_file, header);
    offset_node = offset_node->next;
    i++;
  }

  free(offsets);
  return BAMDB_SUCCESS;
}
