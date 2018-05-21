#include <ck_fifo.h>
#include <ck_pr.h>

#include <inttypes.h>
#include <limits.h>
#include <lmdb.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "bam_api.h"
#include "bam_lmdb.h"
#include "bamdb_status.h"

#include "htslib/bgzf.h"
#include "htslib/hts.h"

#define LMDB_POSTFIX "_lmdb"
#define LMDB_COMMIT_FREQ 500000
#define WORK_BUFFER_SIZE 65536
/* Making this huge because so we don't have to resize while running*/
#define LMDB_INIT_MAPSIZE 100000000000
#define DESERIALIZE_THREAD_COUNT 1
#define MAX_DESERIALIZE_QUEUE_SIZE 16384
#define MAX_WRITE_QUEUE_SIZE 16384
#define MAX_PATH_CHARS 2048

typedef struct _bam_data {
  bam1_t *bam_row;
  int64_t voffset;
} bam_data_t;

typedef struct write_entry {
  char *key;
  int64_t voffset;
} write_entry_t;

typedef struct writer_q {
  char *key;
  ck_fifo_mpmc_t *write_q;
} writer_q_t;

typedef struct _deserialize_thread_data {
  bam_hdr_t *header;
  size_t num_keys;
  writer_q_t **write_queues;
} deserialize_thread_data_t;

typedef struct _writer_thread_data {
  ck_fifo_mpmc_t *write_q;
  char *key_name;
  char *db_path;
} writer_thread_data_t;

ck_fifo_mpmc_t *deserialize_q;

int reader_running;
int deserialize_running;
int write_queue_size;
int deserialize_queue_size;

static char *get_default_dbname(const char *filename) {
  char *db_name;
  char *dot = strrchr(filename, '.');

  /* This needs to be freed by the caller */
  db_name = malloc(dot - filename + sizeof(LMDB_POSTFIX) + 1);
  strncpy(db_name, filename, dot - filename);
  strcpy(db_name + (dot - filename), LMDB_POSTFIX);

  return db_name;
}

static int commit_transaction(MDB_txn *txn) {
  int rc;

  rc = mdb_txn_commit(txn);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error commiting transaction: %s\n", mdb_strerror(rc));
    return 1;
  }

  return 0;
}

static void *deserialize_func(void *arg) {
  deserialize_thread_data_t *data = (deserialize_thread_data_t *)arg;
  char *work_buffer = malloc(WORK_BUFFER_SIZE);
  char *buffer_pos = work_buffer;

  ck_fifo_mpmc_entry_t *garbage;
  bam_data_t *deserialize_entry;

  while (ck_pr_load_int(&reader_running) ||
         CK_FIFO_MPMC_ISEMPTY(deserialize_q) == false) {
    while (ck_fifo_mpmc_trydequeue(deserialize_q, &deserialize_entry,
                                   &garbage) == true) {
      buffer_pos = work_buffer;

      for (size_t i = 0; i < data->num_keys; i++) {
        write_entry_t *w_entry = malloc(sizeof(write_entry_t));
        ck_fifo_mpmc_entry_t *fifo_entry = malloc(sizeof(ck_fifo_mpmc_entry_t));
        char *target_key = data->write_queues[i]->key;

        if (strncmp(target_key, "QNAME", 5) == 0) {
          w_entry->key = strdup(bam_get_qname(deserialize_entry->bam_row));
        } else {
          w_entry->key = strdup(
              bam_str_key(deserialize_entry->bam_row, target_key, buffer_pos));
        }

        w_entry->voffset = deserialize_entry->voffset;

        ck_fifo_mpmc_enqueue(data->write_queues[i]->write_q, fifo_entry,
                             w_entry);
        ck_pr_inc_int(&write_queue_size);

        while (ck_pr_load_int(&write_queue_size) > MAX_WRITE_QUEUE_SIZE) {
          usleep(100);
        }
      }

      ck_pr_dec_int(&deserialize_queue_size);
      bam_destroy1(deserialize_entry->bam_row);
      free(garbage);
      free(deserialize_entry);
    }

    usleep(100);
  }

  ck_pr_dec_int(&deserialize_running);
  pthread_exit(NULL);
}

static MDB_env *set_up_lmdb_env(const char *full_db_path, bool read_only) {
  int rc;
  int flags;
  MDB_env *env = NULL;

  rc = mdb_env_create(&env);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error creating env: %s\n", mdb_strerror(rc));
    return NULL;
  }

  rc = mdb_env_set_maxdbs(env, 1);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error setting maxdbs: %s\n", mdb_strerror(rc));
    return NULL;
  }

  /* We don't ever read while writing */
  flags = MDB_NOLOCK;
  if (read_only) {
    flags = flags | MDB_RDONLY;
  }

  rc = mdb_env_open(env, full_db_path, flags, 0664);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error opening env: %s\n", mdb_strerror(rc));
    return NULL;
  }

  rc = mdb_env_set_mapsize(env, LMDB_INIT_MAPSIZE);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error setting map size: %s\n", mdb_strerror(rc));
    return NULL;
  }

  return env;
}

static void *writer_func(void *arg) {
  writer_thread_data_t *data = (writer_thread_data_t *)arg;

  MDB_env *env = NULL;
  MDB_val key, val;
  MDB_txn *txn = NULL;
  MDB_dbi dbi;
  MDB_cursor *cur = NULL;
  char target_path[MAX_PATH_CHARS];
  uint64_t n = 0;
  int rc;

  write_entry_t *entry;
  ck_fifo_mpmc_entry_t *garbage;

  snprintf(target_path, MAX_PATH_CHARS, "%s/%s", data->db_path, data->key_name);
  mkdir(target_path, 0777);
  env = set_up_lmdb_env(target_path, false);
  if (env == NULL) {
    return NULL;
  }

  rc = mdb_txn_begin(env, NULL, 0, &txn);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error starting transaction: %s\n", mdb_strerror(rc));
    return NULL;
  }

  rc = mdb_dbi_open(txn, NULL, MDB_DUPSORT | MDB_CREATE | MDB_DUPFIXED, &dbi);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error opening database: %s\n", mdb_strerror(rc));
    return NULL;
  }

  rc = mdb_cursor_open(txn, dbi, &cur);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error getting cursor: %s\n", mdb_strerror(rc));
    return NULL;
  }

  while (ck_pr_load_int(&deserialize_running) ||
         CK_FIFO_MPMC_ISEMPTY(data->write_q) == false) {
    while (ck_fifo_mpmc_trydequeue(data->write_q, &entry, &garbage) == true) {
      free(garbage);

      val.mv_size = sizeof(int64_t);
      val.mv_data = &entry->voffset;

      /* insert voffset under bx */
      key.mv_size = strlen(entry->key);
      key.mv_data = entry->key;

      rc = mdb_cursor_put(cur, &key, &val, 0);

      if (rc != MDB_SUCCESS) {
        fprintf(stderr, "Error inserting data: %s\n", mdb_strerror(rc));
        return NULL;
      }

      free(entry->key);
      free(entry);

      n++;
      ck_pr_dec_int(&write_queue_size);
      /* commit every so often for safety */
      if (n % LMDB_COMMIT_FREQ == 0) {
        mdb_cursor_close(cur);
        commit_transaction(txn);
        printf("%" PRIu64
               " records written. Deserialize queue: %d Write queue: %d\n",
               n, ck_pr_load_int(&deserialize_queue_size),
               ck_pr_load_int(&write_queue_size));
        rc = mdb_txn_begin(env, NULL, 0, &txn);
        if (rc != MDB_SUCCESS) {
          fprintf(stderr, "Error starting transaction: %s\n", mdb_strerror(rc));
          return NULL;
        }

        rc = mdb_cursor_open(txn, dbi, &cur);
        if (rc != MDB_SUCCESS) {
          fprintf(stderr, "Error getting cursor: %s\n", mdb_strerror(rc));
          return NULL;
        }
      }
    }

    usleep(100);
  }

  mdb_cursor_close(cur);
  commit_transaction(txn);
  mdb_env_sync(env, 1);
  mdb_dbi_close(env, dbi);
  mdb_env_close(env);

  pthread_exit(NULL);
}

static writer_q_t *init_writer_q(char *key) {
  if (strlen(key) != 2 && strncmp("QNAME", key, 5) != 0) {
    fprintf(stderr, "Target indices must be QNAME or a two letter string key");
    return NULL;
  }

  writer_q_t *new_queue = malloc(sizeof(writer_q_t));
  new_queue->key = strndup(key, 5);
  new_queue->write_q = calloc(1, sizeof(ck_fifo_mpmc_t));
  ck_fifo_mpmc_init(new_queue->write_q, malloc(sizeof(ck_fifo_mpmc_entry_t)));

  return new_queue;
}

int convert_to_lmdb(samFile *input_file, char *db_path,
                    indices_t *target_indices) {
  int rc;
  int r = 0;
  int ret = 0;
  bam_hdr_t *header = NULL;
  bool default_db_path = false;

  deserialize_thread_data_t deserialize_thread_args;

  size_t total_indices = target_indices->num_key_indices +
                         (target_indices->includes_qname ? 1 : 0);
  /* One writer thread per index */
  size_t n_threads = DESERIALIZE_THREAD_COUNT + total_indices;
  pthread_t *threads = calloc(n_threads, sizeof(pthread_t));

  deserialize_thread_args.num_keys = total_indices;
  deserialize_thread_args.write_queues =
      calloc(total_indices, sizeof(writer_q_t));

  for (size_t i = 0; i < target_indices->num_key_indices; i++) {
    writer_q_t *new_queue = init_writer_q(target_indices->key_indices[i]);

    if (new_queue == NULL) {
      ret = 1;
      goto exit;
    }

    deserialize_thread_args.write_queues[i] = new_queue;
  }

  if (target_indices->includes_qname) {
    deserialize_thread_args.write_queues[target_indices->num_key_indices] =
        init_writer_q("QNAME");
  }

  deserialize_q = calloc(1, sizeof(ck_fifo_mpmc_t));
  ck_fifo_mpmc_init(deserialize_q, malloc(sizeof(ck_fifo_mpmc_entry_t)));

  reader_running = 1;
  deserialize_running = DESERIALIZE_THREAD_COUNT;
  write_queue_size = 0;
  deserialize_queue_size = 0;

  if (db_path == NULL) {
    db_path = get_default_dbname(input_file->fn);
    default_db_path = true;
  }
  mkdir(db_path, 0777);
  fprintf(stdout,
          "Attempting to convert bam file %s into lmdb database at path %s\n",
          input_file->fn, db_path);

  header = sam_hdr_read(input_file);
  if (header == NULL) {
    fprintf(stderr, "Unable to read the header from %s\n", input_file->fn);
    ret = 1;
    goto exit;
  }

  for (size_t i = 0; i < DESERIALIZE_THREAD_COUNT; i++) {
    deserialize_thread_args.header = header;
    rc = pthread_create(&threads[i], NULL, deserialize_func,
                        &deserialize_thread_args);
    if (rc != 0) {
      fprintf(stderr,
              "Received non-zero return code when launching deserialize "
              "thread: %d\n",
              rc);
      ret = 1;
      goto exit;
    }
  }

  /* Init writer threads, one per key */
  for (size_t i = 0; i < deserialize_thread_args.num_keys; i++) {
    // TODO: free these after we're done
    writer_thread_data_t *new_writer_args =
        malloc(sizeof(writer_thread_data_t));
    new_writer_args->write_q = deserialize_thread_args.write_queues[i]->write_q;
    new_writer_args->key_name = deserialize_thread_args.write_queues[i]->key;
    new_writer_args->db_path = db_path;

    rc = pthread_create(&threads[DESERIALIZE_THREAD_COUNT + i], NULL,
                        writer_func, new_writer_args);
    if (rc != 0) {
      fprintf(
          stderr,
          "Received non-zero return code when launching writer thread: %d\n",
          rc);
      ret = 1;
      goto exit;
    }
  }

  /* This thread serves as the reader thread */
  while (r >= 0) {
    bam_data_t *entry = malloc(sizeof(bam_data_t));
    ck_fifo_mpmc_entry_t *fifo_entry = malloc(sizeof(ck_fifo_mpmc_entry_t));
    entry->bam_row = bam_init1();

    while (ck_pr_load_int(&deserialize_queue_size) >
           MAX_DESERIALIZE_QUEUE_SIZE) {
      usleep(100);
    }

    entry->voffset = bgzf_tell(input_file->fp.bgzf);
    r = sam_read1(input_file, header, entry->bam_row);
    if (r >= 0) {
      ck_fifo_mpmc_enqueue(deserialize_q, fifo_entry, entry);
      ck_pr_inc_int(&deserialize_queue_size);
    }
  }
  ck_pr_dec_int(&reader_running);

  if (r < -1) {
    fprintf(stderr, "Attempting to process truncated file.\n");
    ret = -1;
    goto exit;
  }
  for (size_t i = 0; i < DESERIALIZE_THREAD_COUNT; i++) {
    pthread_join(threads[i], NULL);
  }

  /* Wait for writers */
  pthread_join(threads[DESERIALIZE_THREAD_COUNT], NULL);
exit:
  free(threads);
  if (default_db_path) {
    free(db_path);
  }

  return ret;
}

static int open_ro_handle(MDB_env *output, MDB_dbi *dbi, MDB_cursor **cur,
                          MDB_txn *txn, const char *key_name,
                          const char *db_path) {
  char target_path[MAX_PATH_CHARS];
  int rc = 0;

  snprintf(target_path, MAX_PATH_CHARS, "%s/%s", db_path, key_name);
  output = set_up_lmdb_env(target_path, true);
  if (output == NULL) {
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
