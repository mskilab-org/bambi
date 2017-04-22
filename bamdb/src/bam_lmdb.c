#include <lmdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "bam_lmdb.h"
#include "bam_api.h"

#include "htslib/bgzf.h"
#include "htslib/hts.h"

#define LMDB_POSTFIX "_lmdb"
#define LMDB_COMMIT_FREQ 100000
#define LMDB_INIT_MAPSIZE 1024 * 1024 * 1024
#define WORK_BUFFER_SIZE 65536

enum lmdb_keys {
  LMDB_QNAME = 0,
  LMDB_BX,
  LMDB_MAX
};

static const char *lmdb_key_names[LMDB_MAX] = {
  "qname",
  "bx"
};


static char *
get_default_dbname(char *filename)
{
  char *db_name;
  char *dot = strrchr(filename, '.');

  db_name = (char *)malloc(dot - filename + sizeof(LMDB_POSTFIX) + 1); /* Leak this */
  strncpy(db_name, filename, dot - filename);
  strcpy(db_name + (dot - filename), LMDB_POSTFIX);

  return db_name;
}

static int
start_transaction(MDB_env *env, MDB_dbi *dbi, MDB_txn **txn)
{
  int rc;

  rc = mdb_txn_begin(env, NULL, 0, txn);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error beginning transaction: %s\n", mdb_strerror(rc));
    return 1;
  }

  for (int i = 0; i < LMDB_MAX; i++) {
    rc = mdb_dbi_open(*txn, lmdb_key_names[i], MDB_DUPSORT | MDB_CREATE, &dbi[i]);
    if (rc != MDB_SUCCESS) {
      fprintf(stderr, "Error opening database: %s\n", mdb_strerror(rc));
      return NULL;
    }
  }

  return 0;
}


static int
commit_transaction(MDB_txn *txn)
{
  int rc;

  rc = mdb_txn_commit(txn);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error commiting transaction: %s\n", mdb_strerror(rc));
    return 1;
  }

  return 0;
}


static int
put_entry(MDB_env *env, MDB_txn *txn, MDB_dbi dbi, MDB_val *key, MDB_val *data, int *mapsize)
{
  int rc = mdb_put(txn, dbi, key, data, 0);

  while (rc != MDB_SUCCESS) {
    if (rc == MDB_MAP_FULL) {
      *mapsize *= 2;
      mdb_env_set_mapsize(env, *mapsize);
    } else {
      fprintf(stderr, "Error putting entry: %s\n", mdb_strerror(rc));
      return 1;
    }

    rc = mdb_put(txn, dbi, key, data, 0);
  }

  return 0;
}

int
convert_to_lmdb(samFile *input_file, char *db_name, int max_rows)
{
  MDB_env *env;
  MDB_dbi dbi[LMDB_MAX];
  MDB_val key, data;
  MDB_txn *txn = NULL;
  int mapsize = LMDB_INIT_MAPSIZE;
  int rc, r;
  int ret = 0;

  bam_hdr_t *header = NULL;
  bam1_t *bam_row;
  char *work_buffer = (char *)malloc(WORK_BUFFER_SIZE);
  char *buffer_pos = work_buffer;

  if (db_name == NULL) {
    db_name = get_default_dbname(input_file->fn);
  }
  printf("Attempting to convert bam file %s into lmdb database at path %s\n", input_file->fn, db_name);

  mkdir(db_name, 0777);

  rc = mdb_env_create(&env);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error creating env: %s\n", mdb_strerror(rc));
    return 1;
  }

  rc = mdb_env_set_maxdbs(env, LMDB_MAX);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error setting maxdbs: %s\n", mdb_strerror(rc));
    return 1;
  }

  rc = mdb_env_set_mapsize(env, mapsize);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error setting map size: %s\n", mdb_strerror(rc));
    return 1;
  }

  rc = mdb_env_open(env, db_name, 0, 0664);
  if (rc != MDB_SUCCESS) {
    fprintf(stderr, "Error opening env: %s\n", mdb_strerror(rc));
    return 1;
  }

  if (start_transaction(env, dbi, &txn) != 0) {
    ret = 1;
    goto exit;
  }

  header = sam_hdr_read(input_file);
  if (header == NULL) {
    fprintf(stderr, "Unable to read the header from %s\n", input_file->fn);
    ret = 1;
    goto exit;
  }

  bam_row = bam_init1();
  uint32_t n = 0;
  while ((r = sam_read1(input_file, header, bam_row)) >= 0) {
    buffer_pos = work_buffer;

    int64_t voffset = bgzf_tell(input_file->fp.bgzf);
    char *qname = bam_get_qname(bam_row);
    char *bx = bam_bx_str(bam_row, buffer_pos);

    data.mv_size = sizeof(int64_t);
    data.mv_data = &voffset;

    /* insert voffset under qname */
    key.mv_size = strlen(qname);
    key.mv_data = qname;

    if (put_entry(env, txn, dbi[LMDB_QNAME], &key, &data, &mapsize) != 0) {
      ret = 1;
      goto exit;
    }

    /* insert voffset under bx */
    key.mv_size = strlen(bx);
    key.mv_data = bx;

    if (put_entry(env, txn, dbi[LMDB_BX], &key, &data, &mapsize) != 0) {
      ret = 1;
      goto exit;
    }

    if (++n % 100000 == 0) printf("%u rows inserted\n", n);

    if (max_rows > 0 && n >= max_rows) {
      break;
    }

    /* commit every so often for safety */
    if (n % LMDB_COMMIT_FREQ == 0) {
      if (commit_transaction(txn) != 0) {
	ret = 1;
	goto exit;
      }

      if (start_transaction(env, dbi, &txn) != 0) {
	ret = 1;
	goto exit;
      }
    }
  }

  if (r < -1) {
    fprintf(stderr, "Attempting to process truncated file.\n");
    ret = 1;
    goto exit;
  }

  commit_transaction(txn);

 exit:
  for (int i = 0; i < LMDB_MAX; i++) {
    mdb_dbi_close(env, dbi[i]);
  }
  mdb_env_close(env);

  return ret;
}
