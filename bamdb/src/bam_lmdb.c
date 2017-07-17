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

#include "bam_lmdb.h"
#include "bam_api.h"

#include "htslib/bgzf.h"
#include "htslib/hts.h"

#define LMDB_POSTFIX "_lmdb"
#define LMDB_COMMIT_FREQ 500000
#define WORK_BUFFER_SIZE 65536
#define LMDB_INIT_MAPSIZE 100000000000 /* Making this huge because we don't want to resize*/

#define DESERIALIZE_THREAD_COUNT 2
#define MAX_DESERIALIZE_QUEUE_SIZE 16384
#define MAX_WRITE_QUEUE_SIZE 16384

enum lmdb_keys {
	LMDB_QNAME = 0,
	LMDB_BX,
	LMDB_MAX
};

typedef struct _bam_data {
	bam1_t *bam_row;
	int64_t voffset;
} bam_data_t;

typedef struct write_entry {
	char *qname;
	char *bx;
	int64_t voffset;
} write_entry_t;

typedef struct _deserialize_thread_data {
	bam_hdr_t *header;
} deserialize_thread_data_t;

typedef struct _writer_thread_data {
	MDB_env *env;
	ck_fifo_mpmc_t *write_q;
} writer_thread_data_t;

static const char *lmdb_key_names[LMDB_MAX] = {
	"qname",
	"bx"
};

ck_fifo_mpmc_t *deserialize_q;
ck_fifo_mpmc_t *write_q;

int reader_running;
int deserialize_running;
int write_queue_size;
int deserialize_queue_size;

static char *
get_default_dbname(const char *filename)
{
	char *db_name;
	char *dot = strrchr(filename, '.');

	/* XXX: Leaking db_name */
	db_name = malloc(dot - filename + sizeof(LMDB_POSTFIX) + 1);
	strncpy(db_name, filename, dot - filename);
	strcpy(db_name + (dot - filename), LMDB_POSTFIX);

	return db_name;
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


static void *
deserialize_func(void *arg)
{
  (void) arg;    
	char *work_buffer = malloc(WORK_BUFFER_SIZE);
	char *buffer_pos = work_buffer;

	ck_fifo_mpmc_entry_t *garbage;
	bam_data_t *deserialize_entry;

	while (ck_pr_load_int(&reader_running) || CK_FIFO_MPMC_ISEMPTY(deserialize_q) == false) {
		while (ck_fifo_mpmc_trydequeue(deserialize_q, &deserialize_entry, &garbage) == true) {
			buffer_pos = work_buffer;

			write_entry_t *w_entry = malloc(sizeof(write_entry_t));
			ck_fifo_mpmc_entry_t *fifo_entry = malloc(sizeof(ck_fifo_mpmc_entry_t));
			w_entry->qname = strdup(bam_get_qname(deserialize_entry->bam_row));
			w_entry->bx = strdup(bam_bx_str(deserialize_entry->bam_row, buffer_pos));
			w_entry->voffset = deserialize_entry->voffset;

			ck_fifo_mpmc_enqueue(write_q, fifo_entry, w_entry);
			ck_pr_inc_int(&write_queue_size);
			ck_pr_dec_int(&deserialize_queue_size);

			while (ck_pr_load_int(&write_queue_size) > MAX_WRITE_QUEUE_SIZE) {
				usleep(100);
			}

			bam_destroy1(deserialize_entry->bam_row);
			free(garbage);
			free(deserialize_entry);
		}

		usleep(100);
	}

	ck_pr_dec_int(&deserialize_running);
	pthread_exit(NULL);
}


static void *
writer_func(void *arg)
{
	writer_thread_data_t *data = (writer_thread_data_t *)arg;

	MDB_val key, val;
	MDB_txn *txn = NULL;
	MDB_dbi dbi[LMDB_MAX];
	MDB_cursor *bx_cur;
	uint64_t n = 0;
	int rc;

	write_entry_t *entry;
	ck_fifo_mpmc_entry_t *garbage;

	uint64_t mapsize = LMDB_INIT_MAPSIZE;
	rc = mdb_env_set_mapsize(data->env, mapsize);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error setting map size: %s\n", mdb_strerror(rc));
	}

	rc = mdb_txn_begin(data->env, NULL, 0, &txn);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error starting transaction: %s\n", mdb_strerror(rc));
		return NULL;
	}

	rc = mdb_dbi_open(txn, lmdb_key_names[LMDB_BX], MDB_DUPSORT | MDB_CREATE | MDB_DUPFIXED, &dbi[LMDB_BX]);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error opening database: %s\n", mdb_strerror(rc));
		return NULL;
	}

	rc = mdb_cursor_open(txn, dbi[LMDB_BX], &bx_cur);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error getting cursor: %s\n", mdb_strerror(rc));
		return NULL;
	}

	while (ck_pr_load_int(&deserialize_running) || CK_FIFO_MPMC_ISEMPTY(deserialize_q) == false) {
		while (ck_fifo_mpmc_trydequeue(write_q, &entry, &garbage) == true) {
			free(garbage);

			val.mv_size = sizeof(int64_t);
			val.mv_data = &entry->voffset;

			/* insert voffset under qname */
			/*
			key.mv_size = strlen(entry->qname);
			key.mv_data = entry->qname;
			put_entry(data->env, txn, dbi[LMDB_QNAME], &key, &val, &mapsize);
			*/

			/* insert voffset under bx */
			key.mv_size = strlen(entry->bx);
			key.mv_data = entry->bx;
			rc = mdb_cursor_put(bx_cur, &key, &val, 0);
			if (rc != MDB_SUCCESS) {
				fprintf(stderr, "Error inserting data: %s\n", mdb_strerror(rc));
				return NULL;
			}

			free(entry->qname);
			free(entry->bx);
			free(entry);

			n++;
			ck_pr_dec_int(&write_queue_size);
			/* commit every so often for safety */
			if (n % LMDB_COMMIT_FREQ == 0) {
				mdb_cursor_close(bx_cur);
				commit_transaction(txn);
				printf("%" PRIu64 " records written. Deserialize queue: %d Write queue: %d\n", n, ck_pr_load_int(&deserialize_queue_size), ck_pr_load_int(&write_queue_size));
				rc = mdb_txn_begin(data->env, NULL, 0, &txn);
				if (rc != MDB_SUCCESS) {
					fprintf(stderr, "Error starting transaction: %s\n", mdb_strerror(rc));
					return NULL;
				}

				rc = mdb_dbi_open(txn, lmdb_key_names[LMDB_BX], MDB_DUPSORT | MDB_DUPFIXED, &dbi[LMDB_BX]);
				if (rc != MDB_SUCCESS) {
					fprintf(stderr, "Error opening database: %s\n", mdb_strerror(rc));
					return NULL;
				}

				rc = mdb_cursor_open(txn, dbi[LMDB_BX], &bx_cur);
				if (rc != MDB_SUCCESS) {
					fprintf(stderr, "Error getting cursor: %s\n", mdb_strerror(rc));
					return NULL;
				}
			}
		}

		usleep(100);
	}

	mdb_cursor_close(bx_cur);
	commit_transaction(txn);
	mdb_env_sync(data->env, 1);
	for (int i = 0; i < LMDB_MAX; i++) {
		mdb_dbi_close(data->env, dbi[i]);
	}

	struct MDB_stat stats;
	mdb_env_stat(data->env, &stats);

	mdb_env_close(data->env);

	pthread_exit(NULL);
}

int
convert_to_lmdb(samFile *input_file, char *db_name)
{
	MDB_env *env;
	int rc;
	int r = 0;
	int ret = 0;
	bam_hdr_t *header = NULL;

	writer_thread_data_t writer_thread_args;
	pthread_t threads[DESERIALIZE_THREAD_COUNT + 1];

	write_q = calloc(1, sizeof(ck_fifo_mpmc_t));
	ck_fifo_mpmc_init(write_q, malloc(sizeof(ck_fifo_mpmc_entry_t)));

	deserialize_q = calloc(1, sizeof(ck_fifo_mpmc_t));
	ck_fifo_mpmc_init(deserialize_q, malloc(sizeof(ck_fifo_mpmc_entry_t)));

	reader_running = 1;
	deserialize_running = DESERIALIZE_THREAD_COUNT;
	write_queue_size = 0;
	deserialize_queue_size = 0;

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

	// rc = mdb_env_open(env, db_name, MDB_WRITEMAP | MDB_MAPASYNC | MDB_NOLOCK, 0664);
	rc = mdb_env_open(env, db_name, MDB_NOLOCK, 0664);
	// rc = mdb_env_open(env, db_name, MDB_WRITEMAP | MDB_NOLOCK, 0664);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error opening env: %s\n", mdb_strerror(rc));
		return 1;
	}

	header = sam_hdr_read(input_file);
	if (header == NULL) {
		fprintf(stderr, "Unable to read the header from %s\n", input_file->fn);
		ret = 1;
		goto exit;
	}

	for (size_t i = 0; i < DESERIALIZE_THREAD_COUNT; i++) {
		deserialize_thread_data_t deserialize_thread_args;
		deserialize_thread_args.header = header;
		rc = pthread_create(&threads[i], NULL, deserialize_func, &deserialize_thread_args);
		if (rc != 0) {
			fprintf(stderr, "Received non-zero return code when launching deserialize thread: %d\n", rc);
			ret = 1;
			goto exit;
		}
	}

	writer_thread_args.env = env;
	rc = pthread_create(&threads[DESERIALIZE_THREAD_COUNT], NULL, writer_func, &writer_thread_args);
	if (rc != 0) {
		fprintf(stderr, "Received non-zero return code when launching writer thread: %d\n", rc);
		ret = 1;
		goto exit;
	}

	while (r >= 0) {
		bam_data_t *entry = malloc(sizeof(bam_data_t));
		ck_fifo_mpmc_entry_t *fifo_entry = malloc(sizeof(ck_fifo_mpmc_entry_t));
		entry->bam_row = bam_init1();

		while (ck_pr_load_int(&deserialize_queue_size) > MAX_DESERIALIZE_QUEUE_SIZE) {
			usleep(100);
		}

		r = sam_read1(input_file, header, entry->bam_row);
		entry->voffset = bgzf_tell(input_file->fp.bgzf);
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

	/* Wait for writer */
	pthread_join(threads[DESERIALIZE_THREAD_COUNT], NULL);
exit:
	return ret;
}


int
get_offsets(offset_list_t *offset_list, const char *lmdb_db_name, const char *bx)
{
	MDB_env *env;
	MDB_dbi dbi[LMDB_MAX];
	MDB_cursor *cur;
	MDB_txn *txn;
	MDB_val key, data;
	int rc;
	int rows = 0;

	rc = mdb_env_create(&env);
	rc = mdb_env_set_maxdbs(env, LMDB_MAX);
	rc = mdb_env_open(env, lmdb_db_name, MDB_RDONLY, 0644);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error opening env: %s\n", mdb_strerror(rc));
		return -1;
	}

	rc = mdb_txn_begin(env, NULL, MDB_RDONLY, &txn);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error beginning transaction: %s\n", mdb_strerror(rc));
		return -1;
	}

	rc = mdb_dbi_open(txn, lmdb_key_names[LMDB_BX],
			MDB_DUPSORT | MDB_DUPFIXED, &dbi[LMDB_BX]);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error opening database %s\n", mdb_strerror(rc));
		return -1;
	}

	key.mv_size = strlen(bx);
	key.mv_data = strdup(bx);

	rc = mdb_cursor_open(txn, dbi[LMDB_BX], &cur);
	if (rc != MDB_SUCCESS) {
		fprintf(stderr, "Error getting cursor: %s\n", mdb_strerror(rc));
		return -1;
	}

	if ((rc = mdb_cursor_get(cur, &key, &data, MDB_SET)) != MDB_SUCCESS) {
		fprintf(stderr, "Error getting key: %s\n", mdb_strerror(rc));
		return -1;
	}
	if ((rc = mdb_cursor_get(cur, &key, &data, MDB_FIRST_DUP)) == 0) {
		offset_node_t *new_node = calloc(1, sizeof(offset_node_t));
		new_node->offset = *(int64_t *)data.mv_data;

		offset_list->head = new_node;
		offset_list->tail = new_node;
        rows++;
		
		while ((rc = mdb_cursor_get(cur, &key, &data, MDB_NEXT_DUP)) == 0) {
			offset_node_t *new_node = calloc(1, sizeof(offset_node_t));
			new_node->offset = *(int64_t *)data.mv_data;

			offset_list->tail->next = new_node;
			offset_list->tail = new_node;
			rows++;
		}
	}

	mdb_cursor_close(cur);
	mdb_txn_abort(txn);
	return rows;
}


bam_row_set_t *
get_bx_rows(char *input_file_name, char *db_path, char *bx)
{
	samFile *input_file = 0;
	bam_hdr_t *header = NULL;
	int n_rows = 0;
	offset_list_t *offset_list = NULL;
	offset_node_t *offset_node;
	bam_row_set_t *row_set = NULL;
	int i = 0;

	if ((input_file = sam_open(input_file_name, "r")) == 0) {
		fprintf(stderr, "Unable to open file %s\n", input_file_name);
		return NULL;
	}

	header = sam_hdr_read(input_file);
	if (header == NULL) {
		fprintf(stderr, "Unable to read the header from %s\n", input_file->fn);
		return NULL;
	}

	offset_list = calloc(1, sizeof(offset_list_t));
	n_rows = get_offsets(offset_list, db_path, bx);
	if (n_rows <= 0) {
		return NULL;
	}

	row_set = malloc(sizeof(bam_row_set_t));
	row_set->n_entries = n_rows;
	row_set->rows = malloc(n_rows * sizeof(bam_row_set_t));

	offset_node = offset_list->head;
	while (offset_node != NULL) {
		/* TODO: make sure we don't overrun row_set */
		row_set->rows[i] = get_bam_row(offset_node->offset, input_file, header);
		offset_node = offset_node->next;
		i++;
	}

	return row_set;
}
