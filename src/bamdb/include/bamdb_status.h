#ifndef BAMDBSTATUS_H
#define BAMDBSTATUS_H

/** Successful result */
#define BAMDB_SUCCESS 0
/** Error reading sequence file */
#define BAMDB_SEQUENCE_FILE_ERROR -16001
/** Error deserializing a sequence record */
#define BAMDB_DESERIALIZE_ERROR -16002
/** Error at the database level */
#define BAMDB_DB_ERROR -16102

#endif
