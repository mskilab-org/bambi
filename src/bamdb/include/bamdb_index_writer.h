/* Only build if we want write functionality so we can avoid the CK dependency
 */
#ifdef BUILD_BAMDB_WRITER
/**
 * @file bamdb_index_writer.h
 * @brief Functions for creating a bamdb index
 *
 */
#ifndef BAMDB_INDEX_WRITER_H
#define BAMDB_INDEX_WRITER_H

#include "bamdb.h"

/** @brief Generate an lmdb based index for a given bam file
 *
 * This function will spawn a dedicated deserialization thread as well as one
 * thread per desired index in order to write in parallel and maximize potential
 * throughput. This can be highly memory and disk IO intensive, so it is
 * advisible to only run this on a machine with no other active workloads.
 *
 * @param[in] input_file The path of the bam file to index
 * @param[in] db_path Optional path of the generated index; a default path
 * based on the input filename will be used if this is NULL
 * @param[in] target_indices Struct containing the desired fields to be indexed
 * @return 0 on success or a non-zero error value on failure
 */
int generate_lmdb_index(samFile *input_file, char *db_path,
                        bamdb_indices_t *target_indices);

#endif
#endif
