/* API for converting between BAM and SQLite */

#include "sqlite3.h"

#include "htslib/sam.h"
#include "bam_api.h"


#ifdef __cplusplus
extern "C" {
    #endif
  void print_bam_sqlite();
  int convert_to_sqlite(samFile *input_file, char *db_name, int max_rows);
  int get_offsets(offset_list_t *offset_list, const char *sqlite_db_name, const char *bx);

  #ifdef __cplusplus
}
#endif
