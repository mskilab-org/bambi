
#include <lmdb.h>
#include "htslib/sam.h"

#include "include/bam_api.h"

//#include "htslib/sam.h"
//#include "htslib/bgzf.h"

//#include "include/bam_sqlite.h"
//#include "include/bam_lmdb.h"
//#include "include/bam_api.h"

#ifdef __cplusplus
extern "C" {
    #endif

  int convert_to_lmdb(samFile *input_file, char *db_path, int max_rows);
  //static int read_file(samFile *input_file, offset_list_t *offset_list);

    #ifdef __cplusplus
}
#endif
