
#include <lmdb.h>

#include "htslib/sam.h"


int convert_to_lmdb(samFile *input_file, char *db_path, int max_rows);

