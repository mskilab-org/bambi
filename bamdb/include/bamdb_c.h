
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "htslib/sam.h"
#include "htslib/bgzf.h"

#include "include/bam_api.h"

#ifdef __cplusplus
extern "C" {
    #endif


  void print_it();
  static void get_bam_tags(const bam1_t *row, char *buffer);
  static int print_bam_row(const bam1_t *row, const bam_hdr_t *header, char *work_buffer);
  int read_file(samFile *input_file, offset_list_t *offset_list);

    #ifdef __cplusplus
}
#endif
