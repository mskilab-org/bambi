#define _GNU_SOURCE

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "htslib/sam.h"
#include "htslib/bgzf.h"

#include "include/bamdb.h"
#include "include/bam_sqlite.h"
#include "include/bam_lmdb.h"
#include "include/bam_api.h"
#include "include/bamdb_c.h"

/* I really hope we don't have sequences longer than this */
#define WORK_BUFFER_SIZE 65536

#define get_int_chars(i) ((i == 0) ? 1 : floor(log10(abs(i))) + 1)


void print_it(){
  printf("print it boys\n");
}

static void
get_bam_tags(const bam1_t *row, char *buffer)
{
  /* Output as TAG:TYPE:VALUE
   * TAG is two characters
   * TYPE is a single character
   * Example: QT:Z:AAFFFKKK */
  uint8_t *aux;
  uint8_t key[2];
  uint8_t type, sub_type;
  size_t buffer_pos = 0;
  uint32_t arr_size;
  char *dummy = 0;

  aux = bam_get_aux(row);
  while (aux+4 <= row->data + row->l_data) {
    key[0] = aux[0];
    key[1] = aux[1];
    sprintf(buffer + buffer_pos, "%c%c:", key[0], key[1]);
    buffer_pos += 3;

    type = aux[2];
    aux += 3;

    /* TODO: add error handling for values that don't conform to type */
    switch(type) {
    case 'A': /* Printable character */
      sprintf(buffer + buffer_pos, "A:%c", *aux);
      buffer_pos += 3;
      aux++;
      break;
    case 'C': /*Signed integer */
      sprintf(buffer + buffer_pos, "i:%d", *aux);
      buffer_pos += 2 + get_int_chars(*aux);
      aux++;
      break;
    case 'c':
      sprintf(buffer + buffer_pos, "i:%" PRId8, *(int8_t*)aux);
      buffer_pos += 2 + get_int_chars(*aux);
      aux++;
      break;
    case 'S':
      sprintf(buffer + buffer_pos, "i:%" PRIu16, *(uint16_t*)aux);
      buffer_pos += 2 + get_int_chars(*aux);
      aux += 2;
      break;
    case 's':
      sprintf(buffer + buffer_pos, "i:%" PRId16, *(int16_t*)aux);
      buffer_pos += 2 + get_int_chars(*aux);
      aux += 2;
      break;
    case 'I':
      sprintf(buffer + buffer_pos, "i:%" PRIu32, *(uint32_t*)aux);
      buffer_pos += 2 + get_int_chars(*aux);
      aux += 4;
      break;
    case 'i':
      sprintf(buffer + buffer_pos, "i:%" PRId32, *(int32_t*)aux);
      buffer_pos += 2 + get_int_chars(*aux);
      aux += 4;
      break;
    case 'f': /* Single precision floating point */
      sprintf(buffer + buffer_pos, "f:%g", *(float*)aux);
      /* Figure out how many chars the fp takes as a string */
      buffer_pos += 2 + snprintf(dummy, 0, "%g", *(float*)aux);
      aux += 4;
      break;
    case 'd':
      /* Double precision floating point. This does appear to be in the BAM spec,
       * I'm copying from samtools which does provide for this */
      sprintf(buffer + buffer_pos, "d:%g", *(float*)aux);
      /* Figure out how many chars the fp takes as a string */
      buffer_pos += 2 + snprintf(dummy, 0, "%g", *(float*)aux);
      aux += 4;
      break;
    case 'Z': /* Printable string */
    case 'H': /* Byte array */
      sprintf(buffer + buffer_pos, "%c:", type);
      buffer_pos += 2;
      while (aux < row->data + row->l_data && *aux) {
	sprintf(buffer + buffer_pos, "%c", *aux++);
	buffer_pos++;
      }
      aux++;
      break;
    case 'B': /* Integer or numeric array */
      sub_type = *(aux++);
      memcpy(&arr_size, aux, 4);

      sprintf(buffer + buffer_pos, "B:%c", sub_type);
      buffer_pos += 3;
      for (int i = 0; i < arr_size; ++i) {
	sprintf(buffer + buffer_pos, ",");
	buffer_pos++;
	switch (sub_type) {
	case 'c':
	  sprintf(buffer + buffer_pos, "%d", *aux);
	  buffer_pos += get_int_chars(*aux);
	  aux++;
	  break;
	case 'C':
	  sprintf(buffer + buffer_pos, "%" PRId8, *(int8_t*)aux);
	  buffer_pos += get_int_chars(*aux);
	  aux++;
	  break;
	case 'S':
	  sprintf(buffer + buffer_pos, "%" PRIu16, *(uint16_t*)aux);
	  buffer_pos += get_int_chars(*aux);
	  aux += 2;
	  break;
	case 's':
	  sprintf(buffer + buffer_pos, "%" PRId16, *(int16_t*)aux);
	  buffer_pos += get_int_chars(*aux);
	  aux += 2;
	  break;
	case 'I':
	  sprintf(buffer + buffer_pos, "i:%" PRIu32, *(uint32_t*)aux);
	  buffer_pos += 2 + get_int_chars(*aux);
	  aux += 4;
	  break;
	case 'i':
	  sprintf(buffer + buffer_pos, "i:%" PRId32, *(int32_t*)aux);
	  buffer_pos += 2 + get_int_chars(*aux);
	  aux += 4;
	  break;
	case 'f': /* Single precision floating point */
	  sprintf(buffer + buffer_pos, "f:%g", *(float*)aux);
	  /* Figure out how many chars the fp takes as a string */
	  buffer_pos += 2 + snprintf(dummy, 0, "%g", *(float*)aux);
	  aux += 4;
	  break;
	}
      }
      break;
    }

    sprintf(buffer + buffer_pos, "\t");
    buffer_pos++;
  }
}

static int
print_bam_row(const bam1_t *row, const bam_hdr_t *header, char *work_buffer, FILE *ifp)
{
  static uint32_t rows = 0;
  char *temp;

  //printf("Row %u:\n", rows);
  fprintf(ifp, "Row %u:\n", rows);

  //printf("\tQNAME: %s\n", bam_get_qname(row));
  //printf("\tFLAG: %u\n", row->core.flag);
  //printf("\tRNAME: %s\n", bam_get_rname(row, header));
  //printf("\tPOS: %d\n", row->core.pos);
  //printf("\tMAPQ: %u\n", row->core.qual);

  fprintf(ifp, "\tQNAME: %s\n", bam_get_qname(row));
  fprintf(ifp, "\tFLAG: %u\n", row->core.flag);
  fprintf(ifp, "\tRNAME: %s\n", bam_get_rname(row, header));
  fprintf(ifp, "\tPOS: %d\n", row->core.pos);
  fprintf(ifp, "\tMAPQ: %u\n", row->core.qual);
  
  temp = work_buffer;
  //printf("\tCIGAR: %s\n", bam_cigar_str(row, work_buffer));
  fprintf(ifp, "\tCIGAR: %s\n", bam_cigar_str(row, work_buffer));
  work_buffer = temp;

  //printf("\tRNEXT: %s\n", bam_get_rnext(row, header));
  //printf("\tPNEXT: %d\n", row->core.mpos + 1);
  //printf("\tTLEN: %d\n", row->core.isize);
  fprintf(ifp, "\tRNEXT: %s\n", bam_get_rnext(row, header));
  fprintf(ifp, "\tPNEXT: %d\n", row->core.mpos + 1);
  fprintf(ifp, "\tTLEN: %d\n", row->core.isize);
  
  temp = work_buffer;
  //printf("\tSEQ: %s\n", bam_seq_str(row, work_buffer));
  fprintf(ifp, "\tSEQ: %s\n", bam_seq_str(row, work_buffer));
  work_buffer = temp;

  temp = work_buffer;
  //printf("\tQUAL: %s\n", bam_qual_str(row, work_buffer));
  fprintf(ifp, "\tQUAL: %s\n", bam_qual_str(row, work_buffer));
  work_buffer = temp;

  temp = work_buffer;
  //printf("\tBX: %s\n", bam_bx_str(row, work_buffer));
  fprintf(ifp, "\tBX: %s\n", bam_bx_str(row, work_buffer));
  work_buffer = temp;

  /* TAGs */
  get_bam_tags(row, work_buffer);
  //printf("\tTAGs: %s\n", work_buffer);
  fprintf(ifp, "\tTAGs: %s\n", work_buffer);

  rows++;
  return 0;
}
 
int read_file(samFile *input_file, offset_list_t *offset_list)
{

  bam_hdr_t *header = NULL;
  bam1_t *bam_row;
  char *work_buffer = NULL;
  static int r = 0;
  int rc = 0;
  int64_t src = 0;
  offset_node_t *offset_node;

  header = sam_hdr_read(input_file);

  if (header == NULL) {
    fprintf(stderr, "Unable to read the header from %s\n", input_file->fn);
    rc = 1;
    return 1;
  }

  work_buffer = malloc(WORK_BUFFER_SIZE);
  if (offset_list != NULL) {
    offset_node = offset_list->head;
    // open file that will store BAM records that share barcode.
    FILE *ifp;
    printf("HAROOOOOOOOOOOOOO\n");
    ifp = fopen("queryResults.txt", "w");

    if (ifp == NULL){
      fprintf(stderr, "Can't open the input file\n");
      exit(1);
    }
    
    while (offset_node != NULL) {
      src = bgzf_seek(input_file->fp.bgzf, offset_node->offset, SEEK_SET);
      if (src != 0) {
	fprintf(stderr, "Error seeking to file offset\n");
	rc = 1;
	goto exit;
      }

      r = sam_read1(input_file, header, bam_row);
      print_bam_row(bam_row, header, work_buffer, ifp);
      offset_node = offset_node->next;
    }
  } else {
    while ((r = sam_read1(input_file, header, bam_row)) >= 0) { // read one alignment from `in'
      FILE *ifp;
      ifp = open("queryResults.txt", "w");
      print_bam_row(bam_row, header, work_buffer, ifp);
    }
    if (r < -1) {
      fprintf(stderr, "Attempting to process truncated file.\n");
      rc = 1;
      goto exit;
    }
  }


 exit:
  free(work_buffer);
  bam_destroy1(bam_row);
  return rc;
}

int main(int argc, char* argv[]){
  int rc = 0;
  int c;
  samFile *input_file = 0;
  
}
