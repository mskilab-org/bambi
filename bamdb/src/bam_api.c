#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bam_api.h"

#define get_int_chars(i) ((i == 0) ? 1 : floor(log10(abs(i))) + 1)

const char *
bam_get_rname(const bam1_t *row, const bam_hdr_t *header)
{
	const char *rname;
	if (row->core.tid >= 0) {
		rname = header->target_name[row->core.tid];
	} else {
		rname = "*";
	}

	return rname;
}


const char *
bam_get_rnext(const bam1_t *row, const bam_hdr_t *header)
{
	char *rnext;
	if (row->core.mtid < 0) {
		rnext = "*";
	} else if (row->core.mtid == row->core.tid) {
		/* "=" to save space, could also use full reference name */
		rnext = "=";
	} else {
		rnext = header->target_name[row->core.mtid];
	}

	return rnext;
}


char *
bam_cigar_str(const bam1_t *row, char *work_buffer)
{
	char *ret = work_buffer;
	/* CIGAR is an array of uint32s. First 4 bits are the CIGAR opration
	 * and the following 28 bits are the number of repeats of the op */
	if (row->core.n_cigar > 0) {
		uint32_t *cigar = bam_get_cigar(row);
		size_t cigar_pos = 0;

		for (int i = 0; i < row->core.n_cigar; ++i) {
			int num_ops = bam_cigar_oplen(cigar[i]);
			char op = bam_cigar_opchr(cigar[i]);
			sprintf(work_buffer + cigar_pos, "%d%c", num_ops, op);
			/* Need enough space for the opcount and the actual opchar
			 * an example would be something like 55F */
			cigar_pos += get_int_chars(num_ops) + 1;
		}

		work_buffer[cigar_pos] = '\0';
		work_buffer += cigar_pos + 1;
	} else {
		strcpy(work_buffer, "*\0");
		work_buffer += 2;
	}

	return ret;
}


char *
bam_seq_str(const bam1_t *row, char *work_buffer)
{
	char *ret = work_buffer;
	int32_t l_qseq = row->core.l_qseq;
	if (l_qseq > 0) {
		uint8_t *seq = bam_get_seq(row);
		for (int i = 0; i < l_qseq; ++i) {
			/* Need to unpack base letter from its 4 bit encoding */
			work_buffer[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)];
		}
		work_buffer[l_qseq] = '\0';
		work_buffer += l_qseq + 1;
	} else {
		strcpy(work_buffer, "*\0");
		work_buffer += 2;
	}

	return ret;
}


char *
bam_qual_str(const bam1_t *row, char *work_buffer)
{
	char *ret = work_buffer;
	int32_t l_qseq = row->core.l_qseq;
	uint8_t *qual = bam_get_qual(row);
	if (qual[0] == 0xff || l_qseq <= 0) {
		/* No quality data, indicate with a '*' */
		strcpy(work_buffer, "*\0");
		work_buffer += 2;
	} else {
		for (int i = 0; i < l_qseq; ++i) {
			/* ASCII of base QUALity plus 33 */
			work_buffer[i] = qual[i] + 33;
		}
		work_buffer[l_qseq] = '\0';
		work_buffer += l_qseq + 1;
	}

	return ret;
}


char *
bam_bx_str(const bam1_t *row, char *work_buffer)
{
	/* Hoping that BXZ doesn't appear for any other purposes in the aux field */
	uint8_t *aux;
	char *ret = work_buffer;
	uint8_t *bx_pos = 0;

	aux = bam_get_aux(row);

	while (aux+4 <= row->data + row->l_data) {
		if (aux[0] == 'B' && aux[1] == 'X' && aux[2] == 'Z') {
			bx_pos = aux + 3;
			break;
		}
		aux++;
	}

	if (bx_pos != NULL) {
		ret = work_buffer;
		while (bx_pos < row->data + row->l_data && *bx_pos) {
			sprintf(work_buffer, "%c", *bx_pos++);
			work_buffer++;
		}
	} else {
		work_buffer[0] = '*';
		work_buffer[1] = '\0';
		work_buffer += 2;
	}

	return ret;
}
