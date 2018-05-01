#include <string>
#include <vector>

#include <Rcpp.h>

extern "C" {
    #include "bam_lmdb.h"
    bam_row_set_t *get_bam_rows(const char *input_file_name, const char *db_path, const char *key_type, const char *key_value);
}

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame
query_bam_index(CharacterVector bam_file_name, CharacterVector bam_index_path,
  CharacterVector key_type, CharacterVector key_value)
{
    bam_row_set_t *bam_rows = NULL;

    /* Series of casts get from R strings to C strings */
    std::string file_name_str = Rcpp::as<std::string>(bam_file_name);
    std::string index_path_str = Rcpp::as<std::string>(bam_index_path);
    std::string key_type_str = Rcpp::as<std::string>(key_type);
    std::string key_val_str = Rcpp::as<std::string>(key_value);
    const char* c_file_name = file_name_str.c_str();
    const char* c_index_path = index_path_str.c_str();
    const char* c_key_type = key_type_str.c_str();
    const char* c_key_val = key_val_str.c_str();

    bam_rows = get_bam_rows(c_file_name, c_index_path, c_key_type, c_key_val);

    NumericVector w = NumericVector::create(9,10);
    
    std::vector<std::string> qname;
    std::vector<int> flag;
    std::vector<std::string> rname;
    std::vector<int> pos;
    std::vector<int> mapq;
    std::vector<std::string> cigar;
    std::vector<std::string> rnext;
    std::vector<int> pnext;
    std::vector<int> tlen;
    std::vector<std::string> seq;
    std::vector<std::string> qual;

    for (size_t i = 0; i < bam_rows->n_entries; i++) {
        qname.push_back(std::string(bam_rows->rows[i]->qname));
        flag.push_back(bam_rows->rows[i]->flag);
        rname.push_back(std::string(bam_rows->rows[i]->rname));
        pos.push_back(bam_rows->rows[i]->pos);
        mapq.push_back(bam_rows->rows[i]->mapq);
        cigar.push_back(std::string(bam_rows->rows[i]->cigar));
        rnext.push_back(std::string(bam_rows->rows[i]->rnext));
        pnext.push_back(bam_rows->rows[i]->pnext);
        tlen.push_back(bam_rows->rows[i]->tlen);
        seq.push_back(std::string(bam_rows->rows[i]->seq));
        qual.push_back(std::string(bam_rows->rows[i]->qual));
    }

    StringVector qname_r(qname.begin(), qname.end());
    NumericVector flag_r(flag.begin(), flag.end());
    StringVector rname_r(rname.begin(), rname.end());
    NumericVector pos_r(pos.begin(), pos.end());
    NumericVector mapq_r(mapq.begin(), mapq.end());
    StringVector cigar_r(cigar.begin(), cigar.end());
    StringVector rnext_r(rnext.begin(), rnext.end());
    NumericVector pnext_r(pnext.begin(), pnext.end());
    NumericVector tlen_r(tlen.begin(), tlen.end());
    StringVector seq_r(seq.begin(), seq.end());
    StringVector qual_r(qual.begin(), qual.end());

    DataFrame ret = DataFrame::create(Named("qname") = qname_r,
                                      Named("flag") = flag_r,
                                      Named("rname") = rname_r,
                                      Named("pos") = pos_r,
                                      Named("mapq") = mapq_r,
                                      Named("cigar") = cigar_r,
                                      Named("rnext") = rnext_r,
                                      Named("pnext") = pnext_r,
                                      Named("tlen") = tlen_r,
                                      Named("seq") = seq_r,
                                      Named("qual") = qual_r);

    return ret;
}


