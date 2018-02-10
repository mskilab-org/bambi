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

    /* Casting to get from R strings to C strings */
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
    
    std::vector<int> pos;
    std::vector<std::string> qname;
    std::vector<std::string> seq;

    for (size_t i = 0; i < bam_rows->n_entries; i++) {
        pos.push_back(bam_rows->rows[i]->pos);
        qname.push_back(std::string(bam_rows->rows[i]->qname));
        seq.push_back(std::string(bam_rows->rows[i]->seq));
    }

    NumericVector pos_r( pos.begin(), pos.end() );
    StringVector qname_r( qname.begin(), qname.end() );
    StringVector seq_r( seq.begin(), seq.end() );

    DataFrame ret = DataFrame::create(Named("pos") = pos_r,
                                      Named("qname") = qname_r,
                                      Named("seq") = seq_r);

    return ret;
}


