#include <string>
#include <vector>

#include <Rcpp.h>

extern "C" {
#include "bam_lmdb.h"
}

using namespace Rcpp;

class StringAux {
 public:
  StringAux(const char* key1, const int num_elm) {
    data = StringVector(num_elm);
    key[0] = key1[0];
    key[1] = key1[1];
    key[2] = 0;
  }

  void add_aux(const char* val, const size_t size, const int index) {
    data[index] = std::string(val, size);
  }

  void add_na(int index) { data[index] = NA_STRING; }

  StringVector get_vector() { return data; }

  char* get_key() { return key; }

  bool compare_key(const char* key1) const {
    return (key1[0] == key[0] && key[1] == key1[1]);
  }

 private:
  char key[3];
  StringVector data;
};

class NumericAux {
 public:
  NumericAux(const char* key1, const char type1, const int num_elm) {
    data = NumericVector(num_elm);
    key[0] = key1[0];
    key[1] = key1[1];
    key[2] = 0;
    type = type1;
  }

  void add_aux(const double val, const int index) { data[index] = val; }

  NumericVector get_vector() { return data; }

  char* get_key() { return key; }

  char get_type() { return type; }

  bool compare_key(const char* key1) const {
    return (key1[0] == key[0] && key[1] == key1[1]);
  }

 private:
  char key[3];
  char type;
  NumericVector data;
};

class IntAux {
 public:
  IntAux(const char* key1, const char type1, const int num_elm) {
    data = IntegerVector(num_elm);
    key[0] = key1[0];
    key[1] = key1[1];
    key[2] = 0;
    type = type1;
  }

  void add_aux(const int val, const int index) { data[index] = val; }

  IntegerVector get_vector() { return data; }

  char* get_key() { return key; }

  char get_type() { return type; }

  bool compare_key(const char* key1) const {
    return (key1[0] == key[0] && key[1] == key1[1]);
  }

 private:
  char key[3];
  char type;
  IntegerVector data;
};

// [[Rcpp::export]]
DataFrame query_bam_index(CharacterVector bam_file_name,
                          CharacterVector bam_index_path,
                          CharacterVector key_type, CharacterVector key_value) {
  bam_row_set_t* bam_rows = NULL;

  /* Series of casts get from R strings to C strings */
  std::string file_name_str = Rcpp::as<std::string>(bam_file_name);
  std::string index_path_str = Rcpp::as<std::string>(bam_index_path);
  std::string key_type_str = Rcpp::as<std::string>(key_type);
  std::string key_val_str = Rcpp::as<std::string>(key_value);
  const char* c_file_name = file_name_str.c_str();
  const char* c_index_path = index_path_str.c_str();
  const char* c_key_type = key_type_str.c_str();
  const char* c_key_val = key_val_str.c_str();
  bam_aux_header_t* avail_aux_tag = NULL;
  std::list<StringAux> stringAuxes;
  std::list<NumericAux> numericAuxes;
  std::list<IntAux> intAuxes;
  int n_cols = 11;
  int rc = 0;

  rc = get_bam_rows(&bam_rows, c_file_name, c_index_path, c_key_type, c_key_val);
  avail_aux_tag = bam_rows->aux_tags.head;
  while (avail_aux_tag) {
    n_cols++;
    switch (avail_aux_tag->key.type) {
      case 'C':
      case 'c':
      case 'S':
      case 's':
      case 'I':
      case 'i':
        intAuxes.push_back(IntAux(avail_aux_tag->key.key,
                                  avail_aux_tag->key.type,
                                  bam_rows->num_entries));
        break;
      case 'A':
      case 'Z':
      case 'H':
        stringAuxes.push_back(
            StringAux(avail_aux_tag->key.key, bam_rows->num_entries));
        break;
      case 'd':
      case 'f':
        numericAuxes.push_back(NumericAux(avail_aux_tag->key.key,
                                          avail_aux_tag->key.type,
                                          bam_rows->num_entries));
        break;
    }

    avail_aux_tag = avail_aux_tag->next;
  }

  StringVector qname(bam_rows->num_entries);
  NumericVector flag(bam_rows->num_entries);
  StringVector rname(bam_rows->num_entries);
  NumericVector pos(bam_rows->num_entries);
  NumericVector mapq(bam_rows->num_entries);
  StringVector cigar(bam_rows->num_entries);
  StringVector rnext(bam_rows->num_entries);
  NumericVector pnext(bam_rows->num_entries);
  NumericVector tlen(bam_rows->num_entries);
  StringVector seq(bam_rows->num_entries);
  StringVector qual(bam_rows->num_entries);

  aux_elm_t* aux_tag = NULL;

  for (size_t i = 0; i < bam_rows->num_entries; i++) {
    qname[i] = std::string(bam_rows->rows[i]->qname);
    flag[i] = bam_rows->rows[i]->flag;
    rname[i] = std::string(bam_rows->rows[i]->rname);
    pos[i] = bam_rows->rows[i]->pos;
    mapq[i] = bam_rows->rows[i]->mapq;
    cigar[i] = std::string(bam_rows->rows[i]->cigar);
    rnext[i] = std::string(bam_rows->rows[i]->rnext);
    pnext[i] = bam_rows->rows[i]->pnext;
    tlen[i] = bam_rows->rows[i]->tlen;
    seq[i] = std::string(bam_rows->rows[i]->seq);
    qual[i] = std::string(bam_rows->rows[i]->qual);
    for (auto& it : intAuxes) {
      bool is_present = false;
      aux_tag = bam_rows->rows[i]->aux_list.head;

      while (aux_tag) {
        if (it.compare_key(aux_tag->key.key)) {
          is_present = true;
          switch (it.get_type()) {
            case 'C':
              it.add_aux(*static_cast<const uint8_t*>(aux_tag->val), i);
              break;
            case 'c':
              it.add_aux(*static_cast<const int8_t*>(aux_tag->val), i);
              break;
            case 'S':
              it.add_aux(*static_cast<const uint16_t*>(aux_tag->val), i);
              break;
            case 's':
              it.add_aux(*static_cast<const int16_t*>(aux_tag->val), i);
              break;
            case 'I':
              it.add_aux(*static_cast<const uint32_t*>(aux_tag->val), i);
              break;
            case 'i':
              it.add_aux(*static_cast<const int32_t*>(aux_tag->val), i);
              break;
            default:
              it.add_aux(*static_cast<const int*>(aux_tag->val), i);
          }
        }

        aux_tag = aux_tag->next;
      }

      if (!is_present) {
        it.add_aux(NA_INTEGER, i);
      }
    }

    for (auto& it : stringAuxes) {
      bool is_present = false;
      aux_tag = bam_rows->rows[i]->aux_list.head;

      while (aux_tag) {
        if (it.compare_key(aux_tag->key.key)) {
          is_present = true;
          it.add_aux(static_cast<const char*>(aux_tag->val), aux_tag->val_size,
                     i);
        }

        aux_tag = aux_tag->next;
      }

      if (!is_present) {
        it.add_na(i);
      }
    }

    for (auto& it : numericAuxes) {
      bool is_present = false;
      aux_tag = bam_rows->rows[i]->aux_list.head;

      while (aux_tag) {
        if (it.compare_key(aux_tag->key.key)) {
          is_present = true;
          if (it.get_type() == 'd') {
            it.add_aux(*static_cast<const double*>(aux_tag->val), i);
          } else {
            it.add_aux(*static_cast<const float*>(aux_tag->val), i);
          }
        }

        aux_tag = aux_tag->next;
      }

      if (!is_present) {
        it.add_aux(NA_REAL, i);
      }
    }
  }

  Rcpp::List colList(n_cols);
  std::vector<std::string> col_names;
  col_names.push_back("qname");
  colList[0] = qname;
  col_names.push_back("flag");
  colList[1] = flag;
  col_names.push_back("rname");
  colList[2] = rname;
  col_names.push_back("pos");
  colList[3] = pos;
  col_names.push_back("mapq");
  colList[4] = mapq;
  col_names.push_back("cigar");
  colList[5] = cigar;
  col_names.push_back("rnext");
  colList[6] = rnext;
  col_names.push_back("pnext");
  colList[7] = pnext;
  col_names.push_back("tlen");
  colList[8] = tlen;
  col_names.push_back("seq");
  colList[9] = seq;
  col_names.push_back("qual");
  colList[10] = qual;

  int j = 11;
  for (auto& it : intAuxes) {
    colList[j] = it.get_vector();
    col_names.push_back(it.get_key());
    j++;
  }

  for (auto& it : numericAuxes) {
    colList[j] = it.get_vector();
    col_names.push_back(it.get_key());
    j++;
  }

  for (auto& it : stringAuxes) {
    colList[j] = it.get_vector();
    col_names.push_back(it.get_key());
    j++;
  }

  colList.attr("names") = wrap(col_names);
  colList.attr("class") = "data.frame";
  colList.attr("row.names") =
      IntegerVector::create(NA_INTEGER, bam_rows->num_entries);

  free(bam_rows);
  return colList;
}
