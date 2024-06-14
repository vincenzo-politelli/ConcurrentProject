#include "alignment.hpp"

SWAlignment_sequential::SWAlignment_sequential(
    char *seq_1, char *seq_2, 
    unsigned int len_seq_1, unsigned int len_seq_2,
    int match, int mismatch, int gap_penalty,
    int is_global) {

    this->seq_1 = seq_1; this->len_seq_1 = len_seq_1;
    this->seq_2 = seq_2; this->len_seq_2 = len_seq_2;
    
    this->match = match;
    this->mismatch = mismatch;
    this->gap_penalty = gap_penalty;

    this->is_global = is_global;

    H = new std::vector<std::vector<int> >(len_seq_1+1, std::vector<int>(len_seq_2+1, 0));
    T = new std::vector<std::vector<int> >( len_seq_1+1, std::vector<int>(len_seq_2+1, -1));

    if (is_global == 0) {
        for (unsigned int i = 1; i < len_seq_1 + 1; ++i) {
            (*T)[i][0] = 2;
        }
        for (unsigned int j = 1; j < len_seq_2 + 1; ++j) {
            (*T)[0][j] = 3;
        }
    }

    this->align_1 = new char[len_seq_1 + len_seq_2 + 2];
    this->align_2 = new char[len_seq_1 + len_seq_2 + 2];

    this->len_align_1 = 0;
    this->len_align_2 = 0;
}

SWAlignment_sequential::~SWAlignment_sequential() {
        delete H;
        delete T;
        delete align_1;
        delete align_2;
}

void SWAlignment_sequential::fill_DP_matrix() {

    for (unsigned int i = 1; i < len_seq_1+1; i++){
        for (unsigned int j = 1; j < len_seq_2+1; j++){
            simple_score(i, j);
        }
    }
}
