#include "alignment.hpp"

void SWAlignment::simple_score(unsigned int i, unsigned int j) {

    int match_mismatch = (seq_1[i-1] == seq_2[j-1]) ? match : mismatch;
    
    int type_1 = (*H)[i-1][j-1] + match_mismatch;
    int type_2 = (*H)[i-1][j] - gap_penalty;
    int type_3 = (*H)[i][j-1] - gap_penalty;

    int score;
    int traceback_path;

    if (is_global == 1) {
        score = 0;
        if (type_1 > score) {
            score = type_1;
        }
        if (type_2 > score) {
            score = type_2;
        }
        if (type_3 > score) {
            score = type_3;
        }
    } 
    
    else {
        score = type_1;
        if (type_2 > score) {
            score = type_2;
        }
        if (type_3 > score) {
            score = type_3;
        }
    }

    if (score == 0 && is_global == 1) {
        traceback_path = -1;
    } 
    else if (score == type_1) {
        traceback_path = 1;
    } 
    else if (score == type_2) {
        traceback_path = 2;
    } 
    else if (score == type_3) {
        traceback_path = 3;
    } 

    (*H)[i][j] = score;
    (*T)[i][j] = traceback_path;
}

void SWAlignment::run_algo() {

    fill_DP_matrix();

    unsigned int begin_align_1;
    unsigned int begin_align_2;

    if (is_global == 1) {
        int best_score = 0;
        for (unsigned int i = 0; i < len_seq_1 + 1; i++) {
            for (unsigned int j = 0; j < len_seq_2 + 1; j++) {
                if ((*H)[i][j] > best_score){
                    best_score = (*H)[i][j];
                    begin_align_1 = i; begin_align_2 = j;
                }
            }
        }
    } 
    else {
        begin_align_1 = len_seq_1;
        begin_align_2 = len_seq_2;
    }

    traceback(begin_align_1, begin_align_2);
}

void SWAlignment::traceback(unsigned int i, unsigned int j){
    while (1) {

        if (is_global == 1 && (*H)[i][j] == 0)
            break;

        if ((*T)[i][j] == 1){
            align_1[len_align_1] = seq_1[i-1];
            align_2[len_align_2] = seq_2[j-1];
            i--;
            j--;
        }
        else if ((*T)[i][j] == 2) {
            align_1[len_align_1] = seq_1[i-1];
            align_2[len_align_2] = '-';
            i--;
        } 
        else if ((*T)[i][j] == 3){
            align_1[len_align_1] = '-';
            align_2[len_align_2] = seq_2[j-1];
            j--;
        }
        else {
            break;
        }

        len_align_1++;
        len_align_2++;
    }

    for (unsigned int i = 0, j = len_align_1 - 1; i < len_align_1/2; i++, j--) {     
        char tmp = align_1[i];  
        align_1[i] = align_1[j];  
        align_1[j] = tmp;  
    }

    for (unsigned int i = 0, j = len_align_2 - 1; i < len_align_2/2; i++, j--) {     
        char tmp = align_2[i];  
        align_2[i] = align_2[j];  
        align_2[j] = tmp;  
    }

}
