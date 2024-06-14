#include "sequence_alignment_sequential.h"
#include "alignment.hpp"
#include "alignment.cpp"
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "test.h"


int main() {
    std::string* seq_1 = new std::string("GACTTAC");
    std::string* seq_2 = new std::string("CGTGAATTCAT");

    int match = 5;
    int mismatch = -3;
    int gap_penalty = 4;
    unsigned int num_threads = 2;
    unsigned int num_rows_block = 2;
    unsigned int num_cols_block = 2;
    ProblemType prob_type = ProblemType::LOCAL;

    Alignment al(seq_1, seq_2, match, mismatch, gap_penalty, prob_type, num_threads, num_rows_block, num_cols_block);
    al.run(Threading::MULTI);
    
    delete seq_1;
    delete seq_2;
}

// int main() {
//     // printf("Testing sequential algorithm\n");
//     // test1();
//     // test2();
//     // test3();
//     // test4();
//     // test5();
//     // write();
//     // return 0;

// }
