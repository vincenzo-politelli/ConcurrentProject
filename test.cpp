#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "sequence_alignment_sequential.h"
#include "sequence_alignment.hpp"
#include "test.h"
using namespace std;

void test1(){
    string seq1 = "ATGTCGA";
    string seq2 = "AGAATCTA";
    int match_score = 5;
    int mismatch_cost = -3;
    int gap_cost = -4;

    SequenceAlignment_Sequential aligner(seq1, seq2, match_score, mismatch_cost, gap_cost);
    aligner.alignSequences();

    cout << "Alignment 1: " << aligner.align1 << endl;
    cout << "Alignment 2: " << aligner.align2 << endl;


    if (aligner.align1=="ATG--TCGA" && aligner.align2=="A-GAATCTA"  ){
        printf("TEST1 PASSED\n\n");
    }
    else{
        printf("TEST1 FAILED\n");

    }

}


void test2(){
    string seq1 = "GACTTAC";;
    string seq2 = "CGTGAATTCAT";
    int match_score = 5;
    int mismatch_cost = -3;
    int gap_cost = -4;

    SequenceAlignment_Sequential aligner(seq1, seq2, match_score, mismatch_cost, gap_cost);
    aligner.alignSequences();

    cout << "Alignment 1: " << aligner.align1 << endl;
    cout << "Alignment 2: " << aligner.align2 << endl;


    if (aligner.align1=="---GACTT-AC"&& aligner.align2=="CGTGAATTCAT"){
        printf("TEST2 PASSED\n\n");
    }
    else{
        printf("TEST2 FAILED\n");
    }
}

void test3(){
    char* sequence_A = "GACTTAC";
    char* sequence_B = "CGTGAATTCAT";
    int n = strlen(sequence_A);
    int m = strlen(sequence_B);
    int gap_penalty = 4;
    int match_score = 5;
    int mismatch_score = -3;
    unsigned int num_threads = 1, block_size_x = 1, block_size_y = 1;
    ProbType at = ProbType::LOCAL_ALIGNMENT;  

    SequenceAlignment aligner(sequence_A, sequence_B, n, m, 
        num_threads, block_size_x, block_size_y, 
        gap_penalty, match_score, mismatch_score, at
    );
    printf("pass\n");
    aligner.alignment();
    printf("Testing parallel\n");
    
    for(int i = 0; i < aligner.len_align1; i++) {
        putchar(aligner.align1[i]);
    } 
    printf("\n");
    // if (seq1=="---GACTT-AC" && seq2=="CGTGAATTCAT"){
    //     printf("TEST3 PASSED\n");
    // }
    // else{
    //     printf("TEST3 FAILED\n");
    //     printf("seq1=%s, seq2=%s",aligner.align1,aligner.align2);
    // }
}


void test4(){
    string seq1 = "TGTTACGG";
    string seq2 = "GGTTGACTA";
    int match_score = 3;
    int mismatch_cost = -3;
    int gap_cost = -2;

    SequenceAlignment_Sequential aligner(seq1, seq2, match_score, mismatch_cost, gap_cost);
    aligner.alignSequences_local();

    cout << "Alignment 1: " << aligner.align1 << endl;
    cout << "Alignment 2: " << aligner.align2 << endl;

    if (aligner.align1=="GTT-AC"&& aligner.align2=="GTTGAC"){
        printf("TEST4 PASSED\n\n");
    }
    else{
        printf("TEST4 FAILED\n");
    }

}









