#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "sequence_alignment_sequential.h"
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
        printf("TEST1 PASSED\n");
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


    if (aligner.align1=="---GACTT-AC"&& aligner.align2=="CGTGAATTCAT"  ){
        printf("TEST2 PASSED\n");
    }
    else{
        printf("TEST2 FAILED\n");

    }

}

