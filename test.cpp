#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include "sequence_alignment_sequential.h"
#include "sequence_alignment_parallel.h"
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
    string seq1 = "GACTTAC";;
    string seq2 = "CGTGAATTCAT";
    int match_score = 5;
    int mismatch_cost = -3;
    int gap_cost = -4;

    SequenceAlignment_Parallel aligner(seq1, seq2, match_score, mismatch_cost, gap_cost);
    aligner.alignSequences();

    cout << "Alignment 1: " << aligner.align1 << endl;
    cout << "Alignment 2: " << aligner.align2 << endl;


    if (aligner.align1=="---GACTT-AC"&& aligner.align2=="CGTGAATTCAT"){
        printf("TEST3 PASSED\n\n");
    }
    else{
        printf("TEST3 FAILED\n");
    }
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

string generateRandomSequence(int length) {
    const string bases = "ACGT";
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, bases.size() - 1);

    string sequence;
    for (int i = 0; i < length; ++i) {
        sequence += bases[dis(gen)];
    }
    return sequence;
}


void test5(){
    int numTests = 100;
    int sequenceLength = 10000;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis_block(1, sequenceLength / 2);

    for (int i = 0; i < numTests; ++i) {
        // Generate random sequences
        string seq1 = generateRandomSequence(sequenceLength);
        string seq2 = generateRandomSequence(sequenceLength-5);
        int match_score = 5;
        int mismatch_cost = -3;
        int gap_cost = -4;

        int block_size_x = dis_block(gen);
        int block_size_y = dis_block(gen);
        //printf("block_size_x=%i ,block_size_y=%i\n",block_size_x,block_size_y);



        SequenceAlignment_Sequential sequentialAligner(seq1, seq2, match_score, mismatch_cost, gap_cost);
        auto startSeq = chrono::high_resolution_clock::now();
        sequentialAligner.alignSequences();
        auto endSeq = chrono::high_resolution_clock::now();
        auto sequentialTime = chrono::duration_cast<chrono::milliseconds>(endSeq - startSeq).count();

        SequenceAlignment_Parallel parallelAligner(seq1, seq2, match_score, mismatch_cost, gap_cost);
        parallelAligner.block_size_x = block_size_x;
        parallelAligner.block_size_y = block_size_y;

        auto startPar = chrono::high_resolution_clock::now();
        parallelAligner.alignSequences();
        auto endPar = chrono::high_resolution_clock::now();
        auto parallelTime = chrono::duration_cast<chrono::milliseconds>(endPar - startPar).count();

        // Validate results
        assert(sequentialAligner.align1 == parallelAligner.align1);
        assert(sequentialAligner.align2 == parallelAligner.align2);

        // Print results
        cout << "Test " << i + 1 << " passed. "
             << "Sequential time: " << sequentialTime << " ms, "
             << "Parallel time: " << parallelTime << " ms." << endl;
    }   
}









