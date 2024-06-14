#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <mutex>

using namespace std;
#pragma once


class SequenceAlignment_Parallel {
public:

    string seq1, seq2;
    int gap_cost, match_score, mismatch_cost;
    vector<vector<int>> dp;
    string align1, align2;
    mutex mtx;
    condition_variable cv;
    vector<vector<bool>> block_done;

    unsigned long nb_threads;
    int block_size_x = 1;  
    int block_size_y = 1; 

    SequenceAlignment_Parallel(const string& s1, const string& s2, int match, int mismatch, int gap)
        : seq1(s1), seq2(s2), gap_cost(gap), match_score(match), mismatch_cost(mismatch) {}

    void initializeDPTable();
    void fillDPTable();
    void traceback();
    void processBlock(int startRow, int endRow, int startCol, int endCol, int block_x, int block_y);
    bool isBlockReady(int block_x, int block_y);
    void fillDPBlock(int startRow, int endRow, int startCol, int endCol, int block_x, int block_y);

    void alignSequences() {
        initializeDPTable();
        fillDPTable();
        traceback();
    }
    
};
