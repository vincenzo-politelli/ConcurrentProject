#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;
#pragma once

class SequenceAlignment_Sequential {
public:
    string seq1, seq2;
    int match_score, mismatch_cost, gap_cost;
    vector<vector<int>> dp;
    string align1, align2;

    SequenceAlignment_Sequential(const string& s1, const string& s2, int match, int mismatch, int gap)
        : seq1(s1), seq2(s2), match_score(match), mismatch_cost(mismatch), gap_cost(gap) {}

    // Function to perform sequence alignment
    void alignSequences() {
        initializeDPTable();
        fillDPTable();
        traceback();
    }

    void initializeDPTable();
    
    void fillDPTable();

    void traceback();
};