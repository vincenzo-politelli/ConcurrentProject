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
    int match_score, mismatch_cost, gap_cost;
    vector<vector<int>> dp;
    string align1, align2;
    mutex lock;


    int max_score; //local alignement
    int max_i; //local alignement
    int max_j; //local alignement



};
