#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include "sequence_alignment_sequential.h"
using namespace std;

void SequenceAlignment_Sequential::initializeDPTable() {
        int m = seq1.length();
        int n = seq2.length();
        dp = vector<vector<int>>(m + 1, vector<int>(n + 1, 0));

        for (int i = 0; i <= m; ++i) {
            dp[i][0] = i * gap_cost;
        }
        for (int j = 0; j <= n; ++j) {
            dp[0][j] = j * gap_cost;
        }
    }

void SequenceAlignment_Sequential::fillDPTable(){    
        int m = seq1.length();
        int n = seq2.length();
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                int match = dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_cost);
                int delete_op = dp[i - 1][j] + gap_cost;
                int insert_op = dp[i][j - 1] + gap_cost;
                dp[i][j] = max({match, delete_op, insert_op});
            }
        }
    // cout << "Completed DP Table:" << endl;
    // for (const auto& row : dp) {
    //     for (int val : row) {
    //         cout << val << " ";
    //     }
    //     cout << endl;
    // }
}

void SequenceAlignment_Sequential::traceback(){
        int i = seq1.length();
        int j = seq2.length();
        align1 = "";
        align2 = "";
        
        while (i > 0 && j > 0) {
            int score_current = dp[i][j];
            int score_diagonal = dp[i - 1][j - 1];
            int score_up = dp[i][j - 1];
            int score_left = dp[i - 1][j];

            if (score_current == score_diagonal + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_cost)) {
                align1 += seq1[i - 1];
                align2 += seq2[j - 1];
                --i;
                --j;
            } else if (score_current == score_left + gap_cost) {
                align1 += seq1[i - 1];
                align2 += '-';
                --i;
            } else if (score_current == score_up + gap_cost) {
                align1 += '-';
                align2 += seq2[j - 1];
                --j;
            }
        }

        while (i > 0) {
            align1 += seq1[i - 1];
            align2 += '-';
            --i;
        }
        while (j > 0) {
            align1 += '-';
            align2 += seq2[j - 1];
            --j;
        }

        std::reverse(this->align1.begin(), this->align1.end());
        std::reverse(this->align2.begin(), this->align2.end());

}


void SequenceAlignment_Sequential::initializeDPTable_local() {
        int m = seq1.length();
        int n = seq2.length();
        dp = vector<vector<int>>(m + 1, vector<int>(n + 1, 0));
    }

void SequenceAlignment_Sequential::fillDPTable_local(){    
    int m = seq1.length();
    int n = seq2.length();
    max_score = 0;
    max_i = 0;
    max_j = 0;
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match = dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_cost);
            int delete_op = dp[i - 1][j] + gap_cost;
            int insert_op = dp[i][j - 1] + gap_cost;
            dp[i][j] = max({0, match, delete_op, insert_op});
            
            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
}

void SequenceAlignment_Sequential::traceback_local(){
    int i = max_i;
    int j = max_j;
    align1 = "";
    align2 = "";

    while (i > 0 && j > 0 && dp[i][j] > 0) {
        int score_current = dp[i][j];
        int score_diagonal = dp[i - 1][j - 1];
        int score_up = dp[i][j - 1];
        int score_left = dp[i - 1][j];

        if (score_current == score_diagonal + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_cost)) {
            align1 += seq1[i - 1];
            align2 += seq2[j - 1];
            --i;
            --j;
        } else if (score_current == score_left + gap_cost) {
            align1 += seq1[i - 1];
            align2 += '-';
            --i;
        } else if (score_current == score_up + gap_cost) {
            align1 += '-';
            align2 += seq2[j - 1];
            --j;
        }
    }

    std::reverse(align1.begin(), align1.end());
    std::reverse(align2.begin(), align2.end());
}

