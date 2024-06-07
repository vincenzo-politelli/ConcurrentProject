#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <thread>
#include <mutex>
#include "sequence_alignment_parallel.h"
using namespace std;

void SequenceAlignment_Parallel::initializeDPTable() {
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

void SequenceAlignment_Parallel::fillDPBlock(int startRow, int endRow, int startCol, int endCol, int block_x, int block_y) {
    for (int i = startRow; i <= endRow; ++i) {
        for (int j = startCol; j <= endCol; ++j) {
            int match = dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_cost);
            int delete_op = dp[i - 1][j] + gap_cost;
            int insert_op = dp[i][j - 1] + gap_cost;
            dp[i][j] = max({match, delete_op, insert_op});
            cout << "i: " << i << ", j: " << j << endl;
            cout << "match: " << match << ", delete_op: " << delete_op << ", insert_op: " << insert_op << endl;
        }
    }

    //unique_lock<mutex> lock(mtx);
    
    block_done[block_x][block_y] = true;

    cv.notify_all();
}

bool SequenceAlignment_Parallel::isBlockReady(int block_x, int block_y) {
    if (block_x == 0 && block_y == 0) return true;
    if (block_x == 0) return block_done[block_x][block_y - 1];
    if (block_y == 0) return block_done[block_x - 1][block_y];
    return block_done[block_x - 1][block_y] && block_done[block_x][block_y - 1];
}

void SequenceAlignment_Parallel::fillDPTable() {
    int m = seq1.length();
    int n = seq2.length();
    int num_blocks_x = (m + this->block_size_x - 1) / this->block_size_x;
    int num_blocks_y = (n + this->block_size_y - 1) / this->block_size_y;
    block_done = vector<vector<bool>>(num_blocks_x, vector<bool>(num_blocks_y, false));
    vector<thread> threads;

    for (int block_x = 0; block_x < num_blocks_x; ++block_x) {
        for (int block_y = 0; block_y < num_blocks_y; ++block_y) {
            int startRow = block_x * block_size_x + 1;
            int endRow = min((block_x + 1) * block_size_x, m);
            int startCol = block_y * block_size_y + 1;
            int endCol = min((block_y + 1) * block_size_y, n);

            threads.emplace_back([=]() {
                unique_lock<mutex> lock(mtx);
                cv.wait(lock, [=] { return isBlockReady(block_x, block_y); });
                lock.unlock();

                fillDPBlock(startRow, endRow, startCol, endCol, block_x, block_y);
            });
        }
    }

    for (auto& th : threads) {
        th.join();
    }
}


void SequenceAlignment_Parallel::traceback() {
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















