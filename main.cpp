//#include "sequence_alignment_sequential.cpp"
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

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

private:
    void initializeDPTable() {
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
    
    void fillDPTable() {
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
    }

    void traceback() {
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
        reverse(align1.begin(), align1.end());
        reverse(align2.begin(), align2.end());
    }
};

// Main function to demonstrate the class usage
int main() {
    string seq1 = "GACTTAC";
    string seq2 = "CGTGAATTCAT";
    int match_score = 5;
    int mismatch_cost = -3;
    int gap_cost = -4;

    SequenceAlignment_Sequential aligner(seq1, seq2, match_score, mismatch_cost, gap_cost);
    aligner.alignSequences();

    cout << "Alignment 1: " << aligner.align1 << endl;
    cout << "Alignment 2: " << aligner.align2 << endl;

    return 0;
}
