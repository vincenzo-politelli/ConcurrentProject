//#include "sequence_alignment_sequential.cpp"
#include <string>
#include <iostream>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

// Function to calculate the T and prev values for a single row
void calculateTandPrev(const string& A, int rowStart, vector<int>& T, vector<int>& prev) {
    int n = A.size();
    T[0] = 0;
    prev[0] = 0;

    for (int j = 1; j < n; j++) {
        int score = max(T[j - 1] - 1, prev[j - 1] - 1);
        T[j] = max(score, static_cast<int>(T[j - 1] + (A[rowStart + j] == A[j])));
        prev[j] = max(T[j - 1], prev[j - 1]);
    }
}

// Function to calculate the T^R and next values for a single row
void calculateTRandNext(const string& A, int rowStart, vector<int>& TR, vector<int>& next) {
    int n = A.size();
    TR[n - 1] = 0;
    next[n - 1] = 0;

    for (int j = n - 2; j >= 0; j--) {
        int score = max(TR[j + 1] - 1, next[j + 1] - 1);
        TR[j] = max(score, static_cast<int>(TR[j + 1] + (A[rowStart + j] == A[j])));
        next[j] = max(TR[j + 1], next[j + 1]);
    }
}

// Function to compute the opt values for a single row
void computeOptValues(const vector<int>& T, const vector<int>& TR, vector<int>& opt) {
    int n = T.size();
    for (int j = 0; j < n; j++) {
        opt[j] = T[j] + TR[j];
    }
}

// DecomposeIteration function
void DecomposeIteration(const string& A, int rowUp, int rowDown, vector<int>& opt, int& rMid, int& pMid) {
    int rowMid = (rowUp + rowDown) / 2;

    vector<int> T(A.size(), 0);
    vector<int> prev(A.size(), 0);
    vector<int> TR(A.size(), 0);
    vector<int> next(A.size(), 0);

    // Initialize T values for rowUp
    calculateTandPrev(A, rowUp, T, prev);

    // Calculate T and prev values for rows rowUp + 1 through rowMid
    for (int row = rowUp + 1; row <= rowMid; row++) {
        calculateTandPrev(A, row, T, prev);
    }

    // Initialize T values for rowDown
    calculateTRandNext(A, rowDown, TR, next);

    // Calculate T^R and next values for rows rowDown - 1 through rowMid
    for (int row = rowDown - 1; row >= rowMid; row--) {
        calculateTRandNext(A, row, TR, next);
    }

    // Compute opt values for rowMid
    computeOptValues(T, TR, opt);

    // Find the leftmost maximal value of opt
    int maxVal = *max_element(opt.begin(), opt.end());
    auto it = find(opt.begin(), opt.end(), maxVal);
    rMid = it - opt.begin();
    pMid = 0; // Single processor case

    // Recursively call DecomposeIteration on the row boundaries
    if (rowMid > rowUp) {
        DecomposeIteration(A, rowUp, rowMid, opt, rMid, pMid);
    }
    if (rowMid < rowDown) {
        DecomposeIteration(A, rowMid + 1, rowDown, opt, rMid, pMid);
    }
}

int main() {
    string A = "ACGTACGT";
    int rowUp = 0;
    int rowDown = A.size() - 1;
    vector<int> opt(A.size(), 0);
    int rMid, pMid;

    DecomposeIteration(A, rowUp, rowDown, opt, rMid, pMid);

    cout << "Optimal alignment scores: ";
    for (int val : opt) {
        cout << val << " ";
    }
    cout << endl;

    return 0;
}