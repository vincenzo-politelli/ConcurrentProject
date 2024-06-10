#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
#pragma once


enum Direction {
    MATCH,
    INSERT,
    DELETE,
    INVALID
};

enum ProbType {
    LOCAL_ALIGNMENT,
    GLOBAL_ALIGNMENT
};

struct Block {
    unsigned int start_x, start_y;
};

inline int max_vector(std::vector<int> ints) {
    if (ints.size() == 0) {
        return std::numeric_limits<int>::min();
    }

    int max = ints[0];
    for (auto it = ints.begin(); it != ints.end(); ++it) {
        if (*it > max) {
            max = *it;
        }
    }

    return max;
}


class SequenceAlignment {
public:
    char *seq_1;
    char *seq_2;
    unsigned int len_seq_1, len_seq_2;
    int gap_cost, match_score, mismatch_cost;

    char *align1;
    char *align2;
    unsigned int len_align1, len_align2;
    ProbType prob_type;

    std::vector<std::vector<int> > *H_mat;
    std::vector<std::vector<Direction> > *traceback_mat;

    unsigned int num_threads, block_size_x, block_size_y;
    unsigned int num_bocks_1, num_bocks_2;

    std::vector<std::thread> threads;
    std::mutex *mutexes;
    std::condition_variable *cond_vars;

    std::atomic_int num_threads_finished;
    std::atomic_int phase;
    int num_phases;

    SequenceAlignment(
        char *seq_1, char *seq_2, 
        unsigned int len_seq_1, unsigned int len_seq_2,
        unsigned int num_threads, unsigned int block_size_x, unsigned int block_size_y,
        int gap_cost, int match_score, int mismatch_cost, ProbType prob_type = ProbType::LOCAL_ALIGNMENT
    );

    ~SequenceAlignment();

    int simple_score(char a_i, char b_j);

    void max_score(unsigned int &k, unsigned int &l);

    void compute_score_matrix();

    void processor_compute(unsigned int processor_id);

    void first_cells(unsigned int processor_id, std::vector<Block> &blocks, unsigned int phase);

    void score_cell(unsigned int i, unsigned int j);

    void traceback(unsigned int i, unsigned int j);

    void alignment();
};

