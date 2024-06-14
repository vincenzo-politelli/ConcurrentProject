#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <atomic>
#include <thread>
#include <vector>
#include <mutex>

struct Block { unsigned int start_x; 
               unsigned int start_y; 
};

class SWAlignment {
public:
    char *seq_1, *seq_2;
    unsigned int len_seq_1, len_seq_2;
    char *align_1, *align_2;
    unsigned int len_align_1, len_align_2;
    
    int match, mismatch, gap_penalty;
    int is_global;

    std::vector<std::vector<int> > *H;
    std::vector<std::vector<int> > *T;

    virtual void fill_DP_matrix() {};
    void simple_score(unsigned int, unsigned int);
    void traceback(unsigned int, unsigned int);
    void run_algo();
};


class SWAlignment_sequential : public SWAlignment {
public:

    SWAlignment_sequential(
        char *seq_1, char *seq_2, 
        unsigned int len_seq_1, unsigned int len_seq_2, 
        int match, int mismatch, int gap_penalty,
        int is_global = 1 );
    ~SWAlignment_sequential();
    void fill_DP_matrix();
};


class SWAlignment_parallel : public SWAlignment {
public:
    unsigned int num_threads, block_size_x, block_size_y;
    unsigned int num_vertical_blocks, num_horizontal_blocks;

    std::vector<std::thread> threads;
    std::atomic_int threads_done;
    std::atomic_int phase;
    int num_phases;

    std::mutex *mutex_array;
    std::condition_variable *condition_vars_array;

    SWAlignment_parallel(
        char *seq_1, char *seq_2, 
        unsigned int len_seq_1, unsigned int len_seq_2,
        unsigned int num_threads, 
        unsigned int block_size_x, unsigned int block_size_y,
        int match, int mismatch, int gap_penalty,
        int is_global = 1
    );

    ~SWAlignment_parallel();

    void fill_DP_matrix();
    void local_compute(unsigned int id);
    void next_phase_prep(unsigned int id, unsigned int curr_phase);
    void start_cell(unsigned int id, std::vector<Block> &blocks, unsigned int phase);
};

int ceil_ratio(int int_1, int int_2) {
    if (int_1 % int_2 == 0)
        return int_1 / int_2;
    else
        return int_1 / int_2 + 1;
}

#endif
