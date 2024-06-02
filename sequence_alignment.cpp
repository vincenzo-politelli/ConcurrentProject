#include "sequence_alignment.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;


SequenceAlignment::SequenceAlignment(
    char *seq_1, char *seq_2, 
    unsigned int len_seq_1, unsigned int len_seq_2,
    unsigned int num_threads, unsigned int block_size_x, unsigned int block_size_y,
    int gap_cost, int match_score, int mismatch_cost, ProbType prob_type) {
    
    this->seq_1 = seq_1;
    this->seq_2 = seq_2;

    this->len_seq_1 = len_seq_1;
    this->len_seq_2 = len_seq_2;
    this->num_threads = num_threads;

    mutexes = new std::mutex[num_threads];
    cond_vars = new std::condition_variable[num_threads];

    this->block_size_x = block_size_x;
    this->block_size_y = block_size_y;
    this->gap_cost = gap_cost;

    this->match_score = match_score;
    this->mismatch_cost = mismatch_cost;
    this->prob_type = prob_type;

    H_mat = new std::vector<std::vector<int>>(len_seq_1+1, std::vector<int>(len_seq_2+1, 0));
    
    traceback_mat = new std::vector<std::vector<Direction>>(
        len_seq_1+1, 
        std::vector<Direction>(len_seq_2+1, Direction::INVALID)
    );

    if (prob_type == ProbType::GLOBAL_ALIGNMENT) {
        for (unsigned int i = 1; i < len_seq_1 + 1; ++i) {
            (*traceback_mat)[i][0] = Direction::INSERT;
        }
        for (unsigned int j = 1; j < len_seq_2 + 1; ++j) {
            (*traceback_mat)[0][j] = Direction::DELETE;
        }
    }

    std::atomic_init(&num_threads_finished, 0);
    std::atomic_init(&phase, 1);

    if (len_seq_1 % block_size_x == 0)
        num_bocks_1 = len_seq_1 / block_size_x;
    else
        num_bocks_1 = (len_seq_1 / block_size_x) + 1;

    if (len_seq_2 % block_size_y == 0)
        num_bocks_2 = len_seq_2 / block_size_y;
    else
        num_bocks_2 = (len_seq_2 / block_size_y) + 1;

    num_phases = num_bocks_1 + num_bocks_2 - 1; // to be computed based on m, n, block_size

    this->align1 = new char[len_seq_1 + len_seq_2 + 2];
    this->align2 = new char[len_seq_1 + len_seq_2 + 2];
}

SequenceAlignment::~SequenceAlignment() {
    delete H_mat;
    delete traceback_mat;
    delete mutexes;
    delete align1;
    delete align2;
}

void SequenceAlignment::compute_score_matrix() {
    for (unsigned int i = 0; i < num_threads; ++i) {
        threads.push_back(std::thread(&SequenceAlignment::processor_compute, this, i+1));
    }

    for (unsigned int i = 0; i < num_threads; ++i) {
        if (threads[i].joinable())
            threads[i].join();
    }
}

void SequenceAlignment::processor_compute(unsigned int processor_id) {
    int local_phase = 1;
    while (local_phase <= num_phases) {
        std::vector<Block> blocks_to_compute;
        first_cells(processor_id, blocks_to_compute, local_phase);

        for (auto b = blocks_to_compute.begin(); b != blocks_to_compute.end(); ++b) {
            Block block = *b;
            for (unsigned int i = block.start_x; i < std::min(block.start_x + block_size_x, len_seq_1 + 1); ++i) {
                for (unsigned int j = block.start_y; j < std::min(block.start_y + block_size_y, len_seq_2 + 1); ++j) {
                    score_cell(i, j);
                }
            }
        }

        local_phase++;
        unsigned int fetched = num_threads_finished.fetch_add(1);
        if (fetched == num_threads - 1) {
            int current_val = fetched + 1;
            while (!num_threads_finished.compare_exchange_weak(current_val, 0));
            phase.fetch_add(1);
            
            for (unsigned int i = 0; i < num_threads; i++) {
                cond_vars[i].notify_one();
            }
            continue;
        }
        std::unique_lock<std::mutex> lk(mutexes[processor_id - 1]);
        while (phase.load() != local_phase)
            cond_vars[processor_id - 1].wait(lk);
    }
}

void SequenceAlignment::first_cells(unsigned int processor_id, std::vector<Block> &blocks, unsigned int phase) {

    if (processor_id > num_threads){
        std::cout << "Invalid : processor_id is greater than the number of threads.\n";
        return;
    }

    int i = 1 + (processor_id - 1) * block_size_x;
    int j = 1 + (phase - 1 - processor_id + 1) * block_size_y;

    while (i < (int) len_seq_1 + 1 && j >= 1) {
        Block new_block;
        new_block.start_x = i;
        new_block.start_y = j;
        blocks.push_back(new_block);

        i += num_threads * block_size_x;
        j -= num_threads * block_size_y;
    }

    if (i < (int) len_seq_1 + 1 && j >= 1) {
        Block new_block;
        new_block.start_x = i;
        new_block.start_y = j;
        blocks.push_back(new_block);
    }
}

int SequenceAlignment::simple_score(char a_i, char b_j) {
    if (a_i == b_j) {
        return match_score;
    } else {
        return mismatch_cost;
    }
}

void SequenceAlignment::max_score(unsigned int &k, unsigned int &l){
    int max = 0;
    for (unsigned int i = 0; i < len_seq_1 + 1; i++) {
        for (unsigned int j = 0; j < len_seq_2 + 1; j++) {
            if ((*H_mat)[i][j] > max){
                max = (*H_mat)[i][j];
                k = i;
                l = j;
            }
        }
    }
}

void SequenceAlignment::score_cell(unsigned int i, unsigned int j) {

    int match_vs_mismatch = (*H_mat)[i-1][j-1] + simple_score(seq_1[i-1], seq_2[j-1]);
    int insert = (*H_mat)[i-1][j] - gap_cost;
    int deletee = (*H_mat)[i][j-1] - gap_cost;

    int score_ij;
    Direction dir;

    if (prob_type == ProbType::LOCAL_ALIGNMENT) {
        std::vector<int> max_v;
        max_v.push_back(0);
        max_v.push_back(match_vs_mismatch);
        max_v.push_back(deletee);
        max_v.push_back(insert);
        score_ij = max_vector(max_v);
    } else {
        std::vector<int> max_v;
        max_v.push_back(match_vs_mismatch);
        max_v.push_back(deletee);
        max_v.push_back(insert);
        score_ij = max_vector(max_v);
    }

    if (score_ij == 0 && prob_type == ProbType::LOCAL_ALIGNMENT) {
        dir = Direction::INVALID;
    } else if (score_ij == match_vs_mismatch) {
        dir = Direction::MATCH;
    } else if (score_ij == deletee) {
        dir = Direction::DELETE;
    } else if (score_ij == insert) {
        dir = Direction::INSERT;
    } else {
        dir = Direction::INVALID;
    }

    (*H_mat)[i][j] = score_ij;
    (*traceback_mat)[i][j] = dir;
}

void SequenceAlignment::traceback(unsigned int i, unsigned int j){
    while (true) {
        if (prob_type == ProbType::LOCAL_ALIGNMENT && (*H_mat)[i][j] == 0)
            break;

        if ((*traceback_mat)[i][j] == Direction::MATCH){
            align1[len_align1] = seq_1[i-1];
            align2[len_align2] = seq_2[j-1];
            i--;
            j--;
        }
        else if ((*traceback_mat)[i][j] == Direction::DELETE){
            align1[len_align1] = '-';
            align2[len_align2] = seq_2[j-1];
            j--;
        }
        else if ((*traceback_mat)[i][j] == Direction::INSERT) {
            align1[len_align1] = seq_1[i-1];
            align2[len_align2] = '-';
            i--;
        } else {
            break;
        }
        len_align1++;
        len_align2++;
    }

    //std::reverse(align1, align1+len_align1);
    //std::reverse(align2, align2+len_align2);
}

void SequenceAlignment::alignment() {
    compute_score_matrix();
    unsigned int align1_start, align2_start;
    if (prob_type == ProbType::LOCAL_ALIGNMENT) {
        max_score(align1_start, align2_start);
    } else {
        align1_start = len_seq_1;
        align2_start = len_seq_2;
    }
    traceback(align1_start, align2_start);
}



