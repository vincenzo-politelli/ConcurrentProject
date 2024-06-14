#include "alignment.hpp"

SWAlignment_parallel::SWAlignment_parallel(
    char *seq_1, char *seq_2, 
    unsigned int len_seq_1, unsigned int len_seq_2,
    unsigned int num_threads,
    unsigned int block_size_x, unsigned int block_size_y,
    int match, int mismatch, int gap_penalty,
    int is_global) {
    
    this->seq_1 = seq_1; this->len_seq_1 = len_seq_1;
    this->seq_2 = seq_2; this->len_seq_2 = len_seq_2;

    this->block_size_x = block_size_x;
    this->block_size_y = block_size_y;

    this->match = match;
    this->mismatch = mismatch;
    this->gap_penalty = gap_penalty;

    this->num_threads = num_threads;
    this->is_global = is_global;

    this->align_1 = new char[1 + len_seq_1 + len_seq_2 + 1];
    this->align_2 = new char[1 + len_seq_1 + len_seq_2 + 1];
    std::cout<<"hiii"<<std::endl;

    mutex_array = new std::mutex[num_threads];
    condition_vars_array = new std::condition_variable[num_threads];

    H = new std::vector<std::vector<int>>(len_seq_1 + 1, std::vector<int>(len_seq_2 + 1, 0));
    T = new std::vector<std::vector<int>>(len_seq_1 + 1, std::vector<int>(len_seq_2 + 1, -1));

    if (is_global == 0) {
        for (unsigned int i = 1; i < len_seq_1 + 1; ++i) {
            (*T)[i][0] = 2;
        }
        for (unsigned int j = 1; j < len_seq_2 + 1; ++j) {
            (*T)[0][j] = 3;
        }
    }
    std::cout<<"hello"<<std::endl;
    
    num_vertical_blocks = ceil_ratio(len_seq_1, block_size_x);
    num_horizontal_blocks = ceil_ratio(len_seq_2, block_size_y);

    num_phases = num_vertical_blocks + num_horizontal_blocks - 1;
    std::cout<<"hello2" <<std::endl;

    std::atomic_init(&threads_done, 0);
    std::atomic_init(&phase, 1);
}

SWAlignment_parallel::~SWAlignment_parallel() {
    
    delete align_1;
    delete align_2;
    delete H;
    delete T;
    delete mutex_array;
    delete condition_vars_array;

}

void SWAlignment_parallel::fill_DP_matrix() {
    for (unsigned int i = 0; i < num_threads; ++i) {
        threads.push_back(std::thread(&SWAlignment_parallel::local_compute, this, i+1));
    }

    for (unsigned int i = 0; i < num_threads; ++i) {
        threads[i].join();
    }
}

void SWAlignment_parallel::local_compute(unsigned int id) {
    unsigned int curr_phase = 1;

    while (num_phases >= curr_phase) {
        
        std::vector<Block> to_compute;
        start_cell(id, to_compute, curr_phase);

        for (auto blockk = to_compute.begin(); blockk != to_compute.end(); ++blockk) {
            Block block = *blockk;

            unsigned int min_1 = std::min(block.start_x + block_size_x, len_seq_1 + 1);
            unsigned int min_2 = std::min(block.start_y + block_size_y, len_seq_2 + 1);

            for (unsigned int i = block.start_x; i < min_1; ++i) {
                for (unsigned int j = block.start_y; j < min_2; ++j) {
                    simple_score(i, j);
                }
            }
        }

        curr_phase++;
        next_phase_prep(id, curr_phase);
    }
}

void SWAlignment_parallel::start_cell(unsigned int id, std::vector<Block> &blocks, unsigned int phase) {

    int i = block_size_x * (id - 1) + 1;
    int j = block_size_y * (phase - id) + 1;

    while (i < len_seq_1 + 1 && j >= 1) {
        Block new_block;

        new_block.start_x = i;
        new_block.start_y = j;

        blocks.push_back(new_block);

        i += num_threads * block_size_x;
        j -= num_threads * block_size_y;
    }

    if (i < len_seq_1 + 1 && j >= 1) {
        Block new_block;

        new_block.start_x = i;
        new_block.start_y = j;

        blocks.push_back(new_block);
    }
}

void SWAlignment_parallel::next_phase_prep(unsigned int id, unsigned int curr_phase) {
    unsigned int fetch_thread = threads_done.fetch_add(1);

    if (fetch_thread == num_threads - 1) {

        int current_val = fetch_thread + 1;
        while (!threads_done.compare_exchange_weak(current_val, 0));
        phase.fetch_add(1);
        
        for (unsigned int i = 0; i < num_threads; i++) {
            condition_vars_array[i].notify_one();
        }

        return;
    }

    std::unique_lock<std::mutex> lockk(mutex_array[id - 1]);
    while (phase.load() != curr_phase) {
        condition_vars_array[id - 1].wait(lockk);
    }
}

