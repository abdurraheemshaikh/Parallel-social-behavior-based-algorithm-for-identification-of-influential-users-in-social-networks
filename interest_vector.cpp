#include <random>
#include <iostream>
#include <omp.h>
#include "interest_vector.h"

void generateInterestVectors(std::vector<Node>& graph) {
    std::cout << "[REPORT] Generating " << D << "-dimensional interest vectors\n";

    // Set number of threads explicitly
    omp_set_num_threads(4);

    // Log thread count
    #pragma omp parallel
    {
        #pragma omp single
        {
            std::cout << "[REPORT] OpenMP parallelization enabled with "
                      << omp_get_num_threads() << " threads\n";
        }
    }

    // Parallelize interest vector generation
    #pragma omp parallel for
    for (size_t u = 0; u < graph.size(); ++u) {
        int thread_id = omp_get_thread_num();
        if (u % 100000 == 0) {
            std::cout << "[REPORT] Thread " << thread_id
                      << " processing node " << u << "\n";
        }

        // Thread-local RNG with unique seed for reproducibility
        std::mt19937 thread_gen(42 + static_cast<unsigned int>(u));
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        auto& vec = graph[u].interest;
        vec.resize(D);

        float sum = 0.0f;
        for (int i = 0; i < D; ++i) {
            vec[i] = dist(thread_gen);
            sum += vec[i];
        }

        for (int i = 0; i < D; ++i) {
            vec[i] /= sum;
        }
    }

    std::cout << "[REPORT] Interest vectors assigned\n";
}

