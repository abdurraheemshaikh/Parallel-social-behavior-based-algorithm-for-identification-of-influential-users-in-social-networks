#include "graph_loader.cpp"
#include "interest_vector.cpp"
#include "partition.cpp"
#include "influence.cpp"
#include "seed_selection.cpp"
#include "node.h"
#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <algorithm>
#include <mpi.h>

void logBarrier(int rank) {
    MPI_Barrier(MPI_COMM_WORLD);
}

void logRoot(int rank, const std::string& msg) {
    if (rank == 0) std::cout << "[REPORT] " << msg << std::endl;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    double start_time = MPI_Wtime();

    logRoot(mpi_rank, "MPI initialized with " + std::to_string(mpi_size) + " processes");

    constexpr uint32_t MAX_TEST_NODES = 10000;
    constexpr size_t N = MAX_TEST_NODES;
    logRoot(mpi_rank, "[LOADING] Initializing graph with " + std::to_string(N) + " nodes");

    std::vector<Node> graph(N);

    // Step 1: Interest Vector Generation
    if (mpi_rank == 0) {
        logRoot(mpi_rank, "[LOADING] Generating interest vectors...");
        generateInterestVectors(graph);
        logRoot(mpi_rank, "Interest vector generation complete.");
    }
    logBarrier(mpi_rank);

    // Step 2: Load Graph
    if (mpi_rank == 0) {
        logRoot(mpi_rank, "[LOADING] Loading graph layers...");
        loadGraph(graph);
        logRoot(mpi_rank, "Graph load complete. Ready for SCC/CAC partitioning.");
    }
    logBarrier(mpi_rank);

    // Step 3: Partition Graph
    if (mpi_rank == 0) {
        logRoot(mpi_rank, "[LOADING] Computing SCC/CAC partitioning...");
        computeSCC_CAC(graph);
        logRoot(mpi_rank, "Partitioning complete.");

        int maxCID = 0;
        for (const auto& node : graph)
            if (node.compID > maxCID) maxCID = node.compID;

        logRoot(mpi_rank, "Found " + std::to_string(maxCID + 1) + " components.");
    }
    logBarrier(mpi_rank);

    // Step 4: Influence Power Calculation
    logRoot(mpi_rank, "[LOADING] Computing influence power...");
    std::array<double, NUM_LAYERS> alpha = {0.0, 0.50, 0.15, 0.35}; // RETWEET, MENTION, REPLY
    double damping = 0.85;
    std::vector<double> IP = computeInfluencePower(graph, alpha, damping);
    logRoot(mpi_rank, "Influence power computation complete.");
    logBarrier(mpi_rank);

    // Step 5: Seed Candidate Selection
    std::vector<int> I_star;
    if (mpi_rank == 0) {
        logRoot(mpi_rank, "[LOADING] Selecting seed candidates (Algorithm 6)...");
        I_star = selectSeedCandidates(graph, IP);
        logRoot(mpi_rank, "Seed candidate selection complete.");
    }
    logBarrier(mpi_rank);

    // Step 6: Final Seed Selection
    std::vector<int> INF;
    if (mpi_rank == 0) {
        logRoot(mpi_rank, "[LOADING] Final influential seed selection (Algorithm 7)...");
        INF = selectFinalSeeds(graph, IP, I_star);

        logRoot(mpi_rank, "Top 10 final seed nodes:");
        for (int i = 0; i < std::min(10, (int)INF.size()); ++i) {
            std::cout << "  Node " << INF[i]
                      << " → IP = " << std::setprecision(6) << IP[INF[i]] << "\n";
        }
    }
    logBarrier(mpi_rank);

    // End timer
    double end_time = MPI_Wtime();
    if (mpi_rank == 0) {
        logRoot(mpi_rank, "All stages complete.");
        std::cout << "[TIME] Total execution time: " 
                  << std::fixed << std::setprecision(2)
                  << (end_time - start_time) << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}

