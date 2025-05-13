#include <omp.h>
#include <mpi.h>
#include <unordered_set>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include "influence.h"
static float jaccard(const std::vector<float>& a, const std::vector<float>& b) {
    float dot = 0, sum_a = 0, sum_b = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        dot += std::min(a[i], b[i]);
        sum_a += std::max(a[i], b[i]);
    }
    return sum_a > 0 ? dot / sum_a : 0.0f;
}

std::vector<double>
computeInfluencePower(const std::vector<Node>& graph,
                      const std::array<double, NUM_LAYERS>& alpha,
                      double d)
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get the rank of the current MPI process
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get the total number of MPI processes

    const int N = graph.size();
    std::vector<double> F(N), IP(N, 0.0), IPnew(N, 0.0);

    // 1) Compute normalized follower count F(u) = |Followers(u)| / N
    std::vector<int> followers(N, 0);
    #pragma omp parallel for
    for (int u = 0; u < N; ++u) {
        for (auto& e : graph[u].out[FOLLOW]) {
            if (e.to < N) {
                #pragma omp atomic
                followers[e.to]++;
            }
        }
    }
    #pragma omp parallel for
    for (int u = 0; u < N; ++u)
        F[u] = (double)followers[u] / N;

    // 2) Find max level from CAC/SCC partitioning
    int maxLevel = 0;
    for (auto& node : graph)
        maxLevel = std::max(maxLevel, node.level);

    // Broadcast maxLevel to all MPI processes
    MPI_Bcast(&maxLevel, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 3) Group nodes by level
    std::vector<std::vector<int>> byLevel(maxLevel + 1);
    for (int u = 0; u < N; ++u)
        byLevel[graph[u].level].push_back(u);

    // 4) Influence propagation level by level
    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        if (mpi_rank == 0) {
            std::cout << "[REPORT] Level " << lvl << ": " << byLevel[lvl].size() << " nodes\n";
        }

        // Distribute nodes at this level across MPI processes
        std::vector<int> local_nodes;
        for (size_t i = mpi_rank; i < byLevel[lvl].size(); i += mpi_size) {
            local_nodes.push_back(byLevel[lvl][i]);
        }

        // Compute influence power for local nodes
        std::vector<double> local_IPnew(N, 0.0);
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < local_nodes.size(); ++i) {
            int u = local_nodes[i];
            double sum = 0.0;

            // Loop through potential followers (v -> u) in FOLLOW layer
            for (int v = 0; v < N; ++v) {
                for (auto& fe : graph[v].out[FOLLOW]) {
                    if (fe.to != u) continue;

                    // Compute Ci(v, u)
                    float Ci = jaccard(graph[v].interest, graph[u].interest);

                    // Compute Nax for each layer (retweet, reply, mention)
                    std::array<int, NUM_LAYERS> Nax = {0};
                    for (auto& re : graph[v].out[RETWEET])
                        if (re.to == u) Nax[RETWEET]++;
                    for (auto& ce : graph[v].out[REPLY])
                        if (ce.to == u) Nax[REPLY]++;
                    for (auto& me : graph[v].out[MENTION])
                        if (me.to == u) Nax[MENTION]++;

                    // Compute ψ(v,u)
                    double psi = 0.0;
                    for (int l = 1; l < NUM_LAYERS; ++l)
                        psi += alpha[l] * Ci * Nax[l];

                    int deg_v = std::max(1, (int)graph[v].out[FOLLOW].size());
                    sum += psi * (IP[v] / deg_v);
                }
            }

            // Final influence power update
            local_IPnew[u] = (1.0 - d) * F[u] + d * sum;
        }

        // Gather results from all MPI processes
        MPI_Allreduce(local_IPnew.data(), IPnew.data(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Commit updates
        #pragma omp parallel for
        for (int u = 0; u < N; ++u) {
            if (graph[u].level == lvl)
                IP[u] = IPnew[u];
        }
    }

    return IP;
}
