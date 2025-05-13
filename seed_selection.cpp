#include "seed_selection.h"
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>
#include <omp.h>

// ------------------------
// Algorithm 6: Candidate Selection
// ------------------------
std::vector<int> selectSeedCandidates(const std::vector<Node>& graph, const std::vector<double>& IP) {
    const int N = graph.size();
    std::unordered_map<int, std::vector<double>> byLevel;

    // Step 1: Group IPs by node level
    #pragma omp parallel
    {
        std::unordered_map<int, std::vector<double>> local_byLevel;

        #pragma omp for nowait
        for (int v = 0; v < N; ++v) {
            local_byLevel[graph[v].level].push_back(IP[v]);
        }

        #pragma omp critical
        {
            for (const auto& [level, vec] : local_byLevel) {
                byLevel[level].insert(byLevel[level].end(), vec.begin(), vec.end());
            }
        }
    }

    // Step 2: Compute IL (average IP for each level)
    std::unordered_map<int, double> IL;

    for (const auto& [level, values] : byLevel) {
        double sum = 0.0;

        #pragma omp parallel for reduction(+:sum)
        for (size_t i = 0; i < values.size(); ++i) {
            sum += values[i];
        }

        IL[level] = sum / std::max(1.0, static_cast<double>(values.size()));
    }

    // Step 3: Select candidates
    std::vector<int> candidates;

    #pragma omp parallel
    {
        std::vector<int> local_candidates;

        #pragma omp for nowait schedule(dynamic)
        for (int v = 0; v < N; ++v) {
            int level = graph[v].level;
            double avg_IP = IL[level];
            double next_IP = IL.count(level + 1) ? IL[level + 1] : 0.0;

            if ((avg_IP - next_IP) > next_IP && IP[v] > avg_IP) {
                local_candidates.push_back(v);
            }
        }

        #pragma omp critical
        candidates.insert(candidates.end(), local_candidates.begin(), local_candidates.end());
    }

    std::cout << "[REPORT] Algorithm 6 selected " << candidates.size() << " seed candidates\n";
    return candidates;
}

// ------------------------
// Algorithm 7: Final Seed Selection
// ------------------------
std::vector<int> selectFinalSeeds(const std::vector<Node>& graph,
                                  const std::vector<double>& IP,
                                  const std::vector<int>& Istar) {
    std::vector<int> INF;
    std::set<int> Iset(Istar.begin(), Istar.end());
    std::unordered_map<int, std::vector<int>> trees;
    std::unordered_map<int, int> treeSize;

    // Step 1: Build Influence-BFS Trees
    #pragma omp parallel
    {
        std::unordered_map<int, std::vector<int>> local_trees;
        std::unordered_map<int, int> local_treeSize;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < Istar.size(); ++i) {
            int root = Istar[i];
            std::queue<int> q;
            std::set<int> visited;

            q.push(root);
            visited.insert(root);

            while (!q.empty()) {
                int u = q.front();
                q.pop();

                local_trees[root].push_back(u);

                for (const auto& edge : graph[u].out[FOLLOW]) {
                    if (Iset.count(edge.to) && visited.insert(edge.to).second) {
                        q.push(edge.to);
                    }
                }
            }

            local_treeSize[root] = local_trees[root].size();
        }

        #pragma omp critical
        {
            for (const auto& [root, nodes] : local_trees) {
                trees[root].insert(trees[root].end(), nodes.begin(), nodes.end());
            }
            for (const auto& [root, size] : local_treeSize) {
                treeSize[root] = size;
            }
        }
    }

    // Step 2: Iteratively select seeds from largest influence trees
    std::set<int> remaining(Istar.begin(), Istar.end());

    while (!remaining.empty()) {
        std::vector<int> active_roots;

        for (const auto& [root, _] : trees) {
            if (remaining.count(root)) {
                active_roots.push_back(root);
            }
        }

        int umax = -1, max_size = -1;

        #pragma omp parallel
        {
            int local_umax = -1, local_max = -1;

            #pragma omp for nowait
            for (size_t i = 0; i < active_roots.size(); ++i) {
                int root = active_roots[i];
                int size = trees[root].size();

                if (size > local_max) {
                    local_umax = root;
                    local_max = size;
                }
            }

            #pragma omp critical
            {
                if (local_max > max_size) {
                    umax = local_umax;
                    max_size = local_max;
                }
            }
        }

        // Step 3: Find BLACK path: intersection of T_umaX and I*
        std::vector<int> BLACK;
        for (int v : trees[umax]) {
            if (remaining.count(v)) {
                BLACK.push_back(v);
            }
        }

        // Step 4: Select node with minimum IP in BLACK
        int vmin = BLACK[0];
        double min_ip = IP[vmin];

        #pragma omp parallel
        {
            int local_vmin = vmin;
            double local_min_ip = min_ip;

            #pragma omp for nowait
            for (size_t i = 1; i < BLACK.size(); ++i) {
                int v = BLACK[i];
                if (IP[v] < local_min_ip) {
                    local_vmin = v;
                    local_min_ip = IP[v];
                }
            }

            #pragma omp critical
            {
                if (local_min_ip < min_ip) {
                    vmin = local_vmin;
                    min_ip = local_min_ip;
                }
            }
        }

        INF.push_back(vmin);

        for (int v : BLACK) {
            remaining.erase(v);
        }
        remaining.erase(vmin);
    }

    std::cout << "[REPORT] Algorithm 7 selected " << INF.size() << " final seeds\n";
    return INF;
}

