#include <stack>
#include <queue>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <omp.h>
#include "partition.h"

void computeSCC_CAC(std::vector<Node>& graph) {
    const int N = static_cast<int>(graph.size());

    // Build FOLLOW and reverse adjacency lists
    std::vector<std::vector<int>> adj(N), radj(N);
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < N; ++u) {
        for (const auto& e : graph[u].out[FOLLOW]) {
            int v = e.to;
            if (v >= 0 && v < N) {
                #pragma omp critical
                {
                    adj[u].push_back(v);
                    radj[v].push_back(u);
                }
            }
        }
    }

    // First DFS pass: compute finish order
    std::vector<char> seen(N, 0);
    std::vector<int> order;
    order.reserve(N);

    #pragma omp parallel
    {
        std::vector<int> local_order;
        local_order.reserve(N);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < N; ++i) {
            if (!seen[i]) {
                std::stack<int> local_dfs;
                local_dfs.push(i << 1); // even=enter, odd=exit
                while (!local_dfs.empty()) {
                    int x = local_dfs.top();
                    local_dfs.pop();
                    int v = x >> 1;
                    bool exiting = x & 1;
                    if (exiting) {
                        local_order.push_back(v);
                    } else if (!seen[v]) {
                        seen[v] = 1;
                        local_dfs.push((v << 1) | 1);
                        for (int w : adj[v]) {
                            if (!seen[w]) local_dfs.push(w << 1);
                        }
                    }
                }
            }
        }

        #pragma omp critical
        order.insert(order.end(), local_order.begin(), local_order.end());
    }

    // Second DFS pass: identify SCCs on reversed graph
    std::vector<int> compID(N, -1);
    int componentCount = 0;
    std::stack<int> dfs;
    for (int i = N - 1; i >= 0; --i) {
        int v = order[i];
        if (compID[v] >= 0) continue;

        dfs.push(v);
        compID[v] = componentCount;
        std::vector<int> members;

        while (!dfs.empty()) {
            int u = dfs.top(); dfs.pop();
            members.push_back(u);
            for (int w : radj[u]) {
                if (compID[w] < 0) {
                    compID[w] = componentCount;
                    dfs.push(w);
                }
            }
        }

        componentCount++;
    }

    // Set node attributes: compID, SCC type, level
    std::vector<int> compSize(componentCount, 0);
    #pragma omp parallel for
    for (int v = 0; v < N; ++v) {
        graph[v].compID = compID[v];
        graph[v].type = Node::SCC;
        graph[v].level = 0;

        #pragma omp atomic
        compSize[compID[v]]++;
    }

    // Mark CACs (singleton components without self-loops)
    #pragma omp parallel for
    for (int v = 0; v < N; ++v) {
        if (compSize[compID[v]] == 1) {
            graph[v].type = Node::CAC;
        }
    }

    // Build component-level adjacency (DAG)
    std::vector<std::vector<int>> cadj(componentCount);
    #pragma omp parallel for
    for (int u = 0; u < N; ++u) {
        int cu = compID[u];
        for (int v : adj[u]) {
            int cv = compID[v];
            if (cu != cv) {
                #pragma omp critical
                cadj[cu].push_back(cv);
            }
        }
    }

    #pragma omp parallel for
    for (auto& nbrs : cadj) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    // Compute topological levels via BFS on component DAG
    std::vector<int> indeg(componentCount, 0);
    for (int u = 0; u < componentCount; ++u)
        for (int v : cadj[u]) indeg[v]++;

    std::queue<int> q;
    for (int u = 0; u < componentCount; ++u)
        if (indeg[u] == 0) q.push(u);

    std::vector<int> clevel(componentCount, 0);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : cadj[u]) {
            clevel[v] = std::max(clevel[v], clevel[u] + 1);
            if (--indeg[v] == 0) q.push(v);
        }
    }

    // Assign node-level topological level
    #pragma omp parallel for
    for (int v = 0; v < N; ++v) {
        graph[v].level = clevel[compID[v]];
    }

    // CAC-merge: merge singleton CACs with one-level-down CAC neighbors
    std::vector<std::vector<int>> membersOf(componentCount);
    for (int v = 0; v < N; ++v)
        membersOf[compID[v]].push_back(v);

    for (int c = 0; c < componentCount; ++c) {
        if (membersOf[c].size() != 1) continue;
        int u = membersOf[c][0];
        if (graph[u].type != Node::CAC) continue;

        int lvl = graph[u].level;
        for (const auto& e : graph[u].out[FOLLOW]) {
            int v = e.to;
            if (v < 0 || v >= N) continue;
            int cv = compID[v];
            if (graph[v].type == Node::CAC && clevel[cv] == lvl - 1) {
                for (int w : membersOf[cv]) {
                    compID[w] = c;
                    graph[w].compID = c;
                }
                break; // Merge once
            }
        }
    }

    std::cout << "[REPORT] computeSCC_CAC: " << componentCount << " initial SCCs, CAC-merged singleton chains\n";
}

