#ifdef __cplusplus
extern "C" {
#endif

#include <metis.h>

#ifdef __cplusplus
}
#endif

#include <vector>
#include <iostream>

// Now you're safe to define C++ functions
void partitionWithMETIS(std::vector<Node>& graph, int num_parts) {
    idx_t n = graph.size();
    std::vector<idx_t> xadj(n + 1, 0);
    std::vector<idx_t> adjncy;

    for (int u = 0; u < n; ++u) {
        xadj[u + 1] = xadj[u] + graph[u].out[FOLLOW].size();
        for (const auto& edge : graph[u].out[FOLLOW]) {
            adjncy.push_back(edge.to);
        }
    }

    std::vector<idx_t> part(n);
    idx_t objval;
    METIS_PartGraphKway(&n, xadj.data(), adjncy.data(), nullptr, nullptr, nullptr, nullptr,
                        &num_parts, nullptr, nullptr, nullptr, &objval, part.data());

    for (int u = 0; u < n; ++u)
        graph[u].level = part[u];
}

