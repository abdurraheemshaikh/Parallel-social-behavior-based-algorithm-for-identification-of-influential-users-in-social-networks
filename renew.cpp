#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <mpi.h>
#include <metis.h>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <queue>
#include <map>
using namespace std;

class Graph {
public:
    vector<vector<int>> adj;
    vector<vector<float>> weights;
    unordered_map<int, int> nodeToIndex;
    vector<int> indexToNode;

    void addEdge(int uOrig, int vOrig, float weight) {
        if (nodeToIndex.find(uOrig) == nodeToIndex.end()) {
            int index = nodeToIndex.size();
            nodeToIndex[uOrig] = index;
            indexToNode.push_back(uOrig);
            adj.emplace_back();
            weights.emplace_back();
        }
        if (nodeToIndex.find(vOrig) == nodeToIndex.end()) {
            int index = nodeToIndex.size();
            nodeToIndex[vOrig] = index;
            indexToNode.push_back(vOrig);
            adj.emplace_back();
            weights.emplace_back();
        }

        int u = nodeToIndex[uOrig];
        int v = nodeToIndex[vOrig];

        adj[u].push_back(v);
        adj[v].push_back(u);
        weights[u].push_back(weight);
        weights[v].push_back(weight);
    }

    void readGraphFromFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            exit(EXIT_FAILURE);
        }

        int u, v;
        float weight;
        while (file >> u >> v >> weight) {
            addEdge(u, v, weight);
        }
        file.close();
    }

    void broadcastGraph(int rank, MPI_Comm comm) {
        int numVertices = adj.size();
        MPI_Bcast(&numVertices, 1, MPI_INT, 0, comm);

        if (rank != 0) {
            adj.resize(numVertices);
            weights.resize(numVertices);
        }

        vector<int> adjSizes(numVertices);
        vector<int> flatAdj;
        vector<float> flatWeights;
        vector<int> nodeLabels;

        if (rank == 0) {
            for (int i = 0; i < numVertices; ++i) {
                adjSizes[i] = adj[i].size();
                flatAdj.insert(flatAdj.end(), adj[i].begin(), adj[i].end());
                flatWeights.insert(flatWeights.end(), weights[i].begin(), weights[i].end());
            }
            nodeLabels = indexToNode;
        }

        MPI_Bcast(adjSizes.data(), numVertices, MPI_INT, 0, comm);

        int totalEdges = 0;
        if (rank == 0) totalEdges = flatAdj.size();
        else {
            for (int sz : adjSizes) totalEdges += sz;
            flatAdj.resize(totalEdges);
            flatWeights.resize(totalEdges);
        }

        MPI_Bcast(flatAdj.data(), totalEdges, MPI_INT, 0, comm);
        MPI_Bcast(flatWeights.data(), totalEdges, MPI_FLOAT, 0, comm);

        if (rank != 0) {
            for (int i = 0, idx = 0; i < numVertices; ++i) {
                adj[i].resize(adjSizes[i]);
                weights[i].resize(adjSizes[i]);
                for (int j = 0; j < adjSizes[i]; ++j) {
                    adj[i][j] = flatAdj[idx];
                    weights[i][j] = flatWeights[idx];
                    idx++;
                }
            }
        }

        nodeLabels.resize(numVertices);
        MPI_Bcast(nodeLabels.data(), numVertices, MPI_INT, 0, comm);
        indexToNode = nodeLabels;
        nodeToIndex.clear();
        for (int i = 0; i < numVertices; ++i) {
            nodeToIndex[indexToNode[i]] = i;
        }
    }

    int size() const { return adj.size(); }

    const vector<int>& operator[](int idx) const { return adj[idx]; }
    vector<int>& operator[](int idx) { return adj[idx]; }

    const vector<float>& getWeights(int idx) const { return weights[idx]; }
    vector<float>& getWeights(int idx) { return weights[idx]; }
};

void partitionGraphMETIS(Graph& g, idx_t numParts, vector<idx_t>& part, int rank) {
    idx_t nvtxs = g.size();
    vector<idx_t> xadj(nvtxs + 1, 0);
    vector<idx_t> adjncy;

    for (idx_t i = 0; i < nvtxs; ++i) {
        xadj[i + 1] = xadj[i] + g[i].size();
        for (int neighbor : g[i]) adjncy.push_back(neighbor);
    }

    part.resize(nvtxs);
    idx_t objval;
    idx_t ncon = 1;

    int result = METIS_PartGraphKway(
        &nvtxs, &ncon,
        xadj.data(), adjncy.data(),
        NULL, NULL, NULL,
        &numParts, NULL, NULL,
        NULL, &objval,
        part.data()
    );

    if (result != METIS_OK) {
        cerr << "[Rank " << rank << "] METIS failed with code: " << result << endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

int computeLevels(Graph& g, const vector<int>& startNodes, const vector<int>& part, int rank, const string& label) {
    unordered_map<int, int> levels;
    unordered_set<int> visited;
    queue<pair<int, int>> q;
    int maxLevel = 0;

    for (int node : startNodes) {
        int idx = g.nodeToIndex[node];
        if (visited.count(idx)) continue;
        q.push({idx, 1});
        visited.insert(idx);
        levels[g.indexToNode[idx]] = 1;

        while (!q.empty()) {
            auto [cur, lvl] = q.front(); q.pop();
            maxLevel = max(maxLevel, lvl);
            for (int nei : g[cur]) {
                if (part[nei] == rank && !visited.count(nei)) {
                    visited.insert(nei);
                    levels[g.indexToNode[nei]] = lvl + 1;
                    q.push({nei, lvl + 1});
                }
            }
        }
    }

    if (!levels.empty()) {
        cout << "[Rank " << rank << "] Max level reached in " << label << ": " << maxLevel << "\n";
    }

    return maxLevel;
}

void computeSACandCAC(Graph &g, const vector<int> &part, int rank) {
    vector<int> sacNodes;
    vector<int> cacNodes;

    for (int u = 0; u < g.adj.size(); ++u) {
        if (part[u] != rank) continue;
        for (int v : g.adj[u]) {
            if (part[v] != rank) {
                sacNodes.push_back(g.indexToNode[u]);
                break;
            }
        }
    }

    for (int u = 0; u < g.adj.size(); ++u) {
        if (part[u] != rank) continue;
        for (int v : g.adj[u]) {
            if (part[v] != rank) {
                cacNodes.push_back(g.indexToNode[u]);
                break;
            }
        }
    }

    if (!sacNodes.empty()) {
        int sacMaxLevel = computeLevels(g, sacNodes, part, rank, "SAC");
        cout << "[Rank " << rank << "] SAC exists with max level " << sacMaxLevel << "\n";
    } else {
        cout << "[Rank " << rank << "] SAC does not exist.\n";
    }

    if (!cacNodes.empty()) {
        int cacMaxLevel = computeLevels(g, cacNodes, part, rank, "CAC");
        cout << "[Rank " << rank << "] CAC exists. Max level: " << cacMaxLevel << "\n";
    } else {
        cout << "[Rank " << rank << "] No CAC found.\n";
    }
}

void computePageRank(Graph& g, vector<float>& pagerank, float alpha = 0.85, int maxIter = 100, float tol = 1e-6) {
    int n = g.size();
    pagerank.assign(n, 1.0f / n);
    vector<float> newRank(n, 0);

    for (int iter = 0; iter < maxIter; ++iter) {
        float maxDiff = 0.0f;
        #pragma omp parallel num_threads(2)
        #pragma omp parallel for reduction(max:maxDiff)
        for (int i = 0; i < n; ++i) {
            newRank[i] = (1 - alpha) / n;
            for (int j : g[i]) {
                newRank[i] += alpha * pagerank[j] / g[j].size();
            }
            maxDiff = max(maxDiff, fabs(newRank[i] - pagerank[i]));
        }

        pagerank = newRank;
        if (maxDiff < tol) break;
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    Graph g;
    if (rank == 0) {
        g.readGraphFromFile("graph.txt");
    }

    vector<idx_t> part;
    if (rank == 0) {
        partitionGraphMETIS(g, numProcs, part, rank);
    }

    g.broadcastGraph(rank, MPI_COMM_WORLD);

    int numNodes = g.size();
    if (rank != 0) part.resize(numNodes);
    MPI_Bcast(part.data(), numNodes, MPI_INT, 0, MPI_COMM_WORLD);

    int localCount = 0;
    for (int i = 0; i < part.size(); ++i) {
        if (part[i] == rank) ++localCount;
    }

    cout << "[Rank " << rank << "] Total nodes received: " << localCount << endl;

    if (rank == 0) {
        cout << "Total number of nodes in the graph: " << numNodes << endl;
    }

    computeSACandCAC(g, part, rank);

    vector<float> pagerank(g.size());
    computePageRank(g, pagerank);

    if (rank == 0) {
        vector<pair<int, float>> pagerankIndex(g.size());
        for (int i = 0; i < g.size(); ++i) {
            pagerankIndex[i] = {i, pagerank[i]};
        }
        sort(pagerankIndex.rbegin(), pagerankIndex.rend(), [](const pair<int, float>& a, const pair<int, float>& b) {
            return a.second < b.second;
        });

        int K = 2;
        cout << "Top " << K << " nodes by PageRank:" << endl;
        for (int i = 0; i < K && i < pagerankIndex.size(); ++i) {
            cout << "Node " << g.indexToNode[pagerankIndex[i].first] << " with rank: " << pagerankIndex[i].second << endl;
        }
    }

    MPI_Finalize();
    return 0;
}
