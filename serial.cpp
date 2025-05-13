#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <queue>        // For std::queue
#include <utility>      // For std::pair

using namespace std;

class Graph {
public:
    vector<vector<int>> adj;       // Adjacency list
    unordered_map<int, string> indexToNode;
    unordered_map<string, int> nodeToIndex;
    vector<vector<float>> weights;

    int size() const { return adj.size(); }

    const vector<int>& operator[](int idx) const { return adj[idx]; }
    vector<int>& operator[](int idx) { return adj[idx]; }

    const vector<float>& getWeights(int idx) const { return weights[idx]; }
    vector<float>& getWeights(int idx) { return weights[idx]; }

 void readGraphFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        return;
    }

    string node1, node2;
    float weight;
    
    // Read the edges and weights from the file
    while (file >> node1 >> node2 >> weight) {
        int u = -1, v = -1;

        // If the node1 or node2 has not been seen before, add it
        if (nodeToIndex.find(node1) == nodeToIndex.end()) {
            u = nodeToIndex.size();  // Assign an index to this node
            nodeToIndex[node1] = u;
            indexToNode[u] = node1;
        } else {
            u = nodeToIndex[node1];
        }

        if (nodeToIndex.find(node2) == nodeToIndex.end()) {
            v = nodeToIndex.size();
            nodeToIndex[node2] = v;
            indexToNode[v] = node2;
        } else {
            v = nodeToIndex[node2];
        }

        // Ensure the adjacency list and weight matrix have space for the new nodes
        while (adj.size() <= max(u, v)) {
            adj.push_back(vector<int>());
            weights.push_back(vector<float>());
        }

        // Add the edge in the adjacency list (undirected)
        adj[u].push_back(v);
        adj[v].push_back(u);

        // Add the weight in the weights matrix
        weights[u].push_back(weight);
        weights[v].push_back(weight);
    }

    cout << "Graph read successfully." << endl;
    cout << "Total nodes: " << nodeToIndex.size() << endl;
    cout << "Total edges: " << file.gcount() << endl;
}



};

void computeSACandCAC(Graph &g) {
    vector<string> sacNodes;  // Changed to vector<string> to store node names
    vector<string> cacNodes;  // Changed to vector<string> to store node names

    for (int u = 0; u < g.adj.size(); ++u) {
        for (int v : g.adj[u]) {
            if (u != v) {
                sacNodes.push_back(g.indexToNode[u]);  // Now pushing back node names
                break;
            }
        }
    }

    for (int u = 0; u < g.adj.size(); ++u) {
        for (int v : g.adj[u]) {
            if (u != v) {
                cacNodes.push_back(g.indexToNode[u]);  // Now pushing back node names
                break;
            }
        }
    }

    auto computeLevels = [&](const vector<string>& startNodes, const string& label) -> int {
        unordered_map<int, int> levels;
        unordered_set<int> visited;
        queue<pair<int, int>> q;
        int maxLevel = 0;

        for (const string& node : startNodes) {
            int idx = g.nodeToIndex[node];  // Correct usage: accessing map with string key
            if (visited.count(idx)) continue;
            q.push({idx, 1});
            visited.insert(idx);
           levels[idx] = 1;



            while (!q.empty()) {
                auto [cur, lvl] = q.front(); q.pop();
                maxLevel = max(maxLevel, lvl);
                for (int nei : g[cur]) {
                    if (!visited.count(nei)) {
                        visited.insert(nei);
                        levels[nei] = lvl + 1;

                        q.push({nei, lvl + 1});
                    }
                }
            }
        }

        if (!levels.empty()) {
            cout << "[Serial] Max level reached in " << label << ": " << maxLevel << "\n";
        }

        return maxLevel;
    };

    if (!sacNodes.empty()) {
        int sacMaxLevel = computeLevels(sacNodes, "SAC");
        cout << "[Serial] SAC exists with max level " << sacMaxLevel << "\n";
    } else {
        cout << "[Serial] SAC does not exist.\n";
    }

    if (!cacNodes.empty()) {
        int cacMaxLevel = computeLevels(cacNodes, "CAC");
        cout << "[Serial] CAC exists. Max level: " << cacMaxLevel << "\n";
    } else {
        cout << "[Serial] No CAC found.\n";
    }
}

void computePageRank(Graph& g, vector<float>& pagerank, float alpha = 0.85, int maxIter = 100, float tol = 1e-6) {
    int n = g.size();
    pagerank.assign(n, 1.0f / n);
    vector<float> newRank(n, 0);

    for (int iter = 0; iter < maxIter; ++iter) {
        float maxDiff = 0.0;
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

int main() {
    Graph g;
    g.readGraphFromFile("graph.txt");

    cout << "Total number of nodes in the graph: " << g.size() << endl;

    computeSACandCAC(g);

    vector<float> pagerank(g.size());
    computePageRank(g, pagerank);

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
        cout << g.indexToNode[pagerankIndex[i].first] << ": " << pagerankIndex[i].second << endl;
    }

    return 0;
}

