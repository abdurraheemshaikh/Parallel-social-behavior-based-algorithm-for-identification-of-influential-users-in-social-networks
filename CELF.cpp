#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <ctime>
#include <sstream>  // For istringstream

using namespace std;

using Graph = vector<vector<int>>;

// Monte Carlo simulation of Independent Cascade (IC) model
int simulate_IC(const Graph& G, const unordered_set<int>& seed_set, double p = 0.01, int mc = 100) {
    cout<<"Monte Carlos simulation....."<<endl;
    int total_spread = 0;
    mt19937 gen(time(0));
    uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < mc; ++i) {
        unordered_set<int> activated(seed_set.begin(), seed_set.end());
        queue<int> q;
        for (int seed : seed_set) q.push(seed);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : G[u]) {
                if (activated.find(v) == activated.end() && dis(gen) < p) {
                    activated.insert(v);
                    q.push(v);
                }
            }
        }
        total_spread += activated.size();
    }

    return total_spread / mc;
}

struct NodeGain {
    int node;
    int gain;
    int last_updated;

    bool operator<(const NodeGain& other) const {
        return gain < other.gain;
    }
};

// CELF++ algorithm
vector<int> CELFpp(const Graph& G, int k, double p = 0.01, int mc = 100) {
    vector<NodeGain> queue;
    unordered_set<int> seed_set;
    vector<int> result;
    cout<<"CELF ++ initiating......"<<endl;
    // Initial marginal gains
    for (int u = 0; u < G.size(); ++u) {
        int spread = simulate_IC(G, {u}, p, mc);
        queue.push_back({u, spread, 0});
    }

    sort(queue.begin(), queue.end(), [](const NodeGain& a, const NodeGain& b) {
        return a.gain > b.gain;
    });

    int iter = 1;
    while (result.size() < k) {
        NodeGain top = queue.front();
        pop_heap(queue.begin(), queue.end()); queue.pop_back();
    
        if (top.last_updated == result.size()) {
            result.push_back(top.node);
            seed_set.insert(top.node);
        } else {
            // Manually compute union of seed_set and {top.node}
            unordered_set<int> temp = seed_set;
            temp.insert(top.node);
            
            int new_gain = simulate_IC(G, temp, p, mc) - simulate_IC(G, seed_set, p, mc);
            queue.push_back({top.node, new_gain, (int)result.size()});
            push_heap(queue.begin(), queue.end());
        }
    
        iter++;
    }
    

    return result;
}

// METIS parser
Graph parse_metis(const string& filename) {
    ifstream infile(filename);
    int n, m;
    infile >> n >> m;
    Graph G(n);

    for (int u = 0; u < n; ++u) {
        string line;
        getline(infile, line);
        if (line.empty()) {
            getline(infile, line);
        }

        istringstream iss(line);
        int v;
        while (iss >> v) {
            if (v - 1 != u)
                G[u].push_back(v - 1);  // Convert to 0-indexed
        }
    }

    return G;
}

int main() {
    string filename = "global.graph";
    int k = 10;  // number of influential nodes to find
    double propagation_prob = 0.01;
    int mc_simulations = 100;

    Graph G = parse_metis(filename);

    vector<int> influential_nodes = CELFpp(G, k, propagation_prob, mc_simulations);

    cout << "Top " << k << " influential nodes:\n";
    for (int u : influential_nodes) {
        cout << u + 1 << " "<<endl; // Output 1-indexed for consistency
    }
    cout << endl;

    return 0;
}
