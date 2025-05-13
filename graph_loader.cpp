#include "graph_loader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <omp.h>
#include <thread>  // For multi-threading
#include <chrono>  // For time management

constexpr uint32_t MAX_TEST_NODES = 10000;

void showLoadingSpinner() {
    const char spinner[] = {'|', '/', '-', '\\'};
    int spinner_index = 0;
    while (true) {
        std::cout << "\r[LOADING..] " << spinner[spinner_index] << std::flush;
        spinner_index = (spinner_index + 1) % 4;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));  // Spinner speed
    }
}

void loadGraph(std::vector<Node>& graph) {
    for (int layer = 0; layer < NUM_LAYERS; ++layer) {
        std::cout << "[REPORT] Loading layer " << layer_names[layer]
                  << " from \"" << files[layer] << "\"\n";

        // Start the spinner only for the first layer (layer == 0)
        if (layer == 0) {
            std::thread spinnerThread(showLoadingSpinner);  // Create a separate thread for the spinner
            spinnerThread.detach();  // Let the spinner run independently
        }

        std::ifstream infile(files[layer]);
        if (!infile.is_open()) {
            std::cerr << "[ERROR] Cannot open " << files[layer] << "\n";
            exit(1);
        }

        std::vector<std::string> lines;
        std::string line;
        while (std::getline(infile, line)) {
            if (!line.empty() && line[0] != '#') {
                lines.push_back(line);
            }
        }

        uint64_t edgeCount = 0;

        #pragma omp parallel for schedule(dynamic) reduction(+:edgeCount)
        for (size_t i = 0; i < lines.size(); ++i) {
            std::istringstream ss(lines[i]);
            uint32_t u, v;
            float w = 1.0f;

            ss >> u >> v;
            if (!(ss >> w)) w = 1.0f;

            if (u >= MAX_TEST_NODES || v >= MAX_TEST_NODES) continue;

            #pragma omp critical
            graph[u].out[layer].emplace_back(Edge{v, w, 0});

            if (++edgeCount % 1000000 == 0) {
                #pragma omp critical
                std::cout << "[REPORT]   " << edgeCount
                          << " edges loaded in " << layer_names[layer] << "\n";
            }
        }

        // Stop spinner after loading the first layer (layer == 0)
        if (layer == 0) {
            std::cout << "\r[LOADING] Done.              " << std::endl;  // Clear spinner line
        }

        std::cout << "[REPORT] Completed layer " << layer_names[layer]
                  << ": total edges = " << edgeCount << "\n";
    }
}

