#pragma once
#include <vector>
#include "node.h"

/// Selects seed candidates I* using level-wise thresholds (Algorithm 6)
std::vector<int> selectSeedCandidates(const std::vector<Node>& graph, const std::vector<double>& IP);

/// Final seed selection using Influence-BFS ranking (Algorithm 7)
std::vector<int> selectFinalSeeds(const std::vector<Node>& graph, const std::vector<double>& IP, const std::vector<int>& candidates);