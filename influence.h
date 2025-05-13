#pragma once
#include "node.h"
#include <vector>
#include <array>

/// Compute the influence power IP[u] for each node u ∈ V
/// Inputs:
///   graph         - graph with FOLLOW + other action layers
///   alpha_factors - weights for [RETWEET, REPLY, MENTION] (must sum to 1)
///   d             - damping factor (e.g., 0.85)
///
/// Returns:
///   vector<double> where IP[u] is the influence of node u
std::vector<double>
computeInfluencePower(const std::vector<Node>& graph,
                      const std::array<double, NUM_LAYERS>& alpha_factors,
                      double d);