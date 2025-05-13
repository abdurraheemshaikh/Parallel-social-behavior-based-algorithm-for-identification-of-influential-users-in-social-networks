#pragma once

#include <vector>
#include "node.h"    // your Node definition with out[], index, lowlink, compID, level, type, onStack

/// Compute SCC + CAC partitioning on the FOLLOW-graph stored in `graph`.
/// After this call:
///   - graph[v].compID is a 0-based component ID
///   - graph[v].type   is SCC or CAC
///   - graph[v].level  is its depth in the component DAG
void computeSCC_CAC(std::vector<Node>& graph);