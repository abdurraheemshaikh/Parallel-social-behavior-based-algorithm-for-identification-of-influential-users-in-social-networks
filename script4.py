from collections import defaultdict

input_file = 'endorsement_weights.txt'
output_file = 'undirected_graph.txt'

# Dictionary to store weights of undirected edges
edge_weights = defaultdict(float)

with open(input_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 3:
            continue  # skip malformed lines
        node1, node2, weight = parts
        weight = float(weight)
        # Sort the node pair to ensure undirected representation
        key = tuple(sorted((node1, node2)))
        edge_weights[key] += weight

# Write the undirected edges with summed weights
with open(output_file, 'w') as f:
    for (u, v), weight in edge_weights.items():
        f.write(f"{u} {v} {weight:.6f}\n")

