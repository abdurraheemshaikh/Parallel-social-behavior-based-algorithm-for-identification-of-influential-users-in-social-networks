import os

# Jaccard similarity as proxy for common interest
def jaccard(set1, set2):
    if not set1 or not set2:
        return 0.0
    return len(set1 & set2) / len(set1 | set2)

# Read features: node -> set of feature indices (with 1's)
def read_features(file_path):
    features = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            node = parts[0]
            vec = set(i for i, v in enumerate(parts[1:]) if v == '1')
            features[node] = vec
    return features

# Read edges from twitter_combined
def read_edges(file_path):
    edges = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                edges.append((parts[0], parts[1]))
    return edges

# Paths
edges_file = 'twitter_combined.txt'
feat_file = 'new_feat.txt'
output_file = 'endorsement_weights.txt'

# Load data
features = read_features(feat_file)
edges = read_edges(edges_file)

# Compute endorsement weight for each edge
with open(output_file, 'w') as out:
    for ux, uy in edges:
        if ux in features and uy in features:
            common_interest = jaccard(features[ux], features[uy])
            N_ai = 1  # assuming 1 interaction occurred
            N_py = len(features[uy]) if features[uy] else 1  # avoid division by 0
            endorsement = (common_interest * N_ai) / N_py
            out.write(f"{ux} {uy} {endorsement:.6f}\n")

