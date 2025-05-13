# File path to your directed graph
input_file = 'endorsement_weights.txt'

# Store edges in a dictionary for quick lookup
edges = {}

# Read all edges and store both directions
with open(input_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 3:
            continue  # Skip malformed lines
        u, v, w = parts
        edges[(u, v)] = w

# Display mutual (bidirectional) edges
printed = set()

print("node1 node2 weight")

for (u, v), w1 in edges.items():
    if (v, u) in edges and (v, u) not in printed:
        w2 = edges[(v, u)]
        print(f"{u} {v} {w1}")
        print(f"{v} {u} {w2}")
        printed.add((u, v))
        printed.add((v, u))

