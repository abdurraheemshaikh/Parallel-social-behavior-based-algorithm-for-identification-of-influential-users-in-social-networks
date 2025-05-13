def analyze_graph(file_path):
    unique_nodes = set()
    edge_count = 0

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 2:
                node1, node2 = parts[0], parts[1]
                unique_nodes.add(node1)
                unique_nodes.add(node2)
                edge_count += 1
            else:
                print(f"Skipping invalid line: {line.strip()}")

    print(f"Number of unique nodes (vertices): {len(unique_nodes)}")
    print(f"Total number of edges: {edge_count}")
    print(f"Total number of vertices: {len(unique_nodes)}")

# Example usage
file_path = 'graph.txt'  # Replace with your actual file path
analyze_graph(file_path)

