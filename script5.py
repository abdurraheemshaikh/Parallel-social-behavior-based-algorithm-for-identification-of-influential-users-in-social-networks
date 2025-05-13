def search_edge(file_path, node_a, node_b):
    found = False
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            u, v = parts[0], parts[1]
            if (u == node_a and v == node_b) or (u == node_b and v == node_a):
                print(f"Found: {line.strip()}")
                found = True
    if not found:
        print("No such edge found.")

# Example usage
if __name__ == "__main__":
    file_path = 'undirected_graph.txt'  # change if using a different file
    node1 = input("Enter first node: ").strip()
    node2 = input("Enter second node: ").strip()
    search_edge(file_path, node1, node2)

