# sort_graph_by_node1.py
input_file = "undirected_graph.txt"
output_file = "sorted_graph.txt"

with open(input_file, "r") as f:
    lines = f.readlines()

# Sort by the first node (column 0)
lines.sort(key=lambda line: int(line.split()[0]))  # or just `line.split()[0]` if nodes are strings

with open(output_file, "w") as f:
    f.writelines(lines)

print(f"Graph sorted by node1 and saved to '{output_file}'.")

