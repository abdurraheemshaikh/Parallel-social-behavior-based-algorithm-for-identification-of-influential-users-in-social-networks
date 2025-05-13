#!/bin/bash

# Input and output filenames
INPUT_FILE="undirected_graph.txt"
OUTPUT_FILE="undirected_graph_10mb.txt"

# Copy the first 500KB (1 kilobyte = 1024 bytes, 500KB = 512000 bytes)
dd if="$INPUT_FILE" of="$OUTPUT_FILE" bs=10M count=1 status=progress

echo "Copied first KB of $INPUT_FILE into $OUTPUT_FILE"

