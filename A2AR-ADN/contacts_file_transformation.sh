#!/bin/bash

# Check if input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input_file> [output_file]"
    exit 1
fi

# Assign input and output file names
INPUT_FILE="$1"
OUTPUT_FILE="${2:-final_contacts.tsv}"  # Default output filename

# Process the file using awk
awk '
/freeSelLabel/ {
    sub(/^freeSelLabel /, "", $0);   # Remove "freeSelLabel"
    split($0, arr, "==");            # Split using "=="
    aa_resid = arr[1];               # Keep only the part before "=="
    next;
}
/^[0-9]+ [0-9]+$/ {
    print aa_resid "\t" $1 "\t" $2;  # Use tabs explicitly
}' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Final contacts saved to $OUTPUT_FILE"

