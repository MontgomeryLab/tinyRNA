#!/bin/bash
# Counts unique sequences in an input fastq file and prints them to the specified output file as pseudo-json
# Format: {"sequence": count, ...}

if (( $# != 2 )); then
	echo "Usage: $(basename -- $0) input-file output-file"
	exit 0
fi

echo "{" > "$2"
awk 'NR % 2 == 0 && NR % 4 != 0 {counts[$0]++} END {for (seq in counts) print "\t\""seq"\"", ":", counts[seq]","}' "$1" >> "$2"
echo "}" >> "$2"
