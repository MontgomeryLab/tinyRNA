#!/bin/bash
# Counts unique sequences in an input fastq file, then prints each unique sequence with a corresponding count.
# Output is a json object

if (( $# != 2 )); then
	echo "Usage: $(basename -- $0) input-file output-file"
	exit 0
fi

awk '
  BEGIN { print "{" }
  NR % 2 == 0 && NR % 4 != 0 {
    if (!($0 in counts)) idx_to_seq[n++] = $0
    counts[$0]++
  }
  END {
    for (i = 0; i < n-1; i++) {
      seq = idx_to_seq[i]
      printf "\t\"%s\": %d,\n",seq,counts[seq]
    }
    # No trailing comma for the last sequence
    last_seq = idx_to_seq[n-1]
    printf "\t\"%s\": %d\n}",last_seq,counts[last_seq]
  }' "$1" >> "$2"