#!/bin/bash
# Counts unique sequences in an input fastq file, then prints each unique sequence with a corresponding count.
# Output is a json object

if (( $# != 2 )); then
	echo "Usage: $(basename -- $0) input-file output-file"
	exit 0
fi

awk '
  BEGIN {print "{"; n = 0; i = 0}
  NR % 2 == 0 && NR % 4 != 0 {if (!($0 in counts)) n++; counts[$0]++}
  END {
    for (seq in counts) {
      printf "\t\"%s\": %d",seq,counts[seq]
      if (i++ < n-1) print ","
    }
  }
  END {print "\n}"}' "$1" >> "$2"
