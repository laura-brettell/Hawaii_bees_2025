#!/bin/bash

# run from Hawaii_old
# Paths to directory
COUNTS_DIR="./mapping_results/counts"




echo "🛠️ Merging all counts into final table..."
# Merge all individual counts into a single matrix
csv_files=($COUNTS_DIR/*_counts.csv)
paste -d, "${csv_files[@]}" | \
awk -F, '{printf "%s",$1; for (i=2; i<=NF; i+=2) printf ",%s", $i; print ""}' \
> virus_counts_matrix.csv

echo "🎉 Pipeline complete! Final counts table: virus_counts_matrix.csv"
