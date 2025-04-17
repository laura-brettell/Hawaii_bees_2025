#!/bin/bash

# Path to mapping results
MAPPING_DIR="./mapping_results"
COUNTS_DIR="./counts_filtered"

# Create output directory if it doesn't exist
mkdir -p $COUNTS_DIR

# Loop through each sample's BAM file
for bam in $MAPPING_DIR/*.sorted.bam; do
    sample_name=$(basename $bam .sorted.bam)

    # Calculate total mapped reads
    total_reads=$(samtools idxstats $bam | awk '{sum += $3} END {print sum}')

    # Filter out viruses with < 0.01% of total reads
    samtools idxstats $bam | \
    awk -v total="$total_reads" -v sample="$sample_name" \
    'BEGIN {print "virus,"sample} {if ($3 / total >= 0.0001) print $1","$3}' \
    > $COUNTS_DIR/${sample_name}_counts_filtered.csv

    echo "✅ Filtered counts for $sample_name"
done

# Merge all filtered counts into final matrix
csv_files=($COUNTS_DIR/*_counts_filtered.csv)

# Get list of all unique viruses (header row)
cut -d, -f1 $COUNTS_DIR/*_counts_filtered.csv | sort | uniq > virus_list.txt

# Start creating the final merged file
echo -n "virus" > virus_counts_filtered_matrix.csv

# Add sample names as header
for file in $COUNTS_DIR/*_counts_filtered.csv; do
    sample_name=$(basename "$file" _counts_filtered.csv)
    echo -n ",$sample_name" >> virus_counts_filtered_matrix.csv
done
echo "" >> virus_counts_filtered_matrix.csv  # New line after header

# Loop through each virus and extract its counts from all samples
while read -r virus; do
    echo -n "$virus" >> virus_counts_filtered_matrix.csv
    for file in $COUNTS_DIR/*_counts_filtered.csv; do
        count=$(grep "^$virus," "$file" | cut -d, -f2)
        echo -n ",${count:-0}" >> virus_counts_filtered_matrix.csv  # Default to 0 if missing
    done
    echo "" >> virus_counts_filtered_matrix.csv  # New line for next virus
done < virus_list.txt


echo "🎉 Final filtered counts matrix: virus_counts_filtered_matrix.csv"
