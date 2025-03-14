#!/bin/bash

# Path to your data directory
DATA_DIR="./data_2012"
OUTPUT_DIR="./trimmed_reads"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through each sample directory
for sample_dir in $DATA_DIR/*; do
    # Extract sample name (e.g. HB_S7, HB_S12, etc.)
    sample_name=$(basename $sample_dir)

    # Find R1 and R2 files
    R1=$(ls $sample_dir/*R1*.fastq.gz)
    R2=$(ls $sample_dir/*R2*.fastq.gz)

    # Run fastp for each sample
    fastp \
        --in1 $R1 \
        --in2 $R2 \
        --out1 $OUTPUT_DIR/${sample_name}_r1_trim.fq.gz \
        --out2 $OUTPUT_DIR/${sample_name}_r2_trim.fq.gz \
        --report_title "${sample_name}_trim_report" \
        --html ${sample_name}_report.html

    echo "Finished trimming $sample_name"
done
