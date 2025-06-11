#!/bin/bash

# Paths to directories and Bowtie index
DATA_DIR="./data_2015_16"
TRIMMED_DIR="./trimmed_reads"
MAPPING_DIR="../../mapping_results"
COUNTS_DIR="../../counts"
BOWTIE_INDEX="../virus_db"

# Create output directories if they don't exist
mkdir -p $TRIMMED_DIR $MAPPING_DIR $COUNTS_DIR

echo "🚀 Starting trimming with fastp..."
# Loop through each sample and run fastp
for sample_dir in $DATA_DIR/*; do
    sample_name=$(basename $sample_dir)
    
    R1=$(ls $sample_dir/*R1*.fastq.gz)
    R2=$(ls $sample_dir/*R2*.fastq.gz)

    fastp --in1 $R1 --in2 $R2 \
          --out1 $TRIMMED_DIR/${sample_name}_r1_trim.fq.gz \
          --out2 $TRIMMED_DIR/${sample_name}_r2_trim.fq.gz \
          --report_title "${sample_name}_trim_report" \
          --html ${sample_name}_report.html
    
    echo "✅ Trimmed $sample_name"
done

echo "🏹 Starting Bowtie2 mapping..."
# Loop through trimmed reads and map with Bowtie2
for r1 in $TRIMMED_DIR/*_r1_trim.fq.gz; do
    sample_name=$(basename $r1 _r1_trim.fq.gz)
    r2="$TRIMMED_DIR/${sample_name}_r2_trim.fq.gz"

    bowtie2 -x $BOWTIE_INDEX \
            -1 $r1 \
            -2 $r2 \
            --threads 4 \
            -S $MAPPING_DIR/${sample_name}_map.sam

    echo "✅ Mapped $sample_name"
done

echo "📈 Extracting read counts..."
# Extract read counts for each sample and save as CSV
for sam in $MAPPING_DIR/*.sam; do
    sample_name=$(basename $sam _map.sam)

    samtools view -bS $sam | samtools sort -o $MAPPING_DIR/${sample_name}.sorted.bam
    samtools index $MAPPING_DIR/${sample_name}.sorted.bam

    samtools idxstats $MAPPING_DIR/${sample_name}.sorted.bam | cut -f 1,3 | \
    awk -v sample="$sample_name" 'BEGIN{print "virus,"sample} {print $1","$2}' \
    > $COUNTS_DIR/${sample_name}_counts.csv

    echo "✅ Extracted counts for $sample_name"
done

echo "🛠️ Merging all counts into final table..."
# Merge all individual counts into a single matrix
csv_files=($COUNTS_DIR/*_counts.csv)
paste -d, "${csv_files[@]}" | \
awk -F, '{printf "%s",$1; for (i=2; i<=NF; i+=2) printf ",%s", $i; print ""}' \
> virus_counts_matrix.csv

echo "🎉 Pipeline complete! Final counts table: virus_counts_matrix.csv"
