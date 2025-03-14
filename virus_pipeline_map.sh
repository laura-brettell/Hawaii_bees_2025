#!/bin/bash

# run from Hawaii_old
# Paths to directories and Bowtie index
DATA_DIR="./raw_data/data_2015_16"
TRIMMED_DIR="./raw_data/trimmed_reads"
MAPPING_DIR="./mapping_results"
COUNTS_DIR="./mapping_results/counts"
#BOWTIE_INDEX="../virus_db"

# Create output directories if they don't exist
#mkdir -p $TRIMMED_DIR $MAPPING_DIR $COUNTS_DIR

mkdir -p $COUNTS_DIR



echo "🏹 Starting Bowtie2 mapping..."
# Loop through trimmed reads and map with Bowtie2
for r1 in $TRIMMED_DIR/*_r1_trim.fq.gz; do
    sample_name=$(basename $r1 _r1_trim.fq.gz)
    r2="$TRIMMED_DIR/${sample_name}_r2_trim.fq.gz"

    /usr/bin/bowtie -x ./virus_db/virus_db \
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
