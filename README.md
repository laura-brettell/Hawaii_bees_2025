# Hawaii_bees_2025
This is a project looking at virome composition, deformed wing virus prevalence, diversity and recombination 15 years after varroa became established on two Hawaiian islands.

This repository contains the scripts to analyse the data and generate the figures. The data comprises samples collected in 2024/2025 which have used long read sequencing (Oxford Nanopore) and samples collected between 2009 and 2016, which have used short read sequencing (Illumina).





## generating summaries

### data preparation

The 2024-25 data were initially processed in Epi2me - The fastq files for each sample were mapped against a custom database using minimap2 to generate abundance tables (tsv) and alignment files (bam and bam.bai).

The 2009-26 data were obtained as fastq format and trimmed with fastp, followed by mapping to a custom database (virus names and accessions in paper supp Table S2) using bowtie2 using 'virus_pipeline_map.sh'. Then, low abundance reads were filtered using 'filter_low_abundance.sh'. The output 'virus_counts_filtered_matrix.csv' was then manually edited to:

-  make the accessions the same format as the long read csv
-  add a column for virus names
-  add three adddiitonal rows containig total_virus (sum of all virus reads/sample), 'total_trim_reads' (the total number of trimmed reads for each sample - from the fastp reports), and  'total_pair' (total trimmed, paired reads, calculated as 0.5* total_trim_reads)
  
and save as 'virus_counts_filtered_matrix_edit_3_No_MiMo.csv'.


### Additional data filtering and generation of plots

To organise and filter the data and generate the virome composition plots, use the R script 'virome_comp_all_years.Rmd'. The required abundance tables and metadata files can be found in thhe folder 'virome_comp_files'.


## DWV recombination

### generate DWV assemblies

Firstly, the 2024-35 data.

1. Extract the DWV reads from bam file generated using minimap2 in Epi2me, eg:

```
/usr/bin/samtools view -b ./bams_BI/Garnets_Hilo1.reference.bam "NC_004830.2" "NC_006494.1" > ./DWV_reads/GA1_DWV_reads.bam

```

2. convert bam files to fasta eg

```
/usr/bin/samtools fasta ./DWV_reads/DE1_DWV_reads.bam > ./DWV_reads/DE1_DWV_reads.fasta

```

3. assemble DWV reads with Canu

```
mkdir DWV_assemblies

conda activate canu-env

canu -p DE2_DWV  \
-d ./  \
-nanopore  \
maxInputCoverage=2000  \
corOutCoverage=all  \
corMinCoverage=5  
corMhapSensitivity=high  \
minOverlapLength=50  \
minReadLength=100  \
genomesize=10000  \  
useGrid=false  /
../../DWV_reads/DE2_DWV_reads.fasta

```

Next, the long read data.....












