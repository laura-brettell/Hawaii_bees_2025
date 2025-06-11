# Hawaii_bees_2025
This is a project looking at virome composition, deformed wing virus prevalence, diversity and recombination 15 years after varroa became established on two Hawaiian islands.

This repository contains the scripts to analyse the data and generate the figures. The data comprises samples collected in 2024/2025 which have used long read sequencing (Oxford Nanopore) and samples collected between 2009 and 2016, which have used short read sequencing (Illumina).





## generating summaries

### data preparation

The 2024-25 data were initially processed in Epi2me - The fastq files for each sample were mapped against a custom database using minimap2 to generate abundance tables (tsv) and alignment files (bam and bam.bai).

The 2009-26 data were obtained as fastq format and trimmed with fastp, followed by mapping to a custom database (virus names and accessions in paper supp Table S2) using bowtie2 using 'virus_pipeline_map.sh'. Then, low abundance reads were filtered using 'filter_low_abundance.sh'. The output 'virus_counts_filtered_matrix.csv' was then manually edited to:
-  make the accessions the same format as the long read csv
-  add acolumn for virus names
-  add three adddiitonal rows containig total_virus (sum of all virus reads/sample), 'total_trim_reads' (the total number of trimmed reads for each sample was calculated using './extract_trimmed_reads.sh' - or maybe I just got these numbers from the fastp reports actually - to check), and  'total_pair' (total trimmed, paired reads, calculated as 0.5* total_trim_reads)
-  and saved as 'virus_counts_filtered_matrix_edit_3_No_MiMo.csv'.


### Additional data filtering and generation of plots

To organise and filter the data and generate the virome composition plots, use the R script 'virome_comp_all_years.Rmd'. The required abundance tables and metadata files can be found in thhe folder 'virome_comp_files'.









