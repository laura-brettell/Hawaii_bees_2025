# Hawaii_bees_2025
This is a project looking at virome composition, deformed wing virus prevalence, diversity and recombination 15 years after varroa became established on two Hawaiian islands.

This repository contains the scripts to analyse the data and generate the figures in the paper draft. The data comprises samples collected in 2024/2025 which have used long read sequencing (Oxford Nanopore) and samples collected between 2009 and 2016, which have used short read sequencing (Illumina).





## generating summaries

### data preparation

The 2024-25 data were initially processed in Epi2me - The fastq files for each sample were mapped against a custom database using minimap2 to generate abundance tables (tsv) and alignment files (bam and bam.bai).

The 2009-26 data were obtained as fastq format and trimmed with fastp, followed by mapping to a custom database (virus names and accessions in paper supplementary materials) using bowtie2 using 'virus_pipeline_map.sh'. Then, low abundance reads were filtered using 'filter_low_abundance.sh'. The output 'virus_counts_filtered_matrix.csv' was then manually edited to:

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

Next, the older, short read data...


1. Extract DWV data in bam format (remember even though it is not explicitly called, this needs bam.bai files to also be in the folder) from the bowtie mapping eg

```
/usr/bin/samtools view -b ./mapping_results/bams/bams_2009/HB_S5.sorted.bam "NC_004830.2" "NC_006494.1" > ./DWV_reads/DWV_reads_2009/HB_S5_DWV_reads.bam
```

2. sort the bams, eg

```
/usr/bin/samtools view -b ./mapping_results/bams/bams_2009/HB_S6.sorted.bam "NC_004830.2" "NC_006494.1" > ./DWV_reads/DWV_reads_2009/HB_S6_DWV_reads.bam
```

3. Sorted bam to paired fastq, eg

```
/usr/bin/bedtools bamtofastq -i ./DWV_reads/DWV_reads_2012/HB_S17_DWV_reads_sorted.bam -fq ./DWV_reads/DWV_reads_2012/HB_S17_DWV_reads_R1.fastq -fq2 ./DWV_reads/DWV_reads_2012/HB_S17_DWV_reads_R2.fastq
```

4. assembly with SPAdes v3.12.0, eg

```
/usr/bin/spades.py -1 ./DWV_reads/DWV_reads_2012/HB_S17_DWV_reads_R1.fastq -2 ./DWV_reads/DWV_reads_2012/HB_S17_DWV_reads_R2.fastq -o ./DWV_assemblies/2012_assemblies/HB_S17_DWV_assembly
```

Then, import contigs to Geneious Prime, as well as individual reads, and inspect contigs and reads that span putative recombination breakpoints.

 




### generate coverage plots

Show coverage of reads to DWV-A and DWV-B ref genomes

Use the script 'Hawaii_dwv_recomb_all_2.Rmd' - this takes bam files from epi2me for the LR data, and those bam files for the short read data generated above. It also caluculates for each sample the dominant variant every 50bp, then when there is a switch (indicative of a recombinant) it records the position into a dataframe and indicates the points on the plots


### make recombination breakpoints diagram

Use the script 'Hawaii2025_dwv_recom_breaks_fig_edit.Rmd'to create schematics of DWV genome organisation for DWV-positive samples with full length coverage. This takes the spreadsheet 'recomb_breakpoints_3.csv' - created with the breakpoints calculated in the previous script, with additional breakpoints manually estimates and curated from the cov plots as well as geneious data. 








