RNAseq_Pipeline_Reference_Method
================================
This Repository contains scripts and a walk-through for analyzing pooled (single-end) RNA-seq data via 
- mapping to a reference transcriptome using Salmon (slurm/bash scripts calling the relevant functions to preprocess and map)

- performing differential gene expression analysis of single stage and time-series RNA-seq data that also touches upon p-value combination for multiple time-points (R workflow using Bioconductor packages)

- custom plotting including base R,ggplot2 and pheatmap

## Demultiplexing using Fastq-multx:
NGS samples may be multiplexed prior to sequencing and depending on the sequencing core/ service (academic or commercial), you may get back the raw reads per sample or you may have to demultiplex this yourself. There are several tools for demultiplexing for example Barcode splitter included in the FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage)

But for this (and other) project(s) I have used fastq-multx (Erik Aronesty (2013). TOBioiJ : "Comparison of Sequencing Utility Programs", DOI:10.2174/1875036201307010001) and can be obtained installed via (bio)conda. And for a small pilot dataset I haven't seen a qualitative difference between the two programs in terms of output and final result.

```
cd /scratch/src/RNAseq/raw/

# for single-end reads the usage is something like this
fastq-multx -l Barcodes.txt Index.fastq.gz Reads.fastq.gz -o n/a -o %.Read1.fastq.gz 

# for paired-end reads the usage is something like this (dual indexed as well)
fastq-multx -l Barcodes.txt Index1.fastq.gz Index2.fastq.gz Read1.fastq.gz Read2.fastq.gz -o n/a -o n/a -o %.Read1.fastq.gz -o %.Read2.fastq.gz

# The -l: Determine barcodes from any read, using BCFIL as a master list
# Use fastq-multx --help to check out other functionalities and their usage
# Tip1: If, for dual-indexed reads, the demultiplexing looks wonky and you've done everything right, try changing the order of the Index files. 
# Tip2: In the Barcodes.txt file you can provide a unique identifier for each index so that when split, the file names reflect this unique Id. I like to give it the form "UniqueId_Rep#_Run#_Lane#"

```
## Adapter trimming
Once the samples have been split into files with a name "UniqueId_Rep_Run_Lane.fastq.gz", the step is to perform adapter trimming. I have extensively used Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic) and Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). Both of these tools are easy to use and once again for a simple case, I haven't seen a qualitative difference between the two. For this project I have used Trim Galore

```
cd /scratch/src/RNAseq/raw/

# For single end reads
trim_galore filename.fastq.gz -o /scratch/src/RNAseq/Trimd_data/

# For paired end reads
trim_galore --paired filename_Read1.fastq.gz filename_Read2.fastq.gz -o /scratch/src/RNAseq/Trimd_data/

```

