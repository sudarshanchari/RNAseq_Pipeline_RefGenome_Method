RNAseq_Pipeline_Reference_Method
================================
This Repository contains scripts for analyzing pooled RNA-seq data via 
- mapping to a reference transcriptome using Salmon (slurm/bash scripts calling the relevant functions to preprocess and map)

- performing differential gene expression analysis of single stage and time-series RNA-seq data that also touches upon p-value combination for multiple time-points (R workflow using Bioconductor packages)

- custom plotting including base R,ggplot2 and pheatmap

Demultiplexing using Fastq-multx:
---------------------------------
Typically samples of ngs are multiplexed prior to sequencing and depending on the sequencing core/ service (academic or commercial)
