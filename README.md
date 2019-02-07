RNAseq_Pipeline_Reference_Method
================================
This Repository contains scripts and a walk-through for analyzing pooled (single-end) RNA-seq data via 
- mapping to a reference transcriptome using Salmon (slurm/bash scripts calling the relevant functions to preprocess and map)

- performing differential gene expression analysis of single stage and time-series RNA-seq data that also touches upon p-value combination for multiple time-points (R workflow using Bioconductor packages)

- custom plotting including base R,ggplot2 and pheatmap

## Environment
- I like doing much of the sequencing related projects in a conda environment (https://conda.io/projects/conda/en/latest/user-guide/getting-started.html). This is really useful especially on university HPCs (and everywhere else really!) for reasons... (ok fine you can see https://medium.freecodecamp.org/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c)

### Demultiplexing using Fastq-multx:
NGS samples may be multiplexed prior to sequencing and depending on the sequencing core/ service (academic or commercial), you may get back the raw reads per sample or you may have to demultiplex this yourself. There are several tools for demultiplexing for example Barcode splitter included in the FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage)

But for this (and other) project(s) I have used fastq-multx (Erik Aronesty (2013). TOBioiJ : "Comparison of Sequencing Utility Programs", DOI:10.2174/1875036201307010001) and can be obtained installed via (bio)conda. 

```
# use fastq-multx --help to check out the usage and functionality of this tool

fastq-multx -l Barcodes_Fastq_multx.txt Index2.fastq.gz Index1.fastq.gz Read1.fastq.gz Read2.fastq.gz -o n/a -o n/a -o %.Read1.fastq.gz -o %.Read2.fastq.gz


