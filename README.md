RNA-seq Pipeline using a reference
==================================
This Repository contains a walk-through for analyzing pooled (single-end) RNA-seq performed on Drosophila melanogaster whole embryos with and without specific defects via 
- mapping to a reference transcriptome using Salmon (slurm/bash scripts calling the relevant functions to preprocess and map)

- performing differential gene expression (DE) analysis of single stage and time-series RNA-seq data that also touches upon p-value combination for multiple time-points (R workflow using Bioconductor packages)

- custom plotting including base R,ggplot2 and pheatmap

## Demultiplexing using Fastq-multx:
NGS samples may be multiplexed prior to sequencing and depending on the sequencing core/ service (academic or commercial), you may get back the raw reads per sample or you may have to demultiplex this yourself. There are several tools for demultiplexing for example Barcode splitter included in the FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage)

But for this (and other) project(s) I have used fastq-multx (Erik Aronesty (2013). TOBioiJ : "Comparison of Sequencing Utility Programs", DOI:10.2174/1875036201307010001) and can be obtained installed via (bio)conda. And for a small pilot dataset I haven't seen a qualitative difference between the two programs in terms of output and final result.

```
# source activate ngs
# conda install fastq-multx

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
# source activate ngs
# conda install trim-galore

# For single end reads
trim_galore filename.fastq.gz -o folder_name

# For paired end reads
trim_galore --paired filename_Read1.fastq.gz filename_Read2.fastq.gz -o folder_name

# TrimGalore automatically prints out filename within the output folder based on the input filename

```
## Mapping and quantification
The main objective of this project was to perform a DE analysis between different conditions for known genes and so I used Salmon (mostly default options) to perform quantification that can then be used downstream.

```
# source activate ngs
# conda install salmon

# For using Salmon we first need the transcriptome (r6.19 from flybase) and index it

wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.19_FB2017_06/fasta/dmel-all-transcript-r6.19.fasta.gz

salmon index -t dmel-all-transcript-r6.19.fasta.gz -i dmel-all-transcript-r6.19_salmon.index # really fast and can be done on the command line directly

# Running the actual quantification 
salmon quant -i dmel-all-transcript-r6.19_salmon.index -l A -r filename_trimmed.fq.gz -p 6 -o filename_salmon_quant

# -l A: Automatically detect library type 
# -r: Unmated or single-end reads
# -p 6: multi-threading on 6 cores

```
If you have used ERCC spike-ins then you can obtain the spike-in sequences (https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip) and add it to the transcriptome prior to indexing.

```
cat dmel-all-transcript-r6.19.fasta ERCC92.fasta > dmel-all-transcript-r6.19_ERCC.fasta

```
## Differential Gene Expression Analysis
Following is a typical DESeq2 differential expression analysis pipeline for bulk RNA-seq adapted from http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

```
# Install the required libraries from Bioconductor or relevant repositories and load them
# For example 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

library(DESeq2)
library(tximport)
library(apeglm) 
library(vsn)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(dendsort)
library(sva)
library(magrittr)
library(matrixStats)
library(metaRNASeq)
library(reshape2)
library(goseq)
library(GO.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)

# Read in the gene conversion file 
ttg <- read.table("~/Desktop/RNAseq_output/Abo_Wt_StagedNC_flybase_transcriptome/scripts_samples_ttg/flybase_transcript_to_gene.txt", h=T)

```
Transcript to gene names/ ID's can be added either during mapping or count matrix generation. Alternatively, the mapping can be obtained from the any of the specific model system databases. Here the mappings have been obtained directly from Flybase and used to generate gene-transcript mapping via tximport during DESeq2 SummerizedExperiment creation.

```
sample.nc14 <- read.csv("~/Desktop/samples_cell_cycle_14.csv", header = TRUE)
str(sample.nc14)

# Replicate, Lane and Sequencing batch number treated as factors
sample.nc14$Replicate <- as.factor(sample.nc14$Replicate)  
sample.nc14$Lane <- as.factor(sample.nc14$Lane)
sample.nc14$Batch <- as.factor(sample.nc14$Batch)

# Check the reference levels and set the base level as WT
# I prefer all comparisons to the control or wild-type
levels(sample.nc14$Genotype)
sample.nc14$Genotype <- relevel(sample.nc14$Genotype, ref="WT")


# Create a tximport object
files.nc14 <- file.path("~/Desktop/desktop_01_19/RNAseq_output/Abo_Wt_StagedNC_flybase_transcriptome/data/salmon_data_by_nc", "NC14", sample.nc14$Sample, "quant.sf")		
names(files.nc14) <- paste0("Sample", 1:24)
all(file.exists(files.nc14))
txi.nc14 <- tximport(files.nc14, type = "salmon", tx2gene = ttg)
names(txi.nc14)
rownames(sample.nc14) <- colnames(txi.nc14$counts) # if you have the sample names as row names in the samples.csv then this converts the row names to column names which holds all of the counts


# Create a DESeq object from tximport object
# Specify the design matrix appropriately to include Main effects and interactions where appropriate
ddsTxi.nc14 <- DESeqDataSetFromTximport(txi.nc14, colData = sample.nc14, design = ~ Genotype)

# There are technical replicates i.e. the same library sequenced multiple times and can be collapsed by summing up the reads
# Before collapsing, check whether the technical replicates are different by performing a PCA and DE analysis between them
# In this case the technical replicates are almost identical and hence collapsed
ddsColl.nc14 <- collapseReplicates(ddsTxi.nc14, ddsTxi.nc14$Sample_by_Batch_Lane, ddsTxi.nc14$Sample_by_Batch_Lane)
colData(ddsColl.nc14)
keep.nc14 <- rowSums(counts(ddsColl.nc14) >1) >= 2 # this keeps genes with counts > 1 in at least 2 of the 6 samples
ddsColl.nc14 <- ddsColl.nc14[keep.nc14,]
ddsColl.nc14 # 7993 genes


# Run the DESeq2 normalization and analysis
deseq.nc14 <- DESeq(ddsColl.nc14)
hist(normalizationFactors(deseq.nc14)) # This should be centered around 1
resultsNames(deseq.nc14)

# Extract normalized counts and save it as a dataframe
counts.nc14 <- as.data.frame(counts(deseq.nc14, normalized=T))
counts.nc14 <- cbind(Gene = rownames(counts.nc14), counts.nc14)


# Extract the appropriate results and perform exploratory analysis
res.abo.nc14 <- results(deseq.nc14,name="Genotype_Mutant_vs_WT",alpha=0.05) # alpha value should be set at adjusted p-value threshold that you call significance at
plotMA(res.abo.nc14, ylim=c(-4,4)) # MA plot

res.abo.nc14.ape <- lfcShrink(deseq.nc14,coef=2,res= res.abo.nc14, type="apeglm") # log fold change shrinkage using apeglm 
df.res.abo.nc14.ape <- as.data.frame(res.abo.nc14.ape)
df.res.abo.nc14.ape <- cbind(Gene = rownames(df.res.abo.nc14.ape),df.res.abo.nc14.ape)
plotMA(res.abo.nc14.ape, ylim=c(-4,4))


res.abo.nc14.ord <- res.abo.nc14.ape[order(res.abo.nc14.ape $padj),]
df.nc14 <- as.data.frame(res.abo.nc14.ord)
df.nc14 <- cbind(Gene = rownames(df.nc14),df.nc14)

dim(df.nc14) # total number of genes
length(which(df.nc14$log2FoldChange>0)) # genes up-regulated
length(which(df.nc14$log2FoldChange<0)) # genes down-regulated


# Filter based on an FDR adjusted p-value threshold < 0.05
res.abo.nc14.sig1 <- subset(res.abo.nc14.ord, padj < 0.05)
df.nc14.sig1 <- as.data.frame(res.abo.nc14.sig1)
df.nc14.sig1 <- cbind(Gene = rownames(df.nc14.sig1),df.nc14.sig1)

dim(df.nc14.sig1)
length(which(df.nc14.sig1$log2FoldChange>0)) 
length(which(df.nc14.sig1$log2FoldChange<0))
df.nc14.sig1.up <- subset(df.nc14.sig1, log2FoldChange>0)
df.nc14.sig1.dn <- subset(df.nc14.sig1, log2FoldChange<0)

# Filter based on a log2 fold-change threshold of 1 i.e. a 2 fold change
res.abo.nc14.l2f <- subset(res.abo.nc14.ord, abs(log2FoldChange) > 1)
df.nc14.l2f <- as.data.frame(res.abo.nc14.l2f)
df.nc14.l2f <- cbind(Gene = rownames(df.nc14.l2f),df.nc14.l2f)

dim(df.nc14.l2f) 
length(which(df.nc14.l2f$log2FoldChange>1)) 
length(which(df.nc14.l2f$log2FoldChange<-1))




