# DESeq2-Differential-Expression-Analysis
Differential expression analysis of RNA seq. data using DESeq2 in R with Salmon gene count estimates, including data preprocessing, statistical analysis, visualization and result interpretation.
.............................................................................................................................................................................................
## Overview
This repository provides a complete workflow for bulk RNA-seq differential expression analysis using the DESeq2 package in R. The workflow starts from Salmon quantification outputs (quant.sf), which are imported into R using tximport to generate gene-level count estimates for differential expression analysis.

The Salmon quantification outputs can be generated using nf-core/rnaseq, Galaxy RNA-seq workflows, standalone Salmon or any other RNA-seq analysis pipeline that produces standard Salmon output files.

This repository demonstrates the downstream analysis required to identify differentially expressed genes (DEGs), visualize transcriptomic differences between experimental conditions and generate publication-quality figures.
.............................................................................................................................................................................................
## Input Files
1. Salmon Quantification Outputs
   Each RNA-seq sample should contain a Salmon quantification file: quant.sf
   These files may be generated using: nf-core/rnaseq, Galaxy RNA-seq workflows, Standalone Salmon or Any RNA-seq analysis pipeline that produces standard Salmon quantification outputs
2. Transcript-to-Gene Mapping File
   A transcript-to-gene mapping file (tx2gene) is required by tximport to convert transcript-level abundance estimates into gene-level count estimates suitable for DESeq2. This mapping file     is typically generated from the reference genome annotation (GTF/GFF3) and contains two columns: Transcript ID and the corresponding Gene ID.
3. Sample Metadata
   The sample metadata file should contain the sample names and their corresponding experimental conditions.
   Example:
     Sample       Condition
     Sample_01  	Control
     Sample_02   	Control
     Sample_03  	Drought
.............................................................................................................................................................................................
  ## Analysis Workflow

The workflow consists of the following steps:

1.Import Salmon quantification outputs using tximport.
2.Convert transcript-level abundance estimates into gene-level count estimates.
  (The tximport package reads the quant.sf files from all samples and combines them into a single object. A transcript-to-gene mapping table called tx2gene, is also supplied.)
3.Import sample metadata.
4.Create the DESeq2 dataset.
5.Estimate size factors for normalization.
6.Estimate gene-wise dispersions.
7.Fit the negative binomial model.
  (produces a log2 fold change).
8.Perform differential expression analysis.
9.Adjust p-values using the Benjamini–Hochberg method.
10.Identify significantly upregulated and downregulated genes.
11.Export result tables and publication-quality figures.
.............................................................................................................................................................................................
## Important packages

| Package             | Purpose                                                                                                            |
| ------------------- | ------------------------------------------------------------------------------------------------------------------ |
| **tximport**        | Imports Salmon `quant.sf` files and converts transcript-level abundance estimates into gene-level count estimates. |
| **DESeq2**          | Performs normalization, dispersion estimation, statistical testing, and differential expression analysis.          |
| **readr**           | Reads input files such as sample metadata and transcript-to-gene (`tx2gene`) mapping files.                        |
| **GenomicFeatures** | Generates the transcript-to-gene (`tx2gene`) mapping from the reference GTF annotation file.                       |
| **ggplot2**         | Creates PCA plots, MA plots, and other publication-quality visualizations.                                         |
| **EnhancedVolcano** | Generates volcano plots of differentially expressed genes.                                                         |
| **pheatmap**        | Creates heatmaps of normalized gene expression.                                                                    |
| **RColorBrewer**    | Provides color palettes for heatmaps and other plots.                                                              |
| **dplyr**           | Data manipulation and filtering.                                                                                   |
| **tibble**          | Creates and manages modern data frames.                                                                            |
| **tidyr**           | Reshapes and tidies data.                                                                                          |
| **stringr**         | Handles string operations such as sample name manipulation.                                                        |

.............................................................................................................................................................................................
## Output files
DESeq2_All_Results.csv
DESeq2_Significant_DEGs.csv
DESeq2_Upregulated_Genes.csv
DESeq2_Downregulated_Genes.csv
DESeq2_Normalized_Counts.csv

PCA.png
Volcano_Plot.png
MA_Plot.png
Top_DEGs_Heatmap.png
Sample_Distance_Heatmap.png
.............................................................................................................................................................................................
> **Note:** DESeq2 can perform differential expression analysis using gene-level count matrices generated by various RNA-seq quantification methods, including STAR gene counts-featureCounts, HTSeq-count and Salmon. This repository uses **Salmon quantification outputs together with tximport**, which summarizes transcript-level abundance estimates into gene-level count estimates for DESeq2. Salmon-based quantification is widely adopted in modern RNA-seq workflows because it provides transcript-aware abundance estimation, particularly for genes with multiple transcript isoforms. However, STAR-derived gene counts are also a valid and widely accepted input for DESeq2.
