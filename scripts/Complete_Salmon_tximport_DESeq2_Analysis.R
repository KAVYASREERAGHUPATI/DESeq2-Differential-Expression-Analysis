# ============================================================
# COMPLETE SALMON + GTF + TXIMPORT + DESEQ2 ANALYSIS
# ============================================================
# This script:
# 1. Installs required packages
# 2. Imports sample metadata
# 3. Generates tx2gene mapping directly from the GTF file
# 4. Imports Salmon quant.sf files using tximport
# 5. Creates the DESeq2 dataset
# 6. Filters low-count genes
# 7. Runs differential-expression analysis
# 8. Identifies significant, upregulated and downregulated genes
# 9. Exports all result tables
#
# Plotting is intentionally not included.
# ============================================================


# ============================================================
# 1. INSTALL REQUIRED PACKAGES
# ============================================================

cran_packages <- c(
  "readr",
  "dplyr",
  "tibble"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioconductor_packages <- c(
  "tximport",
  "DESeq2",
  "rtracklayer"
)

for (pkg in bioconductor_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tximport)
  library(DESeq2)
  library(rtracklayer)
})


# ============================================================
# 3. USER SETTINGS
# ============================================================

# Change this path to your main project folder.

project_folder <- "C:/Users/YourName/Desktop/DESeq2_Project"


# Sample metadata file
metadata_file <- file.path(
  project_folder,
  "data",
  "metadata",
  "sample_metadata.csv"
)


# Reference GTF annotation file
gtf_file <- file.path(
  project_folder,
  "data",
  "reference",
  "Oryza_sativa.IRGSP-1.0.62.gtf.gz"
)


# Main Salmon folder
salmon_folder <- file.path(
  project_folder,
  "data",
  "salmon"
)


# Output folder
output_folder <- file.path(
  project_folder,
  "results"
)

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


# Metadata column names
sample_column <- "Sample"
condition_column <- "Condition"


# Experimental comparison
reference_condition <- "Control"
treatment_condition <- "Drought"


# DEG thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1


# Low-count filtering
minimum_count <- 10
minimum_samples <- 2


# ============================================================
# 4. CHECK INPUT FILES AND FOLDERS
# ============================================================

if (!file.exists(metadata_file)) {
  stop(
    "Metadata file was not found:\n",
    metadata_file
  )
}

if (!file.exists(gtf_file)) {
  stop(
    "GTF annotation file was not found:\n",
    gtf_file
  )
}

if (!dir.exists(salmon_folder)) {
  stop(
    "Salmon folder was not found:\n",
    salmon_folder
  )
}


# ============================================================
# 5. IMPORT SAMPLE METADATA
# ============================================================

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)


required_metadata_columns <- c(
  sample_column,
  condition_column
)


missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)


if (length(missing_metadata_columns) > 0) {
  stop(
    "The metadata file is missing the following column(s): ",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}


metadata[[sample_column]] <- as.character(
  metadata[[sample_column]]
)

metadata[[condition_column]] <- as.character(
  metadata[[condition_column]]
)


# Check duplicate sample names
if (anyDuplicated(metadata[[sample_column]]) > 0) {
  duplicated_samples <- unique(
    metadata[[sample_column]][
      duplicated(metadata[[sample_column]])
    ]
  )

  stop(
    "Duplicate sample names were found: ",
    paste(
      duplicated_samples,
      collapse = ", "
    )
  )
}


# Check whether both conditions are present
available_conditions <- unique(
  metadata[[condition_column]]
)


if (!reference_condition %in% available_conditions) {
  stop(
    "Reference condition was not found in metadata: ",
    reference_condition
  )
}


if (!treatment_condition %in% available_conditions) {
  stop(
    "Treatment condition was not found in metadata: ",
    treatment_condition
  )
}


# ============================================================
# 6. GENERATE tx2gene MAPPING FROM THE GTF FILE
# ============================================================

message("Importing the GTF annotation file...")


gtf_annotation <- rtracklayer::import(
  gtf_file
)


gtf_data <- as.data.frame(
  gtf_annotation
)


required_gtf_columns <- c(
  "transcript_id",
  "gene_id"
)


missing_gtf_columns <- setdiff(
  required_gtf_columns,
  colnames(gtf_data)
)


if (length(missing_gtf_columns) > 0) {
  stop(
    "The GTF file does not contain the following required attributes: ",
    paste(
      missing_gtf_columns,
      collapse = ", "
    )
  )
}


# Prefer transcript rows when present
if (
  "type" %in% colnames(gtf_data) &&
  any(gtf_data$type == "transcript", na.rm = TRUE)
) {

  tx2gene <- gtf_data |>
    filter(
      type == "transcript"
    ) |>
    transmute(
      TXNAME = as.character(transcript_id),
      GENEID = as.character(gene_id)
    )

} else {

  # Use all rows containing transcript_id and gene_id
  tx2gene <- gtf_data |>
    transmute(
      TXNAME = as.character(transcript_id),
      GENEID = as.character(gene_id)
    )
}


# Remove missing and duplicate mappings
tx2gene <- tx2gene |>
  filter(
    !is.na(TXNAME),
    !is.na(GENEID),
    TXNAME != "",
    GENEID != ""
  ) |>
  distinct()


if (nrow(tx2gene) == 0) {
  stop(
    "No transcript-to-gene mappings could be extracted from the GTF file."
  )
}


message(
  "Number of transcript-to-gene mappings generated: ",
  nrow(tx2gene)
)


# Save the generated tx2gene file
write_tsv(
  tx2gene,
  file.path(
    output_folder,
    "Generated_tx2gene_from_GTF.tsv"
  )
)


# ============================================================
# 7. LOCATE SALMON quant.sf FILES
# ============================================================

# Expected structure:
#
# data/salmon/Sample_01/quant.sf
# data/salmon/Sample_02/quant.sf
# data/salmon/Sample_03/quant.sf


quant_files <- file.path(
  salmon_folder,
  metadata[[sample_column]],
  "quant.sf"
)


names(quant_files) <- metadata[[sample_column]]


missing_quant_files <- quant_files[
  !file.exists(quant_files)
]


if (length(missing_quant_files) > 0) {
  stop(
    "The following Salmon quant.sf files were not found:\n",
    paste(
      missing_quant_files,
      collapse = "\n"
    )
  )
}


# ============================================================
# 8. IMPORT SALMON FILES USING TXIMPORT
# ============================================================

message("Importing Salmon quantification files...")


txi <- tximport(
  files = quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  txOut = FALSE,
  countsFromAbundance = "no",
  ignoreTxVersion = TRUE
)


message(
  "Number of genes imported: ",
  nrow(txi$counts)
)

message(
  "Number of samples imported: ",
  ncol(txi$counts)
)


# ============================================================
# 9. PREPARE METADATA FOR DESEQ2
# ============================================================

metadata <- as.data.frame(
  metadata
)


rownames(metadata) <- metadata[[sample_column]]


# Confirm sample order
if (!all(
  colnames(txi$counts) == rownames(metadata)
)) {
  stop(
    "The sample order in tximport does not match the metadata."
  )
}


# Convert condition column to factor
metadata[[condition_column]] <- factor(
  metadata[[condition_column]],
  levels = c(
    reference_condition,
    treatment_condition
  )
)


# ============================================================
# 10. CREATE DESEQ2 DATASET
# ============================================================

design_formula <- as.formula(
  paste(
    "~",
    condition_column
  )
)


dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = metadata,
  design = design_formula
)


# ============================================================
# 11. FILTER VERY LOW-COUNT GENES
# ============================================================

genes_before_filtering <- nrow(
  dds
)


keep_genes <- rowSums(
  counts(dds) >= minimum_count
) >= minimum_samples


dds <- dds[
  keep_genes,
]


genes_after_filtering <- nrow(
  dds
)


message(
  "Genes before filtering: ",
  genes_before_filtering
)

message(
  "Genes after filtering: ",
  genes_after_filtering
)


# ============================================================
# 12. RUN DESEQ2
# ============================================================

# DESeq() performs:
# 1. Size-factor estimation
# 2. Dispersion estimation
# 3. Negative-binomial model fitting
# 4. Statistical testing


dds <- DESeq(
  dds
)


# ============================================================
# 13. EXTRACT DIFFERENTIAL-EXPRESSION RESULTS
# ============================================================

results_deseq2 <- results(
  dds,
  contrast = c(
    condition_column,
    treatment_condition,
    reference_condition
  ),
  alpha = padj_cutoff
)


# Sort by adjusted p-value
results_deseq2 <- results_deseq2[
  order(
    results_deseq2$padj,
    na.last = TRUE
  ),
]


# ============================================================
# 14. CREATE COMPLETE RESULT TABLE
# ============================================================

all_results <- as.data.frame(
  results_deseq2
) |>
  rownames_to_column(
    "Gene_ID"
  ) |>
  mutate(
    Regulation = case_when(

      !is.na(padj) &
        padj < padj_cutoff &
        log2FoldChange > log2fc_cutoff ~
        "Upregulated",

      !is.na(padj) &
        padj < padj_cutoff &
        log2FoldChange < -log2fc_cutoff ~
        "Downregulated",

      TRUE ~
        "Not_significant"
    )
  )


# ============================================================
# 15. SIGNIFICANT DIFFERENTIALLY EXPRESSED GENES
# ============================================================

significant_degs <- all_results |>
  filter(
    !is.na(padj),
    padj < padj_cutoff,
    abs(log2FoldChange) > log2fc_cutoff
  )


# ============================================================
# 16. UPREGULATED GENES
# ============================================================

upregulated_genes <- significant_degs |>
  filter(
    log2FoldChange > log2fc_cutoff
  )


# ============================================================
# 17. DOWNREGULATED GENES
# ============================================================

downregulated_genes <- significant_degs |>
  filter(
    log2FoldChange < -log2fc_cutoff
  )


# ============================================================
# 18. EXTRACT NORMALIZED COUNTS
# ============================================================

normalized_counts <- counts(
  dds,
  normalized = TRUE
) |>
  as.data.frame() |>
  rownames_to_column(
    "Gene_ID"
  )


# ============================================================
# 19. CREATE ANALYSIS SUMMARY
# ============================================================

result_summary <- tibble(
  Category = c(
    "Genes before low-count filtering",
    "Genes after low-count filtering",
    "All tested genes",
    "Significant DEGs",
    "Upregulated genes",
    "Downregulated genes"
  ),

  Count = c(
    genes_before_filtering,
    genes_after_filtering,
    nrow(all_results),
    nrow(significant_degs),
    nrow(upregulated_genes),
    nrow(downregulated_genes)
  )
)


# ============================================================
# 20. EXPORT CSV RESULT FILES
# ============================================================

write_csv(
  all_results,
  file.path(
    output_folder,
    "DESeq2_All_Results.csv"
  )
)


write_csv(
  significant_degs,
  file.path(
    output_folder,
    "DESeq2_Significant_DEGs.csv"
  )
)


write_csv(
  upregulated_genes,
  file.path(
    output_folder,
    "DESeq2_Upregulated_Genes.csv"
  )
)


write_csv(
  downregulated_genes,
  file.path(
    output_folder,
    "DESeq2_Downregulated_Genes.csv"
  )
)


write_csv(
  normalized_counts,
  file.path(
    output_folder,
    "DESeq2_Normalized_Counts.csv"
  )
)


write_csv(
  result_summary,
  file.path(
    output_folder,
    "DESeq2_Result_Summary.csv"
  )
)


# ============================================================
# 21. SAVE R OBJECTS FOR LATER PLOTTING
# ============================================================

saveRDS(
  txi,
  file.path(
    output_folder,
    "tximport_gene_level_object.rds"
  )
)


saveRDS(
  dds,
  file.path(
    output_folder,
    "DESeq2_dataset.rds"
  )
)


saveRDS(
  results_deseq2,
  file.path(
    output_folder,
    "DESeq2_results_object.rds"
  )
)


# ============================================================
# 22. SAVE SESSION INFORMATION
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    output_folder,
    "sessionInfo.txt"
  )
)


# ============================================================
# 23. FINAL SUMMARY
# ============================================================

message(
  "\n============================================"
)

message(
  "DESeq2 analysis completed successfully."
)

message(
  "Comparison: ",
  treatment_condition,
  " versus ",
  reference_condition
)

message(
  "Adjusted p-value cutoff: ",
  padj_cutoff
)

message(
  "Absolute log2 fold-change cutoff: ",
  log2fc_cutoff
)

message(
  "Output folder: ",
  output_folder
)

message(
  "============================================\n"
)


print(
  result_summary
)
