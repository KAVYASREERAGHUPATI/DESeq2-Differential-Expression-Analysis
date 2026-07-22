# ============================================================
# SAMPLE-TO-SAMPLE DISTANCE HEATMAP
# ============================================================

# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install(
    "DESeq2",
    ask = FALSE,
    update = FALSE
  )
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}


# Load packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})


# ============================================================
# USER SETTINGS
# ============================================================

project_folder <- "C:/Users/YourName/Desktop/DESeq2_Project"

results_folder <- file.path(
  project_folder,
  "results"
)

dds_file <- file.path(
  results_folder,
  "DESeq2_dataset.rds"
)

output_file <- file.path(
  results_folder,
  "Sample_Distance_Heatmap.png"
)

condition_column <- "Condition"


# ============================================================
# LOAD DESEQ2 DATASET
# ============================================================

if (!file.exists(dds_file)) {
  stop("DESeq2 dataset was not found:\n", dds_file)
}

dds <- readRDS(
  dds_file
)


# ============================================================
# VARIANCE-STABILIZING TRANSFORMATION
# ============================================================

vsd <- vst(
  dds,
  blind = FALSE
)


# ============================================================
# CALCULATE SAMPLE DISTANCES
# ============================================================

sample_distances <- dist(
  t(
    assay(vsd)
  )
)

sample_distance_matrix <- as.matrix(
  sample_distances
)


# ============================================================
# CREATE SAMPLE LABELS
# ============================================================

sample_labels <- paste(
  colnames(vsd),
  colData(vsd)[[condition_column]],
  sep = " | "
)

rownames(sample_distance_matrix) <- sample_labels
colnames(sample_distance_matrix) <- sample_labels


# ============================================================
# CREATE SAMPLE ANNOTATION
# ============================================================

annotation_data <- data.frame(
  Condition = colData(vsd)[[condition_column]]
)

rownames(annotation_data) <- sample_labels


# ============================================================
# CREATE AND SAVE HEATMAP
# ============================================================

png(
  filename = output_file,
  width = 2800,
  height = 2500,
  res = 300
)

pheatmap(
  sample_distance_matrix,

  annotation_col = annotation_data,

  annotation_row = annotation_data,

  clustering_distance_rows = sample_distances,

  clustering_distance_cols = sample_distances,

  cluster_rows = TRUE,

  cluster_cols = TRUE,

  show_rownames = TRUE,

  show_colnames = TRUE,

  fontsize_row = 8,

  fontsize_col = 8,

  border_color = NA,

  main = "Sample-to-Sample Distance Heatmap"
)

dev.off()

message(
  "Sample-distance heatmap saved:\n",
  output_file
)
