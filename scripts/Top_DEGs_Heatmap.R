# ============================================================
# TOP SIGNIFICANT DEG HEATMAP 
# ============================================================

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

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}


suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
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

results_file <- file.path(
  results_folder,
  "DESeq2_results_object.rds"
)

output_file <- file.path(
  results_folder,
  "Top_DEGs_Heatmap_Red_Blue.png"
)

condition_column <- "Condition"

padj_cutoff <- 0.05
log2fc_cutoff <- 1

number_of_top_genes <- 50


# ============================================================
# LOAD DATA
# ============================================================

if (!file.exists(dds_file)) {
  stop("DESeq2 dataset was not found:\n", dds_file)
}

if (!file.exists(results_file)) {
  stop("DESeq2 results object was not found:\n", results_file)
}

dds <- readRDS(dds_file)

results_deseq2 <- readRDS(results_file)


# ============================================================
# VARIANCE-STABILIZING TRANSFORMATION
# ============================================================

vsd <- vst(
  dds,
  blind = FALSE
)


# ============================================================
# SELECT SIGNIFICANT GENES
# ============================================================

results_table <- as.data.frame(results_deseq2)

results_table$Gene_ID <- rownames(results_table)

significant_genes <- results_table[
  !is.na(results_table$padj) &
    results_table$padj < padj_cutoff &
    abs(results_table$log2FoldChange) > log2fc_cutoff,
]


significant_genes <- significant_genes[
  order(
    significant_genes$padj,
    -abs(significant_genes$log2FoldChange)
  ),
]


if (nrow(significant_genes) == 0) {
  stop(
    "No significant genes were found using the selected thresholds."
  )
}


selected_genes <- head(
  significant_genes$Gene_ID,
  min(
    number_of_top_genes,
    nrow(significant_genes)
  )
)


# ============================================================
# PREPARE HEATMAP MATRIX
# ============================================================

heatmap_matrix <- assay(vsd)[
  selected_genes,
  ,
  drop = FALSE
]


# Convert each gene to a row-wise Z-score
heatmap_matrix <- t(
  scale(
    t(heatmap_matrix)
  )
)


heatmap_matrix[
  is.na(heatmap_matrix)
] <- 0


# ============================================================
# SAMPLE ANNOTATION
# ============================================================

annotation_col <- data.frame(
  Condition = colData(dds)[[condition_column]]
)

rownames(annotation_col) <- colnames(dds)


# Optional annotation colors
annotation_colors <- list(
  Condition = c(
    Control = "#2166AC",
    Drought = "#B2182B"
  )
)


# ============================================================
# BLUE–WHITE–RED COLOR PALETTE
# ============================================================

heatmap_colors <- colorRampPalette(
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)(100)


# Symmetrical breaks keep zero at white
maximum_absolute_value <- max(
  abs(heatmap_matrix),
  na.rm = TRUE
)

heatmap_breaks <- seq(
  -maximum_absolute_value,
  maximum_absolute_value,
  length.out = 101
)


# ============================================================
# CREATE AND SAVE HEATMAP
# ============================================================

png(
  filename = output_file,
  width = 2800,
  height = 3000,
  res = 300
)

pheatmap(
  heatmap_matrix,

  color = heatmap_colors,

  breaks = heatmap_breaks,

  annotation_col = annotation_col,

  annotation_colors = annotation_colors,

  cluster_rows = TRUE,

  cluster_cols = TRUE,

  show_rownames = TRUE,

  show_colnames = TRUE,

  fontsize_row = 6,

  fontsize_col = 8,

  scale = "none",

  border_color = NA,

  main = paste(
    "Top",
    nrow(heatmap_matrix),
    "Significant Differentially Expressed Genes"
  )
)

dev.off()


message(
  "Red-blue DEG heatmap saved:\n",
  output_file
)
