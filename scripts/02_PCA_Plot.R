# ============================================================
# PCA PLOT
# ============================================================

# Install packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}


# Load packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
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
  "PCA.png"
)

condition_column <- "Condition"


# ============================================================
# LOAD DESEQ2 DATASET
# ============================================================

if (!file.exists(dds_file)) {
  stop("DESeq2 dataset was not found:\n", dds_file)
}

dds <- readRDS(dds_file)


# ============================================================
# VARIANCE-STABILIZING TRANSFORMATION
# ============================================================

vsd <- vst(
  dds,
  blind = FALSE
)


# ============================================================
# EXTRACT PCA INFORMATION
# ============================================================

pca_data <- plotPCA(
  vsd,
  intgroup = condition_column,
  returnData = TRUE,
  ntop = min(500, nrow(vsd))
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar"),
  digits = 1
)


# ============================================================
# CREATE PCA PLOT
# ============================================================

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = .data[[condition_column]],
    shape = .data[[condition_column]]
  )
) +
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  xlab(
    paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    )
  ) +
  labs(
    title = "Principal Component Analysis",
    color = "Condition",
    shape = "Condition"
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


# ============================================================
# SAVE PCA PLOT
# ============================================================

ggsave(
  filename = output_file,
  plot = pca_plot,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("PCA plot saved:\n", output_file)
