# ============================================================
# MA PLOT
# ============================================================

# Install package if missing
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


# Load package
suppressPackageStartupMessages({
  library(DESeq2)
})


# ============================================================
# USER SETTINGS
# ============================================================

project_folder <- "C:/Users/YourName/Desktop/DESeq2_Project"

results_folder <- file.path(
  project_folder,
  "results"
)

results_file <- file.path(
  results_folder,
  "DESeq2_results_object.rds"
)

output_file <- file.path(
  results_folder,
  "MA_Plot.png"
)

padj_cutoff <- 0.05
log2fc_cutoff <- 1

reference_condition <- "Control"
treatment_condition <- "Drought"


# ============================================================
# LOAD DESEQ2 RESULTS
# ============================================================

if (!file.exists(results_file)) {
  stop("DESeq2 results object was not found:\n", results_file)
}

results_deseq2 <- readRDS(
  results_file
)


# ============================================================
# CREATE AND SAVE MA PLOT
# ============================================================

png(
  filename = output_file,
  width = 2500,
  height = 2100,
  res = 300
)

plotMA(
  results_deseq2,

  alpha = padj_cutoff,

  ylim = c(-6, 6),

  main = paste(
    "MA Plot:",
    treatment_condition,
    "versus",
    reference_condition
  )
)


# Add fold-change threshold lines
abline(
  h = c(
    -log2fc_cutoff,
    log2fc_cutoff
  ),
  lty = 2,
  lwd = 1.5
)


dev.off()

message("MA plot saved:\n", output_file)
