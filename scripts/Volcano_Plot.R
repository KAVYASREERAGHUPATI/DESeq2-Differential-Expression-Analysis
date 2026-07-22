# ============================================================
# VOLCANO PLOT
# ============================================================

# Install packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install(
    "EnhancedVolcano",
    ask = FALSE,
    update = FALSE
  )
}


# Load package
suppressPackageStartupMessages({
  library(EnhancedVolcano)
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
  "Volcano_Plot.png"
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

volcano_data <- as.data.frame(
  results_deseq2
)

volcano_data$Gene_ID <- rownames(
  volcano_data
)


# Replace missing adjusted p-values with 1 only for plotting
volcano_data$padj_plot <- volcano_data$padj
volcano_data$padj_plot[
  is.na(volcano_data$padj_plot)
] <- 1


# ============================================================
# CREATE AND SAVE VOLCANO PLOT
# ============================================================

png(
  filename = output_file,
  width = 2600,
  height = 2200,
  res = 300
)

EnhancedVolcano(
  volcano_data,

  lab = volcano_data$Gene_ID,

  x = "log2FoldChange",

  y = "padj_plot",

  pCutoff = padj_cutoff,

  FCcutoff = log2fc_cutoff,

  title = paste(
    treatment_condition,
    "versus",
    reference_condition
  ),

  subtitle = paste0(
    "padj < ",
    padj_cutoff,
    " and |log2FC| > ",
    log2fc_cutoff
  ),

  xlab = expression(
    Log[2] ~ "fold change"
  ),

  ylab = expression(
    -Log[10] ~ "adjusted p-value"
  ),

  pointSize = 2.2,

  labSize = 3,

  max.overlaps = 20,

  drawConnectors = TRUE,

  widthConnectors = 0.5,

  legendPosition = "right",

  legendLabSize = 11,

  legendIconSize = 4,

  border = "partial"
)

dev.off()

message("Volcano plot saved:\n", output_file)
