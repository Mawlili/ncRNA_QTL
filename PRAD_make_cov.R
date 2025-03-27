# Load required library
library(dplyr)

# Read each file
setwd("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/PRAD")
gc_metrics <- read.delim("PRAD_gcSummaryMetrics.txt", stringsAsFactors = FALSE)
genome_metrics <- read.delim("PRAD_genomeMetrics.txt", stringsAsFactors = FALSE)
txome_metrics <- read.delim("PRAD_txomeMetrics.txt", stringsAsFactors = FALSE)
metadata <- read.csv("PRAD_metadata.csv", stringsAsFactors = FALSE)
load("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_hidden_cov/PRAD_hcp.RData")

# Merge files step-by-step using sample_ID
merged1 <- merge(gc_metrics, genome_metrics, by = "sample_ID")
merged_all <- merge(merged1, txome_metrics, by = "sample_ID")

# Extract relevant columns
metadata_subset <- metadata %>%
  select(fastq_hash, TCGA_ID, age_at_initial_pathologic_diagnosis)

#formatting and merge
merged_all$fastq_hash <- sub("\\..*$", "", merged_all$sample_ID)
final_merged <- merge(merged_all, metadata_subset, by = "fastq_hash")

# Save merged output
write.table(final_merged, "PRAD_Merged_QC_Metrics_with_age.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#PCA
#remove character
final_merged <- as.data.frame(final_merged)
rownames(final_merged) <- final_merged[[74]]
final_merged_filt <- final_merged[, -c(74)]
cov_num <- final_merged_filt[, -(1:5)]

#merge PCA back
cov_matrix <- as.matrix(cov_num)
cov_matrix <- apply(cov_matrix, 2, as.numeric)
pca_result <- prcomp(cov_matrix)
rownames(pca_result$x) <- rownames(cov_num)
pca_data <- as.data.frame(pca_result$x[, 1:3])  
pca_data$TCGA_ID <- rownames(pca_data)
cov_with_pca <- merge(final_merged, pca_data, by = "TCGA_ID")

#Merge hcp
hcp_z <- as.data.frame(hcp$Z)
hcp_z$TCGA_ID <- rownames(hcp_z)
cov_with_pca_hcp <- merge(cov_with_pca, hcp_z, by = "TCGA_ID")

# Save merged output
write.table(cov_with_pca_hcp, "PRAD_Merged_QC_Metrics_with_age_pca_hcp.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#test different # of hcp
df <- read.delim("PRAD_Merged_QC_Metrics_with_age_pca_hcp.txt", stringsAsFactors = FALSE)

# Find column indices
base_cols <- which(colnames(df) == "PC3")  # get index of PC3

# Define extra variable sets
v_sets <- list(
  V125 = paste0("V", 1:125),
  V150 = paste0("V", 1:150),
  V175 = paste0("V", 1:175),
  V200 = paste0("V", 1:200),
  V250 = paste0("V", 1:250)
)

for (name in names(v_sets)) {
  cols_to_include <- c(1:base_cols, which(colnames(df) %in% v_sets[[name]]))
  cov_subset <- df[, cols_to_include]
  file_name <- paste0("cov_", name, ".tsv")
  write.table(cov_subset, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved:", file_name, "\n")
}
