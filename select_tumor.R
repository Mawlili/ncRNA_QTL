library(data.table)
# Define file paths
input_file <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_BRCA_BED_GENE_LEVEL/TCGA_BRCA_gene_level_log2_lifted_non_coding.bed"
output_file <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_BRCA_BED_GENE_LEVEL/TCGA_BRCA_gene_level_log2_lifted_non_coding_tumor_sample.bed"


bed_data <- fread(input_file, header = TRUE, sep = "\t")

# Extract column names
column_names <- colnames(bed_data)

# Function to check if the fourth number (after the third "-") is < 10
is_selected <- function(name) {
  parts <- unlist(strsplit(name, "-"))
  if (length(parts) >= 4) {
    fourth_number <- as.numeric(gsub("[^0-9]", "", parts[4]))  # Extract numeric part
    return(!is.na(fourth_number) && fourth_number < 10)
  }
  return(FALSE)
}

# Identify selected columns
selected_columns <- column_names[sapply(column_names, is_selected)]

# Keep the first six columns and the selected sample columns
filtered_data <- bed_data[, c(colnames(bed_data)[1:6], selected_columns), drop = FALSE]

write.table(filtered_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

