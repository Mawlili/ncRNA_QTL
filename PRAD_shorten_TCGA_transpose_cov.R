library(data.table)
non_coding <- read.table("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_PRAD_BED_GENE_LEVEL/TCGA_PRAD_gene_level_log2_lifted_non_coding_tumor_sample.bed", header = TRUE, sep = "\t", check.names = FALSE)
coding <- read.table("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_PRAD_BED_GENE_LEVEL/TCGA_PRAD_gene_level_log2_lifted_coding_tumor_sample.bed", header = TRUE, sep = "\t", check.names = FALSE)

shorten_colnames <- function(df) {
  new_colnames <- gsub("^((TCGA-[A-Z0-9]+-[A-Z0-9]+))-.*$", "\\1", colnames(df))  # Extract first three segments
  colnames(df) <- new_colnames
  return(df)
}

non_coding <- shorten_colnames(non_coding)
coding <- shorten_colnames(coding)

write.table(non_coding, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_non_coding.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(coding, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_coding.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#formatting cov
versions <- c(125, 150, 175, 200, 250)
for (v in versions) {
  # Construct file paths
  input_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/PRAD/cov_V", v, ".tsv")
  output_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_cov_V", v, ".txt")
  # short_cov_V125
  # Read and process the file
  cov <- fread(input_path, header = TRUE, sep = "\t", check.names = FALSE)
  cov <- cov[, -c(2, 3, 4, 5, 6)]  # Drop columns 2–6
  cov <- t(cov)                    # Transpose
  write.table(cov, output_path, sep = "\t", row.names = TRUE, col.names = FALSE)
}

for (v in versions) {
  input_path <- paste0("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_cov_V", v, ".txt")
  cov <- fread(input_path, header = TRUE, sep = "\t", check.names = FALSE)
  cov <- shorten_colnames(cov)
  write.table(cov, input_path, sep = "\t", row.names = FALSE, col.names = TRUE)
}
