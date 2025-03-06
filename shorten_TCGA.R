non_coding <- read.table("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_BRCA_BED_GENE_LEVEL/TCGA_BRCA_gene_level_log2_lifted_non_coding_tumor_sample.bed", header = TRUE, sep = "\t", check.names = FALSE)
coding <- read.table("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_BRCA_BED_GENE_LEVEL/TCGA_BRCA_gene_level_log2_lifted_coding_tumor_sample.bed", header = TRUE, sep = "\t", check.names = FALSE)
#cov <-  read.table("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/merged_cov_full_barcode.txt", header = TRUE, sep = "\t", check.names = FALSE)

shorten_colnames <- function(df) {
  new_colnames <- gsub("^((TCGA-[A-Z0-9]+-[A-Z0-9]+))-.*$", "\\1", colnames(df))  # Extract first three segments
  colnames(df) <- new_colnames
  return(df)
}

non_coding <- shorten_colnames(non_coding)
coding <- shorten_colnames(coding)
#cov <- shorten_colnames(cov)

write.table(non_coding, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_non_coding.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(coding, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_coding.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(cov, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_cov.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
