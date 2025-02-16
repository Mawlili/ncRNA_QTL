library(data.table)

samples <- scan("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/vcf/vcf_samples.txt", what = "", sep = "\t", quiet = TRUE)
short_non_coding <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_non_coding.bed")
short_non_coding <- as.data.frame(short_non_coding)
short_coding <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_coding.bed")
short_coding <- as.data.frame(short_coding)
cov <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_cov.txt")
cov <- as.data.frame(cov)

bed_cols <- c("#Chr", "start", "end", "pid", "gid", "strand")
cov_col <- c("SampleID")

overlap_cols_non_coding <- intersect(colnames(short_non_coding), samples)
overlap_cols_coding     <- intersect(colnames(short_coding), samples)
overlap_cov <- intersect(colnames(cov), samples)
subset_non_coding <- short_non_coding[ , c(bed_cols, overlap_cols_non_coding)]
subset_coding     <- short_coding[, c(bed_cols, overlap_cols_coding)]
subset_cov <- cov[,c(cov_col, overlap_cov)]

write.table(subset_non_coding, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_non_coding_overlap.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(subset_coding, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_coding_overlap.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(subset_cov, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_cov_overlap.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
