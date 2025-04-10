library(data.table)



samples <- scan("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/vcf/PRAD/short_sample_names.txt", what = "", sep = "\t", quiet = TRUE)
short_non_coding <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_non_coding.bed")
short_non_coding <- as.data.frame(short_non_coding)
short_coding <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_coding.bed")
short_coding <- as.data.frame(short_coding)
cov <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_cov_V250.txt")
cov <- as.data.frame(cov)

bed_cols <- c("Chr", "start", "end", "pid", "gid", "strand")
cov_col <- c("TCGA_ID")


overlap_cols_non_coding <- intersect(colnames(short_non_coding), samples)
overlap_cols_coding     <- intersect(colnames(short_coding), samples)

non_coding_col <- c(bed_cols, overlap_cols_non_coding)
coding_col <- c(bed_cols, overlap_cols_coding)

subset_non_coding <- short_non_coding[ , non_coding_col]
subset_coding     <- short_coding[, coding_col]

overlap_cov <- intersect(colnames(cov), samples)
subset_cov <- cov[,c(cov_col, overlap_cov)]

colnames(subset_non_coding)[1] <- "#Chr"
colnames(subset_coding)[1] <- "#Chr"

write.table(subset_non_coding, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_non_coding_overlap.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(subset_coding, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_coding_overlap.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(subset_cov, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_cov_overlap_V250.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
