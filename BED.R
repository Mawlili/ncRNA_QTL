library(SummarizedExperiment)
library(dplyr)
library(tximeta)

load("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/se.RData")

expr_data <- assays(se)$counts
granges <- rowRanges(se)
gr_df <- as.data.frame(granges)

bed_data <- cbind(
  chrom = gr_df$seqnames,
  start = gr_df$start,  
  end = gr_df$end,
  phenotype_ID = rownames(gr_df), 
  gene_ID = gr_df$gene_id,
  strand = gr_df$strand,
  expr_data
)

bed_columns <- c("chrom", "start", "end", "phenotype_ID", "gene_ID", "strand")
sample_columns <- colnames(expr_data)
bed_data <- bed_data[, c(bed_columns, sample_columns)]

write.table(
  bed_data,
  file = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BED/qtl_input.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
