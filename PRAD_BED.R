library(SummarizedExperiment)
library(dplyr)
library(tximeta)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(data.table)
library(dplyr)
library(rtracklayer)


#load("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/se_gencodev38.RData")
load("/rsrch5/home/epi/bhattacharya_lab/data/TCGA/PRAD/Processed/tximeta/se_gencodev45.RData")
#txdb <- makeTxDbFromGFF("/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/txome/gencode_v38/gencode.v38.annotation.gtf", format = "gtf")
#se_gene <- summarizeToGene(se)
gse <- stg.gencodev45
expr_data <- assays(gse)$abundance
granges <- rowRanges(gse)
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
bed_df <- as.data.frame(bed_data)
bed_df$strand <- ifelse(bed_df$strand == 2, "-", "+")
colnames(bed_df)[1:6] <- c("Chr", "start", "end", "pid", "gid", "strand")
tpm_threshold <- 0.1
sample_count_threshold <- 0.75 * ( ncol(bed_df) - 6)
filtered_bed_df <- bed_df[rowSums(bed_df >= tpm_threshold) > sample_count_threshold, ]

filtered_bed_df[, 7:ncol(filtered_bed_df)] <- lapply(filtered_bed_df[, 7:ncol(filtered_bed_df)], as.numeric)
log2_normalized_bed_df <- filtered_bed_df
log2_normalized_bed_df[, 7:ncol(filtered_bed_df)] <- log2(filtered_bed_df[, 7:ncol(filtered_bed_df)] + 1)


write.table(
  log2_normalized_bed_df,
  file = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_PRAD_BED_GENE_LEVEL/TCGA_PRAD_gene_level_log2.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

bed_df$start <- as.numeric(bed_df$start)
bed_df$end <- as.numeric(bed_df$end)
gr <- GRanges(
    seqnames = bed_df$Chr,
    ranges = IRanges(start = bed_df$start, end = bed_df$end),
    strand = bed_df$strand
)
mcols(gr)$pid <- bed_df$pid
seqlevelsStyle(gr) <- "UCSC"
chain <- import.chain("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/liftover/hg38ToHg19.over.chain")
lifted <- liftOver(gr, chain)
lifted_df <- as.data.frame(lifted)
lifted_coords <- lifted_df[, c("seqnames", "start", "end")]

merged_df <- merge(bed_df, lifted_df[, c("seqnames", "start", "end", "pid")], 
                   by = "pid", all.x = TRUE)

merged_df$Chr <- merged_df$seqnames
merged_df$start.x <- merged_df$start.y
merged_df$end.x <- merged_df$end.y
  
  merged_df <- merged_df[, !names(merged_df) %in% c("seqnames", "start.y", "end.y")]
  colnames(merged_df)[colnames(merged_df) == "start.x"] <- "start"
  colnames(merged_df)[colnames(merged_df) == "end.x"] <- "end"

merged_df$Chr <- gsub("^chr", "", merged_df$Chr)

write.table(
  merged_df,
  file = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_PRAD_BED_GENE_LEVEL/TCGA_PRAD_gene_level_log2_lifted.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
bed <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_PRAD_BED_GENE_LEVEL/TCGA_PRAD_gene_level_log2_lifted.bed")
bed$pid <- sub("\\..*", "", bed$pid)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
    attributes = c("ensembl_gene_id", "gene_biotype"),
    filters = "ensembl_gene_id",
    values = bed$pid,
    mart = ensembl
)
merged_df <- merge(bed, gene_info, by.x = "pid", by.y = "ensembl_gene_id", all.x = TRUE)
merged_df <- merged_df %>%
  dplyr::select(Chr, start, end, pid, gid, strand, everything())
merged_df <- merged_df %>%
  filter(!is.na(Chr))

protein_coding_df <- subset(merged_df, gene_biotype == "protein_coding")
non_coding_df <- subset(merged_df, gene_biotype != "protein_coding")
protein_coding_df <- protein_coding_df %>%
   dplyr::select(-gene_biotype)
non_coding_df <- non_coding_df %>%
   dplyr::select(-gene_biotype)
write.table(
  protein_coding_df,
  file = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_PRAD_BED_GENE_LEVEL/TCGA_PRAD_gene_level_log2_lifted_coding.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  non_coding_df,
    file = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_PRAD_BED_GENE_LEVEL/TCGA_PRAD_gene_level_log2_lifted_non_coding.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
##linux
(head -n 1 TCGA_PRAD_gene_level_log2_lifted_coding_tumor_sample.bed && tail -n +2 TCGA_PRAD_gene_level_log2_lifted_coding_tumor_sample.bed | sort -k1,1 -k2,2n) > TCGA_PRAD_gene_level_log2_lifted_coding_tumor_sample_sorted.bed

(head -n 1 TCGA_PRAD_gene_level_log2_lifted_non_coding_tumor_sample.bed && tail -n +2 TCGA_PRAD_gene_level_log2_lifted_non_coding_tumor_sample.bed | sort -k1,1 -k2,2n) > TCGA_PRAD_gene_level_log2_lifted_non_coding_tumor_sample_sorted.bed

module load htslib
bgzip TCGA_PRAD_gene_level_log2_lifted_coding_tumor_sample_sorted.bed
bgzip TCGA_PRAD_gene_level_log2_lifted_non_coding_tumor_sample_sorted.bed

tabix -p bed TCGA_PRAD_gene_level_log2_lifted_non_coding_tumor_sample_sorted.bed.gz
tabix -p bed TCGA_PRAD_gene_level_log2_lifted_coding_tumor_sample_sorted.bed.gz
