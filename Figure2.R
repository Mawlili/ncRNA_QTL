#Library
library(ggplot2)
library(dplyr)
library(biomaRt)

#Data processing
data <- readRDS("BRCA_mediation_combined.RDS")
sig_data <- data[data$TME_pval < 0.05, ]
sig_data <- sig_data[!is.na(sig_data$TME), ]
#Get SNP location
ensembl_snp <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

snp_list <- sig_data$snp.snp
snp_info <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                  filters = "snp_filter", values = snp_list, mart = ensembl_snp)
#colnames(snp_info)[3] <- "snp_chrom_start"

merged_data <- left_join(sig_data, snp_info, by = c("snp.snp" = "refsnp_id"))
merged_data$chr_name

#Get coding gene location
ensembl_gene <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = sig_data$pc_gene.pc_gene,
  mart = ensembl_gene
)
colnames(gene_info)[2] <- "pc_gene_chrom_start"

merged_data_gene <- left_join(merged_data, gene_info, by = c("pc_gene.pc_gene" = "ensembl_gene_id"))

#TME normalization
merged_data_gene$TME <- as.numeric(merged_data_gene$TME)
merged_data_gene <- merged_data_gene %>%
  mutate(abs_scaled_TME = abs(TME) / max(abs(TME), na.rm = TRUE))

#Effect direction
merged_data_gene <- merged_data_gene %>%
  mutate(effect_direction = ifelse(TME > 0, "Positive", "Negative"))

#chromosomes
canonical_chromosomes <- c(as.character(1:22), "X", "Y", "MT")
merged_data_gene <- merged_data_gene %>%
  filter(pc_gene_chrom_start %in% canonical_chromosomes, chr_name %in% canonical_chromosomes)

#convert to factor
merged_data_gene <- merged_data_gene %>%
  mutate(pc_gene_chrom_start = factor(pc_gene_chrom_start, levels = canonical_chromosomes),
         chr_name = factor(chr_name, levels = canonical_chromosomes))

#make figure, eSNP (X-axis) position vs. transcription start site (TSS) of pcGene (Y-axis) at FDR-adjusted P < 0.01, sized by absolute scaled TME and colored by direction of effect
p <- ggplot(merged_data_gene, aes(x = chr_name, y = pc_gene_chrom_start)) +
  geom_point(aes(size = abs_scaled_TME, color = effect_direction), alpha = 0.8) +
  
  # color for +/-
  scale_color_manual(values = c("Negative" = "red", "Positive" = "blue")) +
  
  # Size scale for effect magnitude
  scale_size_continuous(range = c(1, 10), name = "Absolute Scaled Effect") +
  
  # Dashed diagonal reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  
  # Labels & themes
  labs(title = "eSNP vs. TSS of pcGene",
       x = "eSNP Position",
       y = "TSS of pcGene",
       color = "Direction",
       size = "Absolute Scaled Effect") +
  
  theme_minimal() +
  theme(legend.position = "bottom")
