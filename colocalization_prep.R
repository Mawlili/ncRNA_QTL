#find_sig_med
data <- readRDS("BRCA_mediation_combined.rds")
data_sig <- data[data$TME_pval < 0.05, ]
pc_sig <- data_sig$pc_gene.pc_gene
pc_sig <- pc_sig[!is.na(pc_sig)]
pc_sig <-unique(pc_sig)
nc_sig <- data_sig[ ,6:ncol(data_sig) ]
nc_sig <- unlist(nc_sig, use.names = FALSE)
nc_sig <- nc_sig[!is.na(nc_sig)]
nc_sig <- unique(nc_sig)
write.table(nc_sig, file = "nc_sig.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(pc_sig, file = "pc_sig.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#filter for sig
library(data.table)
nc_sig <- scan("nc_sig.txt", what = "", quiet = TRUE)
pc_sig <- scan("pc_sig.txt", what = "", quiet = TRUE)
bed_nc <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_non_coding_overlap_sorted_normalized.bed", header = TRUE)
bed_nc_sig <- bed_nc[bed_nc$id %in% nc_sig, ]
bed_pc <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_coding_overlap_sorted_normalized.bed", header = TRUE)
bed_pc_sig <- bed_pc[bed_pc$id %in% pc_sig, ]
fwrite(bed_nc_sig, "nc_sig_mediation_filtered.bed", sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(bed_pc_sig, "pc_sig_mediation_filtered.bed", sep = "\t", quote = FALSE, col.names = TRUE)


#cis_qtl
bgzip -c nc_sig_mediation_filtered.bed > nc_sig_mediation_filtered.bed.gz
tabix -p bed nc_sig_mediation_filtered.bed.gz

bgzip -c pc_sig_mediation_filtered.bed > pc_sig_mediation_filtered.bed.gz
tabix -p bed pc_sig_mediation_filtered.bed.gz


  QTLtools cis --vcf /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/vcf/combined_formatted_vcf_file.vcf.gz \
               --bed /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_mediation_filtered/nc_sig_mediation_filtered.bed.gz \
               --cov /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_cov_overlap.txt \
               --out /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_mediation_filtered/cis_qtl_output/significant_output_cis_qtl_non_protein_coding.txt \
               --nominal 1
               --normal

QTLtools cis --vcf /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/vcf/combined_formatted_vcf_file.vcf.gz \
             --bed /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_mediation_filtered/pc_sig_mediation_filtered.bed.gz \
             --cov /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_cov_overlap.txt \
             --out /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_mediation_filtered/cis_qtl_output/significant_output_cis_qtl_protein_coding.txt \
             --nominal 1
             --normal
