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
bed_nc <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_non_coding_overlap_sorted_normalized.bed", header = TRUE)
bed_nc_sig <- bed_nc[bed_nc$id %in% nc_sig, ]
bed_pc <- fread("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_coding_overlap_sorted_normalized.bed", header = TRUE)
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


#putting them together
library(vcfR)        # fast in‑memory VCF reader
library(data.table) 
vcf_BRCA <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/vcf/combined_formatted_vcf_file.vcf.gz"
vcf_PRAD <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/vcf/PRAD/PRAD_merged_chr_short.vcf.gz"
qtl_files <- list(
  BRCA_protein = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_mediation_filtered/cis_qtl_output/significant_output_cis_qtl_protein_coding.txt",
  BRCA_ncRNA   = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_mediation_filtered/cis_qtl_output/significant_output_cis_qtl_non_protein_coding.txt",
  PRAD_protein = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_mediation/cis_qtl_output/significant_output_cis_qtl_protein_coding.txt",
  PRAD_ncRNA   = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_mediation/cis_qtl_output/significant_output_cis_qtl_non_protein_coding.txt"
)

vcf_to_lookup <- function(vcf_path) {
  v <- read.vcfR(vcf_path, verbose = FALSE)         # ⏱ a minute or two for whole‑genome VCF
  data.table(
    SNP = v@fix[, "ID"],
    REF = v@fix[, "REF"],
    ALT = sub(",.*$", "", v@fix[, "ALT"])  # keep first ALT if multiallelic
  )
}

lookup_BRCA <- vcf_to_lookup(vcf_BRCA);  setkey(lookup_BRCA, SNP)
lookup_PRAD <- vcf_to_lookup(vcf_PRAD);  setkey(lookup_PRAD, SNP)

process_qtl <- function(path, lookup) {
  dt <- fread(path, header = FALSE, showProgress = FALSE)

  ## keep columns 8 1 14 15 12  →  SNP, ncRNA, beta, SE, P_value
  dt <- dt[, .(SNP      = V8,
               ncRNA    = V1,
               beta     = V14,
               SE       = V15,
               P_value  = V12)]

  merge(dt, lookup, by = "SNP", all.x = TRUE)
}

lookup_map <- list(
  BRCA_protein = lookup_BRCA,
  BRCA_ncRNA   = lookup_BRCA,
  PRAD_protein = lookup_PRAD,
  PRAD_ncRNA   = lookup_PRAD
)


for (nm in names(qtl_files)) {
  message("Processing ", nm, " …")
  dt <- process_qtl(qtl_files[[nm]], lookup_map[[nm]])

  ## choose an output name that’s easy to spot:
  ##   e.g. “BRCA_protein_with_REF_ALT.tsv”
  out_name <- paste0(nm, "_with_REF_ALT.tsv")

  fwrite(dt,
         file   = out_name,
         sep    = "\t",
         quote  = FALSE)

  message("  ↳ wrote ", out_name)
}


