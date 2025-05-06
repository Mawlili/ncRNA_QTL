# ▶ Inputs 
med_file   <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_mediation/mediation_paris.txt"
vcf_file   <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/vcf/PRAD/PRAD_merged_chr_short.vcf.gz"
expr_bed   <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_mediation/pc_sig_mediation_filtered.bed.gz"
nc_anno    <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_mediation/ncRNA_annotation.tsv"
cov_file   <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/input/short_cov_overlap_V25.txt" 

# ▶ Runtime settings
window_bp  <- 1e6   # 1 Mb
threads    <- 8
out_root   <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_mediation/trans_qtl_results"  # per‑row subdirs

suppressPackageStartupMessages({
  library(data.table)
})

# ─── 1. Load tables ------------------------------------------------------------
med      <- fread(med_file)
anno     <- fread(nc_anno, col.names = c("chr", "tss", "gene_id"))
expr_dt  <- fread(cmd = sprintf("zcat %s", expr_bed))   # chr start end gene_id samples…
setnames(expr_dt, 4, "gene_id")

# Identify ncRNA columns
nc_cols <- grep("^nc_rna", names(med), value = TRUE)

# ─── 2. Loop over rows ---------------------------------------------------------
for (row_i in seq_len(nrow(med))) {
  cat("▶ Row", row_i, "\n")
  row_dir <- file.path(out_root, sprintf("row_%04d", row_i))
  dir.create(row_dir, showWarnings = FALSE, recursive = TRUE)

  # Extract pc_gene & ncRNAs ----------------------------------------------------
  pc_gene <- med[row_i, pc_gene.pc_gene]
  nc_rnas <- unlist(med[row_i, ..nc_cols], use.names = FALSE)
  nc_rnas <- nc_rnas[!is.na(nc_rnas)]

  # Build BED of ±1 Mb windows --------------------------------------------------
  win_dt <- merge(data.table(gene_id = nc_rnas), anno, by = "gene_id", all.x = TRUE)
  if (any(is.na(win_dt$tss))) {
    cat("   • missing annotation for", sum(is.na(win_dt$tss)), "ncRNA(s) – skipping\n"); next
  }
  win_dt[, `:=`(start = pmax(tss - window_bp, 0), end = tss + window_bp)]
  bed_path <- file.path(row_dir, "regions.bed")
  fwrite(win_dt[, .(chr, start, end)], bed_path, sep = "\t", col.names = FALSE)

  # ── 2a. Slice VCF with plink2 -----------------------------------------------
   vcf_prefix <- file.path(row_dir, "snps")       # will create snps.vcf.gz
vcf_slice  <- paste0(vcf_prefix, ".vcf.gz")

## 1.  Write a ranges.txt file:  chr <tab> start <tab> end
range_file <- file.path(row_dir, "ranges.txt")
fwrite(unique(win_dt[, .(chr, start, end)]), range_file,
       sep = "\t", col.names = FALSE)

## 2.  Call plink2 with --extract range
plink_cmd <- sprintf(
  "plink2 --threads %d --vcf %s --extract range %s --export vcf bgz --out %s",
  threads,
  vcf_file,
  range_file,
  vcf_prefix
)

cat("   • PLINK2 cmd:\n     ", plink_cmd, "\n")
if (system(plink_cmd) != 0) {              # run plink2
  cat("   ✖ plink2 failed → skip row\n"); next
}

system(sprintf("tabix -p vcf %s", vcf_slice)) 

  # ── 2b. Trim expression BED to pc_gene --------------------------------------
 bed_plain  <- file.path(row_dir, "pc_gene.bed")
bed_slice  <- paste0(bed_plain, ".gz")
  
fwrite(
  expr_dt[gene_id == pc_gene],
  file  = bed_plain,
  sep   = "\t",
  quote = FALSE
)
system(sprintf("bgzip -c %s > %s", bed_plain, bed_slice))
system(sprintf("tabix -p bed %s", bed_slice))

  # ── 2c. Run QTLtools ---------------------------------------------------------
  qtl_out <- file.path(row_dir, "qtltools.nominal.txt")
  qtl_cmd <- sprintf(
    "QTLtools trans --vcf %s --bed %s --cov %s --nominal --threshold 1 --normal --out %s",
    vcf_slice, bed_slice, cov_file, qtl_out)
  cat("   • QTLtools cmd:\n     ", qtl_cmd, "\n")
  status <- system(qtl_cmd)
  if (status != 0) { cat("   ✖ QTLtools failed on row", row_i, "\n") }
}

cat("\n✓ All rows processed. Check", out_root, "for per‑row results.\n")

