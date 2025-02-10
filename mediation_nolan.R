# mediation.R
# 2-stage linear regression mediation analysis
# Nolan Cole, UW Biostatistics, 31 January 2024
# libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, readr,
               broom, stringr, future,
               furrr)
source("./mediation_functions.R")

# P = 10000
# num_workers = 50
args <- commandArgs(trailingOnly = TRUE)
P <- args[1] |> as.integer() # number of permutations
num_workers <- args[2] |> as.integer() # number of cores

cat("P: ", P, "\n")
cat("num_workers: ", num_workers, "\n")

# Step 1: Load data
# triplets
triplets <- readRDS("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/out_files/triplets.RDS")
triplet_counts <- readRDS("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/out_files/triplet_counts.RDS")
# covariates
covar <- readr::read_tsv("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/cov_remove_duplicate.txt",
                         col_names = TRUE)[,-1]
# gene expression data
nc_bed <- readr::read_tsv("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_BRCA_BED_GENE_LEVEL/TCGA_BRCA_gene_level_log2_lifted_non_coding_sorted.bed") |> dplyr::select(-pid)
pc_bed <- readr::read_tsv("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/TCGA_BRCA_BED_GENE_LEVEL/TCGA_BRCA_gene_level_log2_lifted_coding_sorted.bed") |> dplyr::select(-pid)
# SNP dosages
snp_dosage <- readRDS("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/out_files/snp_dosage.RDS")

# Step 2: Data engineer
# X_C: Matrix of mean-centered covariates (n samples x k covariates)
covar_ids <- colnames(covar)
X_C <- covar |>
  t() |>
  tidyr::as_tibble() |>
  reframe(across(everything(), ~ . - mean(., na.rm = TRUE))) |>
  as.matrix()

# filter to genes that appear at least once in triplets
# and to samples with covariates
exp_cov_intersect <- intersect(colnames(pc_bed),
                               covar_ids)
pc_intriplet_express <- pc_bed |>
  dplyr::filter(gid %in% unique(triplets$distal_pc)) |>
  select("#Chr", "start", "end", "gid", "strand",
         all_of(exp_cov_intersect))
nc_intriplet_express <- nc_bed |>
  dplyr::filter(gid %in% unique(triplets$local_nc)) |>
  select("#Chr", "start", "end", "gid", "strand",
         all_of(exp_cov_intersect))

# make expression data samples x genes
pc_long <- pc_intriplet_express |>
  select(gid, contains("GTEX")) |>
  tidyr::pivot_longer(contains("GTEX"),
                      names_to = "sample_id",
                      values_to = "exp")
pc_wide <- pc_long |>
  tidyr::pivot_wider(names_from = gid,
                     values_from = exp)

nc_long <- nc_intriplet_express |>
  select(gid, contains("GTEX")) |>
  tidyr::pivot_longer(contains("GTEX"),
                      names_to = "sample_id",
                      values_to = "exp")
nc_wide <- nc_long |>
  tidyr::pivot_wider(names_from = gid,
                     values_from = exp)

# Subset SNP dosage matrix to SNPs of interest
sig_snp_dosage <- snp_dosage[,which(colnames(snp_dosage) %in% triplets$snp_id)]

# Make sure samples by sample_id are aligned across data-sources
nc_pc_aligned <- sum(nc_wide$sample_id == pc_wide$sample_id) == nrow(nc_wide)
snp_dosage_aligned <- sum(nc_wide$sample_id == rownames(sig_snp_dosage)) == nrow(nc_wide)
covar_aligned <- sum(nc_wide$sample_id == covar_ids) == nrow(nc_wide)
if (!(nc_pc_aligned & snp_dosage_aligned & covar_aligned)) stop("data are not aligned.")

# Step 3: Mediation analysis
# perform mediation analysis with triplets:
# {SNP, distal pc gene, local ncRNA(s) }
mediation_res_list <- vector("list", length = nrow(triplet_counts))
for (t in seq_along(1:nrow(triplet_counts))) {
  
  cat("Triplet: ", t, "/", nrow(triplet_counts), "\n")
  
  # Get triplet
  snp <- triplet_counts[t,"snp_id"] |> unlist()
  names(snp) <- "snp"
  pc_gene <- triplet_counts[t,"distal_pc"] |> unlist()
  names(pc_gene) <- "pc_gene"
  nc_rna <- triplets |>
    filter(snp_id == snp, distal_pc == pc_gene) |>
    pull(local_nc)
  names(nc_rna) <- paste0("nc_rna", 1:length(nc_rna))
  
  # X_S: SNP dosage matrix (n samples x p SNPs)
  # Y_G: Gene expression matrix for distal pc gene G (n samples x 1 gene)
  # M: Expression matrix for m local ncRNAs (n samples x m ncRNAs)
  X_S <- tibble(snp = sig_snp_dosage[,snp])
  Y_G <- pc_wide[,pc_gene]
  M <- nc_wide[,colnames(nc_wide) %in% nc_rna]
  
  # Mean center
  Y_G <- Y_G |>
    reframe(across(everything(), ~ . - mean(., na.rm = TRUE))) |>
    as.matrix()
  X_S <- X_S |>
    reframe(across(everything(), ~ . - mean(., na.rm = TRUE))) |>
    as.matrix()
  M <- M |>
    reframe(across(everything(), ~ . - mean(., na.rm = TRUE))) |>
    as.matrix()
  
  # mediation analysis
  tictoc::tic()
  mediation_res_list[[t]] <- mediation_analysis(Y_G = Y_G, X_S = X_S,
                                               M = M, X_C = X_C,
                                               P = P, num_workers = num_workers)
  tictoc::toc()
  
  # store triplet information
  mediation_res_list[[t]] <- c(mediation_res_list[[t]],
                               "snp" = snp,
                               "pc_gene" = pc_gene,
                               nc_rna)
  
  # save as we go
  saveRDS(do.call("bind_rows", mediation_res_list) |> as_tibble(), "mediation_res.RDS")
}
mediation_res <- do.call("bind_rows", mediation_res_list) |> as_tibble()

saveRDS(mediation_res, "mediation_res.RDS")

cat("Finished running. Good workout!\n")
