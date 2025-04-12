#Library
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, readr,
               broom, stringr, future,
               furrr)
#source("./mediation_functions.R")
library(tictoc)
library(data.table)

P = 10000
num_workers = 5

#Function
mediation <- function(Y_G, X_S, M, X_C) {
  
  # Y_G: mean-centered gene expression marix for protein-coding gene G (n samples x 1 gene)
  # X_S: mean-centered SNP dosage matrix (n samples x p SNPs)
  # M: mean-centered expression matrix for m non-coding RNAs (n samples x m ncRNAs)
  # X_S: mean-centered matrix of covariates (n samples x k covariates)
  
  # Step 1: Fit the first regression model to predict Y_G (gene expression)
  # Y_G = X_S * β_S + M * β_M + X_C * β_C + ε_G
  # Where:
  # - Y_G: Gene expression
  # - X_S: SNP dosages
  # - β_S: Effect sizes of SNPs
  # - M: ncRNA expression matrix
  # - β_M: Effect sizes of ncRNAs
  # - X_C: Covariates matrix
  # - ε_G: Random error term
  first_stage <- lm(Y_G ~ X_S + M + X_C) |> broom::tidy()
  beta_M <- first_stage |>
    filter(stringr::str_detect(term, "M")) |>
    pull(estimate)
  beta_S <- first_stage |>
    filter(stringr::str_detect(term, "X_S")) |>
    pull(estimate)
  
  # Step 2: Fit the second regression model to predict the mediator M_j (ncRNA)
  # M_j = X_S * α_Mj + X_C * α_Cj + ε_Mj
  # Where:
  # - M_j: Expression of mediator ncRNA
  # - X_S: SNP dosages
  # - α_Mj: Effect of SNPs on mediator M_j
  # - X_C: Covariates matrix
  # - α_Cj: Effects of covariates on M_j
  # - ε_Mj: Random error term
  second_stage <- lm(M ~ X_S + X_C) |> broom::tidy()
  alpha_M <- second_stage |>
    filter(stringr::str_detect(term, "X_S")) |>
    pull(estimate)
  
  # total mediation effect (TME)
  # \alpha_M \cdot \beta_M
  TME <- sum(alpha_M * beta_M)
  
  # mediation proportion (MP)
  # min{1, TME / (\beta_S + TME) }
  MP <- min(c(1, TME / (beta_S + TME)))
  
  return(c("TME" = TME, "MP" = MP))
}

mediation_perm_test <- function(Y_G, X_S, M, X_C,
                                obs_TME,
                                P = 10000, num_workers = 1) {
  
  # Permutation test for total mediation effect (TME)
  # Y_G: Gene expression matrix for distal pc gene G (n samples x 1 gene)
  # X_S: SNP dosage matrix (n samples x p SNPs)
  # M: Expression matrix for m local ncRNAs (n samples x m ncRNAs)
  # X_C: Matrix of covariates (n samples x k covariates)
  # obs_TME: observed TME
  # P: Number of permutations
  # num_workers: number of workers/cores for parallelization
  
  # Set up parallelization
  plan(multicore, workers = num_workers)
  
  # Use future_map to run the permutations in parallel directly
  perm_results <- furrr::future_map_dfr(1:P, function(i) {
    # Permute Y_G
    indices <- sample(1:nrow(Y_G), replace = FALSE)
    Y_G_perm <- Y_G[indices, ]
    
    # Perform mediation analysis on the permuted data
    mediation_result <- mediation(Y_G = Y_G_perm, X_S = X_S, M = M, X_C = X_C)
    
    # Return the result for this permutation
    return(mediation_result)
  }, .options = furrr_options(seed = 2025))
  
  # H_0: TME = 0 vs H_1: TME \neq 0
  pval <- mean(abs(obs_TME) >= abs(perm_results[,"TME"]))
  
  return(pval)
}

mediation_analysis <- function(Y_G, X_S, M, X_C,
                               obs_TME,
                               P = 10000, num_workers = 1) {
  
  # two-stage mediation and permutation test for total mediation effect (TME)
  # Y_G: Gene expression matrix for distal pc gene G (n samples x 1 gene)
  # X_S: SNP dosage matrix (n samples x p SNPs)
  # M: Expression matrix for m local ncRNAs (n samples x m ncRNAs)
  # X_C: Matrix of covariates (n samples x k covariates)
  # obs_TME: observed TME
  # P: Number of permutations
  # num_workers: number of workers/cores for parallelization
  
  # mediation analysis
  mediation_est <- mediation(Y_G = Y_G, X_S = X_S,
                             M = M, X_C = X_C)
  # permutation test
  TME_pval <- mediation_perm_test(Y_G = Y_G, X_S = X_S,
                                  M = M, X_C = X_C,
                                  obs_TME = mediation_est["TME"],
                                  P = P, num_workers = num_workers)
  
  return(c(mediation_est, "TME_pval" = TME_pval))
}

#load prepared Rdata
triplets <- readRDS("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/out_files/triplets.RDS")
triplet_counts <- readRDS("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/out_files/triplet_counts.RDS")
# covariates
#covar <- readr::read_tsv("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/cov_remove_duplicate.txt",
 #                        col_names = TRUE)[,-1]

covar <- read.table("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/cov_remove_duplicate.txt", 
                 header = TRUE, 
                 sep = "\t", 
                 stringsAsFactors = FALSE, 
                 check.names = FALSE)
covar <- covar[, -1]
covar[102, ] <- ifelse(covar[102, ] == "female", 1,
                            ifelse(covar[102, ] == "male", 2, covar[102, ]))
covar <- data.frame(lapply(covar, as.numeric))

# gene expression data
nc_bed <- readr::read_tsv("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_non_coding_overlap_sorted_normalized.bed") # |> dplyr::select(-pid)
pc_bed <- readr::read_tsv("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/input/short_coding_overlap_sorted_normalized.bed") # |> dplyr::select(-pid)
# SNP dosages
snp_dosage <- readRDS("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/out_files/snp_dosage.RDS")

# Step 2: Data engineer
# X_C: Matrix of mean-centered covariates (n samples x k covariates)
covar_ids <- colnames(covar)

covar_ids <- gsub("\\.", "-", covar_ids)
X_C <- covar |>
  t() |>
  tidyr::as_tibble() |>
  reframe(across(everything(), ~ . - mean(., na.rm = TRUE))) |>
  as.matrix()

# filter to genes that appear at least once in triplets
# and to samples with covariates
exp_cov_intersect <- intersect(colnames(pc_bed),
                               covar_ids)
pc_bed <- pc_bed %>%
  mutate(gid = sub("\\..*", "", gid))
nc_bed <- nc_bed %>%
  mutate(gid = sub("\\..*", "", gid))


pc_intriplet_express <- pc_bed |>
  dplyr::filter(gid %in% unique(triplets$distal_pc)) |>
  dplyr::select("#chr", "start", "end", "gid", "strd",
         all_of(exp_cov_intersect))
nc_intriplet_express <- nc_bed |>
  dplyr::filter(gid %in% unique(triplets$local_nc)) |>
  dplyr::select("#chr", "start", "end", "gid", "strd",
         all_of(exp_cov_intersect))

# make expression data samples x genes
pc_long <- pc_intriplet_express |>
  dplyr::select(gid, contains("TCGA")) |>
  tidyr::pivot_longer(contains("TCGA"),
                      names_to = "sample_id",
                      values_to = "exp")

pc_wide <- pc_long |>
  tidyr::pivot_wider(names_from = gid,
                    values_from = exp,
                      values_fn = mean)

nc_long <- nc_intriplet_express |>
  dplyr::select(gid, contains("TCGA")) |>
  tidyr::pivot_longer(contains("TCGA"),
                      names_to = "sample_id",
                      values_to = "exp")

nc_wide <- nc_long |>
  tidyr::pivot_wider(names_from = gid,
                     values_from = exp,
                        values_fn = mean)

# Subset SNP dosage matrix to SNPs of interest
sig_snp_dosage <- snp_dosage[,which(colnames(snp_dosage) %in% triplets$snp_id)]

# Make sure samples by sample_id are aligned across data-sources
nc_wide <- nc_wide[match(covar_ids, nc_wide$sample_id), ]
pc_wide <- pc_wide[match(covar_ids, pc_wide$sample_id), ]
sig_snp_dosage <- sig_snp_dosage[match(covar_ids, rownames(sig_snp_dosage)), ]

nc_pc_aligned <- sum(nc_wide$sample_id == pc_wide$sample_id) == nrow(nc_wide)
snp_dosage_aligned <- sum(nc_wide$sample_id == rownames(sig_snp_dosage)) == nrow(nc_wide)
covar_aligned <- sum(nc_wide$sample_id == covar_ids) == nrow(nc_wide)
if (!(nc_pc_aligned & snp_dosage_aligned & covar_aligned)) stop("data are not aligned.")



#mediation analysis
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])
total_jobs <- nrow(triplet_counts)
chunk_size <- 1000
start_index <- (job_id - 1) * chunk_size + 1
end_index <- min(job_id * chunk_size, total_jobs)
for (t in  seq(start_index, end_index)) {
  
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
  output_file <- file.path(
  "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_mediation_filtered",
  paste0("mediation_res_", job_id, ".RDS")
)
  
  saveRDS(do.call("bind_rows", mediation_res_list) |> as_tibble(), output_file)
}
mediation_res <- do.call("bind_rows", mediation_res_list) |> as_tibble()

saveRDS(mediation_res, output_file)
