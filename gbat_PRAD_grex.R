# grex.R
# 10-fold CV of genetically-regulated expression
# using SNPs within 1 Mb using elastic net, LASSO, HAL, and SuSiE
# Nolan Cole, UW Biostatistics, 19 February 2025

###########################
## Commandline arguments ##
###########################
# num_workers = 1
args <- commandArgs(trailingOnly = TRUE)

num_workers <- args[1] |> as.integer()
worker_id <- as.integer(args[2])

cat("Running on worker ID:", worker_id, "\n")
cat("Worker is 1/", num_workers, "\n")

set.seed(worker_id)

###########################
## Parameters/filepaths  ##
###########################
bp_window <- 1e6
pred_nc_grex <- paste0(worker_id, "_pred_nc_grex.txt")
data_dir <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat/1_dataprep/data/"
home_dir <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat/"
results_dir <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat/results/"
model_dir <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat/results/best_models/"

#########################
## Libraries/Functions ##
#########################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("isotwas", quietly = TRUE)) devtools::install_github("bhattacharya-a-bt/isotwas")
if (!requireNamespace("hal9001", quietly = TRUE)) devtools::install_github("tlverse/hal9001")
# BiocManager::install(c("GenomicRanges"))
pacman::p_load(isotwas, hal9001,
               dplyr, tidyr, readr,
               broom, stringr, future,
               furrr, magrittr,
               glmnet, susieR, GenomicRanges,
               data.table, tictoc, Rfast,
               caret, rlist)
source(paste0(home_dir, "grex_functions.R"))

##################
## read in data ##
##################
# non-coding
nc_bed <- readr::read_delim(paste0(data_dir, "nc.bed"))
nc_bed_long <- nc_bed |>
  pivot_longer(cols = contains("TCGA"),
               names_to = "sample_id",
               values_to = "expression")
# noncoding<- nc_bed_long |>
#   select(sample_id, expression, gid) |>
#   pivot_wider(names_from = gid, values_from = expression)
# Convert nc_bed to GRanges object
nc_gr <- GRanges(
  seqnames = nc_bed$`#chr`,
  ranges = IRanges(start = nc_bed$start, end = nc_bed$end),
  gid = nc_bed$gid
)

# protein coding
pc_bed <- readr::read_delim(paste0(data_dir, "pc.bed"))
# pivot longer
pc_bed_long <- pc_bed |>
  pivot_longer(cols = contains("TCGA"),
               names_to = "sample_id",
               values_to = "expression")
# Convert pc_bed to GRanges object
pc_gr <- GRanges(
  seqnames = pc_bed$`#chr`,
  ranges = IRanges(start = pc_bed$start, end = pc_bed$end),
  gid = pc_bed$gid
)

# SNP BiomaRt/GRanges
snp_gr <- readRDS(paste0(home_dir, "biomart_snp.RDS"))

# SNP dosages
snp_dosage <- readRDS(paste0(data_dir, "snp_dosage.RDS"))

## Set up which ncRNAs are analyzed on this worker ##
num_ncRNAS <- nrow(nc_bed)
# Compute indices
all_indices <- 1:num_ncRNAS
chunk_size <- ceiling(num_ncRNAS / num_workers)

start_idx <- (worker_id - 1) * chunk_size + 1
end_idx <- min(worker_id * chunk_size, num_ncRNAS)
worker_indices <- all_indices[start_idx:end_idx]

cat("Worker", worker_id, "handles indices", start_idx, "to", end_idx, "\n")

worker_ncrna <- nc_gr$gid[worker_indices]

tictoc::tic()
nc_rna_list <- list() 
for (i in seq_along(1:length(worker_ncrna))) {
  
  ncrna <- worker_ncrna[i]
  
  ## Prep data for specific gene ##
  cat("Gene:", ncrna, "\n")
  
  # filter to pc genes within 1Mb, X
  # Get the GRanges entry for the target gene
  target_gene <- nc_gr[nc_gr$gid == ncrna]
  
  # Define a 1Mb window around the target gene
  target_region <- resize(target_gene, width = bp_window * 2, fix = "center")
  
  # Find overlaps
  hits <- findOverlaps(snp_gr, target_region)
  
  # Extract matching genes
  snps_within_1mb <- snp_gr[queryHits(hits)]
  dosages <- snp_dosage[, colnames(snp_dosage) %in% snps_within_1mb$names] |>
    as_tibble() |>
    mutate("sample_id" = rownames(snp_dosage))
  
  if (nrow(dosages) == 0) {
    cat("No SNPs found for gene: ", ncrna, "\n")
    return(NULL)
  }
  
  #  Prep data for model Y ~ X
  YX <- nc_bed_long |>
    as_tibble() |>
    dplyr::select(sample_id, expression, gid) |>
    filter(gid == ncrna) |>
    inner_join(dosages, by = "sample_id") |>
    dplyr::select(-gid, -sample_id)
  
  Y <- YX |> dplyr::select(expression) |> scale() |> as.matrix()
  X <- YX |> dplyr::select(-expression) |> as.matrix()
  if (anyNA(Y)) {
    cat("  » Expression vector contains NA/NaN – skipping\n")
    next
  }
  if (ncol(X) == 0) next
  
  cat("Dimensions of Y:", dim(Y), "\n")
  cat("Dimensions of X:", dim(X), "\n")
  
  ## Get best model ##
  gene_res <- grex_models(ncrna = ncrna,
                          Y = Y, X = X,
                          dosages = dosages,
                          sample_ids = dosages$sample_id)
  
  # If the result is valid (not NULL), write to the file
  if (!is.null(gene_res)) {
    cat("Saving", ncrna, "(", i, "/", length(worker_ncrna), ")\n")
    # Save predicted gene expression
    data.table::fwrite(gene_res$row,
                       paste0(results_dir, "/", pred_nc_grex),
                       append = TRUE,
                       row.names = FALSE,
                       sep = '\t',
                       quote = FALSE)
    # save model
    saveRDS(gene_res$best_model, paste0(model_dir, "/", ncrna, "_bestmodel.RDS"))
  }
}
tictoc::toc()
