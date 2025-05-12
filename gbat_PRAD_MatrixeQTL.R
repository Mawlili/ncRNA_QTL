# matrixeqtl.R
# Use MatrixEQTL to estimate the association
# between the ncRNA GReX and distal pcGene expression,
# adjusting for the optimized set of covariates from eQTL mapping
# Nolan Cole, UW Biostat, 18 April 2025

###########################
## Parameters/Filepaths  ##
###########################
data_dir <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat/1_dataprep/data/"
home_dir <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat/"
results_dir <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat/results"
model_dir <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat/results/best_models"
biomart_file <- "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_gbat/1_dataprep/data/biomaRt_genome.RData"

# cis
output_file_name <- paste0(results_dir,"/MeQTL_cis_results.txt")
pvOutputThreshold_cis <- 1
cisDist <- 1e6

# trans
output_file_name_tra <- paste0(results_dir, "/MatrixeQTL_trans.txt")
pvOutputThreshold_tra <- 1  # or whatever threshold you want

#########################
## Libraries/Functions ##
#########################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GenomicRanges"))
pacman::p_load(dplyr, tidyr, readr,
               broom, stringr,
               furrr, magrittr,
               GenomicRanges,
               data.table, tictoc,
               MatrixEQTL)

##################
## Read in data ##
##################

## predicted ncRNA expression ##
cat("Reading in GReX\n")
# From gbat.R
# Expression ~ genotypes
grex_files <- list.files(results_dir,
                         pattern = "pred_nc_grex.txt",
                         full.names = TRUE)
grex_list <- list()
for (i in 1:length(grex_files)) {
  grex_list[[i]] <- readr::read_delim(paste0(grex_files[i])) |> 
    dplyr::select(-model) |>
    distinct() |>
    as.data.frame()
  rownames(grex_list[[i]]) <- grex_list[[i]]$gene
  grex_list[[i]]$gene <- NULL
}
grex <- do.call("rbind", grex_list)

# Make sliced data format
nc_grex <- SlicedData$new()
nc_grex$CreateFromMatrix(as.matrix(grex))

# location
load(biomart_file)
# nc_gene_info <- gene_info |> filter(ensembl_gene_id %in% grex$id)
nc_gene_info <- gene_info |> filter(ensembl_gene_id %in% rownames(grex))
nc_genepos <- data.frame(
  # "snp" = grex$id,
  "snp" = rownames(grex),
  "chr" = nc_gene_info$chromosome_name,
  "pos" = nc_gene_info$start_position
)

## pcGene Expression ##
cat("Reading in pcGene Expression\n")
pc_bed <- readr::read_delim(paste0(data_dir, "pc.bed")) |> select(-id)
pc_ex_mat <- pc_bed |>
    select(contains("GTEX")) |>
  as.data.frame()
rownames(pc_ex_mat) <- pc_bed$gid
# Make sliced data format
pc_ex <- SlicedData$new()
pc_ex$CreateFromMatrix(as.matrix(pc_ex_mat))

# location
pc_pos <- pc_bed |>
  select("geneid" = gid,
         "chr" = "#chr",
         "s1" = start,
         "s2" = end)

## covariates ##
cat("Reading in covariates\n")
covar_df <- readr::read_delim(paste0(data_dir, "formatted_covariates.txt")) |>
  as.data.frame()
rownames(covar_df) <- covar_df$id
covar_df <- covar_df |> select(-id)
covariates <- SlicedData$new()
covariates$CreateFromMatrix(as.matrix(covar_df))

# make sure column names align
stopifnot(identical(colnames(grex), colnames(pc_ex_mat)))
stopifnot(identical(colnames(grex), colnames(covar_df)))
covariates <- covariates[, qr(covariates)$rank == 1 | !duplicated(qr.Q(qr(covariates)))]


################
## MatrixeQTL ##
################
cat("Starting Matrix eQTL\n")
me <- Matrix_eQTL_main(
  snps = nc_grex,
  gene = pc_ex,
  cvrt = covariates,
  output_file_name = output_file_name_tra, # ??
  pvOutputThreshold = pvOutputThreshold_tra, # ??
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = output_file_name,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = nc_genepos,
  genepos = pc_pos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

saveRDS(me,
        file = paste0(results_dir, "/MatrixeQTL_results.rds"),
        compress = TRUE)
