# make_snp_dosage.R
# make SNP dosage matrix
# Nolan Cole, UW Biostat, 31 Jan 2025

# make_bed.R
# Format protein coding and non-protein coding .bed files for QTLtools
# https://qtltools.github.io/qtltools/
# Nolan Cole, UW Biostatistics, 13 September 2024

#bigsnpr.bed obtained from running plink make bed from .vcf 
# Install/Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(pacman)
# BiocManager::install(c("DESeq2", "SummarizedExperiment", "bigsnpr", "biomaRt"))
pacman::p_load(dplyr, magrittr, tidyr,
               tibble, ggplot2, janitor,
               DESeq2, SummarizedExperiment,
               bigsnpr, biomaRt, readr)
# read in data
snp_data_file <- bigsnpr::snp_readBed("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/bigsnpr/bigsnpr.bed")
snp_data <- readRDS(snp_data_file)
covar <- readr::read_tsv("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/cov_remove_duplicate.txt",
                         col_names = TRUE)[,-1]
#load(biomart_file)

# Subset to patients we covariates for
covar_ids <- colnames(covar)
sample_intersect_index <- which(snp_data$fam$sample.ID %in% covar_ids)
sample_intersection <- snp_data$fam$sample.ID[sample_intersect_index]
sub_snp_data <- snp_data$genotypes[sample_intersect_index,]

# name SNP locations 
colnames(sub_snp_data) <- snp_data$map$marker.ID
rownames(sub_snp_data) <- sample_intersection

# Save SNP dosage
saveRDS(sub_snp_data,
        file = "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/out_files/snp_dosage.RDS",
        compress = TRUE)
