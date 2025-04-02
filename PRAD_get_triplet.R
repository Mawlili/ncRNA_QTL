# get_triplets.R
# Make list with triplets: (SNP, distal protein RNA, local nc-RNAs)
# Nolan Cole, UW Biostat, 15 January 2025

# libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, readr)

# Read in:
# distal protein-coding results
distal_pc <- readr::read_delim("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/trans_qtl/significant_output_trans_qtl_protein_coding.txt.hits.txt.gz",
                               delim = " ", col_names = FALSE)
colnames(distal_pc) <- c(
  "gene_id",    # Phenotype (gene) ID
  "chr",   # Phenotype chrID
  "start", # Phenotype start
  "snp_id",  # Variant ID
  "variant_chr", # Variant chrID
  "variant_pos", # Variant position
  "p_value",     # Nominal P-value of association
  "dummy",       # Dummy field used in approximated mapping in trans
  "slope"        # Regression slope
)
distal_pc <- distal_pc[,c("gene_id", "snp_id")]

# local non-coding results
local_nc <- readr::read_delim("/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/cis_qtl/cis_qtl_short_cov_overlap_V250.txt",
                              col_names = FALSE)[,c(1, 8)]
# colnames(local_nc) <- c("gene_id", "chr", "start", "end",
#                         "strand", "total_num_tested", "distance",
#                         "snp_id", "chr_snp", "start_snp", "end_snp",
#                         )
colnames(local_nc) <- c("gene_id", "snp_id")

# Subset to SNPs that are associated
# with at least one local nc and distal pc
snp_intersect <- intersect(distal_pc$snp_id, local_nc$snp_id)
local_nc_intersect <- local_nc |>
  filter(snp_id %in% snp_intersect)
# Should see no change in distal_pc
# since we subset the VCF to only be SNPs
# that are associated with local nc RNAs
# but just in case
distal_pc_intersect <- distal_pc |>
  filter(snp_id %in% snp_intersect)

# Get Triplets
# triplet := { single snp, single local pc gene, distal nc gene(s) }
triplets <- distal_pc_intersect |>
  full_join(local_nc_intersect,
            by = "snp_id",
            relationship = "many-to-many") |>
  rename(
    distal_pc = gene_id.x,
    local_nc = gene_id.y
  ) |>
  select(snp_id, distal_pc, local_nc) |>
  group_by(snp_id, distal_pc) |>
  arrange(local_nc) |>
  ungroup()

# Count number of local NCs associated with
# snp-distal PC pairs
# Expect up to hundreds of local ncs in a single triplet
triplet_counts <- triplets |>
  group_by(snp_id, distal_pc) |>
  summarise(local_nc_count = n(),
            .groups = "drop") |>
  arrange(desc(local_nc_count))

# Save triplets
saveRDS(triplets, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/triplets.RDS")
saveRDS(triplet_counts, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_qtl/triplet_counts.RDS")
