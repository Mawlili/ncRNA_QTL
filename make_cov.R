library(data.table)
library(MASS)
library(pracma)
library(dplyr)

#loading data
load("/rsrch5/home/epi/bhattacharya_lab/users/whwu1/temp/hidden_cov/hcp.RData")
cov1 <- fread('/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/merged_cov_with_barcodes.tsv')
sex <- fread('/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/sex_data.tsv')
age <- fread('/rsrch5/home/epi/bhattacharya_lab/data/TCGA/BRCA/age_data.tsv')

#remove character
cov <- as.data.frame(cov)
rownames(cov) <- cov[[75]]
cov <- cov[, -c(74, 75)]
cov <- cov[, -(1:5)]

#PCA
cov_matrix <- as.matrix(cov)
cov_matrix <- apply(cov_matrix, 2, as.numeric)
pca_result <- prcomp(cov_matrix)
rownames(pca_result$x) <- rownames(cov)
pca_data <- as.data.frame(pca_result$x[, 1:3])  

#merge together pca and hcp
cov$SampleID <- rownames(cov)
pca_data$SampleID <- rownames(pca_data)
hcp_z <- as.data.frame(hcp$Z)
hcp_z$SampleID <- rownames(hcp_z)
merged_cov <- merge(cov, pca_data, by = "SampleID")
merged_cov <- merge(merged_cov, hcp_z, by = "SampleID")

#merge sex and age information
merged_cov$SampleID_short <- sub("^((TCGA-[A-Z0-9]+-[A-Z0-9]+)).*$", "\\1", merged_cov$SampleID)
colnames(age) <- c("SampleID_short", "age_at_index")
colnames(sex) <- c("SampleID_short", "SEX")
merged_cov <- merged_cov %>%
  left_join(age, by = "SampleID_short")
merged_cov <- merged_cov %>%
  left_join(sex, by = "SampleID_short")

#export
cov_export <- as.data.frame(t(merged_cov))
cov_export <- cbind(id = rownames(cov_export), cov_export)
write.table(cov_export, "/rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/qtl/merged_cov_full_barcode.txt", sep = "\t", row.names = FALSE, quote = FALSE)
