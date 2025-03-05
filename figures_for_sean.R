library(ggplot2)
library(dplyr)

file_path <- "/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected/ESPRESSO_corrected_classification.txt"

# Set Up
classification_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
myPalette <- c("#6BAED6", "#FC8D59", "#78C679", "#EE6A50", "#969696",
                "#66C2A4", "goldenrod1","darksalmon","#41B6C4")


#Density histogram of distributions of transcript lengths, grouped by structural category

ggplot(classification_data, aes(x = length, fill = structural_category)) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.6, position = "identity") +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  scale_fill_manual(values = myPalette) +
  labs(title = "Density Histogram of Transcript Lengths",
       x = "Transcript Length (bp)", y = "Density")


#Density histogram of distributions of exons per transcript, grouped by structural category
#Density histogram of distributions of exon lengths, grouped by coding potential
#Bar chart of ORF support (whether the isoform has a predicted ORF) by structural category
