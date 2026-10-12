# Principal-Component-Analysis-PCA-
PCA workflow by PLINK
# Principal Component Analysis (PCA)

## PCA Workflow by PLINK

### The Four Main Steps of Genetic PCA

1. **Quality Control (QC):** Remove low-quality SNPs and individuals to make sure the genetic data is reliable.
2. **Genetic Distance Calculation:** Calculate how genetically similar or different the individuals are based on their SNP data.
3. **Dimensionality Reduction (PCA/MDS):** Reduce the large amount of genetic data into a few main components or dimensions that summarize the main genetic differences between individuals.
4. **Visualization:** Create a plot using PC1 and PC2 to see how individuals cluster together and identify genetic relationships between populations or breeds.

## PCA Demo

```r
# PCA demo
# Clear workspace
rm(list = ls())

# Set working directory to your folder path
setwd("/mnt/E04C56E94C56BA4C/GWAS/practics/practices_plink/AdaptMap_three/")

# Run PLINK QC
system("plink --bfile ADAPTmap_genotypeTOP_20160222_full \
       --cow\
       --autosome\
       --geno 0.1 \
       --mind 0.1 \
       --maf 0.05\
       --nonfounders \
       --allow-no-sex\
       --recode \
       --out ADAPTmap_TOP")

######################################
# Principal Component Analysis - PCA #
######################################

## Genetic distances between individuals
system("plink \
       --cow\
       --allow-no-sex \
       --nonfounders \
       --file ADAPTmap_TOP \
       --distance-matrix \
       --out dataForPCA")

## Load data
dist_populations <- read.table("dataForPCA.mdist", header = F)

### Extract breed names
fam <- data.frame(famids = read.table("dataForPCA.mdist.id")[, 1])

### Extract individual names
famInd <- data.frame(IID = read.table("dataForPCA.mdist.id")[, 2])

## Perform PCA using the cmdscale function
# Time intensive step - takes a few minutes with the 4.5K animals
mds_populations <- cmdscale(dist_populations, eig = T, 5)

## Extract the eigen vectors
eigenvec_populations <- cbind(fam, famInd, mds_populations$points)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig) / sum(mds_populations$eig)) * 100, 2)

# Visualize PCA in tidyverse
# Load tidyverse
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

# PCA plot
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `1`, y = `2`, color = famids), show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = "PCA of worldwide goat populations",
       x = paste0("Principal component 1 (", eigen_percent[1], " %)"),
       y = paste0("Principal component 2 (", eigen_percent[2], " %)")) +
  theme_minimal()

ggsave(filename = "pca_plot.png", width = 8, height = 6, dpi = 300)
```
