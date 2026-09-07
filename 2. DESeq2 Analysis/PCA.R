## Fig 1(b) - PCA plot

# Load library
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(ggpubr)
library(svglite)


## Load dataset
# 1. Metadata
meta = read.csv('fig1/metadata.csv', row.names = 1)
meta$Replicate <- factor(meta$Replicate)
meta$Treatment <- factor(meta$Treatment)
meta$Timepoint <- factor(meta$Timepoint)
meta$Treatment <- relevel(meta$Treatment, ref = 'Control')

# 2. Count data
count <- read.delim('fig1/count_table_kallisto_tximport.txt', sep = ' ', header = TRUE)
colnames(count) <- c("Cont_124", "Cont_148", "Cont_172", "Cont_196", "Cont_224", "Cont_248", "Cont_272", "Cont_296", "Cont_324", "Cont_348", "Cont_372", "Cont_396",
                     "LBR124", "LBR148", "LBR172", "LBR196", "LBR224", "LBR248", "LBR272", "LBR296", "LBR324", "LBR348", "LBR372", "LBR396",
                     "PR124", "PR148", "PR172", "PR196", "PR224", "PR248", "PR272", "PR296", "PR324", "PR348", "PR372", "PR396")
# QC
meta <- meta[colnames(count), ] # Rearranging based on count matrix

## DESeq2
count <- round(count)
dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~Treatment)

### Transform to log2(CPM+1) for visualization
counts_mat <- counts(dds)
cpm_mat <- cpm(counts_mat, log = FALSE)

# Log2(CPM + 1)
log2cpm_mat <- log2(cpm_mat + 1)

## Compute PCA
pca <- prcomp(t(log2cpm_mat))

## Extract first two principle components
pca_data <- as.data.frame(pca$x[, 1:2])
pca_data$Sample <- rownames(pca_data) # Add rownames as the Sample column
pca_data$Treatment <- colData(dds)$Treatment # Add treatment as Treatment column
pca_data$Timepoint <- colData(dds)$Timepoint
pca_data

# Calculate variance explained for axis labels
var_explained <- summary(pca)$importance[2, 1:2] * 100
var_explained
summary(pca)

## Data checking
table(pca_data$Treatment, pca_data$Timepoint)

## PCA Plot
pca_treat <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment, shape = Timepoint, fill = Treatment)) +
  geom_point(size = 3, stroke = 1.2) +
  #stat_ellipse(aes(group = Treatment, fill = Treatment, color = Treatment), type = 'norm', alpha = 0.2, geom = 'polygon') +
  scale_shape_manual(values = c(21, 24, 22, 23)) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
  ) +
  theme_pubr(base_size = 16) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_brewer(palette = 'Dark2') +
  scale_fill_brewer(palette = 'Pastel2')

pca_treat
ggsave(pca_treat, file = 'fig1/figs/pca.svg', height = 5, width = 7, dpi = 600)










