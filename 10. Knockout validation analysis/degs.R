# load library
library(tidyverse)
library(DESeq2)
library(edgeR)
library(ggplot2)

# load metadata
metadata <- read.csv('metadata.csv', row.names = 1)
metadata$replicate <- factor(metadata$replicate, levels = c(1, 2, 3)) 
metadata$treatment <- factor(metadata$treatment, levels = c('control', 'heat'))
metadata$treatment <- relevel(metadata$treatment, ref = 'control')
metadata$genotype <- factor(metadata$genotype, levels = c('commonbean', 'teparybean'))
metadata

# load data
commonbean_count <- read.delim('result/commonbean_counts.tsv', header = TRUE, sep = '\t', comment.char = '_', row.names = 1)
teparybean_count <- read.delim('result/teparybean_counts.tsv', header = TRUE, sep = '\t', comment.char = '_', row.names = 1)

# QC
common_meta <- metadata[metadata$genotype == "commonbean", ]
tepary_meta <- metadata[metadata$genotype == "teparybean", ]

all(rownames(common_meta) %in% colnames(commonbean_count))
all(rownames(common_meta) == colnames(commonbean_count))


### commonbean

commonbean_count <- round(commonbean_count)

# Constructing a dds dataset
dds <- DESeqDataSetFromMatrix(countData = commonbean_count, colData = common_meta, design = ~treatment)

# Calculate FPM / CPM
fpmcount <- data.frame(fpm(dds, robust = TRUE))
fpmcount <- fpmcount %>% rownames_to_column(var = 'genes')

head(fpmcount)

# Filtering sum of CPM values less than 1 across treatment
fpmcount <- fpmcount %>%
  filter(rowSums(select(., -genes)) > 1)

write.table(fpmcount, 'result/CPM-commonbean.tsv', row.names = FALSE, quote = FALSE, sep = '\t')

# Filtering original count table with the sumCPM value less than 0
commonbean_count <- commonbean_count %>%
  rownames_to_column(var = 'genes') %>%
  filter(genes %in% fpmcount$genes) %>%
  column_to_rownames(var = 'genes')

# Extracting sample names from metadata for extracting from count table
myControl <- common_meta %>%
  filter(treatment == 'control') %>% {rownames(.)}

myHeat <- common_meta %>%
  filter(treatment == 'heat') %>% {rownames(.)}

# DESeq2
dds = DESeq(dds)

result <- data.frame(results(dds, alpha = 0.05, name = 'treatment_heat_vs_control'))
View(result)

expfpm <- fpmcount %>% filter(rowMeans(select(., -genes)) > 1)

# Filter expressed features for dataset
expressgene <- result %>%
  rownames_to_column(var = 'genes') %>%
  filter(genes %in% expfpm$gene)

View(expressgene)

write.table(expressgene, 'result/commonbean-expressed-genes.tsv', row.names = FALSE, quote = FALSE, sep = '\t')

# DEGs

DEGs <- expressgene %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.0)

View(DEGs)

write.table(DEGs, 'result/commonbean-DEGs.tsv', row.names = FALSE, quote = FALSE, sep = '\t')

### teparybean

teparybean_count <- round(teparybean_count)

# Constructing a dds dataset
dds <- DESeqDataSetFromMatrix(countData = teparybean_count, colData = tepary_meta, design = ~treatment)

# Calculate FPM / CPM
fpmcount <- data.frame(fpm(dds, robust = TRUE))
fpmcount <- fpmcount %>% rownames_to_column(var = 'genes')

head(fpmcount)

# Filtering sum of CPM values less than 1 across treatment
fpmcount <- fpmcount %>%
  filter(rowSums(select(., -genes)) > 1)

write.table(fpmcount, 'result/CPM-teparybean.tsv', row.names = FALSE, quote = FALSE, sep = '\t')

# Filtering original count table with the sumCPM value less than 0
teparybean_count <- teparybean_count %>%
  rownames_to_column(var = 'genes') %>%
  filter(genes %in% fpmcount$genes) %>%
  column_to_rownames(var = 'genes')

# Extracting sample names from metadata for extracting from count table
myControl <- tepary_meta %>%
  filter(treatment == 'control') %>% {rownames(.)}

myHeat <- tepary_meta %>%
  filter(treatment == 'heat') %>% {rownames(.)}

# DESeq2
dds = DESeq(dds)

result <- data.frame(results(dds, alpha = 0.05, name = 'treatment_heat_vs_control'))
View(result)

expfpm <- fpmcount %>% filter(rowMeans(select(., -genes)) > 1)

# Filter expressed features for dataset
expressgene <- result %>%
  rownames_to_column(var = 'genes') %>%
  filter(genes %in% expfpm$gene)

View(expressgene)

write.table(expressgene, 'result/teparybean-expressed-genes.tsv', row.names = FALSE, quote = FALSE, sep = '\t')

# DEGs

DEGs <- expressgene %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.0)

View(DEGs)

write.table(DEGs, 'result/teparybean-DEGs.tsv', row.names = FALSE, quote = FALSE, sep = '\t')
