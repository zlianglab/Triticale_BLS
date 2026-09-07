# Loading necessary packages
library(tidyverse)
library(DESeq2)
library(edgeR)

# Metadata
meta <- read.csv('./meta/metadata.csv', row.names = 1)
meta$Replicate <- factor(meta$Replicate, levels = c(1, 2, 3))
meta$Timepoint <- factor(meta$Timepoint, levels = c(24, 48, 72, 96))
meta$Treatment <- factor(meta$Treatment)
meta$Treatment <- relevel(meta$Treatment, ref = 'Control')

# Raw read count table
count <- read.delim('./result/count_table_kallisto_tximport.txt', sep = ' ', header = TRUE)
colnames(count) <- c("Cont_124", "Cont_148", "Cont_172", "Cont_196", "Cont_224", "Cont_248", "Cont_272", "Cont_296", "Cont_324", "Cont_348", "Cont_372", "Cont_396",
                     "LBR124", "LBR148", "LBR172", "LBR196", "LBR224", "LBR248", "LBR272", "LBR296", "LBR324", "LBR348", "LBR372", "LBR396",
                     "PR124", "PR148", "PR172", "PR196", "PR224", "PR248", "PR272", "PR296", "PR324", "PR348", "PR372", "PR396")

count %>% View()

# QC
meta <- meta[colnames(count),]
all(rownames(meta) %in% colnames(count)) # checking if metadata and count data rowname matches
all(rownames(meta) == colnames(count)) 

count <- round(count) # rounding up the data for DESeq2

# Constructing a DESeq Dataset
dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~Treatment)

# Calculate FPM / CPM
fpmcount <- data.frame(fpm(dds, robust = TRUE))
fpmcount <- fpmcount %>% rownames_to_column(var = 'Feature')

fpmcount %>% View()

# Filtering sum of CPM values less than 0 across treatment
fpmcount <- fpmcount %>%
  filter(rowSums(select(., -Feature)) > 0)

write.table(file = 'result/CPM_total_no_filtering.txt', fpmcount, quote = FALSE, row.names = FALSE)

# Filtering original count table with the sumCPM value less than 0
count <- count %>%
  rownames_to_column(var = 'Feature') %>%
  filter(Feature %in% fpmcount$Feature) %>%
  column_to_rownames(var = 'Feature')

# Extracting sample names from metadata for extracting from count table
myControl <- meta %>%
  filter(Treatment == 'Control') %>% {rownames(.)}

myLBR <- meta %>%
  filter(Treatment == 'LBR') %>% {rownames(.)}

myPR <- meta %>%
  filter(Treatment == 'PR') %>% {rownames(.)}

# Subsetting with LBR and PR
LBR_v_PRfpm <- fpmcount %>%
  select(Feature, all_of(myLBR), all_of(myPR))
View(LBR_v_PRfpm)

#### Getting data for different timepoint

## LBR vs. PR-------------------------------------------------------------------

for (timepoint in unique(meta$Timepoint)) {
 # Organizing metadata
 meta_subset = meta %>% filter(Timepoint == timepoint & Treatment != 'Control')
 meta_subset$Treatment <- factor(meta_subset$Treatment)
 meta_subset$Treatment <- relevel(meta_subset$Treatment, ref = 'LBR')

 # Organizing count data
 sample_subset = rownames(meta_subset)
 count_subset = count[, sample_subset, drop = FALSE]

 # Create DESEq2 Matrix from the dataset
 dds = DESeqDataSetFromMatrix(countData = count_subset, colData = meta_subset, design = ~Treatment)
 dds = DESeq(dds)

 # Perform pairwise differential expression for PR vs LBR
 result <- data.frame(results(dds, alpha = 0.05, name = 'Treatment_PR_vs_LBR'))

 # Calculate FPM / CPM for the count subset
 fpmcount_24 <- data.frame(fpm(dds, robust = TRUE))
 expfpm <- fpmcount_24[rowMeans(fpmcount_24) > 1, ] %>%
  rownames_to_column(var = 'Feature')

 # Filter expressed features for dataset
 expressgene <- result %>%
  rownames_to_column(var = 'Feature') %>%
  filter(Feature %in% expfpm$Feature)

 write.table(file = paste0('./result/PR_v_LBR/Expressed_genes_PR_v_LBR_', timepoint, '.txt'), expressgene, sep = ' ', row.names = FALSE)

 # Filter DEGs from dataset
 DEGs <- expressgene %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.0)

 write.table(file = paste0('./result/PR_v_LBR/DEGs_PR_v_LBR_', timepoint, '.txt'), DEGs, sep = ' ', row.names = FALSE)

} 
 
## Treatment vs. Control--------------------------------------------------------

# Organizing metadata
for (timepoint in unique(meta$Timepoint)) {
  
  meta_subset = meta %>% filter(Timepoint == timepoint & Treatment != 'LBR')
  meta_subset$Treatment <- factor(meta_subset$Treatment)

  meta_subset$Treatment <- relevel(meta_subset$Treatment, ref = 'Control')

 # Organizing count data
 sample_subset = rownames(meta_subset)
 count_subset = count[, sample_subset, drop = FALSE]

 # Create DESEq2 Matrix from the dataset
 dds = DESeqDataSetFromMatrix(countData = count_subset, colData = meta_subset, design = ~Treatment)
 dds = DESeq(dds)

 # Perform pairwise differential expression for Treatment vs Cont
 result <- data.frame(results(dds, alpha = 0.05, name = 'Treatment_PR_vs_Control'))

 # Calculate FPM / CPM for the count subset
 fpmcount_24 <- data.frame(fpm(dds, robust = TRUE))
 expfpm <- fpmcount_24[rowMeans(fpmcount_24) > 1, ] %>%
  rownames_to_column(var = 'Feature')

 # Filter expressed features for dataset
 expressgene <- result %>%
  rownames_to_column(var = 'Feature') %>%
  filter(Feature %in% expfpm$Feature)

 write.table(file = paste0('./result/PR_v_Cont/Expressed_genes_PR_v_Cont_', timepoint, '.txt'), expressgene, sep = ' ', row.names = FALSE)

 # Filter DEGs from dataset
 DEGs <- expressgene %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.0)

 write.table(file = paste0('./result/PR_v_Cont/DEGs_PR_v_Cont_', timepoint, '.txt'), DEGs, sep = ' ', row.names = FALSE)

}
 
# Classify DEGs
classify_deg <- function(df) {
  df %>%
    transmute(
      Feature = Feature,
      Expression = case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange >= 1 ~ 'upregulated',
        !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ 'downregulated',
        TRUE ~ 'expressed' # All the other genes
      )
    )
}

degcount <- classify_deg(DEGs)
print(table(degcount$Expression))

