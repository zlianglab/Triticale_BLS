# Load library
library(tidyverse)

# LBR v Cont Expressed Genes
lbr24 <- read.delim('./data/LBR_v_Cont/Expressed_genes_LBR_v_Cont_24.txt', sep = ' ', header = TRUE) 
lbr48 <- read.delim('./data/LBR_v_Cont/Expressed_genes_LBR_v_Cont_48.txt', sep = ' ', header = TRUE)
lbr72 <- read.delim('./data/LBR_v_Cont/Expressed_genes_LBR_v_Cont_72.txt', sep = ' ', header = TRUE)
lbr96 <- read.delim('./data/LBR_v_Cont/Expressed_genes_LBR_v_Cont_24.txt', sep = ' ', header = TRUE)


# PR v Cont Expressed Genes
pr24 <- read.delim('./data/PR_v_Cont/Expressed_genes_PR_v_Cont_24.txt', sep = ' ', header = TRUE) 
pr48 <- read.delim('./data/PR_v_Cont/Expressed_genes_PR_v_Cont_48.txt', sep = ' ', header = TRUE)
pr72 <- read.delim('./data/PR_v_Cont/Expressed_genes_PR_v_Cont_72.txt', sep = ' ', header = TRUE)
pr96 <- read.delim('./data/PR_v_Cont/Expressed_genes_PR_v_Cont_24.txt', sep = ' ', header = TRUE)

# Combined DEG
expgene <- Reduce(union, list(lbr24, lbr48, lbr72, lbr96, pr24, pr48, pr72, pr96))
expgene %>% View()

expgene <- expgene %>%
  distinct(Feature, .keep_all = TRUE)

background <- expgene %>%
  filter(abs(log2FoldChange) < 1 | pvalue > 0.05)

background %>% View()

bg_list <- background %>%
  select(Feature)

write.table(file = './result/bg_list.txt', bg_list, row.names = FALSE, col.names = FALSE, quote = FALSE)

