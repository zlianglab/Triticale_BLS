# Load libraries
library(tidyverse)

# Load the CPM file generated from DESeq2
cpm <- read.delim('data/CPM_total_kallisto_tximport.txt', sep = ' ', header = TRUE)
View(cpm)

genes <- cpm$Feature %>% as_data_frame()

genes %>% View()
genes %>% as_tibble()
colnames(genes) <- 'Feature'

write.table(genes, './result/cpm_background_list.txt', sep = ' ', row.names = FALSE, col.names = FALSE, quote = FALSE)
