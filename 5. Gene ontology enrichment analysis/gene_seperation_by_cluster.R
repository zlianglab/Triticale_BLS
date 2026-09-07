# Load library
library(tidyverse)

# load dataset

#----------------------------------------LBR------------------------------------
lbr <- read.delim('./data/gene_grouping_by_clusters/LBR-gene-group-by-cluster-16.txt', sep = ' ', header = TRUE)
lbr %>% View()
as_tibble(lbr)
colnames(lbr)

# Loop through each columns and save seperately
for(i in 1:ncol(lbr)) {
  filename <- paste0('./result/lbr/', colnames(lbr)[i], '.txt')
  
  # Omit NA from each column before saving
  temp <- na.omit(lbr[,i])
  
  # Write the column to a file
  write.table(temp, filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#----------------------------------------PR-------------------------------------
pr <- read.delim('./data/gene_grouping_by_clusters/PR-gene-group-by-cluster-16.txt', sep = ' ', header = TRUE)
View(pr)

# Loop through each columns and save seperately
for(i in 1:ncol(pr)) {
  filename <- paste0('./result/pr/', colnames(pr)[i], '.txt')
  
  # Omit NA from each column before saving
  temp <- na.omit(pr[,i])
  
  # Write the column to a file
  write.table(temp, filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#----------------------------------------PR v LBR-------------------------------
pr <- read.delim('./data/Rerun/Cont_v_LBR_PR-gene-group-by-cluster-16.txt', sep = ' ', header = TRUE)
View(pr)

# Loop through each columns and save seperately
for(i in 1:ncol(pr)) {
  filename <- paste0('./result/Rerun/cont_v_lbr_pr_', colnames(pr)[i], '.txt')
  
  # Omit NA from each column before saving
  temp <- na.omit(pr[,i])
  
  # Write the column to a file
  write.table(temp, filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
