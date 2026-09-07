# Load library
library(tidyverse)

## Load data and QC
# Syntenic Table (for 1000bp - gene IDs)
syntenic <- read.delim('data/fig3/Durum_Rye_ABR_unique_triplet.trt.restrict.txt', header = FALSE)
colnames(syntenic) <- c('A', 'B', 'R')

as_tibble(syntenic)
as_tibble(cluster)

# Cluster Table
cluster <- read.delim('./data/fig3/Rerun/Cont-v-LBR-v-PR-cluster-genelist-k-16.tsv', sep = '\t', header = TRUE)

# Join
syntenic_cluster <- syntenic %>%
  mutate(
    A_gene = sub("\\.\\d+$", "", A),
    B_gene = sub("\\.\\d+$", "", B),
    R_gene = sub("\\.\\d+$", "", R)
  ) %>%
  left_join(cluster, by = c('A_gene' = 'Feature')) %>%
  rename(Cluster_A = Cluster) %>%
  left_join(cluster, by = c('B_gene' = 'Feature')) %>%
  rename(Cluster_B = Cluster) %>%
  left_join(cluster, by = c('R_gene' = 'Feature')) %>%
  rename(Cluster_R = Cluster) %>%
  drop_na(Cluster_A, Cluster_B, Cluster_R) %>%
  select(-A_gene, -B_gene , -R_gene)

syntenic_cluster %>% View()

# Define subgenome categories
categories <- list(
  'ABR_same' = list(
    condition = quote(Cluster_A == Cluster_B & Cluster_B == Cluster_R)
  ),
  'ABR_diff' = list(
    condition = quote(Cluster_A != Cluster_B & Cluster_B != Cluster_R & Cluster_A != Cluster_R)
  ),
  'AB_same_R_diff' = list(
    condition = quote(Cluster_A == Cluster_B & Cluster_A != Cluster_R)
  ),
  'AR_same_B_diff' = list(
    condition = quote(Cluster_A == Cluster_R & Cluster_A != Cluster_B)
  ),
  'BR_same_A_diff' = list(
    condition = quote(Cluster_B == Cluster_R & Cluster_B != Cluster_A)
  )
)

# Output directory
outdir <- '/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/30-mafft-rerun/output/'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (category in names(categories)) {
  # Filter data based on category
  filter_data <- syntenic_cluster %>%
    filter(eval(categories[[category]]$condition))
  
  if (nrow(filter_data) > 0) {
    subgenomes <- c('A', 'B', 'R')
    
    for (subgenome in subgenomes) {
      # 1000bp version
      bp1000_list <- filter_data %>%
        select(all_of(subgenome)) %>%
        distinct()
      
      bp1000_file <- file.path(outdir, paste0(category, '_', subgenome, '.txt'))
      write.table(bp1000_list, bp1000_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    cat(category, ":", nrow(filter_data), "triplets\n")
  } 
}








































