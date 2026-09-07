# Load library
library(tidyverse)

# Load dataset
syntenic <- read.delim('./data/fig3/Durum_Rye_ABR_unique_triplet.trt.restrict.txt', header = FALSE)
colnames(syntenic) <- c('A', 'B', 'R')
View(syntenic)
head(syntenic)

syntenic$A <- gsub("\\.\\d+$", "", syntenic$A)
syntenic$B <- gsub("\\.\\d+$", "", syntenic$B)
syntenic$R <- gsub("\\.\\d+$", "", syntenic$R)

cluster <- read.delim('./data/fig3/Rerun/Cont-v-LBR-v-PR-cluster-genelist-k-16.tsv', sep = '\t', header = TRUE)
View(cluster)

# longestTrt <- read.delim('data/fig3/triticale.longestTrt.txt', header = FALSE)
# tail(longestTrt)
# longestTrt$V1 <- gsub('gene:', '', longestTrt$V1)
# longestTrt$V2 <- gsub('transcript:', '', longestTrt$V2)

# Joing cluster with syntenic table
syntenic_cluster <- syntenic %>%
  left_join(cluster, by = c('A' = 'Feature')) %>%
  rename(cluster_A = Cluster) %>%
  left_join(cluster, by = c('B' = 'Feature')) %>%
  rename(cluster_B = Cluster) %>%
  left_join(cluster, by = c('R' = 'Feature')) %>%
  rename(cluster_R = Cluster) %>%
  drop_na(cluster_A, cluster_B, cluster_R)

View(syntenic_cluster)

# ratio <- data.frame()
# 
# for (cluster in 1:16) {
#   df <- cluster %>%
#     filter(Cluster == cluster) %>%
#     mutate(
#       ratio = Cluster %in% 
#     )
# }

# Find out the categories
# A, B, R same
abr.same <- syntenic_cluster %>%
  filter(cluster_A == cluster_B & cluster_B == cluster_R) %>% select(A, B, R)

write.table(abr.same, file = '/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/30-mafft-rerun/genelist/ABR_same.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# A, B, R different
abr.diff <- syntenic_cluster %>%
  filter(cluster_A != cluster_B & cluster_B != cluster_R & cluster_A != cluster_R) %>% select(A, B, R)
write.table(abr.diff, file = '/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/30-mafft-rerun/genelist/ABR_diff.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# A, B same, R different
ab.same.r.diff <- syntenic_cluster %>%
  filter(cluster_A == cluster_B & cluster_B != cluster_R & cluster_A != cluster_R) %>% select(A, B, R)
write.table(ab.same.r.diff, file = '/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/30-mafft-rerun/genelist/AB_same_R_diff.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# A, R same, B different
ar.same.b.diff <- syntenic_cluster %>%
  filter(cluster_A == cluster_R & cluster_B != cluster_A & cluster_B != cluster_R) %>% select(A, B, R)
write.table(ar.same.b.diff, file = '/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/30-mafft-rerun/genelist/AR_same_B_diff.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# B, R same, A different
br.same.a.diff <- syntenic_cluster %>%
  filter(cluster_B == cluster_R & cluster_B != cluster_A & cluster_A != cluster_R) %>% select(A, B, R)
write.table(br.same.a.diff, file = '/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/30-mafft-rerun/genelist/BR_same_A_diff.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# Filtered list
filter_list <- syntenic_cluster %>%
  filter(cluster_A == cluster_B & cluster_B != cluster_R & cluster_A != cluster_R) %>%
  select(A, B, R, cluster_A) %>%
  rename(cluster = cluster_A) %>%
  arrange(cluster)

filter_list %>% View()

write.table(filter_list['R'], '/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/12_mafft/ABR-diff/R.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# Create lookup vector
lookup <- setNames(longestTrt$V2, longestTrt$V1)

# Filter and replace in one step
genelist <- filter_list %>%
  filter(A %in% names(lookup) | B %in% names(lookup) | R %in% names(lookup)) %>%
  mutate(
    A = lookup[A],
    B = lookup[B], 
    R = lookup[R]
  ) %>% select(-cluster)

View(genelist)

write.table(genelist['R'], file = '/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/12_mafft/ABR-diff/R_cds.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)




