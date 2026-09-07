# Load library
library(tidyverse)

# load dataset
motif <- read.delim('/mmfs1/thunder/projects/zhikai.liang/Fahad/triticale/20-motif/tomtom-result/arabidopsis_cluster_motif_info.txt',
                    sep = ' ', header = FALSE)
colnames(motif) <- c("cluster", "motif_sequence", "arabidopsis_tf_id", "tf_family")
motif$cluster <- gsub('cluster', '', motif$cluster)
motif <- motif %>% mutate(cluster = as.numeric(cluster)) %>% arrange(cluster)
head(motif)
View(motif)

write.table(motif, 'result/arabidopsis-TF.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
