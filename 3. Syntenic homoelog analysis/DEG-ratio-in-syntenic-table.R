# load library
library(tidyverse)

# List all directories
samples <- list.files(path = './data/DEGs', full.names = T)
files <- file.path(samples)
names(files) <- str_replace(samples, './data/DEGs/DEGs_', '') %>%
  str_replace('_v_Cont_', '') %>%
  str_replace('.txt', '')

names(files)

# Load dataset
deg_list <- lapply(files, read.delim, sep = ' ', header = TRUE)
triplet_list <- read.delim('./data/Durum_Rye_ABR_unique_triplet.trt.restrict.txt', header = FALSE)
colnames(triplet_list) <- c('A', 'B', 'R')
triplet_list <- triplet_list %>%
  mutate(across(everything(), ~ gsub("\\.\\d+$", "", .x)))

head(deg_list)
head(triplet_list)

triplet_list_unique <- triplet_list %>%
  add_count(A, name = 'A_count') %>%
  add_count(B, name = 'B_count') %>%
  add_count(R, name = 'R_count') %>%
  filter(!(A_count > 1 | B_count > 1 | R_count > 1)) %>%
  select(-A_count, -B_count, -R_count)

combine <- triplet_list_unique %>%
  pivot_longer(cols = everything(), values_to = 'gene') %>%
  select(gene)



lbr24_triad <- length(intersect(combine$gene, deg_list$LBR24$Feature))

deg_count <- data.frame(
  LBR24 = length(deg_list$LBR24$Feature),
  LBR48 = length(deg_list$LBR48$Feature),
  LBR72 = length(deg_list$LBR72$Feature),
  LBR96 = length(deg_list$LBR96$Feature),
  PR24 = length(deg_list$PR24$Feature),
  PR48 = length(deg_list$PR48$Feature),
  PR72 = length(deg_list$PR72$Feature),
  PR96 = length(deg_list$PR96$Feature)
)

deg_count

triad_ratio <- data_frame(
  'LBR24' = (length(intersect(combine$gene, deg_list$LBR24$Feature)) /
    length(deg_list$LBR24$Feature)) * 100,
  'LBR48' = (length(intersect(combine$gene, deg_list$LBR48$Feature)) /
               length(deg_list$LBR48$Feature)) * 100,
  'LBR72' = (length(intersect(combine$gene, deg_list$LBR72$Feature)) /
               length(deg_list$LBR72$Feature)) * 100,
  'LBR96' = (length(intersect(combine$gene, deg_list$LBR96$Feature)) /
               length(deg_list$LBR96$Feature)) * 100,
  'PR24' = (length(intersect(combine$gene, deg_list$PR24$Feature)) /
               length(deg_list$PR24$Feature)) * 100,
  'PR48' = (length(intersect(combine$gene, deg_list$PR48$Feature)) /
               length(deg_list$PR48$Feature)) * 100,
  'PR72' = (length(intersect(combine$gene, deg_list$PR72$Feature)) /
               length(deg_list$PR72$Feature)) * 100,
  'PR96' = (length(intersect(combine$gene, deg_list$PR96$Feature)) /
               length(deg_list$PR96$Feature)) * 100,
  
)

triad_ratio
min(triad_ratio)
max(triad_ratio)
