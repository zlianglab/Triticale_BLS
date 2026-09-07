# Load library
library(tidyverse)
library(gridExtra)
library(svglite)
library(ggtext)

# Load dataset
# TAL Templace
template <- c(1:8)

for (temp in template) {
  assign(
    paste0('p3.tal', temp),
    read.delim(paste0('data/P3_template/Predictions_for_P3-tempTALE', temp, '.tsv'), header = TRUE, sep = '\t')
  )
}

cpm <- read.delim('data/CPM_total_no_filtering.txt', header = TRUE, sep = ' ')

# # PR v Cont
# timepoint <- c('24', '48', '72', '96')
# 
# for (time in timepoint) {
#   assign(
#     paste0("pr.cont", time),
#     read.delim(
#       paste0("data/PR_v_Cont/Expressed_genes_PR_v_Cont_", time, ".txt"),
#       header = TRUE, sep = " "
#     )
#   )
# }
# 
# # PR v LBR
# for (time in timepoint) {
#   assign(
#     paste0("pr.lbr", time),
#     read.delim(
#       paste0("data/PR_v_Cont/Expressed_genes_PR_v_Cont_", time, ".txt"),
#       header = TRUE, sep = " "
#     )
#   )
# }


# QC
p3.tal1$Sequence.ID <- gsub(':.*', '', p3.tal1$Sequence.ID)
p3.tal2$Sequence.ID <- gsub(':.*', '', p3.tal2$Sequence.ID)
p3.tal3$Sequence.ID <- gsub(':.*', '', p3.tal3$Sequence.ID)
p3.tal4$Sequence.ID <- gsub(':.*', '', p3.tal4$Sequence.ID)
p3.tal5$Sequence.ID <- gsub(':.*', '', p3.tal5$Sequence.ID)
p3.tal6$Sequence.ID <- gsub(':.*', '', p3.tal6$Sequence.ID)
p3.tal7$Sequence.ID <- gsub(':.*', '', p3.tal7$Sequence.ID)
p3.tal8$Sequence.ID <- gsub(':.*', '', p3.tal8$Sequence.ID)

p3.tal1 <- p3.tal1 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
p3.tal2 <- p3.tal2 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
p3.tal3 <- p3.tal3 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
p3.tal4 <- p3.tal4 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
p3.tal5 <- p3.tal5 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
p3.tal6 <- p3.tal6 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
p3.tal7 <- p3.tal7 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
p3.tal8 <- p3.tal8 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)

tal_names <- paste0("tal", 1:8)
tp_names  <- c("24 hpi", "48 hpi", "72 hpi", "96 hpi")

# Two aligned matrices: same rows (templates), same cols (timepoints)
logfc_mat <- matrix(NA_real_, 8, 4, dimnames = list(tal_names, tp_names))
pval_mat  <- matrix(NA_real_, 8, 4, dimnames = list(tal_names, tp_names))

safe_t <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || sd(x) == 0) return(NA_real_)
  t.test(x, mu = 0)$p.value
}
safe_mean <- function(x) mean(x[is.finite(x)])  # matches safe_t's filtering

for (temp in 1:8) {
  tal <- get(paste0("p3.tal", temp))
  
  P3 <- cpm %>%
    filter(Feature %in% tal$Sequence.ID) %>%
    mutate(
      Cont24 = (Cont_124 + Cont_224 + Cont_324) / 3,
      Cont48 = (Cont_148 + Cont_248 + Cont_348) / 3,
      Cont72 = (Cont_172 + Cont_272 + Cont_372) / 3,
      Cont96 = (Cont_196 + Cont_296 + Cont_396) / 3,
      PR24 = (PR124 + PR224 + PR324) / 3,
      PR48 = (PR148 + PR248 + PR348) / 3,
      PR72 = (PR172 + PR272 + PR372) / 3,
      PR96 = (PR196 + PR296 + PR396) / 3
    ) %>%
    mutate(
      H24 = log2((PR24 + 1) / (Cont24 + 1)),
      H48 = log2((PR48 + 1) / (Cont48 + 1)),
      H72 = log2((PR72 + 1) / (Cont72 + 1)),
      H96 = log2((PR96 + 1) / (Cont96 + 1))
    ) %>%
    select(Feature, H24, H48, H72, H96)
  
  H <- list(P3$H24, P3$H48, P3$H72, P3$H96)
  
  cat("template", temp, "(n =", nrow(P3), "genes)\n")
  
  logfc_mat[temp, ] <- vapply(H, safe_mean, numeric(1))
  pval_mat[temp, ]  <- vapply(H, safe_t,    numeric(1))
  
}

logfc_mat
pval_mat

write.table(logfc_mat, 'result/fig4c/PR_v_Cont_logFC.tsv', quote = FALSE, row.names = FALSE, sep = '\t')
