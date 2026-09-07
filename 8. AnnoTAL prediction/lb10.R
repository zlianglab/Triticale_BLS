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
    paste0('lb10.tal', temp),
    read.delim(paste0('data/LB10_template/Predictions_for_LB10-tempTALE', temp, '.tsv'), header = TRUE, sep = '\t')
  )
}

cpm <- read.delim('data/CPM_total_no_filtering.txt', header = TRUE, sep = ' ')

head(lb10.tal1)

# QC
lb10.tal1$Sequence.ID <- gsub(':.*', '', lb10.tal1$Sequence.ID)
lb10.tal2$Sequence.ID <- gsub(':.*', '', lb10.tal2$Sequence.ID)
lb10.tal3$Sequence.ID <- gsub(':.*', '', lb10.tal3$Sequence.ID)
lb10.tal4$Sequence.ID <- gsub(':.*', '', lb10.tal4$Sequence.ID)
lb10.tal5$Sequence.ID <- gsub(':.*', '', lb10.tal5$Sequence.ID)
lb10.tal6$Sequence.ID <- gsub(':.*', '', lb10.tal6$Sequence.ID)
lb10.tal7$Sequence.ID <- gsub(':.*', '', lb10.tal7$Sequence.ID)
lb10.tal8$Sequence.ID <- gsub(':.*', '', lb10.tal8$Sequence.ID)

lb10.tal1 <- lb10.tal1 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
lb10.tal2 <- lb10.tal2 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
lb10.tal3 <- lb10.tal3 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
lb10.tal4 <- lb10.tal4 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
lb10.tal5 <- lb10.tal5 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
lb10.tal6 <- lb10.tal6 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
lb10.tal7 <- lb10.tal7 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)
lb10.tal8 <- lb10.tal8 %>% arrange(desc(Score)) %>% select(Sequence.ID, Score)

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

safe_mean <- function(x) mean(x[is.finite(x)])

for (temp in 1:8) {
  tal <- get(paste0("lb10.tal", temp))
  
  P3 <- cpm %>%
    filter(Feature %in% tal$Sequence.ID) %>%
    mutate(
      LBR24 = (LBR124 + LBR224 + LBR324) / 3,
      LBR48 = (LBR148 + LBR248 + LBR348) / 3,
      LBR72 = (LBR172 + LBR272 + LBR372) / 3,
      LBR96 = (LBR196 + LBR296 + LBR396) / 3,
      # Cont24 = (Cont_124 + Cont_224 + Cont_324) / 3,
      # Cont48 = (Cont_148 + Cont_248 + Cont_348) / 3,
      # Cont72 = (Cont_172 + Cont_272 + Cont_372) / 3,
      # Cont96 = (Cont_196 + Cont_296 + Cont_396) / 3
      PR24 = (PR124 + PR224 + PR324) / 3,
      PR48 = (PR148 + PR248 + PR348) / 3,
      PR72 = (PR172 + PR272 + PR372) / 3,
      PR96 = (PR196 + PR296 + PR396) / 3
    ) %>%
    mutate(
      H24 = log2((LBR24 + 1) / (PR24 + 1)),
      H48 = log2((LBR48 + 1) / (PR48 + 1)),
      H72 = log2((LBR72 + 1) / (PR72 + 1)),
      H96 = log2((LBR96 + 1) / (PR96 + 1))
    ) %>%
    select(Feature, H24, H48, H72, H96)
  
  H <- list(P3$H24, P3$H48, P3$H72, P3$H96)
  
  logfc_mat[temp, ] <- vapply(H, safe_mean, numeric(1))
  pval_mat[temp, ]  <- vapply(H, safe_t,    numeric(1))
  
  cat("template", temp, "(n =", nrow(P3), "genes)\n")
}

logfc_mat
pval_mat


write.table(pval_mat, 'result/fig4c/LBR/LBR_v_PR_pval.tsv', quote = FALSE, row.names = FALSE, sep = '\t')
