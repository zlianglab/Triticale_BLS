# Load library
library(tidyverse)
library(broom)
library(emmeans)

# Load dataset
df.p3.l1.8 <- read.delim('data/round3/8hpi/p3-l1.txt', header = FALSE, sep = '\t')
df.p3.l2.8 <- read.delim('data/round3/8hpi/p3-l2.txt', header = FALSE, sep = '\t')
df.ko.l1.8 <- read.delim('data/round3/8hpi/ko-l1.txt', header = FALSE, sep = '\t')
df.ko.l2.8 <- read.delim('data/round3/8hpi/ko-l2.txt', header = FALSE, sep = '\t')
df.p3.l1.16 <- read.delim('data/round3/16hpi/p3-l1.txt', header = FALSE, sep = '\t')
df.p3.l2.16 <- read.delim('data/round3/16hpi/p3-l2.txt', header = FALSE, sep = '\t')
df.ko.l1.16 <- read.delim('data/round3/16hpi/ko-l1.txt', header = FALSE, sep = '\t')
df.ko.l2.16 <- read.delim('data/round3/16hpi/ko-l2.txt', header = FALSE, sep = '\t')
df.p3.l1.24 <- read.delim('data/round3/24hpi/p3-l1.txt', header = FALSE, sep = '\t')
df.p3.l2.24 <- read.delim('data/round3/24hpi/p3-l2.txt', header = FALSE, sep = '\t')
df.ko.l1.24 <- read.delim('data/round3/24hpi/ko-l1.txt', header = FALSE, sep = '\t')
df.ko.l2.24 <- read.delim('data/round3/24hpi/ko-l2.txt', header = FALSE, sep = '\t')
df.p3.l1.48 <- read.delim('data/round3/48hpi/p3-l1.txt', header = FALSE, sep = '\t')
df.p3.l2.48 <- read.delim('data/round3/48hpi/p3-l2.txt', header = FALSE, sep = '\t')
df.ko.l1.48 <- read.delim('data/round3/48hpi/ko-l1.txt', header = FALSE, sep = '\t')
df.ko.l2.48 <- read.delim('data/round3/48hpi/ko-l2.txt', header = FALSE, sep = '\t')
df.p3.l1.72 <- read.delim('data/round3/72hpi/p3-l1.txt', header = FALSE, sep = '\t')
df.p3.l2.72 <- read.delim('data/round3/72hpi/p3-l2.txt', header = FALSE, sep = '\t')
df.ko.l1.72 <- read.delim('data/round3/72hpi/ko-l1.txt', header = FALSE, sep = '\t')
df.ko.l2.72 <- read.delim('data/round3/72hpi/ko-l2.txt', header = FALSE, sep = '\t')

process_df <- function(df) {
  df %>%
    select(-V1, -V2, -V6, -V7) %>%
    setNames(c('Rep', 'Metric', 'Value')) %>%
    filter(Metric %in% c("AriIdx", "ChlIdx", "FvFmMax", "NDVI", "NPQ_5.0", "NPQ_10.0")) %>%
    mutate(Value = as.numeric(Value))
}

p3.l1.8 <- process_df(df.p3.l1.8)
p3.l2.8 <- process_df(df.p3.l2.8)
ko.l1.8 <- process_df(df.ko.l1.8)
ko.l2.8 <- process_df(df.ko.l2.8)
p3.l1.16 <- process_df(df.p3.l1.16)
p3.l2.16 <- process_df(df.p3.l2.16)
ko.l1.16 <- process_df(df.ko.l1.16)
ko.l2.16 <- process_df(df.ko.l2.16)
p3.l1.24 <- process_df(df.p3.l1.24)
p3.l2.24 <- process_df(df.p3.l2.24)
ko.l1.24 <- process_df(df.ko.l1.24)
ko.l2.24 <- process_df(df.ko.l2.24)
p3.l1.48 <- process_df(df.p3.l1.48)
p3.l2.48 <- process_df(df.p3.l2.48)
ko.l1.48 <- process_df(df.ko.l1.48)
ko.l2.48 <- process_df(df.ko.l2.48)
p3.l1.72 <- process_df(df.p3.l1.72)
p3.l2.72 <- process_df(df.p3.l2.72)
ko.l1.72 <- process_df(df.ko.l1.72)
ko.l2.72 <- process_df(df.ko.l2.72)

# Combine datapoints into dataframe
grid <- expand.grid(
  treatment = c('p3', 'ko'),
  leaf = c('l1', 'l2'),
  time = c(8, 16, 24, 48, 72),
  stringsAsFactors = FALSE
)

grid

grid$obj <- with(grid, paste(treatment, leaf, time, sep = '.'))

# dat <- pmap_dfr(grid, function(treatment, leaf, time, obj) {
#   get(obj) %>%
#     mutate(treatment = treatment, leaf = leaf, time = time)
# }
#   ) %>%
#   mutate(
#     position = case_when(
#       leaf == "l1" & treatment == "p3" ~ "top",
#       leaf == "l1" & treatment == "ko" ~ "bottom",
#       leaf == "l2" & treatment == "ko" ~ "top",
#       leaf == "l2" & treatment == "p3" ~ "bottom"
#     ),
#     treatment = factor(treatment, levels = c('p3', 'ko')),
#     leaf = factor(leaf),
#     time = factor(time, levels = c(8, 16, 24, 48, 72))
#   ) %>%
#   select(-leaf)

dat <- pmap_dfr(grid, function(treatment, leaf, time, obj) {
  get(obj) %>%
    mutate(treatment = treatment, leaf = leaf, time = time)
}) %>%
  mutate(
    position = case_when(
      leaf == "l1" & treatment == "p3" ~ "top",
      leaf == "l1" & treatment == "ko" ~ "bottom",
      leaf == "l2" & treatment == "ko" ~ "top",
      leaf == "l2" & treatment == "p3" ~ "bottom"
    ),
    Rep = if_else(position == "bottom", Rep + 6, Rep),
    treatment = factor(treatment, levels = c('p3', 'ko')),
    time = factor(time, levels = c(8, 16, 24, 48, 72))
  ) %>%
  select(-leaf, -position)

nrow(dat)
View(dat)

dat <- dat %>% arrange(treatment, time, Metric, Rep)

## Fit one linear model per metric
## metric × timepoint: y ~ position + treatment + position:treatment
# fits <- dat %>%
#   group_by(Metric, time) %>%
#   group_split() %>%
#   set_names(map_chr(., ~ paste(unique(.x$Metric), unique(.x$time), sep = "_t"))) %>%
#   map(~ lm(Value ~ position + treatment + position:treatment, data = .x))

fits <- dat %>%
  group_by(Metric, time) %>%
  group_split() %>%
  set_names(map_chr(., ~ paste(unique(.x$Metric), unique(.x$time), sep = "_t"))) %>%
  map(~ lm(Value ~ treatment, data = .x))

View(fits)

## F-test for each terms
anova.summary <- map_dfr(fits, ~ tidy(anova(.x)), .id = "Metric_time") %>%
  separate(Metric_time, into = c("Metric", "time"), sep = "_t") %>%
  mutate(time = as.integer(time))

print(anova.summary, n = Inf)

## P3 vs KO at each position per timepoint
# contrast.position <- map_dfr(
#   fits,
#   ~ as_tibble(pairs(emmeans(.x, ~ treatment | position))),
#   .id = "Metric_time"
# ) %>%
#   separate(Metric_time, into = c("Metric", "time"), sep = "_t") %>%
#   mutate(time = as.integer(time),
#          p.adj_BH = p.adjust(p.value, method = "BH"),
#          sig = case_when(
#            p.adj_BH < 0.001 ~ "***",
#            p.adj_BH < 0.01  ~ "**",
#            p.adj_BH < 0.05  ~ "*",
#            TRUE             ~ "ns"
#          )) %>%
#   arrange(Metric, time, position)

contrast <- map_dfr(
  fits,
  ~ as_tibble(pairs(emmeans(.x, ~ treatment))),
  .id = "Metric_time"
) %>%
  separate(Metric_time, into = c("Metric", "time"), sep = "_t") %>%
  mutate(time = as.integer(time),
         sig = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE             ~ "ns"
         )) %>%
  arrange(Metric, time)

print(contrast, n = Inf)

write.table(anova.summary, file = 'result/mean/3-point/anova-summary.csv', sep = ',', row.names = FALSE)
write.table(contrast, file = 'result/mean/3-point/p3-ko-cotrast.csv', sep = ',', row.names = FALSE)




















