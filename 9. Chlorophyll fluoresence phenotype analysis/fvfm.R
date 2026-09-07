# Load library
library(tidyverse)
library(ggpubr)
library(rstatix)
library(svglite)

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

View(df.p3.l1.8)

process_df <- function(df) {
  df %>%
    select(-V1, -V2, -V6, -V7) %>%
    setNames(c('Rep', 'Metric', 'Value')) %>%
    filter(Metric %in% c("FvFmMax"))
}

process_df <- function(df) {
  df %>%
    select(-V1, -V2, -V6, -V7) %>%
    setNames(c('Rep', 'Metric', 'Value')) %>%
    filter(Metric %in% c("NPQ_5.0"))
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

combine.fv.fm <- bind_rows(
  p3.l1.8 %>% mutate(Group = 'P3', Time = '8hpi'),
  p3.l2.8 %>% mutate(Group = 'P3', Rep = Rep + 6, Time = '8hpi'),
  p3.l1.16 %>% mutate(Group = 'P3', Time = '16hpi'),
  p3.l2.16 %>% mutate(Group = 'P3', Rep = Rep + 6, Time = '16hpi'),
  p3.l1.24 %>% mutate(Group = 'P3', Time = '24hpi'),
  p3.l2.24 %>% mutate(Group = 'P3', Rep = Rep + 6, Time = '24hpi'),
  p3.l1.48 %>% mutate(Group = 'P3', Time = '48hpi'),
  p3.l2.48 %>% mutate(Group = 'P3', Rep = Rep + 6, Time = '48hpi'),
  p3.l1.72 %>% mutate(Group = 'P3', Time = '72hpi'),
  p3.l2.72 %>% mutate(Group = 'P3', Rep = Rep + 6, Time = '72hpi'),
  ko.l1.8 %>% mutate(Group = 'KO', Time = '8hpi'),
  ko.l2.8 %>% mutate(Group = 'KO', Rep = Rep + 6, Time = '8hpi'),
  ko.l1.16 %>% mutate(Group = 'KO', Time = '16hpi'),
  ko.l2.16 %>% mutate(Group = 'KO', Rep = Rep + 6, Time = '16hpi'),
  ko.l1.24 %>% mutate(Group = 'KO', Time = '24hpi'),
  ko.l2.24 %>% mutate(Group = 'KO', Rep = Rep + 6, Time = '24hpi'),
  ko.l1.48 %>% mutate(Group = 'KO', Time = '48hpi'),
  ko.l2.48 %>% mutate(Group = 'KO', Rep = Rep + 6, Time = '48hpi'),
  ko.l1.72 %>% mutate(Group = 'KO', Time = '72hpi'),
  ko.l2.72 %>% mutate(Group = 'KO', Rep = Rep + 6, Time = '72hpi'),
)

combine.fv.fm

combine.fv.fm <- combine.fv.fm %>%
  mutate(
    Value = as.numeric(Value),
    Time = factor(Time, levels = c('8hpi', '16hpi', '24hpi', '48hpi', '72hpi')),
    Group = factor(Group, levels = c('P3', 'KO'))
  )

## T-test
stat.test <- combine.fv.fm %>%
  group_by(Time) %>%
  t_test(Value ~ Group) %>%
  add_significance("p") %>%
  add_xy_position(x = "Group")

stat.test

p <- ggplot(combine.fv.fm, aes(Group, Value)) +
  geom_boxplot(aes(fill = Group), width = 0.6, outlier.shape = NA) +
  # geom_jitter(width = 0.12, size = 1.2, alpha = 1.0) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, binwidth = 0.001) +
  facet_wrap(~ Time, nrow = 1) +
  stat_pvalue_manual(stat.test, label = "p.signif",
                     tip.length = 0.01, bracket.size = 0.4) +
  scale_fill_manual(values = c(P3 = "#4C72B0", KO = "#e76f51")) +
  labs(x = NULL, y = expression(F[v]/F[m])) +
  # labs(x = NULL, y = 'NPQ_5.0') +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0.2), axis.ticks = element_line(linewidth = 0.2),
        strip.background = element_blank(),
        
        strip.text = element_text(size = 12))
p

ggsave(p, file = 'result/final-figure/mean-fvfmboxplot.png', height = 6, width = 6, dpi = 300)
