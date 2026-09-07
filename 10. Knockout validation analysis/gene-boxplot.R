# Load library
library(tidyverse)
library(ggplot2)
library(svglite)
library(rstatix)
library(ggpubr)

# Load dataset
cpm.24 <- read.delim('result/CPM-table-24.tsv', header = TRUE, sep = '\t')
cpm.72 <- read.delim('result/CPM-table-72.tsv', header = TRUE, sep = '\t')
exp.24 <- read.delim('result/expressed-genes-24.tsv', header = TRUE, sep = '\t')
exp.72 <- read.delim('result/expressed-genes-72.tsv', header = TRUE, sep = '\t')

gene.24 <- cpm.24 %>%
  filter(Feature %in% c('TRITD3Av1G249280', 'TRITD1Av1G150360', 'TRITD2Bv1G161480'))
gene.24

gene.72 <- cpm.72 %>%
  filter(Feature %in% c('TRITD3Av1G249280', 'TRITD1Av1G150360', 'TRITD2Bv1G161480'))

gene.72

exp.24 <- exp.24 %>%
  filter(Feature %in% c('TRITD3Av1G249280', 'TRITD1Av1G150360', 'TRITD2Bv1G161480')) %>%
  select(Feature, padj)

head(exp.24)

exp.72 <- exp.72 %>%
  filter(Feature %in% c('TRITD3Av1G249280', 'TRITD1Av1G150360', 'TRITD2Bv1G161480')) %>%
  select(Feature, padj) %>%
  mutate(padj = replace_na(padj, 1))

head(exp.72)

colnames(gene.72) <- colnames(gene.24)

gene <- rbind(gene.24 %>% mutate(time = 24),
              gene.72 %>% mutate(time = 72))
gene

gene <- gene %>%
  left_join(exp.24, by = c('Feature' = 'Feature')) %>%
  left_join(exp.72, by = c('Feature' = 'Feature')) %>%
  mutate(padj.24 = padj.x, padj.72 = padj.y) %>%
  select(-padj.x, -padj.y)

# Plot
df <- gene |>
  pivot_longer(cols = matches("^(KO|P3)_"),
               names_to = c("strain", "rep"), names_sep = "_",
               values_to = "tpm") |>
  mutate(
    strain  = factor(strain, levels = c("P3", "KO")),          # P3 = control, reference
    Feature = factor(Feature, levels = c("TRITD3Av1G249280",   # 3Av first
                                         "TRITD2Bv1G161480",
                                         "TRITD1Av1G150360"))   # 1Av third
  )

# significance table: pick padj matching each panel's timepoint, map to stars
sig <- df |>
  distinct(Feature, time, padj.24, padj.72) |>
  mutate(
    padj     = if_else(time == 24, padj.24, padj.72),
    group1   = "P3", group2 = "KO",
    p.signif = cut(padj, breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
                   labels = c("****", "***", "**", "*", "ns"))
  ) |>
  left_join(df |> group_by(Feature, time) |> summarise(ymax = max(tpm), .groups = "drop"),
            by = c("Feature", "time")) |>
  mutate(y.position = ymax * 1.06)

p <- ggplot(df, aes(strain, tpm)) +
  geom_boxplot(aes(fill = strain), width = 0.62, colour = "grey25",
               linewidth = 0.5, fatten = 1.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.09, height = 0, seed = 1),
             size = 1.9, colour = "black") +
  ggpubr::stat_pvalue_manual(sig, label = "p.signif",
                             tip.length = 0.01, bracket.size = 0.4,
                             hide.ns = FALSE) +
  facet_wrap(vars(time, Feature), scales = "free_y", ncol = 3,
             labeller = labeller(time = function(x) paste0(x, " hpi"),
                                 .multi_line = TRUE)) +
  scale_fill_manual(values = c(P3 = "#4472B0", KO = "#E8734A"), guide = "none") +
  scale_y_continuous(n.breaks = 4,
                     expand = expansion(mult = c(0.08, 0.18))) +   # was 0.12; more headroom for brackets
  labs(x = NULL, y = "Expression (meanCPM)") +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 13, colour = "black", margin = margin(b = 8)),
    axis.text        = element_text(colour = "black", size = 13),
    axis.title.y     = element_text(size = 14, margin = margin(r = 8)),
    axis.line        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(4, "pt"),
    panel.spacing    = unit(1.6, "lines"),
    plot.margin      = margin(10, 10, 5, 5)
  )
p

ggsave(p, file = 'result/fig5/boxplot.svg', height = 6, width = 10, dpi = 300)
