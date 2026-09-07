# Load library
library(tidyverse)
library(ggplot2)

# 
p3.ko.deg <- read.delim('result/DEGs.tsv', header = TRUE, sep = '\t')$Feature
p3.ko.cpm <- read.delim('result/CPM-table.tsv', header = TRUE, sep = '\t')

head(p3.ko.cpm)
head(p3.ko.deg)

deg.cpm <- p3.ko.cpm %>%
  filter(Feature %in% p3.ko.deg)

head(deg.cpm)

deg.long <- deg.cpm %>%
  pivot_longer(-Feature, names_to = "Sample", values_to = "CPM") %>%
  mutate(Group   = str_remove(Sample, "_[0-9]+$"),
         Feature = factor(Feature, levels = p3.ko.deg)) %>%
  mutate(Group = factor(Group, levels = c('P3', 'KO')))

View(deg.long)
head(deg.cpm)

p <- ggplot(deg.long, aes(Group, CPM, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.85) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(KO = "#D55E00", P3 = "#0072B2")) +
  labs(x = NULL, y = "Expression (CPM)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.text       = element_text(size = 8, face = "italic"),
        panel.grid.minor = element_blank())

p <- ggplot(deg.long, aes(Group, CPM, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9,
               linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.4, shape = 21,
              fill = "black", alpha = 1.0) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(P3 = "#357266", KO = "#a30015")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "Expression (CPM)") +
  theme_classic() +
  theme(legend.position  = "none",
        strip.text       = element_text(size = 5, hjust = 0.5),
        strip.background  = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_text(size = 6),
        axis.text.y      = element_text(size = 6),
        axis.title.y = element_blank(),
        panel.spacing    = unit(0.6, "lines"))
p

ggsave(p, file = 'result/boxplot.svg', height = 6, width = 4, dpi = 300)
