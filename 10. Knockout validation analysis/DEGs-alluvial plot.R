# Load library
library(tidyverse)
library(ggplot2)
library(ggalluvial)
library(svglite)

# Load dataset
deg.24 <- read.delim('result/DEGs.tsv', header = TRUE, sep = '\t')
deg.72 <- read.delim('result/DEGs-72.tsv', header = TRUE, sep = '\t')

head(deg.24)
head(deg.72)

deg <- rbind(
  deg.24 %>% mutate(Time = 24),
  deg.72 %>% mutate(Time = 72)
)

lvls <- c("Upregulated", "Downregulated", "Non-DEG")
pal  <- c(Upregulated = "#2a9d8f", Downregulated = "#f4a261", "Non-DEG" = "#264653")

s24 <- deg.24 %>% transmute(Feature, t24 = if_else(log2FoldChange > 1, "Upregulated", "Downregulated"))
s72 <- deg.72 %>% transmute(Feature, t72 = if_else(log2FoldChange > 1, "Upregulated", "Downregulated"))

head(s24)
table(s24$t24)
table(s72$t72)

flow <- full_join(s24, s72, by = "Feature") %>%
  mutate(
    t24 = factor(replace_na(t24, "Non-DEG"), lvls),
    t72 = factor(replace_na(t72, "Non-DEG"), lvls)
  ) %>%
  count(t24, t72)

p <- ggplot(flow, aes(axis1 = t24, axis2 = t72, y = n)) +
  geom_alluvium(aes(fill = t24), width = 0.22, alpha = 0.3,
                color = NA, decreasing = NA, discern = FALSE) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.22,
               color = "#ffffff", linewidth = 0, decreasing = NA, discern = FALSE) +
  # geom_text(stat = "stratum", decreasing = NA,
  #           aes(label = paste0(after_stat(stratum), "\nn = ", after_stat(count))),
  #           color = "white", size = 3.8, lineheight = 0.9) +
  scale_x_discrete(limits = c("24", "72"), expand = c(0.06, 0.06)) +
  scale_fill_manual(values = pal, name = "") +
  labs(y = "Number of genes") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid       = element_blank(),
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    axis.text.y      = element_blank(),
    axis.text.x      = element_text(size = 14, color = "grey10"),
    legend.position  = "top",
    plot.margin      = margin(12, 20, 12, 20)
  )
p

ggsave(p, file = 'result/fig5/alluvial.svg', height = 6, width = 4, dpi = 300)
