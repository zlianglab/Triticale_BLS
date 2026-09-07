# Load library
library(tidyverse)
library(purrr)
library(ggpubr)

# Load directories
samples <- list.files(path = 'data/Rerun', full.names = T)
file <- file.path(samples)
names(file) <- str_replace(samples, 'data/Rerun/', '') %>%
  str_replace('.txt', '')


names(file)

# Load dataset
subgenome <- lapply(file, read.delim, sep = ' ', header = FALSE)

combined <- map_dfr(names(subgenome), ~{
  data.frame(
    value = subgenome[[.x]]$V2,
    dataset = .x,
    stringsAsFactors = TRUE
  )
})

combined %>% View()  

combined_long <- combined %>%
  mutate(
    region_type = case_when(
      grepl("1000bp", dataset) ~ "1000bp",
      grepl("CDS", dataset) ~ "CDS",
      TRUE ~ "Other"
    ),
    comparison_type = case_when(
      grepl("AB_same_R_diff", dataset) ~ "AB_same_R_diff",
      grepl("AR_same_B_diff", dataset) ~ "AR_same_B_diff", 
      grepl("BR_same_A_diff", dataset) ~ "BR_same_A_diff",
      grepl("ABR_diff", dataset) ~ "ABR_diff",
      grepl("ABR_same", dataset) ~ "ABR_same",
      TRUE ~ "Other"
    )
  )%>% select(-dataset)

combined_long$region_type <- factor(combined_long$region_type, levels = c('CDS', '1000bp'))
combined_long$comparison_type <- factor(combined_long$comparison_type, levels = c('ABR_same', 'ABR_diff', 'AB_same_R_diff', 'AR_same_B_diff', 'BR_same_A_diff'))

combined_long %>% View()

as_tibble(combined_long)

# Combined boxplot
p <- ggplot(combined_long, aes(x = comparison_type, y = value, fill = region_type)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("1000bp" = "#E69F00", "CDS" = "#56B4E9")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = "Similarity Score",
    fill = "Region Type"
  )

p

p <- ggplot(combined_long, aes(x = comparison_type, y = value, fill = comparison_type)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~ region_type, ncol = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("comparison1", "comparison2")),
                     label = "p.signif") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.position = "none",
    strip.text = element_text(size = 8, face = "bold", color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  labs(
    x = NULL, 
    y = "Similarity Score"
  ) +
  scale_fill_brewer(type = "qual", palette = "Dark2")

## Wilcox test

p <- ggplot(combined_long, aes(x = comparison_type, y = value, fill = comparison_type)) +
  geom_boxplot(alpha = 0.8) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",
                     label.y.npc = 0.9) +  # Position as fraction of panel height
  facet_wrap(~ region_type, ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  labs(x = NULL, y = "Similarity Score") +
  scale_fill_brewer(type = "qual", palette = "Dark2")

p

ggsave(p, file = 'result/mafft_boxplot.png', height = 5, width = 4, dpi = 300)
