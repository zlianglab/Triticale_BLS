# Load library
library(tidyverse)
library(purrr)
library(ggpubr)

# Load directories
samples <- list.files(path = 'fig3/mafft_result', full.names = T)
file <- file.path(samples)
names(file) <- str_replace(samples, 'fig3/mafft_result/', '') %>%
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

combined_long <- combined %>%
  mutate(
    region_type = case_when(
      grepl("1000bp", dataset) ~ "1000bp",
      grepl("CDS", dataset) ~ "CDS",
      TRUE ~ "Other"
    ),
    comparison_type = case_when(
      grepl("AB-same-R-diff", dataset) ~ "AB-same-R-diff",
      grepl("AR-same-B-diff", dataset) ~ "AR-same-B-diff", 
      grepl("BR-same-A-diff", dataset) ~ "BR-same-A-diff",
      grepl("ABR-diff", dataset) ~ "ABR-diff",
      grepl("ABR-same", dataset) ~ "ABR-same",
      TRUE ~ "Other"
    )
  )%>% select(-dataset)

combined_long$region_type <- factor(combined_long$region_type, levels = c('CDS', '1000bp'))
combined_long$comparison_type <- factor(combined_long$comparison_type, levels = c('ABR-same', 'ABR-diff', 'AB-same-R-diff', 'AR-same-B-diff', 'BR-same-A-diff'))

# Combined boxplot
p <- ggplot(combined_long, aes(x = comparison_type, y = value, fill = comparison_type)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~ region_type, ncol = 2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.position = "none",
    strip.text = element_text(size = 8, color = "black"),
    strip.background = element_rect(fill = "white", color = "white"),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  labs(
    x = NULL, 
    y = "Similarity Score"
  ) +
  scale_fill_brewer(type = "qual", palette = "Dark2")

p

ggsave(p, file = 'fig3/mafft_boxplot.svg', height = 5, width = 6, dpi = 600)

## T-Test

# Create group
combined_long$group <- case_when(
  combined_long$comparison_type == "ABR-same" ~ "ABR_same",
  combined_long$comparison_type == "ABR-diff" ~ "ABR_diff", 
  combined_long$comparison_type == "AB-same-R-diff" ~ "AB_same_R_diff",
  combined_long$comparison_type == "AR-same-B-diff" ~ "AR_same_B_diff",
  combined_long$comparison_type == "BR-same-A-diff" ~ "BR_same_A_diff"
)

groups <- c("ABR_same", "ABR_diff", "AB_same_R_diff", "AR_same_B_diff", "BR_same_A_diff")
regions <- c("CDS", "1000bp")

# Initialize empty results dataframe with proper column types
results <- data.frame(
  region = character(),
  g1 = character(),
  g2 = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Run t-tests
for(region in regions) {
  for(i in 1:(length(groups)-1)) {
    for(j in (i+1):length(groups)) {
      
      # Get data for each group
      data1 <- combined_long$value[combined_long$region_type == region & combined_long$group == groups[i]]
      data2 <- combined_long$value[combined_long$region_type == region & combined_long$group == groups[j]]
      
      # Perform t-test if both groups have data
      if(length(data1) > 1 & length(data2) > 1) {
        test_result <- t.test(data1, data2)
        
        # Add to results
        results <- rbind(results, data.frame(
          region = region,
          g1 = groups[i],
          g2 = groups[j],
          p_value = test_result$p.value,
          stringsAsFactors = FALSE
        ))
        
        # Print progress
        cat(sprintf("%s: %s vs %s, p = %.6f\n", region, groups[i], groups[j], test_result$p.value))
      }
    }
  }
}

# Add multiple testing correction
results$p_adj <- p.adjust(results$p_value, method = "BH")
results$significance <- case_when(
  results$p_adj < 0.001 ~ "***",
  results$p_adj < 0.01 ~ "**", 
  results$p_adj < 0.05 ~ "*",
  results$p_adj < 0.1 ~ ".",
  TRUE ~ "ns"
)

results <- results %>%
  group_by(region) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

results$significance <- case_when(
  results$p_adj < 0.001 ~ "***",
  results$p_adj < 0.01  ~ "**",
  results$p_adj < 0.05  ~ "*",
  results$p_adj < 0.1   ~ ".",
  TRUE ~ "ns"
)


View(results)













