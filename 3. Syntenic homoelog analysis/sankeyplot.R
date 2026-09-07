# Load library
library(tidyverse)
library(ggplot2)
library(networkD3)
library(ggalluvial)
library(svglite)

# Load dataset
lbr <- read.delim('fig2/LBR.txt', sep = ' ', header = TRUE)
pr <- read.delim('fig2/PR.txt', sep = ' ', header = TRUE)



head(lbr)

lbr_stable <- lbr %>%
  filter(pattern.24 == pattern.48 & pattern.48 == pattern.72 & pattern.24 == pattern.72)
pr_stable <- pr %>%
  filter(pattern.24 == pattern.48 & pattern.48 == pattern.72 & pattern.24 == pattern.72)

(nrow(lbr_stable) / nrow(lbr)) * 100
(nrow(pr_stable) / nrow(pr)) * 100

lbr <- lbr %>%
  mutate(treatment = 'LBR')
pr <- pr %>%
  mutate(treatment = 'PR')

# Sankey data
lbr_sankey <- lbr %>%
  select(triad_id, treatment, pattern.24, pattern.48, pattern.72, pattern.96) %>%
  pivot_longer(
    cols = starts_with('pattern.'),
    names_to = 'timepoint',
    values_to = 'pattern'
  ) %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c('pattern.24', 'pattern.48', 'pattern.72', 'pattern.96'),
                       labels = c('24h', '48h', '72h', '96h'))
  )


lbr_sankey %>% View()

pr_sankey <- pr %>%
  select(triad_id, treatment, pattern.24, pattern.48, pattern.72, pattern.96) %>%
  pivot_longer(
    cols = starts_with('pattern.'),
    names_to = 'timepoint',
    values_to = 'pattern'
  ) %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c('pattern.24', 'pattern.48', 'pattern.72', 'pattern.96'),
                       labels = c('24h', '48h', '72h', '96h'))
  )

combined <- bind_rows(
  lbr_sankey %>% mutate(triad_id = paste0(triad_id, "_LBR")),
  pr_sankey %>% mutate(triad_id = paste0(triad_id, "_PR"))
) %>%
  mutate(treatment = factor(treatment, levels = c("LBR", "PR"))) %>%
  mutate(pattern = factor(pattern, levels = c(
    'Balanced', 'A-dominant', 'B-dominant', 'R-dominant', 'A-suppressed', 'B-suppressed', 'R-suppressed'
  )))

combined %>% View()

# Comparison between fixed timepoint in LBR and PR
sankey_dataset <- function(timepoint) {
  lbr_timepoint <- paste0('LBR_', timepoint)
  pr_timepoint <- paste0('PR_', timepoint)
  
  df.lbr <- lbr %>%
    select(triad_id, all_of(paste0('pattern.', timepoint))) %>%
    rename(pattern = paste0('pattern.', timepoint)) %>%
    mutate(
      timepoint = lbr_timepoint,
    ) %>%
    select(triad_id, timepoint, pattern)
  
  df.pr <- pr %>%
    select(triad_id, all_of(paste0('pattern.', timepoint))) %>%
    rename(pattern = paste0('pattern.', timepoint)) %>%
    mutate(
      timepoint = pr_timepoint,
    ) %>%
    select(triad_id, timepoint, pattern)
  
  # Combine
  combine <- bind_rows(df.lbr, df.pr) %>%
    mutate(
      timepoint = factor(timepoint, levels = c(lbr_timepoint, pr_timepoint)),
      pattern = factor(pattern, levels = c(
        'Balanced', 'A-dominant', 'B-dominant', 'R-dominant', 'A-suppressed', 'B-suppressed', 'R-suppressed'
      ))
    )
  
  return(combine)
}

combine24 <- sankey_dataset(24)
combine48 <- sankey_dataset(48)
combine72 <- sankey_dataset(72)
combine96 <- sankey_dataset(96)


# Define color scheme
pattern_colors <- c(
  "Balanced" = "#8DD3C7",
  "A-dominant" = "#FB8072",
  "B-dominant" = "#BEBADA",
  "R-dominant" = "#80B1D3",
  "A-suppressed" = "#FDB462",
  "B-suppressed" = "#B3DE69",
  "R-suppressed" = "#FCCDE5",
  "None" = "#E0E0E0"
)

# Sankey
sankey <- ggplot(combine96, aes(x = timepoint, stratum = pattern, alluvium = triad_id, fill = pattern)) +
  geom_flow(stat = 'alluvium',
            aes(fill = pattern),
            lode.guidance = 'frontback',
            width = 0.3,
            alpha = 0.4,
            curve_type = 'sigmoid' 
  ) +
  geom_stratum(
    aes(fill = pattern),
    width = 0.3,
    alpha = 0.9,
    color = 'black',
    size = 0.5
  ) +
  scale_fill_manual(values = pattern_colors,
                    breaks = c("Balanced", "A-dominant", "B-dominant", "R-dominant",
                               "A-suppressed", "B-suppressed", "R-suppressed"))+
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank()
  )

sankey

ggsave(sankey, file = 'fig2/supplemental/sankey_LBR_PR_96.png', bg = 'transparent', height = 5, width = 3, dpi = 300)
