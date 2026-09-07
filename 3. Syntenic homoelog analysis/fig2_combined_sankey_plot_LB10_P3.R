# Load library
library(tidyverse)
library(ggplot2)
library(networkD3)
library(ggalluvial)

# Load dataset
lbr <- read.delim('./result/fig2/triad_pattern/LBR.txt', sep = ' ', header = TRUE)
pr <- read.delim('./result/fig2/triad_pattern/PR.txt', sep = ' ', header = TRUE)

lbr %>% View()

as_tibble(lbr)
as_tibble(pr)

table(lbr$pattern.24, lbr$pattern.48)
table(lbr$pattern.48, lbr$pattern.72)
table(lbr$pattern.72, lbr$pattern.96)

table(pr$pattern.24, pr$pattern.48)
table(pr$pattern.48, pr$pattern.72)
table(pr$pattern.72, pr$pattern.96)

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
                       labels = c('24h', '48h'))
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
                       labels = c('24h', '48h'))
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

# Define treatment colors
treatment_colors <- c("LBR" = "#4A90E2", "PR" = "#E24A4A")

# Sankey plot
sankey <- ggplot(combined, aes(x = timepoint, stratum = pattern, alluvium = triad_id, fill = pattern)) +
  geom_flow(stat = 'alluvium',
    aes(fill = pattern),
    lode.guidance = 'frontback',
    width = 0.3,
    alpha = 0.6,
    curve_type = 'sigmoid'
  ) +
  geom_stratum(
    aes(fill = pattern),
    width = 0.3,
    alpha = 0.9,
    color = 'white',
    size = 0.5
  ) +
  scale_fill_manual(values = pattern_colors,
                    breaks = c("Balanced", "A-dominant", "B-dominant", "R-dominant",
                               "A-suppressed", "B-suppressed", "R-suppressed"))+
  facet_wrap(~treatment, ncol = 2) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.text = element_text(size = 12, face = 'bold'),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

  
sankey


# Comparison between fixed timepoint in LBR and PR
# 24 hours
lbr.24 <- lbr %>%
  select(triad_id, pattern.48) %>%
  mutate(
    timepoint = 'LBR_48',
    pattern = pattern.48,
  ) %>%
  select(triad_id, timepoint, pattern)

pr.24 <- pr %>%
  select(triad_id, pattern.48) %>%
  mutate(
    timepoint = 'PR_48',
    pattern = pattern.48,
  ) %>%
  select(triad_id, timepoint, pattern)

# Combine
combine.24 <- bind_rows(lbr.24, pr.24) %>%
  mutate(timepoint = factor(timepoint, levels = c('LBR_48', 'PR_48'))) %>%
  
  mutate(pattern = factor(pattern, levels = c(
    'Balanced', 'A-dominant', 'B-dominant', 'R-dominant', 'A-suppressed', 'B-suppressed', 'R-suppressed'
  )))

combine.24 %>% View()


# # Function
# sankey_dataset <- function(timepoint) {
#   df.lbr <- lbr %>%
#     select(triad_id, paste0('pattern.', timepoint)) %>%
#     mutate(
#       timepoint = paste0('LBR_', timepoint),
#       pattern = paste0('pattern.', timepoint),
#     ) %>%
#     select(triad_id, timepoint, pattern)
#   
#   df.pr <- pr %>%
#     select(triad_id, paste0('pattern.', timepoint)) %>%
#     mutate(
#       timepoint = paste0('LBR_', timepoint),
#       pattern = paste0('pattern.', timepoint),
#     ) %>%
#     select(triad_id, timepoint, pattern)
#   
#   # Combine
#   combine <- bind_rows(df.lbr, df.pr) %>%
#     mutate(timepoint = factor(timepoint, levels = c(paste0('LBR_', timepoint),  paste0('PR_', timepoint)))) %>%
#     
#     mutate(pattern = factor(pattern, levels = c(
#       'Balanced', 'A-dominant', 'B-dominant', 'R-dominant', 'A-suppressed', 'B-suppressed', 'R-suppressed'
#     )))
#   
#   return(combine)
# }

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


# Sankey
sankey <- ggplot(combine24, aes(x = timepoint, stratum = pattern, alluvium = triad_id, fill = pattern)) +
  geom_flow(stat = 'alluvium',
            aes(fill = pattern),
            lode.guidance = 'frontback',
            width = 0.3,
            alpha = 0.6,
            curve_type = 'sigmoid'
  ) +
  geom_stratum(
    aes(fill = pattern),
    width = 0.3,
    alpha = 0.9,
    color = 'white',
    size = 0.5
  ) +
  scale_fill_manual(values = pattern_colors,
                    breaks = c("Balanced", "A-dominant", "B-dominant", "R-dominant",
                               "A-suppressed", "B-suppressed", "R-suppressed"))+
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.text = element_text(size = 12, face = 'bold'),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank()
  )

sankey

ggsave(sankey, file = 'result/fig2/sankey_LBR_PR_24.png', bg = 'transparent', height = 6, width = 4, dpi = 300)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  