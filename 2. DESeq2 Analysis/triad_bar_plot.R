# Load library
library(tidyverse)
library(ggplot2)

# Load data
triad_data <- data.frame(
  treatment = c(rep('LB10', 4), rep('P3', 4)),
  timepoint = c(24, 48, 72, 96, 24, 48, 72, 96),
  up = c(4119, 3968, 3129, 3223, 1843, 2562, 3334, 3640),
  up_triad = c(2026, 1930, 1422, 1574, 903, 1188, 1567, 1757),
  down = c(1299, 1919, 1609, 1104, 396, 1633, 2128, 2222),
  down_triad = c(633, 1064, 818, 548, 166, 916, 1135, 1222)
)

# Reshape for figuring plot
datapoint <- triad_data %>%
  transmute(
    treatment, timepoint, up_triad, down_triad,
    up_other = up - up_triad,
    down_other = down - down_triad
  ) %>%
  pivot_longer(
    cols = c(up_other, up_triad, down_other, down_triad),
    names_to = 'category',
    values_to = 'count'
  ) %>%
  mutate(
    count = ifelse(str_detect(category, 'down'), -count, count),
    timepoint = as.factor(timepoint),
    category = factor(
      category,
      levels = c('down_other', 'down_triad', 'up_other', 'up_triad')
    )
  )

head(datapoint)

# Plotting
category_color <- c(
  "up_other" = "#f5edbe",
  "up_triad" = "#a0c29c",  
  "down_other" = "#fdbd8c",
  "down_triad" = "#e9834a"
)

triad_plot <- ggplot(datapoint, aes(x = timepoint, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.7, color = 'black') +
  geom_text(
    aes(label = abs(count)),
    position = position_stack(vjust = 0.5),
    size = 2.5,
    color = 'black',
    fontface = 'bold'
  ) +
  facet_wrap(~treatment, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    values = category_color,
    labels = c(
      "up_other" = "Upregulated",
      "up_triad" = "Upregulated in Triads",
      "down_other" = "Downregulated",
      "down_triad" = "Downregulated in Triads"
    ),
    name = NULL
  ) +
  labs(
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(size = 14),
    strip.background.x = element_blank(),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(hjust = 1),
    strip.background = element_rect(fill = "#ffffff"),
    legend.position = "top"
  ) +
  geom_hline(yintercept = 0, color = "black")

triad_plot

ggsave(triad_plot, file = 'fig1/figs/triad_plot.svg', height = 5, width = 12, dpi = 600)
