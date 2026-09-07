# Load library
library(tidyverse)
library(ggplot2)

# Load dataset
triad_ratio <- data.frame(
  A_B_R_same = 281,
  A_B_R_diff = 169,
  A_B_same_R_diff = 115,
  A_R_same_B_diff = 102,
  B_R_same_A_diff = 106
)


# Convert to long format
triad_long <- triad_ratio %>%
  pivot_longer(cols = everything(), 
               names_to = "category", 
               values_to = "count") %>%
  mutate(category = case_when(
    category == "A_B_R_same" ~ "A, B, R Same",
    category == "A_B_R_diff" ~ "A, B, R Different", 
    category == "A_B_same_R_diff" ~ "A, B Same; R Different",
    category == "A_R_same_B_diff" ~ "A, R Same; B Different",
    category == "B_R_same_A_diff" ~ "B, R Same; A Different"
  ))

triad_long

# Pie chart

color <- c("#3c8f6f", "#a2bc8c", "#f8ebc3", "#edaf85", "#de6c6c")

piechart <- ggplot(triad_long, aes(x = '', y = count, fill = category)) +
  geom_col() +
  coord_polar('y', start = 0) +
  scale_fill_manual(values = color) +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold', size = 8) +
  theme_void() +
  theme(
    legend.position = 'none',
    legend.direction = 'horizontal',
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = 'bold')
  )

piechart

ggsave(piechart, file = 'fig3/piechart_2.svg', height = 5, width = 5, dpi = 600)
