# Load Libraries
# remotes::install_version("ggplot2", version = '3.5.0')
library(ggplot2)
library(tidyverse)
library(ggtern)
library(svglite)

# Read expressed genes list (CPM)
expressgene <- read.table('fig2/CPM_total_kallisto_tximport.txt', sep = ' ', header = TRUE)


# Make a table with the means of expressed genes on different time period
meangeneexp <- expressgene %>%
  mutate(
    # Control Averages
    Cont24 = (Cont_124 + Cont_224 + Cont_324) / 3,
    Cont48 = (Cont_148 + Cont_248 + Cont_348) / 3,
    Cont72 = (Cont_172 + Cont_272 + Cont_372) / 3,
    Cont96 = (Cont_196 + Cont_296 + Cont_396) / 3,
    
    # LBR Averages
    LBR24 = (LBR124 + LBR224 + LBR324) / 3,
    LBR48 = (LBR148 + LBR248 + LBR348) / 3,
    LBR72 = (LBR172 + LBR272 + LBR372) / 3,
    LBR96 = (LBR196 + LBR296 + LBR396) / 3,
    
    # PR Averages
    PR24 = (PR124 + PR224 + PR324) / 3,
    PR48 = (PR148 + PR248 + PR348) / 3,
    PR72 = (PR172 + PR272 + PR372) / 3,
    PR96 = (PR196 + PR296 + PR396) / 3
  ) %>%
  select(Feature, Cont24, Cont48, Cont72, Cont96, LBR24, LBR48, LBR72, LBR96, PR24, PR48, PR72, PR96)

#--------------------------------Control24Hour----------------------------------

dataset <- meangeneexp %>%
  select(Feature, PR24) %>%
  rename(TPM = PR24)

# Load syntenic gene list
syntenic <- read.delim('fig2/Durum_Rye_ABR_unique_triplet.txt', header = FALSE)
colnames(syntenic) <- c('A', 'B', 'R')


# Add TPM values to each of the syntenic gene table
# A gene
syntenic <- syntenic %>%
  left_join(dataset, by = c('A' = 'Feature')) %>%
  rename(TPM_A = TPM) %>%
  relocate(TPM_A, .after = A)

# B Gene
syntenic <- syntenic %>%
  left_join(dataset, by = c('B' = 'Feature')) %>%
  rename(TPM_B = TPM) %>%
  relocate(TPM_B, .after = B)

# R Gene
syntenic <- syntenic %>%
  left_join(dataset, by = c('R' = 'Feature')) %>%
  rename(TPM_R = TPM) %>%
  relocate(TPM_R, .after = R)

# Filter out rows where sum of CPM < 0.5
syntenic <- syntenic %>%
  filter((TPM_A + TPM_B + TPM_R) > 0.5)


# Calculate relative expression level
relative_expression <- syntenic %>%
  mutate(
    sum_tpm = TPM_A + TPM_B + TPM_R,
    A = (TPM_A / sum_tpm) * 100,
    B = TPM_B / sum_tpm * 100,
    R = TPM_R / sum_tpm * 100
  ) %>%
  select(-sum_tpm, -TPM_A, -TPM_B, -TPM_R)

#-----------------------Expression Level----------------------------------------

relative_expression <- relative_expression %>%
  mutate(
    expression = case_when(
      A > B & A > R & B < 20 & R < 20 ~ "A_dom",
      B > A & B > R & A < 20 & R < 20 ~ "B_dom",
      R > A & R > B & A < 20 & B < 20 ~ "R_dom",
      A < 20 ~ "A_sup",
      B < 20 ~ "B_sup",
      R < 20 ~ "R_sup",
      TRUE ~ "balanced"
    )
  )

table(relative_expression$expression)
head(relative_expression)

test <- relative_expression %>% filter(expression == 'R_dom')
test <- test[4, ]
test

#----------------------------------Ternary Plots--------------------------------

# Color
# group_colors <- c(
#   "Balanced" = "gray60",       # Updated to match new labels
#   "A Suppressed" = "#1B9E77",  
#   "B Suppressed" = "#D95F02",  
#   "R Suppressed" = "#7570B3",  
#   "A Dominant" = "#E7298A",    
#   "B Dominant" = "#66A61E",    
#   "R Dominant" = "#E6AB02"     
# )
# 
# group_colors <- c(
#   "Balanced" = "#8DD3C7",
#   "A Suppressed" = "#FB8072",
#   "B Suppressed" = "#BEBADA",
#   "R Suppressed" = "#80B1D3",
#   "A Dominant" = "#FDB462",
#   "B Dominant" = "#B3DE69",
#   "R Dominant" = "#FCCDE5"
# )

group_colors <- c(
  "Balanced" = "#8DD3C7",
  "A Suppressed" = "#FDB462",
  "B Suppressed" = "#B3DE69",
  "R Suppressed" = "#FCCDE5",
  "A Dominant" = "#FB8072",
  "B Dominant" = "#BEBADA",
  "R Dominant" = "#80B1D3"
)

# Make sure your factor levels and labels are correct
relative_expression$expression <- factor(
  relative_expression$expression,
  levels = c("balanced", "A_sup", "B_sup", "R_sup", "A_dom", "B_dom", "R_dom"),
  labels = c("Balanced", "A Suppressed", "B Suppressed", "R Suppressed", "A Dominant", "B Dominant", "R Dominant")
)

# Create the ternary plot
ternaryplot <- ggtern(data = relative_expression, aes(x = A, y = B, z = R, color = expression)) +
  geom_point(size = 1.5, alpha = 0.3) +
  scale_color_manual(values = group_colors) +  # Removed the undefined 'labels' parameter
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_bw() +
  theme_showarrows() +
  theme_clockwise() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    # tern.axis.ticks.length.major = unit(0.8, "cm"),
    axis.ticks = element_line(linewidth = 1, colour = "black"),
    legend.position = 'none',
    # legend.direction = 'horizontal',
    legend.box = 'vertical',
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, 'cm'),
    legend.spacing.x = unit(0.3, 'cm'),
    legend.margin = margin(b = 10),
    tern.axis.arrow.sep = 0.12,
    tern.axis.arrow = element_line(linewidth = 2, colour = "black")
  )

ternaryplot

ggsave(ternaryplot, file = 'fig2/ternary_plot_new.svg',  height = 7, width = 7, dpi = 600)

## Without Legend
ternaryplot1 <- ggtern(data = test, aes(x = A, y = B, z = R, color = expression)) +
  geom_point(size = 1.5, alpha = 0.3, show.legend = FALSE) +
  scale_color_manual(values = group_colors, labels = labels) +
  theme_bw() +
  theme_showarrows() +
  theme_clockwise() +
  theme(
    axis.text = element_text(size = 18, face = 'bold'),
    tern.axis.ticks.length.major = unit(0.8, "cm"),
    axis.ticks = element_line(linewidth = 1, colour = "black"),
    # tern.axis.arrow.sep = 0.12,
    tern.axis.arrow = element_line(size = 2, colour = "black")
  )

ternaryplot1

## Single point
# 3 arrows: start at each vertex, end at the datapoint
pt <- as.numeric(test[1, c("A", "B", "R")])          # the data point
V  <- rbind(A = c(100, 0, 0),                          # the three vertices
            B = c(0, 100, 0),
            R = c(0, 0, 100))

gap <- 4   # distance from arrow tip to the point (composition units, 0–100). Increase to push further.

ends <- t(apply(V, 1, function(v) {
  d <- v - pt                       # direction from point toward this vertex
  pt + gap * d / sqrt(sum(d^2))     # step back from the point by `gap`
}))

arrows_df <- data.frame(
  A = V[, 1], B = V[, 2], R = V[, 3],
  Aend = ends[, 1], Bend = ends[, 2], Rend = ends[, 3]
)

group_colors <- c("R_dom" = "#c1121f")

ternaryplot1 <- ggtern(data = test, aes(x = A, y = B, z = R, color = expression)) +
  geom_segment(
    data = arrows_df,
    aes(x = A, y = B, z = R, xend = Aend, yend = Bend, zend = Rend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    linewidth = 0.3, colour = "black"
  ) +
  geom_point(size = 2, alpha = 1, show.legend = FALSE) +
  scale_color_manual(values = group_colors, labels = labels) +
  theme_bw() +
  theme_showgrid() +          
  theme_hidetitles() +          # <- removes A / B / R at the corners
  theme_clockwise() +
  theme_hidelabels() +          # <- removes 0 / 20 / 40 ... tick labels
  theme_hideticks() +           # <- removes the tick marks
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_blank(),
    # tern.axis.ticks.length.major = unit(0.8, "cm"),
    axis.ticks = element_line(linewidth = 1, colour = "black"),
    legend.position = 'none',
    # legend.direction = 'horizontal',
    legend.box = 'vertical',
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.8, 'cm'),
    legend.spacing.x = unit(0.3, 'cm'),
    legend.margin = margin(b = 10),
    tern.axis.arrow.sep = 0.12,
    tern.axis.arrow = element_line(linewidth = 0.2, colour = "black")
  )
ternaryplot1


ggsave(ternaryplot1, file = 'fig2/ternary_plot_single.svg',  height = 7, width = 7, dpi = 600)
