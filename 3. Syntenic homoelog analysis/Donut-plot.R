# Load library
library(tidyverse)
library(ggplot2)
library(patchwork)

# Load dataset
cont <- read.delim('fig2/Cont.txt', sep = ' ', header = TRUE)
lbr <- read.delim('fig2/LBR.txt', sep = ' ', header = TRUE)
pr <- read.delim('fig2/PR.txt', sep = ' ', header = TRUE)

# QC
cont <- cont %>% select(A, B, R, pattern.24, pattern.48, pattern.72, pattern.96)
lbr <- lbr %>% select(A, B, R, pattern.24, pattern.48, pattern.72, pattern.96)
pr <- pr %>% select(A, B, R, pattern.24, pattern.48, pattern.72, pattern.96)

head(cont)

# count 
# Stable triads
cont_stable <- cont %>%
  filter(pattern.24 == pattern.48 & pattern.48 == pattern.72 & pattern.72 == pattern.96) %>% nrow()
lbr_stable <- lbr %>%
  filter(pattern.24 == pattern.48 & pattern.48 == pattern.72 & pattern.72 == pattern.96) %>% nrow()
pr_stable <- pr %>%
  filter(pattern.24 == pattern.48 & pattern.48 == pattern.72 & pattern.72 == pattern.96) %>% nrow()

# Dynamic triads
cont_dyn <- nrow(cont) - cont_stable
lbr_dyn <- nrow(lbr) - lbr_stable
pr_dyn <- nrow(pr) - pr_stable

### Control--------------------
# R suppressed
cont_R_sup <- cont %>%
  filter(pattern.24 == 'R-suppressed' | pattern.48 == 'R-suppressed' | pattern.72 == 'R-suppressed' | pattern.96 == 'R-suppressed') %>% pull(R) %>% n_distinct()

# Stable R suppressed
cont_R_sup_stable <- cont %>%
  filter(pattern.24 == 'R-suppressed' & pattern.48 == 'R-suppressed' & pattern.72 == 'R-suppressed' & pattern.96 == 'R-suppressed') %>% pull(R) %>% n_distinct()

# Dynamic R suppressed
cont_R_sup_dyn <- cont_R_sup - cont_R_sup_stable


# R dominant
cont_R_dom <- cont %>%
  filter(pattern.24 == 'R-dominant' | pattern.48 == 'R-dominant' | pattern.72 == 'R-dominant' | pattern.96 == 'R-dominant') %>% pull(R) %>% n_distinct()

# Stable R Dominant
cont_R_dom_stable <- cont %>%
  filter(pattern.24 == 'R-dominant' & pattern.48 == 'R-dominant' & pattern.72 == 'R-dominant' & pattern.96 == 'R-dominant') %>% pull(R) %>% n_distinct()

# Dynamic R Dominant
cont_R_dom_dyn <- cont_R_dom - cont_R_dom_stable

### LBR-----------------
lbr_R_sup <- lbr %>%
  filter(pattern.24 == 'R-suppressed' | pattern.48 == 'R-suppressed' | pattern.72 == 'R-suppressed' | pattern.96 == 'R-suppressed') %>% pull(R) %>% n_distinct()

# Stable R suppressed
lbr_R_sup_stable <- lbr %>%
  filter(pattern.24 == 'R-suppressed' & pattern.48 == 'R-suppressed' & pattern.72 == 'R-suppressed' & pattern.96 == 'R-suppressed') %>% pull(R) %>% n_distinct()

# Dynamic R suppressed
lbr_R_sup_dyn <- lbr_R_sup - lbr_R_sup_stable


# R dominant
lbr_R_dom <- lbr %>%
  filter(pattern.24 == 'R-dominant' | pattern.48 == 'R-dominant' | pattern.72 == 'R-dominant' | pattern.96 == 'R-dominant') %>% pull(R) %>% n_distinct()

# Stable R Dominant
lbr_R_dom_stable <- lbr %>%
  filter(pattern.24 == 'R-dominant' & pattern.48 == 'R-dominant' & pattern.72 == 'R-dominant' & pattern.96 == 'R-dominant') %>% pull(R) %>% n_distinct()

# Dynamic R Dominant
lbr_R_dom_dyn <- lbr_R_dom - lbr_R_dom_stable

### PR-----------------
pr_R_sup <- pr %>%
  filter(pattern.24 == 'R-suppressed' | pattern.48 == 'R-suppressed' | pattern.72 == 'R-suppressed' | pattern.96 == 'R-suppressed') %>% pull(R) %>% n_distinct()

# Stable R suppressed
pr_R_sup_stable <- pr %>%
  filter(pattern.24 == 'R-suppressed' & pattern.48 == 'R-suppressed' & pattern.72 == 'R-suppressed' & pattern.96 == 'R-suppressed') %>% pull(R) %>% n_distinct()

# Dynamic R suppressed
pr_R_sup_dyn <- pr_R_sup - pr_R_sup_stable


# R dominant
pr_R_dom <- pr %>%
  filter(pattern.24 == 'R-dominant' | pattern.48 == 'R-dominant' | pattern.72 == 'R-dominant' | pattern.96 == 'R-dominant') %>% pull(R) %>% n_distinct()

# Stable R Dominant
pr_R_dom_stable <- pr %>%
  filter(pattern.24 == 'R-dominant' & pattern.48 == 'R-dominant' & pattern.72 == 'R-dominant' & pattern.96 == 'R-dominant') %>% pull(R) %>% n_distinct()

# Dynamic R Dominant
pr_R_dom_dyn <- pr_R_dom - pr_R_dom_stable


df <- data.frame(
  'Treatment' = c('Control', 'LBR', 'PR', 'Control-R-sup', 'Control-R-dom', 'LBR-R-sup', 'LBR-R-dom', 'PR-R-sup', 'PR-R-dom'),
  'Total' = c(nrow(cont), nrow(lbr), nrow(pr), cont_R_sup, cont_R_dom, lbr_R_sup, lbr_R_dom, pr_R_sup, pr_R_dom),
  'Stable' = c(cont_stable, lbr_stable, pr_stable, cont_R_sup_stable, cont_R_dom_stable, lbr_R_sup_stable, lbr_R_dom_stable, pr_R_sup_stable, pr_R_dom_stable),
  'Dynamic' = c(cont_dyn, lbr_dyn, pr_dyn, cont_R_sup_dyn, cont_R_dom_dyn, lbr_R_sup_dyn, lbr_R_dom_dyn, pr_R_sup_dyn, pr_R_dom_dyn)
)

# df <- df %>% 
#   mutate(
#     Total = ifelse(Treatment %in% c("Control", "LBR", "PR"), Total * 3, Total),
#     Stable = ifelse(Treatment %in% c("Control", "LBR", "PR"), Stable * 3, Stable),
#     Dynamic = ifelse(Treatment %in% c("Control", "LBR", "PR"), Dynamic * 3, Dynamic)
#   )

print(as_tibble(df), n = 999)

               a### Donut plots
# Colors per ring: outer, middle, inner — stable/dynamic shades per ring
ring_colors <- list(
  "3" = c("Stable" = "#5a189a", "Dynamic" = "#9d4edd"),
  "2" = c("Stable" = "#006d77", "Dynamic" = "#84a98c"),
  "1" = c("Stable" = "#ff9f1c", "Dynamic" = "#ffbf69")
)

groups <- list(
  Control = c("Control-R-dom", "Control-R-sup", "Control"),
  LBR     = c("LBR-R-dom",     "LBR-R-sup",     "LBR"),
  PR      = c("PR-R-dom",      "PR-R-sup",       "PR")
)

make_plot <- function(treatments, title) {
  sub <- df %>%
    filter(Treatment %in% treatments) %>%
    mutate(Treatment = factor(Treatment, levels = treatments),
           ring_id = as.integer(Treatment))
  
  long <- bind_rows(
    sub %>% mutate(category = "Stable",  xmin = 0,      xmax = Stable),
    sub %>% mutate(category = "Dynamic", xmin = Stable, xmax = Total)
  ) %>%
    mutate(
      ymin = ring_id - 0.4,
      ymax = ring_id + 0.4,
      fill_key = paste0(ring_id, "_", category)
    )
  
  fill_vals <- unlist(lapply(names(ring_colors), function(r) {
    cols <- ring_colors[[r]]
    setNames(cols, paste0(r, "_", names(cols)))
  }))
  
  # each ring has its own scale: set xmax per ring to its Total
  ring_maxes <- sub %>% select(ring_id, Total)
  long <- long %>% left_join(ring_maxes, by = "ring_id", suffix = c("", "_max"))
  
  # normalize xmin/xmax per ring to a common 0-1 scale
  long <- long %>%
    mutate(
      xmin_n = xmin / Total_max,
      xmax_n = xmax / Total_max
    )
  
  ggplot(long) +
    geom_rect(aes(xmin = xmin_n, xmax = xmax_n, ymin = ymin, ymax = ymax, fill = fill_key),
              color = "black", linewidth = 0) +
    coord_polar(theta = "x") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(breaks = 1:3, labels = treatments) +
    scale_fill_manual(
      values = fill_vals,
      labels = c("1_Stable" = "R-dom Stable", "1_Dynamic" = "R-dom Dynamic",
                 "2_Stable" = "R-sup Stable", "2_Dynamic" = "R-sup Dynamic",
                 "3_Stable" = "Total Stable",  "3_Dynamic" = "Total Dynamic"),
      name = NULL
    ) +
    labs(title = title) +
    theme_void() +
    theme(
      axis.text.y     = element_text(size = 9, face = "bold"),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "right"
    )
}

p1 <- make_plot(groups$Control, "Control")
p2 <- make_plot(groups$LBR, "LBR")
p3 <- make_plot(groups$PR, "PR")

p <- wrap_plots(p1, p2, p3, nrow = 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.direction = 'horizontal')


ggsave(p, file = 'fig2/donut-plot.svg', height = 8, width = 8, dpi = 300)

# # Donut plots
# max_val <- max(df$Total)
# 
# # Each treatment = 1 ring, split into Stable and Dynamic segments
# rings <- bind_rows(lapply(seq_len(nrow(df)), function(i) {
#   data.frame(
#     treatment = df$Treatment[i],
#     category  = c("Stable", "Dynamic"),
#     xmin      = c(0, df$Stable[i]),
#     xmax      = c(df$Stable[i], df$Total[i]),
#     ring_id   = i
#   )
# })) %>%
#   mutate(
#     ymin = ring_id - 0.4,
#     ymax = ring_id + 0.4,
#     treatment = factor(treatment, levels = c("Control", "LBR", "PR"))
#   )
# 
# p <- ggplot(rings) +
#   geom_rect(aes(xmin = 0, xmax = max_val, ymin = ymin, ymax = ymax),
#             fill = "grey90", color = NA) +
#   geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
#             color = "white", linewidth = 0.3) +
#   # Label: Stable
#   geom_text(data = distinct(rings, treatment, ring_id) %>%
#               left_join(df, by = c("treatment" = "Treatment")),
#             aes(x = Stable / 2, y = ring_id, label = Stable),
#             size = 3.2, fontface = "bold", color = "white") +
#   # Label: Dynamic
#   geom_text(data = distinct(rings, treatment, ring_id) %>%
#               left_join(df, by = c("treatment" = "Treatment")),
#             aes(x = Stable + Dynamic / 2, y = ring_id, label = Dynamic),
#             size = 3.2, fontface = "bold", color = "white") +
#   scale_y_continuous(
#     breaks = 1:3,
#     labels = c("Control", "LBR", "PR")
#   ) +
#   coord_polar(theta = "x", start = 0, clip = "off") +
#   xlim(0, max_val) +
#   ylim(0, 3.6) +
#   scale_fill_manual(values = c("Stable" = "#98c1d9", "Dynamic" = "#ee6c4d"), name = NULL) +
#   theme_void() +
#   theme(
#     legend.position = "top",
#     axis.text.y     = element_text(size = 10),
#     plot.margin     = margin(10, 20, 10, 10),
#     plot.title      = element_text(hjust = 0.5,, size = 13),
#     legend.direction = 'horizontal'
#   )
# p
# 
# 
# ggsave(p, file = 'fig2/donut-plot.svg', height = 4, width = 4, dpi = 300)



