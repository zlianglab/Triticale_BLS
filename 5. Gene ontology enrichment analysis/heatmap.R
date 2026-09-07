# Load library
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(svglite)
library(reshape2)
library(RColorBrewer)

# Load dataset
# dat = read.table("fig3/PostClusterGO.txt",sep='\t',head=T)
dat = read.delim('fig3/rerun/go_collapse.txt', sep = '\t', header = TRUE)
head(dat)


View(dat)

# Multiply gene ratio
dat <- dat %>%
  mutate(Count = Count * 100)

# Matrix formatting
mat <- dat %>%
  select(Cluster, Representative, Count) %>%
  pivot_wider(names_from = Cluster, values_from = Count, values_fill = 0) %>%
  column_to_rownames('Representative') %>%
  as.matrix()

head(mat)

print(mat)

# Normalize per row
mat_norm <- t(apply(mat, 1, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}))

head(mat_norm)

# Reshape for plotting
plot_data <- melt(as.matrix(mat_norm), varnames = c('Function', 'Cluster'), value.name = 'Ratio')
head(plot_data)
plot_data$Function <- str_to_sentence(plot_data$Function)

head(plot_data)


# Create plot
p <- ggplot(plot_data, aes(x = Cluster, y = Function, fill = Ratio)) +
  geom_tile(color = 'white', size = 1) +
  scale_fill_gradientn(colors = colorRampPalette(brewer.pal(8, "Blues"))(25)) +
  scale_x_discrete(expand = c(0,0), limits = 1:16) +
  scale_y_discrete(expand = c(0,0)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.position = 'none',
    # legend.direction = 'vertical',
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.width = unit(0.5, 'cm'),
    legend.key.height = unit(1.5, 'cm'),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 5, 'mm')
  )

p <- ggplot(plot_data, aes(x = Cluster, y = Function, fill = Ratio)) +
  geom_tile(color = 'white', size = 1) +
  scale_fill_gradientn(colors = colorRampPalette(brewer.pal(8, "Blues"))(25)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0), limits = 1:16, position = 'top') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 6),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 6),
    axis.title = element_blank(),
    legend.position = 'none',
    # legend.position = 'right',
    # legend.direction = 'vertical',
    # legend.title = element_blank(),
    # legend.text = element_text(size = 12),
    # legend.key.width = unit(0.5, 'cm'),
    # legend.key.height = unit(1.5, 'cm'),
    panel.grid = element_blank(),
    # plot.margin = margin(5, 5, 5, 5, 'mm')
  )

p

ggsave(p, file = 'fig3/heatmap_rerun.svg', height = 6, width = 12, dpi = 300)









