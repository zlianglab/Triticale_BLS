# Load library
library(tidyverse)
library(ggplot2)

# Load dataset
cpm <- read.delim('./data/CPM_total_kallisto_tximport.txt', sep = ' ', header = TRUE)

# LBR v Cont DEG
lbr24 <- read.delim('./data/LBR_v_Cont/DEGs_LBR_v_Cont_24.txt', sep = ' ', header = TRUE)
lbr48 <- read.delim('./data/LBR_v_Cont/DEGs_LBR_v_Cont_48.txt', sep = ' ', header = TRUE)
lbr72 <- read.delim('./data/LBR_v_Cont/DEGs_LBR_v_Cont_72.txt', sep = ' ', header = TRUE)
lbr96 <- read.delim('./data/LBR_v_Cont/DEGs_LBR_v_Cont_96.txt', sep = ' ', header = TRUE)

lbrdeg24 <- lbr24$Feature
lbrdeg48 <- lbr48$Feature
lbrdeg72 <- lbr72$Feature
lbrdeg96 <- lbr96$Feature

# PR v Cont DEG
pr24 <- read.delim('./data/PR_v_Cont/DEGs_PR_v_Cont_24.txt', sep = ' ', header = TRUE) 
pr48 <- read.delim('./data/PR_v_Cont/DEGs_PR_v_Cont_48.txt', sep = ' ', header = TRUE)
pr72 <- read.delim('./data/PR_v_Cont/DEGs_PR_v_Cont_72.txt', sep = ' ', header = TRUE)
pr96 <- read.delim('./data/PR_v_Cont/DEGs_PR_v_Cont_96.txt', sep = ' ', header = TRUE)

prdeg24 <- pr24$Feature
prdeg48 <- pr48$Feature
prdeg72 <- pr72$Feature
prdeg96 <- pr96$Feature

# Combined DEG
deg <- Reduce(union, list(lbrdeg24, lbrdeg48, lbrdeg72, lbrdeg96, prdeg24, prdeg48, prdeg72, prdeg96))

View(deg)

meangeneexp <- cpm %>%
  filter(Feature %in% deg) %>%
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
  mutate(
    LBR24 = log2((LBR24 + 1) / (Cont24 + 1)),
    LBR48 = log2((LBR48 + 1) / (Cont48 + 1)),
    LBR72 = log2((LBR72 + 1) / (Cont72 + 1)),
    LBR96 = log2((LBR96 + 1) / (Cont96 + 1)),
    
    PR24 = log2((PR24 + 1) / (Cont24 + 1)),
    PR48 = log2((PR48 + 1) / (Cont48 + 1)),
    PR72 = log2((PR72 + 1) / (Cont72 + 1)),
    PR96 = log2((PR96 + 1) / (Cont96 + 1))
  ) %>%
  select(Feature, LBR24, LBR48, LBR72, LBR96, PR24, PR48, PR72, PR96)

meangeneexp %>% View()

# Scaling so that mean = 1 and sd = 0 for each feature across 8 different columns

scale_rows <- function(df, cols) {
  df[cols] <- t(apply(df[cols], 1, function(x) {
    if(sd(x) == 0) rep(0, length(x)) else as.vector(scale(x))
  }))
  return(df)
}

data_scaled <- scale_rows(meangeneexp, c("LBR24", "LBR48", "LBR72", "LBR96", "PR24", "PR48", "PR72", "PR96")) %>%
  select(LBR24, LBR48, LBR72, LBR96, PR24, PR48, PR72, PR96)


write.table(file = './result/Cont_v_LBR_PR/scaled_CPM_DEGs_LBR_PR_combined.txt', data_scaled, sep = ' ', row.names = FALSE)

data_scaled %>% View()

# ------------------Calculate Within Sum of Squares for K values----------------

# Set seed for reporducibility
set.seed(2105)

# Function to calculate Within Cluster Sum of Squares (WSS) for a range of k values
wss <- function(k) {
  kmeans(data_scaled, centers = k, nstart = 10, iter.max = 100)$tot.withinss
}

# Calculate WSS for k = 1 to 20
k_values <- 1:20
wss_values <- sapply(k_values, wss)

# Initialize an empty vector to store proportions
proportions <- numeric(length(wss_values) - 1)

# Loop through wss_values starting from the second value
for (i in 2:length(wss_values)) {
  proportions[i] <-(wss_values[i] - wss_values[i-1]) / wss_values[i - 1]  # Current value / Previous value
}  

## Save base R elbow plot
elbow_df <- data.frame(
  k = k_values,
  wss = wss_values
)

elbow_plot <- ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(
    breaks = elbow_df$k
  ) +
  labs(
    x = 'Number of clusters (k)',
    y = 'Total within-cluster sum of squares',
    title = 'Elbow Method for Optimal k'
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.3),
    plot.title = element_text(hjust = 0.5, size = 15)
  )

elbow_plot

ggsave('./result/Rerun/elbow_plot_for_cluster_k_16.png', elbow_plot, width= 6, height= 6, dpi = 300)

## First and second derivative plot
first_derivative <- diff(wss_values) * (-1)  # absolute drop in WSS
second_derivative <- diff(first_derivative)  # change in WSS drop (inflection)
first_derivative
second_derivative

# Define the corresponding k-values
k_first <- 2:length(wss_values)         # k=2:12
k_second <- 3:length(wss_values)        # k=3:12

# Create data frames
first_derivative_df <- data.frame(
  k = k_first,
  first_derivative = first_derivative
)

second_derivative_df <- data.frame(
  k = k_second,
  second_derivative = second_derivative
)

# Print them to console
print(first_derivative_df)
print(second_derivative_df)

## Save second derivative proportion plot for wss
second_derivative_df

second_derivative_plot <- ggplot(second_derivative_df, aes(x = k, y = second_derivative)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(
    breaks = second_derivative_df$k
  ) +
  labs(
    x = 'Number of clusters (k)',
    y = 'Second Derivative of WSS',
    title = 'Second Derivative of WSS Curve'
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', linewidth = 0.3),
    plot.title = element_text(hjust = 0.5, size = 15)
  )

second_derivative_plot

ggsave('./result/Rerun/second_derovative_plot_for_cluster_k_16.png', second_derivative_plot, width= 6, height= 6, dpi = 300)

# ------------------------------ K-Means Clustering-----------------------------

# Set seed for reporducibility
set.seed(2105)

# Run K-Means clustering
kmeans_result <- kmeans(data_scaled, centers = 16, nstart = 100, iter.max = 10000, algorithm = 'Lloyd')

# Add cluster assignment to the data (k = 16)
data_scaled$Cluster <- kmeans_result$cluster
data_scaled$Feature <- meangeneexp$Feature
data_scaled <- data_scaled %>% select(Feature, LBR24, LBR48, LBR72, LBR96, PR24, PR48, PR72, PR96, Cluster)

data_scaled %>% View()

write.table(data_scaled, file = './result/Rerun/expgene_cpm_cluster_log2FCscale_PR_k_16.txt', sep = ' ', row.names = FALSE)
export <- data_scaled %>% select(Feature, Cluster)
export %>% View()
write.table(export, file = './result/Rerun/Cont-v-LBR-v-PR-cluster-genelist-k-16.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

# Read dataset
export <- read.delim('./result/Rerun/Cont-v-LBR-v-PR-cluster-genelist-k-16.tsv', sep = '\t', header = TRUE)
head(export)

data_scaled$Cluster <- export$Cluster
data_scaled$Feature <- meangeneexp$Feature

data_scaled <- data_scaled %>% select(Feature, LBR24, LBR48, LBR72, LBR96, PR24, PR48, PR72, PR96, Cluster)
data_scaled %>% View()

#### Gene Counts in each clusters
gene_counts <- data_scaled %>%
  group_by(Cluster) %>%
  summarise(n_genes = n())
View(gene_counts)

# Label for the facets
cluster_labels <- setNames(
  paste0('Cluster ', gene_counts$Cluster, ': ', gene_counts$n_genes),
  gene_counts$Cluster
)

# Converting to long format
data_long <- data_scaled %>%
  pivot_longer(
    cols = c(LBR24, LBR48, LBR72, LBR96, PR24, PR48, PR72, PR96),
    names_to = 'Condition_Time',
    values_to = 'Log2FC'
  ) %>%
  separate(Condition_Time, into = c('Condition', 'Time'), sep = '(?<=R)(?=\\d)') %>%
  mutate(
    Time = paste0('H', Time),
  )

View(data_long)

# Create summary statistics
summary <- data_long %>%
  group_by(Cluster, Condition, Time) %>%
  summarise(
    mean = mean(Log2FC, na.rm = TRUE),
    sd = sd(Log2FC, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(Timepoint = Time)

summary

# Create cluster labels (adjust as needed)
cluster_labels <- setNames(
  paste("Cluster", sort(unique(data_scaled$Cluster))),
  sort(unique(data_scaled$Cluster))
)

# Create the plot
clustering_plot <- ggplot(data_long, aes(x = Time, y = Log2FC, group = interaction(Feature, Condition))) +
  # Add ribbons for both conditions
  geom_ribbon(
    data = summary %>% filter(Condition == "LBR"),
    aes(x = Timepoint, ymin = mean - sd, ymax = mean + sd, group = Cluster),
    fill = '#0072BC', alpha = 0.1, inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = summary %>% filter(Condition == "PR"),
    aes(x = Timepoint, ymin = mean - sd, ymax = mean + sd, group = Cluster),
    fill = '#7338A0', alpha = 0.1, inherit.aes = FALSE
  ) +
  # Add mean lines for both conditions
  stat_summary(
    data = data_long %>% filter(Condition == "LBR"),
    aes(group = Cluster),
    fun = mean, geom = 'line', color = '#007191', linewidth = 0.5
  ) +
  stat_summary(
    data = data_long %>% filter(Condition == "PR"),
    aes(group = Cluster),
    fun = mean, geom = 'line', color = '#d31f11', linewidth = 0.5
  ) +
  scale_x_discrete(breaks = c("H24", "H48", "H72", "H96"),
                   labels = c("24", "48", "72", "96")) +
  facet_wrap(~ Cluster, ncol = 16, strip.position = "top", scales = 'free_y',
             labeller = labeller(Cluster = cluster_labels)) +
  theme_void() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    # axis.text.x = element_text(color = "black", size = 10),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.ticks.x = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0, "lines"),
    # strip.text = element_text(size = 10, face = 'bold'),
    strip.text = element_blank(),
    strip.background = element_rect(fill = 'white', color = 'white')
  ) +
  labs(x = '')



clustering_plot


ggsave('./result/Rerun/cluster_rerun.svg', clustering_plot, width = 16, height = 1.5, dpi = 600)
ggsave('./result/Rerun/cluster_long2.png', clustering_plot, width = 6, height = 6, dpi = 600)
# Create the plot with individual gene trajectories only

clustering_plot <- ggplot(data_long, aes(x = Time, y = Log2FC, group = interaction(Feature, Condition), color = Condition)) + 
  # geom_line(alpha = 0.05, size = 0.1) + 
  scale_color_manual(values = c("LBR" = '#007191', "PR" = '#d31f11')) + 
  scale_x_discrete(breaks = c("H24", "H48", "H72", "H96"), labels = c("24", "48", "72", "96")) +
  stat_summary(
    data = data_long %>% filter(Condition == "LBR"),
    aes(group = Cluster),
    fun = mean, geom = 'line', color = '#007191', linewidth = 1.0
  ) +
  stat_summary(
    data = data_long %>% filter(Condition == "PR"),
    aes(group = Cluster),
    fun = mean, geom = 'line', color = '#d31f11', linewidth = 1.0
  ) +
  facet_wrap(~ Cluster, ncol = 4, strip.position = "top", scales = "free_y", 
             labeller = labeller(Cluster = cluster_labels)) + 
  theme_void() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(color = "black"),
    panel.spacing = unit(0, "lines"),
    strip.text = element_text(size = 10, face = 'bold'),
    strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = "none"
  ) + 
  labs(x = '')

# scale_color_manual(values = c("LBR" = '#007191', "PR" = '#d31f11')) +

clustering_plot

ggsave('./result/Rerun/cluster_plot_k_16_Cont_LBR_v_PR_poster_5.png', clustering_plot, width = 8, height = 8, dpi = 300)






