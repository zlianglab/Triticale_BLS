# load library
library(tidyverse)
library(ggplot2)
library(patchwork)
library(svglite)


# Load dataset
# go_data <- data.frame(
#   Cluster = 1:16,
#   GO = c(22, 40, 265, 70, 68, 16, 85, 129, 109, 160, 17, 151, 67, 221, 234, 76),
#   Collapse = c(11, 9, 148, 19, 14, 5, 19, 24, 22, 24, 4, 36, 11, 54, 48, 13)
# )

go_data <- read.csv('fig3/rerun/rep-ratio.csv', header = TRUE)
head(go_data)

# Calculate removed GO terms
go_data <- go_data %>%
  mutate(Removed = GO - Collapse)

# Create a list to store individual plots
plot_list <- list()

# Generate pie chart for each cluster
for(i in 1:16) {
  # Prepare data for this cluster
  pie_data <- data.frame(
    Category = c("Collapsed", "Removed"),
    Count = c(go_data$Collapse[i], go_data$Removed[i])
  ) %>%
    mutate(
      Percentage = Count / sum(Count) * 100,
      Label = paste0(Count, "\n(", round(Percentage, 1), "%)")
    )
  
  # Create pie chart
  p <- ggplot(pie_data, aes(x = "", y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 0) +
    coord_polar("y", start = 0) +
    # geom_text(aes(label = Label), 
    #           position = position_stack(vjust = 0.5),
    #           size = 3, fontface = "bold") +
    scale_fill_manual(values = c("Collapsed" = "#AFD1E6", 
                                 "Removed" = "#519BCB")) +
    theme_void() +
    theme(
      plot.title = element_blank(),
      legend.position = "none"
      # plot.margin = margin(5, 5, 5, 5)
    ) +
    labs(title = paste0("Cluster ", i))
  
  plot_list[[i]] <- p
}

# Arrange all plots in a single row
combined_plot <- wrap_plots(plot_list, nrow = 1)

# Display the plot
print(combined_plot)

 # Save as high-resolution image
ggsave("fig3/GO_piecharts_horizontal_rerun.svg", 
       plot = combined_plot, 
       width = 16, 
       height = 5, 
       dpi = 300)
