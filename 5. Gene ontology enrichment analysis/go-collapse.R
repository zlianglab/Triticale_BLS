# Load library
library(tidyverse)

# Load dataset
rep <- read.csv('data/rerun/representative.csv', header = TRUE)
cluster <- read.csv('data/rerun/go-cluster.csv', header = TRUE)

head(rep)
head(cluster)

df <- rep %>%
  left_join(cluster, by = c('TermID' = 'GO')) %>% arrange(Cluster)
head(df)
View(df)

# Filtering NA and null
df <- df %>% filter(Representative != 'null') %>% filter(Cluster != 'NA')

write.table(df, file = 'result/rerun/cluster-rep-go.tsv', row.names = FALSE, quote = FALSE, sep = '\t')

# Finding number of GO in each cluster each cluster
df %>% count(Cluster)

# Finding number of representatives in each cluster
df %>% count(Cluster, Representative) %>% filter(n >= 3) %>% arrange(Cluster, desc(n)) # Filtering n less than 2 and arranging by number

# Pirvot the data to wider format
wide <- df %>%
  count(Cluster, Representative) %>%
  filter(n >= 3) %>%
  # Ensure all clusters 1-16 are represented
  complete(Cluster = 1:16, Representative, fill = list(n = 0)) %>%
  pivot_wider(names_from = Cluster,
              values_from = n,
              values_fill = 0) %>%
  arrange(Representative)

head(wide)
View(wide)

# Create the "Others" row from filtered out data (n < 2)
others_row <- df %>%
  count(Cluster, Representative) %>%
  filter(n < 3) %>%
  group_by(Cluster) %>%
  summarise(n = sum(n), .groups = 'drop') %>%
  pivot_wider(names_from = Cluster,
              values_from = n,
              values_fill = 0) %>%
  mutate(Representative = "Others", .before = 1)

# Combine the main table with the Others row
wide_combined <- bind_rows(wide, others_row)

print(wide_combined, n = Inf)
print(wide, n = Inf)
print(others_row, n = Inf)

wide_filtered <- wide %>% filter(rowSums(select(., where(is.numeric))) >= 3)

cellular_machinery <- c(
  "Golgi apparatus subcompartment",
  "Golgi organization", 
  "ribosome",
  "stromule", 
  # "light-harvesting complex",
  "external encapsulating structure"
)

wide_filtered <- wide_filtered %>%
  filter(!Representative %in% cellular_machinery)

# Convert to ratios (proportions) across columns
wide_ratios <- wide_filtered %>%
  mutate(
    across(-Representative, ~ .x / rowSums(select(cur_data(), -Representative), na.rm = TRUE))
  ) %>%
  # Replace NaN values (from 0/0) with 0
  mutate(across(-Representative, ~ ifelse(is.nan(.), 0, .)))

print(wide_ratios, n = Inf)

# Convert wide format to long format
long_data <- wide_ratios %>%
  pivot_longer(
    cols = -Representative,
    names_to = "Cluster",
    values_to = "Count"
  ) %>%
  mutate(Cluster = as.numeric(Cluster)) %>%
  arrange(Cluster, Representative)  # Sort by cluster first

View(long_data)

write.table(long_data, 'result/rerun/go_collapse.txt', sep = '\t', quote = FALSE, row.names = FALSE)
