# benchmarks/benchmark_r.R

# Load necessary libraries
library(microbenchmark)
library(adespatial)  # Assuming this is the package with beta.div.comp
library(data.table)  # For reading CSV files
library(tidyverse)
library(bench)
#Read in the sample data
df <- read.csv("~/.julia/dev/MetaCommunityMetrics/data/metacomm_rodent_df.csv")

# Benchmarking R functions
matrix_with_abundance <- df %>%
  filter(Sampling_date_order == 50) %>% 
  select(-Presence) %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

### The binary matrix
matrix_with_presence <- df %>%
  filter(Sampling_date_order == 50) %>% 
  select(-Abundance) %>%
  pivot_wider(names_from = Species, values_from = Presence, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

#### Benchmark the `beta.div.comp` function
beta_diversity_1 <- mark(beta.div.comp(matrix_with_abundance, coef = "J", quant = TRUE),
  iterations = 10000,
  check = TRUE,
  time_unit = "us")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_microseconds <- as.numeric(mean(beta_diversity_1$time[[1]])) * 1e+6
memory_usage_kib <- as.numeric(beta_diversity_1$mem_alloc) / 1024

# Print the results
cat("Execution Time (Microseconds):", execution_time_microseconds, "\n")
cat("Memory Usage (KiB):", memory_usage_kib, "\n")

beta_diversity_2 <- mark(beta.div.comp(matrix_with_abundance, coef = "J", quant = FALSE),
                         iterations = 10000,
                         check = TRUE,
                         time_unit = "us")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_microseconds <- mean(beta_diversity_2$time[[1]])* 1e+6
memory_usage_kib <- as.numeric(beta_diversity_2$mem_alloc) / 1024

# Print the results
cat("Execution Time (Microseconds):", execution_time_microseconds, "\n")
cat("Memory Usage (KiB):", memory_usage_kib, "\n")

beta_diversity_3 <- mark(beta.div.comp(matrix_with_presence, coef = "J", quant = FALSE),
                         iterations = 10000,
                         check = TRUE,
                         time_unit = "us")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_microseconds <- as.numeric(mean(beta_diversity_3$time[[1]])) * 1e+6
memory_usage_kib <- as.numeric(beta_diversity_3$mem_alloc) / 1024

# Print the results
cat("Execution Time (Microseconds):", execution_time_microseconds, "\n")
cat("Memory Usage (KiB):", memory_usage_kib, "\n")
