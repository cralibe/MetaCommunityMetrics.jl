# benchmarks/benchmark_r.R

# Load necessary libraries
library(data.table)  
library(tidyverse)
library(bench)
## R packages that provide equivalent functions of my package  
library(adespatial) 
library(DNCImper)

#Read in the sample data
df <- read.csv("~/.julia/dev/MetaCommunityMetrics/data/metacomm_rodent_df.csv")
comm_df_for_DNCI<-read.csv("~/.julia/dev/MetaCommunityMetrics/benchmarks/benchmark_r/data/DNCI_comm.csv")
grouping_for_DNCI<-read.csv("~/.julia/dev/MetaCommunityMetrics/benchmarks/benchmark_r/data/cluster_list_t1.csv")

#Data Wrangling
community_matrix<-df %>%
  select(-Presence) %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

## The abundance matrix
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

####Benchmark the spatial beta diversity function
spatial_beta_div_1 <- mark(df %>% 
                           group_by(plot,Species) %>%
                           dplyr::summarise(Abundance = sum(Abundance)) %>% 
                           spread(key = Species, value = Abundance, fill = 0) %>% 
                           ungroup() %>% 
                           dplyr::select(-plot) %>% 
                           beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE),
                          iterations = 1000,
                          check = TRUE,
                          time_unit = "us")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_millisecond <- as.numeric(mean(spatial_beta_div_1$time[[1]])) * 1000
memory_usage_mib <- as.numeric(spatial_beta_div_1$mem_alloc) / 1.048576e+6

# Print the results
cat("Execution Time (Milliseconds):", execution_time_millisecond, "\n")
cat("Memory Usage (MiB):", memory_usage_mib, "\n")

spatial_beta_div_2 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Abundance = sum(Abundance)) %>% 
                             spread(key = Species, value = Abundance, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                           iterations = 1000,
                           check = TRUE,
                           time_unit = "us")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_millisecond <- as.numeric(mean(spatial_beta_div_2$time[[1]])) * 1000
memory_usage_mib <- as.numeric(spatial_beta_div_2$mem_alloc) / 1.048576e+6

# Print the results
cat("Execution Time (Milliseconds):", execution_time_millisecond, "\n")
cat("Memory Usage (MiB):", memory_usage_mib, "\n")


spatial_beta_div_3 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                             spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                           iterations = 1000,
                           check = TRUE,
                           time_unit = "ms")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_millisecond <- as.numeric(mean(spatial_beta_div_3$time[[1]])) * 1000
memory_usage_mib <- as.numeric(spatial_beta_div_3$mem_alloc) / 1.048576e+6

# Print the results
cat("Execution Time (Milliseconds):", execution_time_millisecond, "\n")
cat("Memory Usage (MiB):", memory_usage_mib, "\n")


temporal_beta_div_1 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              summarise(Abundance = sum(Abundance)) %>% 
                              spread(key = Species, value = Abundance, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE),
                           iterations = 1000,
                           check = TRUE,
                           time_unit = "ms")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_millisecond <- as.numeric(mean(temporal_beta_div_1$time[[1]])) * 1000
memory_usage_mib <- as.numeric(temporal_beta_div_1$mem_alloc) / 1.048576e+6

# Print the results
cat("Execution Time (Milliseconds):", execution_time_millisecond, "\n")
cat("Memory Usage (MiB):", memory_usage_mib, "\n")

temporal_beta_div_2 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              summarise(Abundance = sum(Abundance)) %>% 
                              spread(key = Species, value = Abundance, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                            iterations = 1000,
                            check = TRUE,
                            time_unit = "ms")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_millisecond <- as.numeric(mean(temporal_beta_div_2$time[[1]])) * 1000
memory_usage_mib <- as.numeric(temporal_beta_div_2$mem_alloc) / 1.048576e+6

# Print the results
cat("Execution Time (Milliseconds):", execution_time_millisecond, "\n")
cat("Memory Usage (MiB):", memory_usage_mib, "\n")

temporal_beta_div_3 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                              spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                            iterations = 1000,
                            check = TRUE,
                            time_unit = "ms")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_millisecond <- as.numeric(mean(temporal_beta_div_3$time[[1]])) * 1000
memory_usage_mib <- as.numeric(temporal_beta_div_3$mem_alloc) / 1.048576e+6

# Print the results
cat("Execution Time (Milliseconds):", execution_time_millisecond, "\n")
cat("Memory Usage (MiB):", memory_usage_mib, "\n")







#### Benchmark the DNCI function
groups <- grouping_for_DNCI$Group

DNCI_multigroup_result <- mark(DNCImper:::DNCI_multigroup(comm_df_for_DNCI, 
                                                       groups, Nperm = 100, 
                                                       symmetrize = FALSE, 
                                                       plotSIMPER = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_millisecond <- as.numeric(mean(temporal_beta_div_3$time[[1]])) * 1000
memory_usage_mib <- as.numeric(temporal_beta_div_3$mem_alloc) / 1.048576e+6

# Print the results
cat("Execution Time (Milliseconds):", execution_time_millisecond, "\n")
cat("Memory Usage (MiB):", memory_usage_mib, "\n")




#### Benchmark the prop_patches function
prop_patches_result <- mark(prop_patches <- df %>% 
                              group_by(Species, plot) %>%
                              dplyr::summarise(mean_abundance = mean(Abundance)) %>% 
                              filter(mean_abundance  > 0) %>% 
                              dplyr::summarise(n_patches = n()) %>% 
                              mutate(prop_patches = n_patches/max(df$plot)) %>% 
                              summarise(mean_prop_patches = mean(prop_patches), min_prop_patches = min(prop_patches), max_prop_patches = max(prop_patches)),
                            iterations = 1000,
                            check = TRUE,
                            time_unit = "ms")

# Convert execution time to microseconds (µs) and memory allocation to kibibytes (KiB)
execution_time_millisecond <- as.numeric(mean(prop_patches_result$time[[1]])) * 1000
memory_usage_mib <- as.numeric(prop_patches_result$mem_alloc) / 1.048576e+6

# Print the results
cat("Execution Time (Milliseconds):", execution_time_millisecond, "\n")
cat("Memory Usage (MiB):", memory_usage_mib, "\n")

#### Benchmark the CV_meta_simple function
# Extract unique values for Species, Sampling_date_order, and plot
species_vals <- unique(df$Species)
date_vals <- unique(df$Sampling_date_order)
plot_vals <- unique(df$plot)

# Create the array with dimensions based on unique Species, Sampling_date_order, and plot
metacomm_tsdata <- array(0, dim = c(length(species_vals), length(date_vals), length(plot_vals)))

# Step 3: Populate the array with values from the data frame
for(i in 1:nrow(df)){
  # Map Species, Sampling_date_order, and plot to their corresponding indices
  species_index <- which(species_vals == df$Species[i])
  date_index <- which(date_vals == df$Sampling_date_order[i])
  plot_index <- which(plot_vals == df$plot[i])
  
  # Assign the value from the data frame (Abundance) to the array
  metacomm_tsdata[species_index, date_index, plot_index] <- df$Abundance[i]
}

var.partition <- function(metacomm_tsdata){
  ## The function "var.partition" performs the partitioning of variability
  ## across hierarchical levesl within a metacommunity.
  ## The input array "metacomm_tsdata" is an N*T*M array. The first dimension represents N species,
  ## the second represents time-series observations of length T, and the third represents M local communities.
  ## The output includes four variability and four synchrony metrics as defined in the main text.
  ## Note that, to be able to handle large metacommunities, this code has avoided calculating all covariance.
  ts_metacom <- apply(metacomm_tsdata,2,sum)
  ts_patch <- apply(metacomm_tsdata,c(2,3),sum)
  ts_species <- apply(metacomm_tsdata,c(1,2),sum)
  sd_metacom <- sd(ts_metacom)
  sd_patch_k <- apply(ts_patch,2,sd)
  sd_species_i <- apply(ts_species,1,sd)
  sd_species_patch_ik <- apply(metacomm_tsdata,c(1,3),sd)
  mean_metacom <- mean(ts_metacom)
  CV_S_L <- sum(sd_species_patch_ik)/mean_metacom
  CV_C_L <- sum(sd_patch_k)/mean_metacom
  CV_S_R <- sum(sd_species_i)/mean_metacom
  CV_C_R <- sd_metacom/mean_metacom

  summary <- c(CV_S_L=CV_S_L, CV_C_L=CV_C_L, CV_S_R=CV_S_R, CV_C_R=CV_C_R)
  return(summary)
}

test=var.partition(metacomm_tsdata)




