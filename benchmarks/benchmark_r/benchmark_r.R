# benchmarks/benchmark_r.R
# Load necessary libraries
library(data.table)  
library(tidyverse)
library(bench)
## R packages that provide equivalent functions of my package  
library(adespatial) 
library(DNCImper)
library(MVNH)

#Read in the sample data
full_df <- read.csv("metacomm_rodent_df.csv",na.strings = c(""), stringsAsFactors = FALSE) #There is a species called "NA", need to handle NAs with caution
comm_full<- read.csv("data_for_testing/comm_full_df.csv")
groups_full<- read.csv("data_for_testing/groups_full_df.csv")

medium_df <- read.csv("data_for_testing/medium_dataset.csv")
comm_medium <- read.csv("data_for_testing/comm_medium_df.csv")
groups_medium <- read.csv("data_for_testing/groups_medium_df.csv")

small_df <- read.csv("data_for_testing/small_dataset.csv")
comm_small <- read.csv("data_for_testing/comm_small_df.csv")
groups_small <-read.csv("data_for_testing/groups_small_df.csv")

#Data Wrangling
#Full Dataset####
df <- full_df%>%
  mutate(Species = ifelse(Species=="NA", "N_A", Species))
## The abundance matrix
matrix_with_abundance <- df %>%
  filter(Sampling_date_order == 50) %>% 
  select(-Presence) %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude, normalized_temperature, normalized_precipitation)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

### The presence/absence matrix
matrix_with_presence <- df %>%
  filter(Sampling_date_order == 50) %>% 
  select(-Abundance) %>%
  pivot_wider(names_from = Species, values_from = Presence, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude, normalized_temperature, normalized_precipitation)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

#### Benchmark the `beta.div.comp` function
beta_diversity_1 <- mark(beta.div.comp(matrix_with_abundance, coef = "J", quant = TRUE),
  iterations = 100,
  check = TRUE,
  time_unit = "ms")

beta_diversity_2 <- mark(beta.div.comp(matrix_with_abundance, coef = "J", quant = FALSE),
                         iterations = 100,
                         check = TRUE,
                         time_unit = "ms")

beta_diversity_3 <- mark(beta.div.comp(matrix_with_presence, coef = "J", quant = FALSE),
                         iterations = 100,
                         check = TRUE,
                         time_unit = "ms")

####Benchmark the spatial beta diversity function
spatial_beta_div_1 <- mark(df %>% 
                           group_by(plot,Species) %>%
                           dplyr::summarise(Abundance = sum(Abundance)) %>% 
                           spread(key = Species, value = Abundance, fill = 0) %>% 
                           ungroup() %>% 
                           dplyr::select(-plot) %>% 
                           beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE),
                          iterations = 100,
                          check = TRUE,
                          time_unit = "ms")

spatial_beta_div_2 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Abundance = sum(Abundance)) %>% 
                             spread(key = Species, value = Abundance, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")

spatial_beta_div_3 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                             spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")
####Benchmark the temproal beta diversity function
temporal_beta_div_1 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              summarise(Abundance = sum(Abundance)) %>% 
                              spread(key = Species, value = Abundance, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")

temporal_beta_div_2 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              summarise(Abundance = sum(Abundance)) %>% 
                              spread(key = Species, value = Abundance, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

temporal_beta_div_3 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                              spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

#### Benchmark the DNCI function
DNCI_multigroup_result <- mark(DNCImper:::DNCI_multigroup(comm_full,
                                                       groups_full$Group, Nperm = 100,
                                                       symmetrize = FALSE,
                                                       plotSIMPER = FALSE,
                                                      parallelComputing = TRUE),
                               iterations = 100,
                               check = FALSE,
                            time_unit = "ms")
saveRDS(DNCI_multigroup_result, "../benchmarks/result/DNCI_full_result.rds")
DNCI_multigroup_result<-readRDS("../benchmarks/result/DNCI_full_result.rds")
#### Benchmark the prop_patches function
prop_patches_result <- mark(df %>% 
                              group_by(Species, plot) %>%
                              dplyr::summarise(mean_abundance = mean(Abundance)) %>% 
                              filter(mean_abundance  > 0) %>% 
                              dplyr::summarise(n_patches = n()) %>% 
                              mutate(prop_patches = n_patches/max(df$plot)) %>% 
                              summarise(mean_prop_patches = mean(prop_patches), min_prop_patches = min(prop_patches), max_prop_patches = max(prop_patches)),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

#### Benchmark the CV_meta function
# The R function “var.partition” to calculate variability and synchrony across hierarchical levels
var.partition <- function(metacomm_tsdata){
  ## The function "var.partition" performs the partitioning of variability
  ## across hierarchical levesl within a metacommunity.
  ## The input array "metacomm_tsdata" is an N*T*M array. The first dimension represents N species,
  ## the second represents time-series observations of length T, and the third represents M local communities.
  ## The output includes four variability and four synchrony metrics as defined in the main text.
  ## Note that, to be able to handle large metacommunities, this code has avoided calculating all covariance.
  ts_metacom <- apply(metacomm_tsdata,2,sum) #summing the total biovolume across all species and patch at every time point
  ts_patch <- apply(metacomm_tsdata,c(2,3),sum) #summing the total biovolume across all species in every patch at every time point
  ts_species <- apply(metacomm_tsdata,c(1,2),sum)#summing the total biovolume across all patches for every species at every time point
  sd_metacom <- sd(ts_metacom)
  sd_patch_k <- apply(ts_patch,2,sd)
  sd_species_i <- apply(ts_species,1,sd)
  sd_species_patch_ik <- apply(metacomm_tsdata,c(1,3),sd)
  mean_metacom <- mean(ts_metacom)
  CV_S_L <- sum(sd_species_patch_ik)/mean_metacom
  CV_C_L <- sum(sd_patch_k)/mean_metacom
  CV_S_R <- sum(sd_species_i)/mean_metacom
  CV_C_R <- sd_metacom/mean_metacom
  #phi_S_L2R <- CV_S_R/CV_S_L
  #phi_C_L2R <- CV_C_R/CV_C_L
  #phi_S2C_L <- CV_C_L/CV_S_L
  #phi_S2C_R <- CV_C_R/CV_S_R
  partition_3level <- c(CV_S_L=CV_S_L, CV_C_L=CV_C_L, CV_S_R=CV_S_R, CV_C_R=CV_C_R)
                        #,
                        #phi_S_L2R=phi_S_L2R, phi_C_L2R=phi_C_L2R, phi_S2C_L=phi_S2C_L,
                        #phi_S2C_R=phi_S2C_R)
  return(partition_3level)
}

#The master funciton to be benchmarked
CV_meta <- function(data, species, time, plot, abundance){
  # Extract unique values for Species, Sampling_date_order, and plot
  species_vals <- unique(species)
  date_vals <- unique(time)
  plot_vals <- unique(plot)
  
  # Create the array with dimensions based on unique Species, Sampling_date_order, and plot
  metacomm_tsdata <- array(0, dim = c(length(species_vals), length(date_vals), length(plot_vals)))
  
  # Populate the array with values from the data frame
  for(i in 1:nrow(data)){
    # Map Species, Sampling_date_order, and plot to their corresponding indices
    species_index <- which(species_vals == species[i])
    date_index <- which(date_vals == time[i])
    plot_index <- which(plot_vals == plot[i])
    
    # Assign the value from the data frame (Abundance) to the array
    metacomm_tsdata[species_index, date_index, plot_index] <- abundance[i]
  }
  result <- var.partition(metacomm_tsdata)
  
  return(result)
}
  


CV_meta_result <- mark(CV_meta(df, df$Species, df$Sampling_date_order, df$plot, df$Abundance),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

#### Benchmark the hypervolume functions
data_1 <- df%>%
  filter(Presence >0)%>%
  filter(Species =="BA")%>%
  select(normalized_temperature, normalized_precipitation)

data_2 <- df%>%
  filter(Presence >0)%>%
  filter(Species =="SH")%>%
  select(normalized_temperature, normalized_precipitation)

hypervolume_det_result <- mark(MVNH_det(data_1, var.names = c("Temperature", "Precipitation")), iterations = 100,
                               check = TRUE,
                               time_unit = "ms")

hypervolume_dis_result <- mark(MVNH_dissimilarity(data_1, data_2, var.names = c("Temperature", "Precipitation")), iterations = 100,
                                check = TRUE,
                                time_unit = "ms")

#Save all the results into a dataframe
test_case_list = c("beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                   "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                   "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                   "DNCI_multigroup_result",
                   "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result")

# Initialize empty data frame
all_time_full <- data.frame()

for (testcase in test_case_list) {
# Get the object with that name
testcase_obj <- get(testcase)
  
# Get times from the object
times <- as.numeric(testcase_obj$time[[1]]) * 1e+3
  
# Create data frame for this test case
data <- data.frame(TestCase = rep(testcase, length(times)),Time = times)
  
  # Append to the full data frame - this is the key fix
  all_time_full <- rbind(all_time_full, data)
}

write.csv(all_time_full, "../benchmarks/result/all_time_full_df_r.csv")


benchmark_result_full_df<-data.frame(TestCase = c("beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                                "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                                "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                                "DNCI_multigroup_result",
                                                "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result"),
                                  Time_minimum = c(as.numeric(min(beta_diversity_1$time[[1]])) * 1e+3,
                                                  as.numeric(min(beta_diversity_2$time[[1]])) * 1e+3,
                                                  as.numeric(min(beta_diversity_3$time[[1]])) * 1e+3,
                                                  as.numeric(min(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                  as.numeric(min(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                  as.numeric(min(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                  as.numeric(min(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                  as.numeric(min(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                  as.numeric(min(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                  as.numeric(min(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                  as.numeric(min(prop_patches_result$time[[1]])) * 1e+3,
                                                  as.numeric(min(CV_meta_result$time[[1]])) * 1e+3,
                                                  as.numeric(min(hypervolume_det_result$time[[1]])) * 1e+3,
                                                  as.numeric(min(hypervolume_dis_result$time[[1]])) * 1e+3),
                                  Time_median = c(as.numeric(median(beta_diversity_1$time[[1]])) * 1e+3,
                                                   as.numeric(median(beta_diversity_2$time[[1]])) * 1e+3,
                                                   as.numeric(median(beta_diversity_3$time[[1]])) * 1e+3,
                                                   as.numeric(median(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                   as.numeric(median(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                   as.numeric(median(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                   as.numeric(median(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                   as.numeric(median(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                   as.numeric(median(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                   as.numeric(median(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                   as.numeric(median(prop_patches_result$time[[1]])) * 1e+3,
                                                   as.numeric(median(CV_meta_result$time[[1]])) * 1e+3,
                                                  as.numeric(median(hypervolume_det_result$time[[1]])) * 1e+3,
                                                  as.numeric(median(hypervolume_dis_result$time[[1]])) * 1e+3),
                                  Time_mean = c(as.numeric(mean(beta_diversity_1$time[[1]])) * 1e+3,
                                                as.numeric(mean(beta_diversity_2$time[[1]])) * 1e+3,
                                                as.numeric(mean(beta_diversity_3$time[[1]])) * 1e+3,
                                                as.numeric(mean(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                as.numeric(mean(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                as.numeric(mean(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                as.numeric(mean(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                as.numeric(mean(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                as.numeric(mean(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                as.numeric(mean(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                as.numeric(mean(prop_patches_result$time[[1]])) * 1e+3,
                                                as.numeric(mean(CV_meta_result$time[[1]])) * 1e+3,
                                                as.numeric(mean(hypervolume_det_result$time[[1]])) * 1e+3,
                                                as.numeric(mean(hypervolume_dis_result$time[[1]])) * 1e+3),
                                  Time_maximum = c(as.numeric(max(beta_diversity_1$time[[1]])) * 1e+3,
                                                as.numeric(max(beta_diversity_2$time[[1]])) * 1e+3,
                                                as.numeric(max(beta_diversity_3$time[[1]])) * 1e+3,
                                                as.numeric(max(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                as.numeric(max(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                as.numeric(max(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                as.numeric(max(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                as.numeric(max(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                as.numeric(max(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                as.numeric(max(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                as.numeric(max(prop_patches_result$time[[1]])) * 1e+3,
                                                as.numeric(max(CV_meta_result$time[[1]])) * 1e+3,
                                                as.numeric(max(hypervolume_det_result$time[[1]])) * 1e+3,
                                                as.numeric(max(hypervolume_dis_result$time[[1]])) * 1e+3),
                                  Time_std = c(as.numeric(sd(beta_diversity_1$time[[1]])) * 1e+3,
                                                   as.numeric(sd(beta_diversity_2$time[[1]])) * 1e+3,
                                                   as.numeric(sd(beta_diversity_3$time[[1]])) * 1e+3,
                                                   as.numeric(sd(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                   as.numeric(sd(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                   as.numeric(sd(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                   as.numeric(sd(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                   as.numeric(sd(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                   as.numeric(sd(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                   as.numeric(sd(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                   as.numeric(sd(prop_patches_result$time[[1]])) * 1e+3,
                                                   as.numeric(sd(CV_meta_result$time[[1]])) * 1e+3,
                                               as.numeric(sd(hypervolume_det_result$time[[1]])) * 1e+3,
                                               as.numeric(sd(hypervolume_dis_result$time[[1]])) * 1e+3),
                                  memory = c(as.numeric(beta_diversity_1$mem_alloc) / 1024^2,
                                               as.numeric(beta_diversity_2$mem_alloc) / 1024^2,
                                               as.numeric(beta_diversity_3$mem_alloc) / 1024^2,
                                               as.numeric(spatial_beta_div_1$mem_alloc) / 1024^2,
                                               as.numeric(spatial_beta_div_2$mem_alloc) / 1024^2,
                                               as.numeric(spatial_beta_div_3$mem_alloc) / 1024^2,
                                               as.numeric(temporal_beta_div_1$mem_alloc) / 1024^2,
                                               as.numeric(temporal_beta_div_2$mem_alloc) / 1024^2,
                                               as.numeric(temporal_beta_div_3$mem_alloc) / 1024^2,
                                               as.numeric(DNCI_multigroup_result$mem_alloc) / 1024^2,
                                               as.numeric(prop_patches_result$mem_alloc) / 1024^2,
                                               as.numeric(CV_meta_result$mem_alloc) / 1024^2,
                                  as.numeric(hypervolume_det_result$mem_alloc) / 1024^2,
                                  as.numeric(hypervolume_dis_result$mem_alloc) / 1024^2))
           
write.csv(benchmark_result_full_df, "../benchmarks/result/benchmark_result_full_df_r.csv")

#Medium Dataset####
df <- medium_df%>%
  mutate(Species = ifelse(Species=="NA", "N_A", Species))
## The abundance matrix
matrix_with_abundance <- df %>%
  filter(Sampling_date_order == 50) %>% 
  select(-Presence) %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude, normalized_temperature, normalized_precipitation)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

### The presence/absence matrix
matrix_with_presence <- df %>%
  filter(Sampling_date_order == 50) %>% 
  select(-Abundance) %>%
  pivot_wider(names_from = Species, values_from = Presence, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude, normalized_temperature, normalized_precipitation)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

#### Benchmark the `beta.div.comp` function
beta_diversity_1 <- mark(beta.div.comp(matrix_with_abundance, coef = "J", quant = TRUE),
                         iterations = 100,
                         check = TRUE,
                         time_unit = "ms")

beta_diversity_2 <- mark(beta.div.comp(matrix_with_abundance, coef = "J", quant = FALSE),
                         iterations = 100,
                         check = TRUE,
                         time_unit = "ms")

beta_diversity_3 <- mark(beta.div.comp(matrix_with_presence, coef = "J", quant = FALSE),
                         iterations = 100,
                         check = TRUE,
                         time_unit = "ms")

####Benchmark the spatial beta diversity function
spatial_beta_div_1 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Abundance = sum(Abundance)) %>% 
                             spread(key = Species, value = Abundance, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")

spatial_beta_div_2 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Abundance = sum(Abundance)) %>% 
                             spread(key = Species, value = Abundance, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")

spatial_beta_div_3 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                             spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")
####Benchmark the temproal beta diversity function
temporal_beta_div_1 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              summarise(Abundance = sum(Abundance)) %>% 
                              spread(key = Species, value = Abundance, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

temporal_beta_div_2 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              summarise(Abundance = sum(Abundance)) %>% 
                              spread(key = Species, value = Abundance, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

temporal_beta_div_3 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                              spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

#### Benchmark the DNCI function
DNCI_multigroup_result <- mark(DNCImper:::DNCI_multigroup(comm_medium,
                                                          groups_medium$Group, Nperm = 100,
                                                          symmetrize = FALSE,
                                                          plotSIMPER = FALSE,
                                                          parallelComputing = TRUE),
                               iterations = 100,
                               check = FALSE,
                               time_unit = "ms")
saveRDS(DNCI_multigroup_result, "../benchmarks/result/DNCI_medium_result.rds")
DNCI_multigroup_result<-readRDS("../benchmarks/result/DNCI_medium_result.rds")
#### Benchmark the prop_patches function
prop_patches_result <- mark(df %>% 
                              group_by(Species, plot) %>%
                              dplyr::summarise(mean_abundance = mean(Abundance)) %>% 
                              filter(mean_abundance  > 0) %>% 
                              dplyr::summarise(n_patches = n()) %>% 
                              mutate(prop_patches = n_patches/max(df$plot)) %>% 
                              summarise(mean_prop_patches = mean(prop_patches), min_prop_patches = min(prop_patches), max_prop_patches = max(prop_patches)),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

CV_meta_result <- mark(CV_meta(df, df$Species, df$Sampling_date_order, df$plot, df$Abundance),
                       iterations = 100,
                       check = TRUE,
                       time_unit = "ms")

#### Benchmark the hypervolume functions
data_1 <- df%>%
  filter(Presence >0)%>%
  filter(Species =="BA")%>%
  select(normalized_temperature, normalized_precipitation)

data_2 <- df%>%
  filter(Presence >0)%>%
  filter(Species =="SH")%>%
  select(normalized_temperature, normalized_precipitation)

hypervolume_det_result <- mark(MVNH_det(data_1, var.names = c("Temperature", "Precipitation")), iterations = 100,
                               check = TRUE,
                               time_unit = "ms")

hypervolume_dis_result <- mark(MVNH_dissimilarity(data_1, data_2, var.names = c("Temperature", "Precipitation")), iterations = 100,
                               check = TRUE,
                               time_unit = "ms")

#Save all the results into a dataframe
test_case_list = c("beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                   "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                   "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                   "DNCI_multigroup_result",
                   "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result")

# Initialize empty data frame
all_time_medium <- data.frame()

for (testcase in test_case_list) {
  # Get the object with that name
  testcase_obj <- get(testcase)
  
  # Get times from the object
  times <- as.numeric(testcase_obj$time[[1]]) * 1e+3
  
  # Create data frame for this test case
  data <- data.frame(TestCase = rep(testcase, length(times)),Time = times)
  
  # Append to the full data frame - this is the key fix
  all_time_medium <- rbind(all_time_medium, data)
}

write.csv(all_time_medium, "../benchmarks/result/all_time_medium_df_r.csv")

benchmark_result_medium_df<-data.frame(TestCase = c("beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                                  "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                                  "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                                  "DNCI_multigroup_result",
                                                  "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result"),
                                     Time_minimum = c(as.numeric(min(beta_diversity_1$time[[1]])) * 1e+3,
                                                      as.numeric(min(beta_diversity_2$time[[1]])) * 1e+3,
                                                      as.numeric(min(beta_diversity_3$time[[1]])) * 1e+3,
                                                      as.numeric(min(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                      as.numeric(min(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                      as.numeric(min(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                      as.numeric(min(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                      as.numeric(min(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                      as.numeric(min(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                      as.numeric(min(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                      as.numeric(min(prop_patches_result$time[[1]])) * 1e+3,
                                                      as.numeric(min(CV_meta_result$time[[1]])) * 1e+3,
                                                      as.numeric(min(hypervolume_det_result$time[[1]])) * 1e+3,
                                                      as.numeric(min(hypervolume_dis_result$time[[1]])) * 1e+3),
                                     Time_median = c(as.numeric(median(beta_diversity_1$time[[1]])) * 1e+3,
                                                     as.numeric(median(beta_diversity_2$time[[1]])) * 1e+3,
                                                     as.numeric(median(beta_diversity_3$time[[1]])) * 1e+3,
                                                     as.numeric(median(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                     as.numeric(median(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                     as.numeric(median(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                     as.numeric(median(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                     as.numeric(median(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                     as.numeric(median(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                     as.numeric(median(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                     as.numeric(median(prop_patches_result$time[[1]])) * 1e+3,
                                                     as.numeric(median(CV_meta_result$time[[1]])) * 1e+3,
                                                     as.numeric(median(hypervolume_det_result$time[[1]])) * 1e+3,
                                                     as.numeric(median(hypervolume_dis_result$time[[1]])) * 1e+3),
                                     Time_mean = c(as.numeric(mean(beta_diversity_1$time[[1]])) * 1e+3,
                                                   as.numeric(mean(beta_diversity_2$time[[1]])) * 1e+3,
                                                   as.numeric(mean(beta_diversity_3$time[[1]])) * 1e+3,
                                                   as.numeric(mean(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                   as.numeric(mean(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                   as.numeric(mean(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                   as.numeric(mean(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                   as.numeric(mean(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                   as.numeric(mean(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                   as.numeric(mean(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                   as.numeric(mean(prop_patches_result$time[[1]])) * 1e+3,
                                                   as.numeric(mean(CV_meta_result$time[[1]])) * 1e+3,
                                                   as.numeric(mean(hypervolume_det_result$time[[1]])) * 1e+3,
                                                   as.numeric(mean(hypervolume_dis_result$time[[1]])) * 1e+3),
                                     Time_maximum = c(as.numeric(max(beta_diversity_1$time[[1]])) * 1e+3,
                                                      as.numeric(max(beta_diversity_2$time[[1]])) * 1e+3,
                                                      as.numeric(max(beta_diversity_3$time[[1]])) * 1e+3,
                                                      as.numeric(max(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                      as.numeric(max(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                      as.numeric(max(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                      as.numeric(max(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                      as.numeric(max(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                      as.numeric(max(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                      as.numeric(max(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                      as.numeric(max(prop_patches_result$time[[1]])) * 1e+3,
                                                      as.numeric(max(CV_meta_result$time[[1]])) * 1e+3,
                                                      as.numeric(max(hypervolume_det_result$time[[1]])) * 1e+3,
                                                      as.numeric(max(hypervolume_dis_result$time[[1]])) * 1e+3),
                                     Time_std = c(as.numeric(sd(beta_diversity_1$time[[1]])) * 1e+3,
                                                  as.numeric(sd(beta_diversity_2$time[[1]])) * 1e+3,
                                                  as.numeric(sd(beta_diversity_3$time[[1]])) * 1e+3,
                                                  as.numeric(sd(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                  as.numeric(sd(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                  as.numeric(sd(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                  as.numeric(sd(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                  as.numeric(sd(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                  as.numeric(sd(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                  as.numeric(sd(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                  as.numeric(sd(prop_patches_result$time[[1]])) * 1e+3,
                                                  as.numeric(sd(CV_meta_result$time[[1]])) * 1e+3,
                                                  as.numeric(sd(hypervolume_det_result$time[[1]])) * 1e+3,
                                                  as.numeric(sd(hypervolume_dis_result$time[[1]])) * 1e+3),
                                     memory = c(as.numeric(beta_diversity_1$mem_alloc) / 1024^2,
                                                as.numeric(beta_diversity_2$mem_alloc) / 1024^2,
                                                as.numeric(beta_diversity_3$mem_alloc) / 1024^2,
                                                as.numeric(spatial_beta_div_1$mem_alloc) / 1024^2,
                                                as.numeric(spatial_beta_div_2$mem_alloc) / 1024^2,
                                                as.numeric(spatial_beta_div_3$mem_alloc) / 1024^2,
                                                as.numeric(temporal_beta_div_1$mem_alloc) / 1024^2,
                                                as.numeric(temporal_beta_div_2$mem_alloc) / 1024^2,
                                                as.numeric(temporal_beta_div_3$mem_alloc) / 1024^2,
                                                as.numeric(DNCI_multigroup_result$mem_alloc) / 1024^2,
                                                as.numeric(prop_patches_result$mem_alloc) / 1024^2,
                                                as.numeric(CV_meta_result$mem_alloc) / 1024^2,
                                                as.numeric(hypervolume_det_result$mem_alloc) / 1024^2,
                                                as.numeric(hypervolume_dis_result$mem_alloc) / 1024^2))



write.csv(benchmark_result_medium_df, "../benchmarks/result/benchmark_result_medium_df_r.csv")

#Small Dataset####
df <- small_df%>%
  mutate(Species = ifelse(Species=="NA", "N_A", Species))
## The abundance matrix
matrix_with_abundance <- df %>%
  filter(Sampling_date_order == 55) %>% 
  select(-Presence) %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude, normalized_temperature, normalized_precipitation)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

### The presence/absence matrix
matrix_with_presence <- df %>%
  filter(Sampling_date_order == 55) %>% 
  select(-Abundance) %>%
  pivot_wider(names_from = Species, values_from = Presence, values_fill=0) %>%
  select(-c(Year, Month, Day, Sampling_date_order, plot, Longitude, Latitude, normalized_temperature, normalized_precipitation)) %>% 
  select(which(colSums(.) !=0)) %>%
  filter(rowSums(.) != 0)

#### Benchmark the `beta.div.comp` function
beta_diversity_1 <- mark(beta.div.comp(matrix_with_abundance, coef = "J", quant = TRUE),
                         iterations = 100,
                         check = TRUE,
                         time_unit = "ms")

beta_diversity_2 <- mark(beta.div.comp(matrix_with_abundance, coef = "J", quant = FALSE),
                         iterations = 100,
                         check = TRUE,
                         time_unit = "ms")

beta_diversity_3 <- mark(beta.div.comp(matrix_with_presence, coef = "J", quant = FALSE),
                         iterations = 100,
                         check = TRUE,
                         time_unit = "ms")

####Benchmark the spatial beta diversity function
spatial_beta_div_1 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Abundance = sum(Abundance)) %>% 
                             spread(key = Species, value = Abundance, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")

spatial_beta_div_2 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Abundance = sum(Abundance)) %>% 
                             spread(key = Species, value = Abundance, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")

spatial_beta_div_3 <- mark(df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                             spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                           iterations = 100,
                           check = TRUE,
                           time_unit = "ms")
####Benchmark the temproal beta diversity function
temporal_beta_div_1 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              summarise(Abundance = sum(Abundance)) %>% 
                              spread(key = Species, value = Abundance, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

temporal_beta_div_2 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              summarise(Abundance = sum(Abundance)) %>% 
                              spread(key = Species, value = Abundance, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

temporal_beta_div_3 <- mark(df %>% 
                              group_by(Sampling_date_order,Species) %>%
                              dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                              spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                              ungroup() %>% 
                              select(-Sampling_date_order) %>% 
                              beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

#### Benchmark the DNCI function
DNCI_multigroup_result <- mark(DNCImper:::DNCI_multigroup(comm_small,
                                                          groups_small$Group, Nperm = 100,
                                                          symmetrize = FALSE,
                                                          plotSIMPER = FALSE,
                                                          parallelComputing = TRUE),
                               iterations = 100,
                               check = FALSE,
                               time_unit = "ms")
saveRDS(DNCI_multigroup_result, "../benchmarks/result/DNCI_small_result.rds")
DNCI_multigroup_result<-readRDS("../benchmarks/result/DNCI_small_result.rds")

#### Benchmark the prop_patches function
prop_patches_result <- mark(df %>% 
                              group_by(Species, plot) %>%
                              dplyr::summarise(mean_abundance = mean(Abundance)) %>% 
                              filter(mean_abundance  > 0) %>% 
                              dplyr::summarise(n_patches = n()) %>% 
                              mutate(prop_patches = n_patches/max(df$plot)) %>% 
                              summarise(mean_prop_patches = mean(prop_patches), min_prop_patches = min(prop_patches), max_prop_patches = max(prop_patches)),
                            iterations = 100,
                            check = TRUE,
                            time_unit = "ms")

CV_meta_result <- mark(CV_meta(df, df$Species, df$Sampling_date_order, df$plot, df$Abundance),
                       iterations = 100,
                       check = TRUE,
                       time_unit = "ms")

#### Benchmark the hypervolume functions
data_1 <- df%>%
  filter(Presence >0)%>%
  filter(Species =="BA")%>%
  select(normalized_temperature, normalized_precipitation)

data_2 <- df%>%
  filter(Presence >0)%>%
  filter(Species =="SH")%>%
  select(normalized_temperature, normalized_precipitation)

hypervolume_det_result <- mark(MVNH_det(data_1, var.names = c("Temperature", "Precipitation")), iterations = 100,
                               check = TRUE,
                               time_unit = "ms")

hypervolume_dis_result <- mark(MVNH_dissimilarity(data_1, data_2, var.names = c("Temperature", "Precipitation")), iterations = 100,
                               check = TRUE,
                               time_unit = "ms")

#Save all the results into a dataframe
test_case_list = c("beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                   "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                   "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                   "DNCI_multigroup_result",
                   "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result")

# Initialize empty data frame
all_time_small <- data.frame()

for (testcase in test_case_list) {
  # Get the object with that name
  testcase_obj <- get(testcase)
  
  # Get times from the object
  times <- as.numeric(testcase_obj$time[[1]]) * 1e+3
  
  # Create data frame for this test case
  data <- data.frame(TestCase = rep(testcase, length(times)),Time = times)
  
  # Append to the full data frame - this is the key fix
  all_time_small <- rbind(all_time_small, data)
}

write.csv(all_time_small, "../benchmarks/result/all_time_small_df_r.csv")

benchmark_result_small_df<-data.frame(TestCase = c("beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                                    "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                                    "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                                    "DNCI_multigroup_result",
                                                    "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result"),
                                       Time_minimum = c(as.numeric(min(beta_diversity_1$time[[1]])) * 1e+3,
                                                        as.numeric(min(beta_diversity_2$time[[1]])) * 1e+3,
                                                        as.numeric(min(beta_diversity_3$time[[1]])) * 1e+3,
                                                        as.numeric(min(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                        as.numeric(min(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                        as.numeric(min(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                        as.numeric(min(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                        as.numeric(min(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                        as.numeric(min(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                        as.numeric(min(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                        as.numeric(min(prop_patches_result$time[[1]])) * 1e+3,
                                                        as.numeric(min(CV_meta_result$time[[1]])) * 1e+3,
                                                        as.numeric(min(hypervolume_det_result$time[[1]])) * 1e+3,
                                                        as.numeric(min(hypervolume_dis_result$time[[1]])) * 1e+3),
                                       Time_median = c(as.numeric(median(beta_diversity_1$time[[1]])) * 1e+3,
                                                       as.numeric(median(beta_diversity_2$time[[1]])) * 1e+3,
                                                       as.numeric(median(beta_diversity_3$time[[1]])) * 1e+3,
                                                       as.numeric(median(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                       as.numeric(median(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                       as.numeric(median(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                       as.numeric(median(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                       as.numeric(median(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                       as.numeric(median(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                       as.numeric(median(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                       as.numeric(median(prop_patches_result$time[[1]])) * 1e+3,
                                                       as.numeric(median(CV_meta_result$time[[1]])) * 1e+3,
                                                       as.numeric(median(hypervolume_det_result$time[[1]])) * 1e+3,
                                                       as.numeric(median(hypervolume_dis_result$time[[1]])) * 1e+3),
                                       Time_mean = c(as.numeric(mean(beta_diversity_1$time[[1]])) * 1e+3,
                                                     as.numeric(mean(beta_diversity_2$time[[1]])) * 1e+3,
                                                     as.numeric(mean(beta_diversity_3$time[[1]])) * 1e+3,
                                                     as.numeric(mean(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                     as.numeric(mean(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                     as.numeric(mean(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                     as.numeric(mean(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                     as.numeric(mean(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                     as.numeric(mean(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                     as.numeric(mean(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                     as.numeric(mean(prop_patches_result$time[[1]])) * 1e+3,
                                                     as.numeric(mean(CV_meta_result$time[[1]])) * 1e+3,
                                                     as.numeric(mean(hypervolume_det_result$time[[1]])) * 1e+3,
                                                     as.numeric(mean(hypervolume_dis_result$time[[1]])) * 1e+3),
                                       Time_maximum = c(as.numeric(max(beta_diversity_1$time[[1]])) * 1e+3,
                                                        as.numeric(max(beta_diversity_2$time[[1]])) * 1e+3,
                                                        as.numeric(max(beta_diversity_3$time[[1]])) * 1e+3,
                                                        as.numeric(max(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                        as.numeric(max(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                        as.numeric(max(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                        as.numeric(max(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                        as.numeric(max(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                        as.numeric(max(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                        as.numeric(max(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                        as.numeric(max(prop_patches_result$time[[1]])) * 1e+3,
                                                        as.numeric(max(CV_meta_result$time[[1]])) * 1e+3,
                                                        as.numeric(max(hypervolume_det_result$time[[1]])) * 1e+3,
                                                        as.numeric(max(hypervolume_dis_result$time[[1]])) * 1e+3),
                                       Time_std = c(as.numeric(sd(beta_diversity_1$time[[1]])) * 1e+3,
                                                    as.numeric(sd(beta_diversity_2$time[[1]])) * 1e+3,
                                                    as.numeric(sd(beta_diversity_3$time[[1]])) * 1e+3,
                                                    as.numeric(sd(spatial_beta_div_1$time[[1]])) * 1e+3,
                                                    as.numeric(sd(spatial_beta_div_2$time[[1]])) * 1e+3,
                                                    as.numeric(sd(spatial_beta_div_3$time[[1]])) * 1e+3,
                                                    as.numeric(sd(temporal_beta_div_1$time[[1]])) * 1e+3,
                                                    as.numeric(sd(temporal_beta_div_2$time[[1]])) * 1e+3,
                                                    as.numeric(sd(temporal_beta_div_3$time[[1]])) * 1e+3,
                                                    as.numeric(sd(DNCI_multigroup_result$time[[1]])) * 1e+3,
                                                    as.numeric(sd(prop_patches_result$time[[1]])) * 1e+3,
                                                    as.numeric(sd(CV_meta_result$time[[1]])) * 1e+3,
                                                    as.numeric(sd(hypervolume_det_result$time[[1]])) * 1e+3,
                                                    as.numeric(sd(hypervolume_dis_result$time[[1]])) * 1e+3),
                                       memory = c(as.numeric(beta_diversity_1$mem_alloc) / 1024^2,
                                                  as.numeric(beta_diversity_2$mem_alloc) / 1024^2,
                                                  as.numeric(beta_diversity_3$mem_alloc) / 1024^2,
                                                  as.numeric(spatial_beta_div_1$mem_alloc) / 1024^2,
                                                  as.numeric(spatial_beta_div_2$mem_alloc) / 1024^2,
                                                  as.numeric(spatial_beta_div_3$mem_alloc) / 1024^2,
                                                  as.numeric(temporal_beta_div_1$mem_alloc) / 1024^2,
                                                  as.numeric(temporal_beta_div_2$mem_alloc) / 1024^2,
                                                  as.numeric(temporal_beta_div_3$mem_alloc) / 1024^2,
                                                  as.numeric(DNCI_multigroup_result$mem_alloc) / 1024^2,
                                                  as.numeric(prop_patches_result$mem_alloc) / 1024^2,
                                                  as.numeric(CV_meta_result$mem_alloc) / 1024^2,
                                                  as.numeric(hypervolume_det_result$mem_alloc) / 1024^2,
                                                  as.numeric(hypervolume_dis_result$mem_alloc) / 1024^2))


write.csv(benchmark_result_small_df, "../benchmarks/result/benchmark_result_small_df_r.csv")


####Additional benchmarking result
DNCI_multigroup_result_p <- mark(DNCImper:::DNCI_multigroup(comm_full,
                                                            groups_full$Group, 
                                                            Nperm = 100,
                                                            symmetrize = FALSE,
                                                            plotSIMPER = FALSE,                                              
                                                        parallelComputing = TRUE), iterations = 100,
                                                                 check = FALSE,
                                                               time_unit = "ms")

saveRDS(DNCI_multigroup_result_p, "../benchmarks/result/DNCI_full_result_with_parallelComputing.rds")





