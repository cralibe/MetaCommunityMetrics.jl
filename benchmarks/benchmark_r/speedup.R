# benchmarks/result.R
#This script is used to calculate the speedup and combine the results from julia and R.
#Load necessary libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(boot)
#Read in the result
#all time samples
all_time_full_julia<-read.csv("result/all_time_full_julia.csv")
all_time_full_r<-read.csv("result/all_time_full_df_r.csv")

all_time_medium_julia<-read.csv("result/all_time_medium_julia.csv")
all_time_medium_r<-read.csv("result/all_time_medium_df_r.csv")

all_time_small_julia<-read.csv("result/all_time_small_julia.csv")
all_time_small_r<-read.csv("result/all_time_small_df_r.csv")

#min,median,mean,max time and memory data
full_result_julia<-read.csv("result/benchmark_result_full_df_julia.csv")
full_result_r<-read.csv("result/benchmark_result_full_df_r.csv")

medium_result_julia<-read.csv("result/benchmark_result_medium_df_julia.csv")
medium_result_r<-read.csv("result/benchmark_result_medium_df_r.csv")

small_result_julia<-read.csv("result/benchmark_result_small_df_julia.csv")
small_result_r<-read.csv("result/benchmark_result_small_df_r.csv")

#Organize the data
full_result<-full_result_julia%>%
  select(TestCase, Time_median)%>%
  rename(Time_median_julia=Time_median)%>%
  left_join(full_result_r)%>%
  select(TestCase, Time_median_julia, Time_median)%>%
  rename(Time_median_r=Time_median)%>%
  filter(!is.na(Time_median_r))%>%
  mutate(Speedup=Time_median_r/Time_median_julia)

medium_result<-medium_result_julia%>%
  select(TestCase, Time_median)%>%
  rename(Time_median_julia=Time_median)%>%
  left_join(medium_result_r)%>%
  select(TestCase, Time_median_julia, Time_median)%>%
  rename(Time_median_r=Time_median)%>%
  filter(!is.na(Time_median_r))%>%
  mutate(Speedup=Time_median_r/Time_median_julia)


small_result<-small_result_julia%>%
  select(TestCase, Time_median)%>%
  rename(Time_median_julia=Time_median)%>%
  left_join(small_result_r)%>%
  select(TestCase, Time_median_julia, Time_median)%>%
  rename(Time_median_r=Time_median)%>%
  filter(!is.na(Time_median_r))%>%
  mutate(Speedup=Time_median_r/Time_median_julia)

write.csv(small_result, "result/small_df_speedup.csv")
write.csv(medium_result, "result/medium_df_speedup.csv")
write.csv(full_result, "result/full_df_speedup.csv")


##Bootstrap 95% confidence interval for the median to deal with outlier in the higher end
combined_data<-rbind(all_time_full_julia%>%mutate(DataSize="Large", Language="Julia"),
                     all_time_medium_julia%>%mutate(DataSize="Medium", Language="Julia"),
                    all_time_small_julia%>%mutate(DataSize="Small", Language="Julia"),
                     all_time_full_r%>%select(-X)%>%mutate(DataSize="Large", Language="R"),
                     all_time_medium_r%>%select(-X)%>%mutate(DataSize="Medium", Language="R"),
                     all_time_small_r%>%select(-X)%>%mutate(DataSize="Small", Language="R"))
                     
# Get all unique combinations of DataSize, Language,TestCase variables
combinations <- combined_data %>%
  distinct(TestCase, DataSize)

# Initialize empty data frame
median_ci_all <- data.frame()

for (combo in 1:nrow(combinations)) {
  # Get the object with that name
  current_combo<-combinations[combo,]
  
  julia_data<- combined_data %>% filter(TestCase == current_combo$TestCase, 
                                  DataSize == current_combo$DataSize, 
                                  Language == "Julia")
  
  r_data<- combined_data %>% filter(TestCase == current_combo$TestCase, 
                                        DataSize == current_combo$DataSize, 
                                        Language == "R")

  # Calculate the ratio of medians
  median_julia <- median(julia_data$Time)
  median_r <- median(r_data$Time)
  speedup_ratio <- median_r / median_julia
  
  # Bootstrap CI for ratio of medians
  n_boots <- 2000
  bootstrap_ratios <- numeric(n_boots)
  
  for (i in 1:n_boots) {
    # Resample independently
    julia_sample <- sample(julia_data$Time, length(julia_data$Time), replace = TRUE)
    r_sample <- sample(r_data$Time, length(r_data$Time), replace = TRUE)
    
    # Calculate ratio of medians for this bootstrap sample
    bootstrap_ratios[i] <- median(r_sample) / median(julia_sample)
  }
  
  # Calculate 95% bootstrap CI
  ci_lower <- quantile(bootstrap_ratios, 0.025)
  ci_upper <- quantile(bootstrap_ratios, 0.975)  
  
  result<-data.frame(TestCase = current_combo$TestCase, 
                     DataSize = current_combo$DataSize, 
                     julia_median=median_julia,
                     r_median=median_r,
                     Speedup_median = speedup_ratio,
                     Lower_CI = ci_lower,
                     Upper_CI = ci_upper)
  
  median_ci_all<-rbind(median_ci_all, result)
}

test_case_order <- c(
  "Beta Diversity (Abundance, quant=true)",
  "Beta Diversity (Abundance, quant=false)",
  "Beta Diversity (Presence, quant=false)",
  "Spatial Beta Diversity (Abundance, quant=true)",
  "Spatial Beta Diversity (Abundance, quant=false)",
  "Spatial Beta Diversity (Presence, quant=false)",
  "Temporal Beta Diversity (Abundance, quant=true)",
  "Temporal Beta Diversity (Abundance, quant=false)",
  "Temporal Beta Diversity (Presence, quant=false)",
  "Dispersal-niche continuum index",
  "Occupied Patches Proportion",
  "Variability Metrics",
  "Hypervolume Estimation",
  "Hypervolume Dissimilarity"
)
median_ci_all_save <- median_ci_all %>%
  mutate(TestCase = case_when(
    TestCase == "beta_diversity_1" ~ "Beta Diversity (Abundance, quant=true)",
    TestCase == "beta_diversity_2" ~ "Beta Diversity (Abundance, quant=false)",
    TestCase == "beta_diversity_3" ~ "Beta Diversity (Presence, quant=false)",
    TestCase == "spatial_beta_div_1" ~ "Spatial Beta Diversity (Abundance, quant=true)",
    TestCase == "spatial_beta_div_2" ~ "Spatial Beta Diversity (Abundance, quant=false)",
    TestCase == "spatial_beta_div_3" ~ "Spatial Beta Diversity (Presence, quant=false)",
    TestCase == "temporal_beta_div_1" ~ "Temporal Beta Diversity (Abundance, quant=true)",
    TestCase == "temporal_beta_div_2" ~ "Temporal Beta Diversity (Abundance, quant=false)",
    TestCase == "temporal_beta_div_3" ~ "Temporal Beta Diversity (Presence, quant=false)",
    TestCase == "DNCI_multigroup_result" ~ "Dispersal-niche continuum index",
    TestCase == "prop_patches_result" ~ "Occupied Patches Proportion",
    TestCase == "CV_meta_result" ~ "Variability Metrics",
    TestCase == "hypervolume_det_result" ~ "Hypervolume Estimation",
    TestCase == "hypervolume_dis_result" ~ "Hypervolume Dissimilarity",
    TRUE ~ TestCase  # Keep original value if no match
  ))%>%
  mutate(TestCase = factor(TestCase, levels = test_case_order)) %>%
  arrange(TestCase)

write.csv(median_ci_all_save, "result/median_ci_all.csv")
##plot
# First, extract the base test name (before the underscore and number)
median_ci_df <- median_ci_all %>%
  mutate(
    # Extract the base test name (everything before the last underscore or before _#)
    TestGroup = case_when(
      grepl("beta_diversity", TestCase) ~ "beta_diversity",
      grepl("spatial_beta_div", TestCase) ~ "spatial_beta_div",
      grepl("temporal_beta_div", TestCase) ~ "temporal_beta_div",
      TRUE ~ TestCase
    ),
    # Extract the number/variant (for assigning shapes)
    TestVariant = case_when(
      grepl("_[123]$", TestCase) ~ as.numeric(sub(".*_([123])$", "\\1", TestCase)),
      TRUE ~ 1
    ),
    # Flag to identify beta diversity-related test cases
    IsBetaDiv = grepl("beta_div", TestCase)
  )

# Create a custom color palette for the test groups
test_group_colors <- c(
  "Beta Diversity" = "purple",
  "Spatial Beta Diversity" = "#56B4E9",
  "Temporal Beta Diversity" = "#009E73",
  "Dispersal-niche continuum index" = "pink",
  "Occupied Patches Proportion" = "orange",
  "Variability Metrics" = "red",
  "Hypervolume Estimation" = "brown",
  "Hypervolume Dissimilarity" = "black"
)
# First, create a mapping of TestCase to TestGroup
test_case_to_group <- c(
  # Beta Diversity
  "beta_diversity_1" = "Beta Diversity",
  "beta_diversity_2" = "Beta Diversity",
  "beta_diversity_3" = "Beta Diversity",
  
  # Spatial Beta Diversity
  "spatial_beta_div_1" = "Spatial Beta Diversity",
  "spatial_beta_div_2" = "Spatial Beta Diversity",
  "spatial_beta_div_3" = "Spatial Beta Diversity",
  
  # Temporal Beta Diversity
  "temporal_beta_div_1" = "Temporal Beta Diversity",
  "temporal_beta_div_2" = "Temporal Beta Diversity",
  "temporal_beta_div_3" = "Temporal Beta Diversity",
  
  # Dispersal-niche continuum index
  "DNCI_multigroup_result" = "Dispersal-niche continuum index",
  
  # Occupied Patches Proportion
  "prop_patches_result" = "Occupied Patches Proportion",
  
  # Variability Metrics
  "CV_meta_result" = "Variability Metrics",
  
  # Hypervolume Estimation
  "hypervolume_det_result" = "Hypervolume Estimation",
  
  # Hypervolume Dissimilarity
  "hypervolume_dis_result" = "Hypervolume Dissimilarity"
)



# Create TestCase order based on the group order
test_case_order <- c(
  # Beta Diversity
  "beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
  
  # Spatial Beta Diversity
  "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
  
  # Temporal Beta Diversity
  "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
  
  # Dispersal-niche continuum index
  "DNCI_multigroup_result",
  
  # Occupied Patches Proportion
  "prop_patches_result",
  
  # Variability Metrics
  "CV_meta_result",
  
  # Hypervolume Estimation
  "hypervolume_det_result",
  
  # Hypervolume Dissimilarity
  "hypervolume_dis_result"
)

# Convert TestGroup and TestCase to factors with the specified order
median_ci_df  <- median_ci_df  %>%
  # Ensure TestGroup is assigned correctly based on TestCase
  mutate(
    TestGroup = test_case_to_group[TestCase],
    # Then convert both to factors with the proper ordering
    TestGroup = factor(TestGroup, levels = test_group_order),
    TestCase = factor(TestCase, levels = test_case_order)
  )
# Now create the plot with the ordered data
p <- ggplot(median_ci_df , 
            aes(x = DataSize, y = Speedup_median, 
                group = TestCase, 
                color = TestGroup)) +  
  # Add horizontal line at y=1 first
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +
  # All lines are now solid
  geom_line(size = 0.7) +
  # Only add points for beta diversity test cases with different shapes
  geom_point(data = median_ci_df  %>% filter(IsBetaDiv), 
             aes(shape = factor(TestVariant)),
             size = 4) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), 
                width = 0.2,  
                alpha = 0.7) +
  scale_color_manual(
    values = test_group_colors, 
    name = "Test Group"
  ) +
  scale_shape_manual(
    values = c(16, 17, 15), 
    name = "Variant", 
    label = c(
      "Abundance, quant=true", 
      "Abundance, quant=false", 
      "Presence, quant=false"
    )
  ) +
  theme_cowplot() +
  scale_x_discrete(limits = c("Small", "Medium", "Large")) +
  labs(y = "Speedup", x = "Data Size")

#save te plot
ggsave("result/speedup.pdf", dpi=300, width = 10, height = 5, bg="white")
#Memory
full_result_memory<-full_result_julia%>%
  select(TestCase, memory)%>%
  rename(memory_julia=memory)%>%
  left_join(full_result_r)%>%
  select(TestCase, memory_julia, memory)%>%
  rename(memory_r=memory)%>%
  filter(!is.na(memory_r))%>%
  mutate(TestCase = case_when(
    TestCase == "beta_diversity_1" ~ "Beta Diversity (Abundance, quant=true)",
    TestCase == "beta_diversity_2" ~ "Beta Diversity (Abundance, quant=false)",
    TestCase == "beta_diversity_3" ~ "Beta Diversity (Presence, quant=false)",
    TestCase == "spatial_beta_div_1" ~ "Spatial Beta Diversity (Abundance, quant=true)",
    TestCase == "spatial_beta_div_2" ~ "Spatial Beta Diversity (Abundance, quant=false)",
    TestCase == "spatial_beta_div_3" ~ "Spatial Beta Diversity (Presence, quant=false)",
    TestCase == "temporal_beta_div_1" ~ "Temporal Beta Diversity (Abundance, quant=true)",
    TestCase == "temporal_beta_div_2" ~ "Temporal Beta Diversity (Abundance, quant=false)",
    TestCase == "temporal_beta_div_3" ~ "Temporal Beta Diversity (Presence, quant=false)",
    TestCase == "DNCI_multigroup_result" ~ "Dispersal-niche continuum index",
    TestCase == "prop_patches_result" ~ "Occupied Patches Proportion",
    TestCase == "CV_meta_result" ~ "Variability Metrics",
    TestCase == "hypervolume_det_result" ~ "Hypervolume Estimation",
    TestCase == "hypervolume_dis_result" ~ "Hypervolume Dissimilarity",
    TRUE ~ TestCase  # Keep original value if no match
  ))%>%
  mutate(TestCase = factor(TestCase, levels = test_case_order)) %>%
  arrange(TestCase)


medium_result_memory<-medium_result_julia%>%
  select(TestCase, memory)%>%
  rename(memory_julia=memory)%>%
  left_join(medium_result_r)%>%
  select(TestCase, memory_julia, memory)%>%
  rename(memory_r=memory)%>%
  filter(!is.na(memory_r))%>%
  mutate(TestCase = case_when(
    TestCase == "beta_diversity_1" ~ "Beta Diversity (Abundance, quant=true)",
    TestCase == "beta_diversity_2" ~ "Beta Diversity (Abundance, quant=false)",
    TestCase == "beta_diversity_3" ~ "Beta Diversity (Presence, quant=false)",
    TestCase == "spatial_beta_div_1" ~ "Spatial Beta Diversity (Abundance, quant=true)",
    TestCase == "spatial_beta_div_2" ~ "Spatial Beta Diversity (Abundance, quant=false)",
    TestCase == "spatial_beta_div_3" ~ "Spatial Beta Diversity (Presence, quant=false)",
    TestCase == "temporal_beta_div_1" ~ "Temporal Beta Diversity (Abundance, quant=true)",
    TestCase == "temporal_beta_div_2" ~ "Temporal Beta Diversity (Abundance, quant=false)",
    TestCase == "temporal_beta_div_3" ~ "Temporal Beta Diversity (Presence, quant=false)",
    TestCase == "DNCI_multigroup_result" ~ "Dispersal-niche continuum index",
    TestCase == "prop_patches_result" ~ "Occupied Patches Proportion",
    TestCase == "CV_meta_result" ~ "Variability Metrics",
    TestCase == "hypervolume_det_result" ~ "Hypervolume Estimation",
    TestCase == "hypervolume_dis_result" ~ "Hypervolume Dissimilarity",
    TRUE ~ TestCase  # Keep original value if no match
  ))%>%
  mutate(TestCase = factor(TestCase, levels = test_case_order)) %>%
  arrange(TestCase)


small_result_memory<-small_result_julia%>%
  select(TestCase, memory)%>%
  rename(memory_julia=memory)%>%
  left_join(small_result_r)%>%
  select(TestCase, memory_julia, memory)%>%
  rename(memory_r=memory)%>%
  filter(!is.na(memory_r))%>%
  mutate(TestCase = case_when(
    TestCase == "beta_diversity_1" ~ "Beta Diversity (Abundance, quant=true)",
    TestCase == "beta_diversity_2" ~ "Beta Diversity (Abundance, quant=false)",
    TestCase == "beta_diversity_3" ~ "Beta Diversity (Presence, quant=false)",
    TestCase == "spatial_beta_div_1" ~ "Spatial Beta Diversity (Abundance, quant=true)",
    TestCase == "spatial_beta_div_2" ~ "Spatial Beta Diversity (Abundance, quant=false)",
    TestCase == "spatial_beta_div_3" ~ "Spatial Beta Diversity (Presence, quant=false)",
    TestCase == "temporal_beta_div_1" ~ "Temporal Beta Diversity (Abundance, quant=true)",
    TestCase == "temporal_beta_div_2" ~ "Temporal Beta Diversity (Abundance, quant=false)",
    TestCase == "temporal_beta_div_3" ~ "Temporal Beta Diversity (Presence, quant=false)",
    TestCase == "DNCI_multigroup_result" ~ "Dispersal-niche continuum index",
    TestCase == "prop_patches_result" ~ "Occupied Patches Proportion",
    TestCase == "CV_meta_result" ~ "Variability Metrics",
    TestCase == "hypervolume_det_result" ~ "Hypervolume Estimation",
    TestCase == "hypervolume_dis_result" ~ "Hypervolume Dissimilarity",
    TRUE ~ TestCase  # Keep original value if no match
  ))%>%
  mutate(TestCase = factor(TestCase, levels = test_case_order)) %>%
  arrange(TestCase)

write.csv(small_result_memory, "result/small_df_memory.csv")
write.csv(medium_result_memory, "result/medium_df_memory.csv")
write.csv(full_result_memory, "result/full_df_memory.csv")


small_result_memory%>%
  mutate(diff=memory_julia-memory_r)%>%
  arrange(diff)

medium_result_memory%>%
  mutate(diff=memory_julia-memory_r)%>%
  arrange(diff)

full_result_memory%>%
  mutate(diff=memory_julia-memory_r)%>%
  arrange(diff)


##Speedup for the DNCI with parallelComputing in R using the full dataset
DNCI_multigroup_result<-read.csv("../benchmarks/result/benchmark_result_full_df_julia.csv")%>%
  filter(TestCase=="DNCI_multigroup_result")
DNCI_multigroup_result_p<-readRDS("../benchmarks/result/DNCI_full_result_with_parallelComputing.rds")


result_df<-data.frame(TestCase = "DNCI_multigroup_result",
                      Speedup = as.numeric(mean(DNCI_multigroup_result_p$time[[1]])) * 1e+3/
                        DNCI_multigroup_result$Time_median,
                      julia_memory = DNCI_multigroup_result$memory,
                      r_memory = as.numeric(DNCI_multigroup_result_p$mem_alloc) / 1024^2)
                      
write.csv(result_df, "result/benchmark_DNCI_with_parallelComputing_full_df.csv")



