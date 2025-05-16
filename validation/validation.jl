# validation/result_from_julia.jl
using Pkg
Pkg.activate("validation")

using MetaCommunityMetrics
using Random
using Pipe: @pipe
using CSV
using DataFrames
using RCall
using JLD2
using Printf
using Statistics

#Load the sample data
df = load_sample_data()

###Data Wrangling for running functions in this package and the equivalent function in R 
##Preparing the matrices for calculating beta diveristy
df = CSV.read("data/metacomm_rodent_df.csv", DataFrame)


#matrix with the species abundance
sample_matrix_abundance = @pipe df |> 
select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Presence)) |> # remove the unused four columns
filter(row -> row[:Sampling_date_order] == 50, _) |> # filter rows where Sampling_date_order == 50
unstack(_, :Species, :Abundance) |> #convert it back to the wide format 
select(_, Not(:Sampling_date_order, :plot)) |> # remove the unused columns
Matrix(_) |> #convert the dataframe to a matrix
_[:, sum(_, dims=1)[1, :] .!= 0] |> # Remove columns with sum of zero
_[sum(_, dims=2)[:, 1] .!= 0,:] # Remove rows with sum of zero


#matrix with the species presence-absence data
sample_matrix_presence = @pipe df |> 
select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Abundance)) |> # remove the unused four columns
filter(row -> row[:Sampling_date_order] == 50, _) |># filter rows where Sampling_date_order == 50
unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
select(_, Not(:Sampling_date_order, :plot)) |> # remove the unused columns
Matrix(_) |> #convert the dataframe to a matrix      
_[:, sum(_, dims=1)[1, :] .!= 0] |> # Remove columns with sum of zero
_[sum(_, dims=2)[:, 1] .!= 0,:] # Remove rows with sum of zero

#Data for the DNCI analysis
total_presence_df=@pipe df|>
filter(row -> row[:Sampling_date_order] == 50, _) |># filter rows where Sampling_date_order == 50
groupby(_,[:Species,:Sampling_date_order])|>
combine(_,:Presence=>sum=>:Total_Presence)

#Remove singletons 
total_richness_df= @pipe df|>
filter(row -> row[:Sampling_date_order] == 50, _) |># filter rows where Sampling_date_order == 50
leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
filter(row -> row[:Total_Presence] > 1, _)|> 
groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
combine(_,:Presence=>sum=>:Total_Richness)

# The groupings at t=50
cluster_result=create_clusters(total_richness_df.Sampling_date_order, 
total_richness_df.Latitude, 
total_richness_df.Longitude, 
total_richness_df.plot,
total_richness_df.Total_Richness) 

#The community presence matrix at t=50 with singletons removed
comm= @pipe df|>
filter(row -> row[:Sampling_date_order] == 50, _) |>
leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
filter(row -> row[:Total_Presence] > 1, _)|> 
select(_, [:plot, :Species, :Presence]) |>
unstack(_, :Species, :Presence, fill=0) |>
sort(_, :plot) |>
select(_, Not(:plot)) 

#The grouping result at t=50
group = @pipe cluster_result[50] |>
    select(_, [:Time, :Patch, :Group]) |>
    rename(_, :Patch => :plot) |>
    sort(_, :plot)


#Environmental Data for hypervolume calculations
sp1_data = @pipe df |> 
filter(row -> row[:Presence] > 0, _) |>
filter(row -> row[:Species] == "BA", _) |>
select(_, [:normalized_temperature, :normalized_precipitation])

sp2_data = @pipe df |> 
filter(row -> row[:Presence] > 0, _) |>
filter(row -> row[:Species] == "SH", _) |>
select(_, [:normalized_temperature, :normalized_precipitation])

## Make data available to R
R"df <- $df"  
R"sample_matrix_abundance <- $sample_matrix_abundance" 
R"sample_matrix_presence <- $sample_matrix_presence"
R"comm <- $comm"
R"groups <- as.factor($group[,'Group'])"
R"sp1_data <- $sp1_data"
R"sp2_data <- $sp2_data"

### Equivalent implementation in R 
## Install R package if not already installed
R"""
install.packages("data.table", repos="http://cran.us.r-project.org")
install.packages("tidyverse", repos="http://cran.us.r-project.org")
install.packages("adespatial", repos="http://cran.us.r-project.org")
install.packages("DNCImper", repos="http://cran.us.r-project.org")
install.packages("devtools", repos="http://cran.us.r-project.org")
library(devtools)
devtools::install_github("Corentin-Gibert-Paleontology/DNCImper")
install.packages("MVNH", repos="http://cran.us.r-project.org")
"""

## Load the R packages
R"""
library(data.table)  
library(tidyverse)
library(adespatial) 
library(DNCImper)
library(MVNH)
"""

## Calculate beta diversity using the function from the adespatial package in R
R"""
# beta.div.comp function from the adespatial package
beta_diversity_1 <- beta.div.comp(sample_matrix_abundance, coef = "J", quant = TRUE)
beta_diversity_2 <- beta.div.comp(sample_matrix_abundance, coef = "J", quant = FALSE)
beta_diversity_3 <- beta.div.comp(sample_matrix_presence, coef = "J", quant = FALSE)

# Calculate the mean spatial beta diversity
#using abundance data, with quant = TRUE
spatial_beta_div_1 <- df %>% 
                           group_by(plot,Species) %>%
                           dplyr::summarise(Abundance = sum(Abundance)) %>% 
                           spread(key = Species, value = Abundance, fill = 0) %>% 
                           ungroup() %>% 
                           dplyr::select(-plot) %>% 
                           beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE)

#using abundance data, with quant = FALSE                          
spatial_beta_div_2 <- df %>% 
                           group_by(plot,Species) %>%
                           dplyr::summarise(Abundance = sum(Abundance)) %>% 
                           spread(key = Species, value = Abundance, fill = 0) %>% 
                           ungroup() %>% 
                           dplyr::select(-plot) %>% 
                           beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE)

#using presence.absence data, with quant = FALSE
spatial_beta_div_3 <- df %>% 
                             group_by(plot,Species) %>%
                             dplyr::summarise(Sum_Presence = sum(Presence)) %>% 
                             spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                             ungroup() %>% 
                             dplyr::select(-plot) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE)

# Calculate the mean temporal beta diversity
#using abundance data, with quant = TRUE
temporal_beta_div_1 <- df %>% 
                             group_by(Sampling_date_order,Species) %>%
                             summarise(Abundance = sum(Abundance)) %>% 
                             spread(key = Species, value = Abundance, fill = 0) %>% 
                             ungroup() %>% 
                             select(-Sampling_date_order) %>% 
                             beta.div.comp(coef = "J", quant = TRUE, save.abc = FALSE)  

#using abundance data, with quant = FALSE
temporal_beta_div_2 <- df %>%
                             group_by(Sampling_date_order,Species) %>%
                             summarise(Abundance = sum(Abundance)) %>% 
                             spread(key = Species, value = Abundance, fill = 0) %>% 
                             ungroup() %>% 
                             select(-Sampling_date_order) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE)

#using presence.absence data, with quant = FALSE
temporal_beta_div_3 <- df %>%
                             group_by(Sampling_date_order,Species) %>%
                             summarise(Sum_Presence = sum(Presence)) %>% 
                             spread(key = Species, value = Sum_Presence, fill = 0) %>% 
                             ungroup() %>% 
                             select(-Sampling_date_order) %>% 
                             beta.div.comp(coef = "J", quant = FALSE, save.abc = FALSE)
"""


## The DNCI analysis
R"""
# Repeat DNCI_multigroup() multiple times
n_reps <- 100
dnci_values_df<-data.frame()
set.seed(123)  # For reproducibility
for (i in 1:n_reps) {
  result <- DNCImper:::DNCI_multigroup(comm, 
                                       groups, Nperm = 1000, 
                                       symmetrize = FALSE, 
                                       plotSIMPER = FALSE)
  # Save the DNCI values for all group pairs
  dnci_values_df<-rbind(dnci_values_df, result[,2:4] )
}
"""

## The Occupied Patches Proportion
R"""
prop_patches_result <- df %>% 
                             group_by(Species, plot) %>%
                             dplyr::summarise(mean_abundance = mean(Abundance)) %>% 
                             filter(mean_abundance  > 0) %>% 
                             dplyr::summarise(n_patches = n()) %>% 
                             mutate(prop_patches = n_patches/max(df$plot)) %>% 
                             summarise(mean_prop_patches = mean(prop_patches), min_prop_patches = min(prop_patches), max_prop_patches = max(prop_patches))

"""


R"""
## The Variability metric
# The R function “var.partition” written by Wang et al. (2019) to calculate variability and synchrony across hierarchical levels
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
  phi_S_L2R <- CV_S_R/CV_S_L
  phi_C_L2R <- CV_C_R/CV_C_L
  phi_S2C_L <- CV_C_L/CV_S_L
  phi_S2C_R <- CV_C_R/CV_S_R
  partition_3level <- c(CV_S_L=CV_S_L, CV_C_L=CV_C_L, CV_S_R=CV_S_R, CV_C_R=CV_C_R,
                        phi_S_L2R=phi_S_L2R, phi_C_L2R=phi_C_L2R, phi_S2C_L=phi_S2C_L,
                        phi_S2C_R=phi_S2C_R)
  return(partition_3level)
}

species_vals <- unique(df$Species)
date_vals <- unique(df$Sampling_date_order)
plot_vals <- unique(df$plot)

# Create the array with dimensions based on unique Species, Sampling_date_order, and plot
metacomm_tsdata <- array(0, dim = c(length(species_vals), length(date_vals), length(plot_vals)))

# Populate the array with values from the data frame
for(i in 1:nrow(df)){
  # Map Species, Sampling_date_order, and plot to their corresponding indices
  species_index <- which(species_vals == df$Species[i])
  date_index <- which(date_vals == df$Sampling_date_order[i])
  plot_index <- which(plot_vals == df$plot[i])
  
  # Assign the value from the data frame (Abundance) to the array
  metacomm_tsdata[species_index, date_index, plot_index] <- df$Abundance[i]
}

CV_result <- var.partition(metacomm_tsdata)

"""

## Hypervolume
R"""
# using functions from the R package MVNH (https://github.com/lvmuyang/MVNH)

MVNH_det_result <- MVNH_det(sp1_data, var.names = c("normalized_temperature", "normalized_precipitation"))

MVNH_dissimilarity <- MVNH_dissimilarity(sp1_data, sp2_data,  var.names = c("normalized_temperature", "normalized_precipitation"))
"""

## Get results back from R to Julia
#Beta Diversity
beta_diversity_1_r = R"beta_diversity_1$part"
beta_diversity_1_r = rcopy(beta_diversity_1_r)

beta_diversity_2_r = R"beta_diversity_2$part"
beta_diversity_2_r = rcopy(beta_diversity_2_r)

beta_diversity_3_r = R"beta_diversity_3$part"
beta_diversity_3_r = rcopy(beta_diversity_3_r)

spatial_beta_div_1_r = R"spatial_beta_div_1$part"
spatial_beta_div_1_r = rcopy(spatial_beta_div_1_r)

spatial_beta_div_2_r = R"spatial_beta_div_2$part"
spatial_beta_div_2_r = rcopy(spatial_beta_div_2_r)

spatial_beta_div_3_r = R"spatial_beta_div_3$part"
spatial_beta_div_3_r = rcopy(spatial_beta_div_3_r)

temporal_beta_div_1_r = R"temporal_beta_div_1$part"
temporal_beta_div_1_r = rcopy(temporal_beta_div_1_r)

temporal_beta_div_2_r = R"temporal_beta_div_2$part"
temporal_beta_div_2_r = rcopy(temporal_beta_div_2_r)

temporal_beta_div_3_r = R"temporal_beta_div_3$part"
temporal_beta_div_3_r = rcopy(temporal_beta_div_3_r)

#DNCI
dnci_values_r = R"dnci_values_df"
dnci_values_r = rcopy(dnci_values_r)
save_object("validation/validation_output/dnci_values_r.jld2", dnci_values_r)
load_object("validation/validation_output/dnci_values_r.jld2")

# Occupied Patches Proportion
prop_patches_result_r = R"prop_patches_result"
prop_patches_result_r = rcopy(prop_patches_result_r)

# Variability metric
CV_result_r = R"CV_result[1:4]"
CV_result_r = rcopy(CV_result_r)

#Hypervolume
MVNH_det_result_r = R"MVNH_det_result"
MVNH_det_result_r = rcopy(MVNH_det_result_r)

MVNH_dissimilarity_r = R"MVNH_dissimilarity"
MVNH_dissimilarity_r = rcopy(MVNH_dissimilarity_r)

###Same caculations in Julia
##Beta Diversity
beta_diversity_1_julia = beta_diversity(sample_matrix_abundance,quant = true)
beta_diversity_2_julia = beta_diversity(sample_matrix_abundance, quant = false)
beta_diversity_3_julia = beta_diversity(sample_matrix_presence, quant = false)

spatial_beta_div_1_julia = spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
spatial_beta_div_2_julia = spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
spatial_beta_div_3_julia = spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)

temporal_beta_div_1_julia = temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
temporal_beta_div_2_julia = temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
temporal_beta_div_3_julia = temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)

##DNCI
comm_mat=Matrix(comm)
n_rep=100
dnci_values_julia=DataFrame( group1= Int[], group2=Int[], DNCI=Float64[])
Random.seed!(123) 
for i in 1:n_rep
    result = DNCI_multigroup(comm_mat, group.Group, 1000; count = false)
    append!(dnci_values_julia, result[:,1:3])
end

save_object("validation/validation_output/dnci_values_julia.jld2", dnci_values_julia)
dnci_values_julia=load_object("validation/validation_output/dnci_values_julia.jld2")

##Occupied Patches Proportion
prop_patches_result_julia = prop_patches(df.Presence, df.Species, df.plot)

##Variability metric
CV_result_julia = CV_meta(df.Abundance, df.Sampling_date_order, df.plot, df.Species)

##Hypervolume
MVNH_det_result_julia = MVNH_det(sp1_data; var_names=["Temperature", "Precipitation"])
MVNH_dissimilarity_result_julia = MVNH_dissimilarity(sp1_data, sp2_data; var_names=["Temperature", "Precipitation"])


##Compare the results (absolute value only)
#Beta Diveristy
beta_diversity_summary_df = DataFrame(
    TestCase = ["Beta Diversity (Abundance, quant=true)", "Beta Diversity (Abundance, quant=false)", "Beta Diversity (Presence, quant=false)",
              "Spatial Beta Diversity (Abundance, quant=true)", "Spatial Beta Diversity (Abundance, quant=false)", "Spatial Beta Diversity (Presence, quant=false)",
              "Temporal Beta Diversity (Abundance, quant=true)", "Temporal Beta Diversity (Abundance, quant=false)", "Temporal Beta Diversity (Presence, quant=false)"],
    BDtotal = [abs.(beta_diversity_1_r[1] - beta_diversity_1_julia.BDtotal[1]), abs.(beta_diversity_2_r[1] - beta_diversity_2_julia.BDtotal[1]), abs.(beta_diversity_3_r[1] - beta_diversity_3_julia.BDtotal[1]),
              abs.(spatial_beta_div_1_r[1] - spatial_beta_div_1_julia.spatial_BDtotal[1]), abs.(spatial_beta_div_2_r[1] - spatial_beta_div_2_julia.spatial_BDtotal[1]), abs.(spatial_beta_div_3_r[1] - spatial_beta_div_3_julia.spatial_BDtotal[1]), 
              abs(temporal_beta_div_1_r[1] - temporal_beta_div_1_julia.temporal_BDtotal[1]), abs.(temporal_beta_div_2_r[1] - temporal_beta_div_2_julia.temporal_BDtotal[1]), abs.(temporal_beta_div_3_r[1] - temporal_beta_div_3_julia.temporal_BDtotal[1])],
    Repl = [abs.(beta_diversity_1_r[2] - beta_diversity_1_julia.Repl[1]), abs.(beta_diversity_2_r[2] - beta_diversity_2_julia.Repl[1]), abs.(beta_diversity_3_r[2] - beta_diversity_3_julia.Repl[1]),
              abs.(spatial_beta_div_1_r[2] - spatial_beta_div_1_julia.spatial_Repl[1]), abs.(spatial_beta_div_2_r[2] - spatial_beta_div_2_julia.spatial_Repl[1]), abs.(spatial_beta_div_3_r[2] - spatial_beta_div_3_julia.spatial_Repl[1]),
              abs.(temporal_beta_div_1_r[2] - temporal_beta_div_1_julia.temporal_Repl[1]), abs.(temporal_beta_div_2_r[2] - temporal_beta_div_2_julia.temporal_Repl[1]), abs.(temporal_beta_div_3_r[2] - temporal_beta_div_3_julia.temporal_Repl[1])],
    RichDif = [abs.(beta_diversity_1_r[3] - beta_diversity_1_julia.RichDif[1]), abs.(beta_diversity_2_r[3] - beta_diversity_2_julia.RichDif[1]), abs.(beta_diversity_3_r[3] - beta_diversity_3_julia.RichDif[1]),
              abs.(spatial_beta_div_1_r[3] - spatial_beta_div_1_julia.spatial_RichDif[1]), abs.(spatial_beta_div_2_r[3] - spatial_beta_div_2_julia.spatial_RichDif[1]), abs.(spatial_beta_div_3_r[3] - spatial_beta_div_3_julia.spatial_RichDif[1]),
              abs.(temporal_beta_div_1_r[3] - temporal_beta_div_1_julia.temporal_RichDif[1]), abs.(temporal_beta_div_2_r[3] - temporal_beta_div_2_julia.temporal_RichDif[1]), abs.(temporal_beta_div_3_r[3] - temporal_beta_div_3_julia.temporal_RichDif[1])],
              )

## DNCI
dnci_values_r.programming_language .= "R"
dnci_values_julia.programming_language .= "Julia"

DNCI_combined_df = vcat(dnci_values_r, dnci_values_julia)

dnci_values_r_1_2 = @pipe dnci_values_r |>
    filter(row -> row[:group1] == "1" && row[:group2] == "2", _)
dnci_values_r_1_3 = @pipe dnci_values_r |>
    filter(row -> row[:group1] == "1" && row[:group2] == "3", _)
dnci_values_r_2_3 = @pipe dnci_values_r |>
    filter(row -> row[:group1] == "2" && row[:group2] == "3", _)

dnci_values_julia_1_2 = @pipe dnci_values_julia |>
    filter(row -> row[:group1] == 1 && row[:group2] == 2, _)
dnci_values_julia_1_3 = @pipe dnci_values_julia |>
    filter(row -> row[:group1] == 1 && row[:group2] == 3, _)
dnci_values_julia_2_3 = @pipe dnci_values_julia |>
    filter(row -> row[:group1] == 2 && row[:group2] == 3, _)

dnci_summary_df = DataFrame(
    group1 = [1,1,2],
    group2 = [2,3,3],
    DNCI_minmum_R = [minimum(dnci_values_r_1_2[:,3]), minimum(dnci_values_r_1_3[:,3]), minimum(dnci_values_r_2_3[:,3])],
    DNCI_mean_R = [mean(dnci_values_r_1_2[:,3]), mean(dnci_values_r_1_3[:,3]), mean(dnci_values_r_2_3[:,3])],
    DNCI_maximum_R = [maximum(dnci_values_r_1_2[:,3]), maximum(dnci_values_r_1_3[:,3]), maximum(dnci_values_r_2_3[:,3])],
    DNCI_minimum_Julia = [minimum(dnci_values_julia_1_2[:,3]), minimum(dnci_values_julia_1_3[:,3]), minimum(dnci_values_julia_2_3[:,3])],
    DNCI_mean_Julia = [mean(dnci_values_julia_1_2[:,3]), mean(dnci_values_julia_1_3[:,3]), mean(dnci_values_julia_2_3[:,3])],
    DNCI_maximum_Julia = [maximum(dnci_values_julia_1_2[:,3]), maximum(dnci_values_julia_1_3[:,3]), maximum(dnci_values_julia_2_3[:,3])],
)

## Occupied Patches Proportion
prop_patches_summary_df = DataFrame(
    min_prop_patches = [abs.(prop_patches_result_r.min_prop_patches[1] .- prop_patches_result_julia.min_prop_patches[1])],
    mean_prop_patches = [abs.(prop_patches_result_r.mean_prop_patches[1] .- prop_patches_result_julia.mean_prop_patches[1])],
    max_prop_patches =  [abs.(prop_patches_result_r.max_prop_patches[1] .- prop_patches_result_julia.max_prop_patches[1])],
)

## Variability metric
CV_summary_df = DataFrame(
    CV_S_L = [abs.(CV_result_r[1] .- CV_result_julia.CV_s_l[1])],
    CV_C_L = [abs.(CV_result_r[2] .- CV_result_julia.CV_c_l[1])],
    CV_S_R = [abs.(CV_result_r[3] .- CV_result_julia.CV_s_r[1])],
    CV_C_R = [abs.(CV_result_r[4] .- CV_result_julia.CV_c_r[1])],
)

## Hypervolume
# Volume
MVNH_det_summary_df = DataFrame(
    total = [abs.(MVNH_det_result_r[1] .- MVNH_det_result_julia.total[1])],
    correlation = [abs.(MVNH_det_result_r[4] .- MVNH_det_result_julia.correlation[1])],
    Temperature = [abs.(MVNH_det_result_r[2] .- MVNH_det_result_julia.Temperature[1])],
    Precipitation = [abs.(MVNH_det_result_r[3] .- MVNH_det_result_julia.Precipitation[1])],
)
# Dissimilarity
MVNH_dissimilarity_summary_df = DataFrame(
    Metric = ["Bhattacharyya_distance", "Mahalanobis_distance", "Determinant_ratio"],
    total = [abs.(MVNH_dissimilarity_r[:Bhattacharyya_distance][1] .- MVNH_dissimilarity_result_julia.total[1]), 
             abs.(MVNH_dissimilarity_r[:Mahalanobis_distance][1] .- MVNH_dissimilarity_result_julia.total[2]), 
             abs.(MVNH_dissimilarity_r[:Determinant_ratio][1] .- MVNH_dissimilarity_result_julia.total[3])],
    correlation = [abs.(MVNH_dissimilarity_r[:Bhattacharyya_distance][4] .- MVNH_dissimilarity_result_julia.correlation[1]), 
    abs.(MVNH_dissimilarity_r[:Mahalanobis_distance][4] .- MVNH_dissimilarity_result_julia.correlation[2]), 
    abs.(MVNH_dissimilarity_r[:Determinant_ratio][4] .- MVNH_dissimilarity_result_julia.correlation[3])],
    Temperature = [abs.(MVNH_dissimilarity_r[:Bhattacharyya_distance][2] .- MVNH_dissimilarity_result_julia.Temperature[1]),
    abs.(MVNH_dissimilarity_r[:Mahalanobis_distance][2] .- MVNH_dissimilarity_result_julia.Temperature[2]), 
    abs.(MVNH_dissimilarity_r[:Determinant_ratio][2] .- MVNH_dissimilarity_result_julia.Temperature[3])],
    Precipitation = [abs.(MVNH_dissimilarity_r[:Bhattacharyya_distance][3] .- MVNH_dissimilarity_result_julia.Precipitation[1]),
    abs.(MVNH_dissimilarity_r[:Mahalanobis_distance][3] .- MVNH_dissimilarity_result_julia.Precipitation[2]),
    abs.(MVNH_dissimilarity_r[:Determinant_ratio][3] .- MVNH_dissimilarity_result_julia.Precipitation[3])],
)

# A function to format all numeric columns to 4 decimal places
function format_numeric_columns(df::DataFrame)
    for col in names(df)
        if eltype(df[!, col]) <: Number
            df[!, col] = [@sprintf("%.4f", val) for val in df[!, col]]
        end
    end
    return df
end

# Format the numeric columns in all dataframes
beta_diversity_summary_df = format_numeric_columns(beta_diversity_summary_df)
dnci_summary_df = format_numeric_columns(dnci_summary_df)
prop_patches_summary_df = format_numeric_columns(prop_patches_summary_df)
CV_summary_df = format_numeric_columns(CV_summary_df)
MVNH_det_summary_df = format_numeric_columns(MVNH_det_summary_df)
MVNH_dissimilarity_summary_df = format_numeric_columns(MVNH_dissimilarity_summary_df)




# Save the results to CSV 
CSV.write("validation/validation_output/beta_diversity_summary_df.csv", beta_diversity_summary_df)
CSV.write("validation/validation_output/dnci_combined_df.csv", DNCI_combined_df)
CSV.write("validation/validation_output/dnci_summary_df.csv", dnci_summary_df)
CSV.write("validation/validation_output/prop_patches_summary_df.csv", prop_patches_summary_df)
CSV.write("validation/validation_output/CV_summary_df.csv", CV_summary_df)
CSV.write("validation/validation_output/MVNH_det_summary_df.csv", MVNH_det_summary_df)
CSV.write("validation/validation_output/MVNH_dissimilarity_summary_df.csv", MVNH_dissimilarity_summary_df)








comm=load_object("/Users/yc2864/Downloads/comm.jld2")
groups=load_object("/Users/yc2864/Downloads/clustering_result_4.jld2")

groups=groups.Group
result = DNCI_multigroup(comm, groups.Group, 1000; count = false)
DNCI_result = Internal.DNCI_ses(paired_x, group_pair, Nperm; count=false)
comm=paired_x
groups=group_pair
Random.seed!(123)
Internal.PerSIMPER(comm, groups, Nperm; count=false)
AnaSimp = Internal.simper(comm, groups)



function svg_plot_clusters(latitude::Vector{Float64}, longitude::Vector{Float64}, group::Union{AbstractVector, String}, output_file="clusters.svg")
    # Get unique cluster IDs and assign numeric identifiers
    unique_clusters = unique(group)
    cluster_map = Dict(cluster => i for (i, cluster) in enumerate(unique_clusters))
    numeric_ids = [cluster_map[cluster] for cluster in group]
    
    # Define a set of distinctive colors
    base_colors = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
        "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5", 
        "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f"
    ]
    
    # Generate colors for each cluster
    colors = if length(unique_clusters) <= length(base_colors)
        base_colors[1:length(unique_clusters)]
    else
        # Generate additional colors if needed
        c = copy(base_colors)
        while length(c) < length(unique_clusters)
            push!(c, "#" * join(rand('0':'9', 6)))
        end
        c
    end
    
    # Calculate bounds with padding
    lon_padding = 0.05 * (maximum(longitude) - minimum(longitude))
    lat_padding = 0.05 * (maximum(latitude) - minimum(latitude))
    
    min_lon = minimum(longitude) - lon_padding
    max_lon = maximum(longitude) + lon_padding
    min_lat = minimum(latitude) - lat_padding
    max_lat = maximum(latitude) + lat_padding
    
    # SVG dimensions
    width = 800
    height = 600
    margin = 100
    
    # Map coordinates to SVG space
    function map_coords(lon, lat)
        x = margin + (lon - min_lon) / (max_lon - min_lon) * (width - 2 * margin)
        # Flip y-axis (SVG has origin at top-left)
        y = height - margin - (lat - min_lat) / (max_lat - min_lat) * (height - 2 * margin)
        return round(x, digits=1), round(y, digits=1)
    end
    
    # Start building SVG content
    svg = """<?xml version="1.0" encoding="UTF-8"?>
    <svg width="$(width)" height="$(height)" xmlns="http://www.w3.org/2000/svg">
    <rect width="100%" height="100%" fill="white"/>
    
    <!-- Title -->
    <text x="$(width/2)" y="25" font-family="Arial" font-size="18" text-anchor="middle" font-weight="bold">Cluster Visualization</text>
    
    <!-- Axes -->
    <line x1="$(margin)" y1="$(height-margin)" x2="$(width-margin)" y2="$(height-margin)" stroke="black" stroke-width="1.5"/>
    <line x1="$(margin)" y1="$(margin)" x2="$(margin)" y2="$(height-margin)" stroke="black" stroke-width="1.5"/>
    
    <!-- Axis labels -->
    <text x="$(width/2)" y="$(height-10)" font-family="Arial" font-size="14" text-anchor="middle">Longitude</text>
    <text x="15" y="$(height/2)" font-family="Arial" font-size="14" text-anchor="middle" transform="rotate(-90, 15, $(height/2))">Latitude</text>
    
    <!-- Data points -->
    """
    
    # Add all data points
    for i in 1:length(latitude)
        x, y = map_coords(longitude[i], latitude[i])
        cluster_id = numeric_ids[i]
        color = colors[cluster_id]
        
        svg *= """
        <circle cx="$(x)" cy="$(y)" r="5" fill="$(color)" stroke="black" stroke-width="0.5" opacity="0.8"/>
        """
    end
    
    # Add legend
    legend_x = width - margin + 15
    legend_y = margin + 20
    
    svg *= """
    <!-- Legend -->
    <text x="$(legend_x)" y="$(legend_y - 15)" font-family="Arial" font-size="12" font-weight="bold">Clusters</text>
    """
    
    for (i, cluster) in enumerate(unique_clusters)
        y_pos = legend_y + (i-1) * 20
        svg *= """
        <rect x="$(legend_x)" y="$(y_pos-10)" width="10" height="10" fill="$(colors[i])" stroke="black" stroke-width="0.5"/>
        <text x="$(legend_x + 20)" y="$(y_pos)" font-family="Arial" font-size="12">$(cluster)</text>
        """
    end
    
    # Close SVG
    svg *= "</svg>"
    
    # Write to file
    open(output_file, "w") do f
        write(f, svg)
    end
    
    println("Cluster visualization saved to $(abspath(output_file))")
    println("Open this file in any web browser or image viewer to see the visualization")
    
    return output_file
end
svg_plot_clusters(groups.Latitude, groups.Longitude, groups.Group)


DNCI_result = DNCI_ses(paired_x, group_pair, Nperm; count)
