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

###Data Wrangling for running functions in this package and the equivalent functions in R 

##Preparing the data for calculating beta diveristy
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
clustering_result = create_groups(df.Sampling_date_order, 
                                            df.Latitude, 
                                            df.Longitude,                                      
                                            df.plot, 
                                            df.Species, 
                                            df.Presence)

# The groupings at t=60
group_df = @pipe df |>
                filter(row -> row[:Sampling_date_order] == 60, _) |>
                select(_, [:plot, :Species, :Presence]) |>
                innerjoin(_, clustering_result[60], on = [:plot => :Site, :Species], makeunique = true)|>
                select(_, [:plot, :Species, :Presence, :Group]) |>
                unstack(_, :Species, :Presence, fill=0)

#The community presence matrix at t=60
comm= @pipe group_df |>
                  select(_, Not([:plot,:Group]))



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
R"groups <- $group_df[,'Group']"
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


## The DNCI analysis in R
R"""
# Repeat DNCI_multigroup() multiple times
n_reps <- 100
dnci_values_df<-data.frame()
set.seed(123)  # For reproducibility
for (i in 1:n_reps) {
  result <- DNCImper:::DNCI_multigroup(comm, 
                                       groups, Nperm = 1000, 
                                       symmetrize = FALSE, 
                                       plotSIMPER = FALSE,
                                       parallelComputing = TRUE)
  # Save the DNCI values for all group pairs
  dnci_values_df<-rbind(dnci_values_df, result[,2:6] )
}
"""

## The Occupied Patches Proportion in R
R"""
prop_patches_result <- df %>% 
                             group_by(Species, plot) %>%
                             dplyr::summarise(mean_abundance = mean(Abundance)) %>% 
                             filter(mean_abundance  > 0) %>% 
                             dplyr::summarise(n_patches = n()) %>% 
                             mutate(prop_patches = n_patches/max(df$plot)) %>% 
                             summarise(mean_prop_patches = mean(prop_patches), min_prop_patches = min(prop_patches), max_prop_patches = max(prop_patches))

"""

## The Variability Metric in R
# The R function “var.partition” written by Wang et al. (2019) to calculate variability and synchrony across hierarchical levels
R"""
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

## Hypervolume Measurements in R
# using functions from the R package MVNH (https://github.com/lvmuyang/MVNH)
R"""
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
save_object("validation/output/dnci_values_r.jld2", dnci_values_r)
dnci_values_r = load_object("validation/output/dnci_values_r.jld2")

#Occupied Patches Proportion
prop_patches_result_r = R"prop_patches_result"
prop_patches_result_r = rcopy(prop_patches_result_r)

#Variability metric
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
dnci_values_julia=DataFrame()
Random.seed!(123) 
for i in 1:n_rep
    result = DNCI_multigroup(comm_mat, group_df.Group, 1000; Nperm_count = false)
    append!(dnci_values_julia, result)
end

save_object("validation/output/dnci_values_julia.jld2", dnci_values_julia)
dnci_values_julia=load_object("validation/output/dnci_values_julia.jld2")

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

DNCI_combined_df = vcat(dnci_values_r, dnci_values_julia[:,(vcat(1:5, 7))])

dnci_values_r_1_2 = @pipe dnci_values_r |>
    filter(row -> row[:group1] == 1 && row[:group2] == 2, _) |> 
    filter(row -> !isnan(row[:DNCI]) && !isinf(row[:DNCI]), _) #remove all NaN and Inf values
dnci_values_r_1_3 = @pipe dnci_values_r |>
    filter(row -> row[:group1] == 1 && row[:group2] == 3, _) |> 
    filter(row -> !isnan(row[:DNCI]) && !isinf(row[:DNCI]), _) #remove all NaN and Inf values
dnci_values_r_2_3 = @pipe dnci_values_r |>
    filter(row -> row[:group1] == 2 && row[:group2] == 3, _) |> 
    filter(row -> !isnan(row[:DNCI]) && !isinf(row[:DNCI]), _) #remove all NaN and Inf values

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
CSV.write("validation/output/beta_diversity_summary_df.csv", beta_diversity_summary_df)
CSV.write("validation/output/dnci_combined_df.csv", DNCI_combined_df)
CSV.write("validation/output/dnci_summary_df.csv", dnci_summary_df)
CSV.write("validation/output/prop_patches_summary_df.csv", prop_patches_summary_df)
CSV.write("validation/output/CV_summary_df.csv", CV_summary_df)
CSV.write("validation/output/MVNH_det_summary_df.csv", MVNH_det_summary_df)
CSV.write("validation/output/MVNH_dissimilarity_summary_df.csv", MVNH_dissimilarity_summary_df)
