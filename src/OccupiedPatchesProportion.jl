using Pkg
Pkg.activate(".")

using Distributions
using Plots
using LinearAlgebra
using GaussianRandomFields
using Distances
using DataFrames
using CSV
using StatsBase
using DataStructures
using ProgressMeter
using Pipe: @pipe
using Combinatorics

###Proportion of patches occupied across species
function prop_patches(dynamic_df)
    species_patches_df= 
        @pipe dynamic_df[:,[:Presence, :Species, :Patch]]|>#select column N, Species, Patch
        groupby(_, [:Patch, :Species])|>
    combine(_, :Presence => sum => :Total_N)|>
    transform(_, :Total_N=> ByRow(x->ifelse(x> 0, 1.0, 0.0))=> :Presence)|>
    groupby(_, [:Species])


    patches_occupied_matrix=hcat(unique(dynamic_df.Species), zeros(Float64, size(unique(dynamic_df.Species), 1))) 
    for group in 1:size(patches_occupied_matrix,1)
        patches_occupied_matrix[group,2]= sum(species_patches_df[group].Presence)/nrow(species_patches_df[group])
    end

    prop_patches_df=DataFrames.DataFrame(mean_prop_patches = mean(patches_occupied_matrix[:,2]),
        min_prop_patches = minimum(patches_occupied_matrix[:,2]),
        max_prop_patches = maximum(patches_occupied_matrix[:,2]))   
end
