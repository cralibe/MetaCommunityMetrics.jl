# src/MetaCommunityMetrics.jl

module MetaCommunityMetrics

using DataFrames
using CSV
using Pipe
using Statistics
using Distances
using Clustering
using Distributions
using LinearAlgebra
using GaussianRandomFields
using StatsBase
using DataStructures
using ProgressMeter
using Combinatorics
using Random


# Include the internal utilities
include("Internal.jl")
using .Internal  


# Public function
include("BetaDiversity.jl")
include("DNCI.jl")
include("NicheOverlapIndex.jl")
include("OccupiedPatchesProportion.jl")
include("VariabilityMetrics.jl")
include("Hypervolume.jl")

# Function to load sample data
function load_sample_data()
    path = joinpath(pkgdir(MetaCommunityMetrics), "data", "metacomm_rodent_df.csv")
    sample_df = CSV.read(path, DataFrame)
    df = select(sample_df, All())   #convert all columns back to regular vectors
    return df
end


export 
    load_sample_data,
    beta_diversity, spatial_beta_div, temporal_beta_div, 
    create_clusters, plot_clusters, DNCI_multigroup, 
    niche_overlap, 
    prop_patches, 
    CV_meta,
    MVNH_det, MVNH_dissimilarity, average_MVNH_det, average_MVNH_dissimilarity
end