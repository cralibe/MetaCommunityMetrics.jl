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
using Plots
using Random
using Pkg.Artifacts


# Include the internal utilities
include("Internal.jl")
using .Internal  


# Public function
include("BetaDiversity.jl")
include("DNCI.jl")
include("NicheOverlapIndex.jl")
include("OccupiedPatchesProportion.jl")
include("VariabilityMetrics.jl")

# Function to load sample data
function load_sample_data()
    artifact_path = artifact"rodent_abundance_data"
    sample_df = CSV.read(joinpath(artifact_path, "rodent_abundance_data.csv"), DataFrame)
    return sample_df
end


export 
    beta_diversity, mean_spatial_beta_div, mean_temporal_beta_div, 
    create_clusters, plot_clusters, DNCI_multigroup, 
    niche_overlap, prop_patches, CV_meta, CV_meta_simple,
    load_sample_data

end