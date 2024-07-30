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

# Include the internal utilities
include("Internal.jl")
using .Internal  


export greet_your_package_name, beta_diversity, mean_spatial_beta_div, mean_temporal_beta_div

# Public function
include("BetaDiversity.jl")
#include("DNCI.jl")
include("functions.jl")



end