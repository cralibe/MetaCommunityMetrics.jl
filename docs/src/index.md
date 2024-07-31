
# MetaCommunityMetrics.jl 
*A collection of tools and utilities for analyzing meta-communities in Julia.*

Welcome to the documentation for MetaCommunityMetrics. Here you can find guides and reference material on how to use the functions.

## Getting Started

###Installation

To install MetaCommunityMetrics, use the following command:

```julia
using Pkg
Pkg.add("MetaCommunityMetrics")
using MetaCommunityMetrics
```

```@meta
CurrentModule = MetaCommunityMetrics
```
##Functions
```@docs
beta_diversity(mat::Matrix; quant::Bool)
mean_spatial_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String}; quant::Bool)
mean_temporal_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String};quant::Bool)
create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int})
plot_clusters(grouped_data::DataFrame)
DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000, count::Bool=true)
```



