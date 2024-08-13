# MetaCommunityMetrics.jl 
*A collection of tools and utilities for analyzing meta-communities in Julia.*

Welcome to the documentation for `MetaCommunityMetrics`. Here you can find guides and reference material on how to use the functions.

## An Overview
This package is a comprehensive toolkit designed to characterize the spatiotemporal structure and dynamics of a metacommunity—a network of communities linked by the dispersal of multiple, interacting species, each with unique niche breadths. It includes functions to calculate a range of specific metrics, which have been previously implemented in R and proven valuable for metacommunity analysis. 

However, they come with high computational costs, especially for large species community datasets. To address this issue, MetaCommunityMetrics.jl was developed in Julia, a programming language known for its efficiency in handling computationally intensive tasks. This implementation significantly improves the efficiency of calculating these metrics, making it a powerful tool for metacommunity analysis. 

These metrics include:
- Beta diversity and its components: total diversity, species replacement (turnover), and richness differences for both presence-absence and abundance data
- Dispersal-niche continuum index to evaluate the degree to which communities are influenced by dispersal processes and niche breadth
- Niche overlap indices to determine the extent of niche sharing among species within the metacommunity
- The proportion of habitat patches occupied by each species
- The variability of community composition across different spatial and temporal scales


## Getting Started

### Installation

To install MetaCommunityMetrics, use the following command:

```julia
using Pkg

Pkg.add("MetaCommunityMetrics")

using MetaCommunityMetrics
```

## Function Documentation

- [Beta Diversity Functions](@ref BetaDiveristy)
- [Dispersal-Niche Continuum Index (DNCI) Functions](@ref DNCI)
- [Niche Overlap Index Function](@ref NicheOverlapIndex)
- [Occupied Patches Proportion Function](@ref OccupiedPatchesProportion)
- [Variability Metrics Functions](@ref VariabilityMetrics)