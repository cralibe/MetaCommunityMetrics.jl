# Dispersal-Niche Continuum Index (DNCI) Functions
The Dispersal-Niche Continuum Index (DNCI) functions in MetaCommunityMetrics are adapted from the R package `DNCImper`: Assembly process identification based on SIMPER analysis. These methods, originally developed by Clarke in 1993 and later refined by Gibert & Escarguel in 2019 and Vilmi, Gibert et al. in 2021, offer powerful tools for identifying the processes underlying species assembly in metacommunities.

## Background
The DNCI functions is built around the Per-SIMPER and DNCI analyses. PerSIMPER, based on the Similarity Percentage (SIMPER) analysis developed by Clarke (1993), assesses the contribution of individual taxa to overall dissimilarity (OAD) between groups of assemblages. PerSIMPER enhances this by comparing empirical SIMPER plots with randomized plots generated through matrix permutation, which helps identify whether niche, dispersal, or both processes are driving community assembly.

The DNCI (Dispersal-Niche Continuum Index) further extends this approach by transforming the qualitative results of PerSIMPER into a quantitative index, providing a straightforward measure of the influence of niche and dispersal processes on community structure.

## Functionality Overview

The DNCI functions in `MetaCommunityMetrics` allow you to analyze the processes driving species assembly within your dataset. By comparing empirical data with randomized permutations, you can determine the extent to which niche and dispersal processes have influenced the structure of metacommunities.

## The Functions 
- create_clusters: Groups sampling locations based on their spatial attributes and species richness, which can then be used to assess DNCI.
- plot_clusters: Visualizes the clusters created, allowing for an intuitive understanding of spatial groupings.
- DNCI_multigroup: Computes the Dispersal-Niche Continuum Index (DNCI) across multiple groups, helping to quantify the relative influence of niche versus dispersal processes.

```@meta
CurrentModule = MetaCommunityMetrics
```
```@docs
create_clusters(time::Vector{Int}, latitude::Vector{Float64}, longitude::Vector{Float64}, patch::Vector{Int})
plot_clusters(grouped_data::DataFrame)
DNCI_multigroup(comm::Matrix, groups::Vector, Nperm::Int=1000, count::Bool=true)
```

