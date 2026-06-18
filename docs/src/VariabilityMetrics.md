# Variability Metrics Functions
```@meta
CurrentModule = MetaCommunityMetrics
```
The variability metrics functions in `MetaCommunityMetrics.jl` are designed to quantify temporal variability in species and community abundances across spatial scales and organizational levels within a metacommunity. By partitioning variability at the species and community levels, both locally and regionally, these metrics implicitly reflect the net outcomes of ecological processes, including biotic interactions, dispersal, and environmental correlation, through the patterns of spatial synchrony and temporal covariance in abundance data. These functions are based on Wang et al. (2019), which provides a framework for partitioning variability at different hierarchical levels within a metacommunity.

## Functionality Overview
In `MetaCommunityMetrics`, the `CV_meta()` function is translated from the R function `var.partition()` in the supplementary material of Wang et al. (2019), published in Ecography (https://doi.org/10.1111/ecog.04290). This function was first translated from R to Julia in August 2024 and is used here for non-commercial scientific research purposes in accordance with Wiley's terms and conditions for use of published content.

The function provides four metrics that are designed to quantify variability at different scales and contexts within a metacommunity:
- `CV_s_l` (local-scale average species variability): the sum of the temporal standard deviations of each species' abundance (or biomass) at each site, divided by the mean total metacommunity abundance (or biomass). It captures species-level temporal variability at the local scale.
- `CV_s_r` (regional-scale average species variability): the sum of the temporal standard deviations of each species' abundance (or biomass) summed across all sites, divided by the mean total metacommunity abundance (or biomass). Because each species' regional abundance is the sum across all sites, this metric implicitly incorporates pairwise temporal covariances between sites for the same species, reflecting spatial synchrony.
- `CV_c_l` (local-scale average community variability): the sum of the temporal standard deviations of total community abundance (or biomass) at each site, divided by the mean total metacommunity abundance (or biomass). Because total community abundance at each site is the sum of all species abundances, by the variance sum law this metric implicitly incorporates pairwise temporal covariances between species within each site, reflecting the net outcome of biotic interactions such as competition and facilitation.
- `CV_c_r` (regional-scale community variability): the temporal standard deviation of total metacommunity abundance (or biomass) across all species and sites, divided by the mean total metacommunity abundance (or biomass). This metric implicitly incorporates all pairwise temporal covariances: between species within sites (biotic interactions), between sites for the same species (spatial synchrony), and between different species at different sites.

## The Function
```@docs
CV_meta
```

## References
- Wang, S., Lamy, T., Hallett, L. M. & Loreau, M. Stability and synchrony across ecological hierarchies in heterogeneous metacommunities: linking theory to data. Ecography 42, 1200-1211 (2019). https://doi.org:https://doi.org/10.1111/ecog.04290

