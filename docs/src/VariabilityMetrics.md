# Variability Metrics Functions
```@meta
CurrentModule = MetaCommunityMetrics
```
The variability metrics functions in `MetaCommunityMetrics.jl` are designed to quantify temporal variability in species and community abundances across spatial scales and organizational levels within a metacommunity. By partitioning variability at the species and community levels, both locally and regionally, these metrics implicitly reflect the net outcomes of ecological processes such as biotic interactions and spatial synchrony. These functions are based on Wang et al. (2019), which provides a framework for partitioning variability at different hierarchical levels within a metacommunity.

## Functionality Overview
In `MetaCommunityMetrics`, the `CV_meta()` function is translated from the R function `var.partition()` in the supplementary material of Wang et al. (2019), published in Ecography (https://doi.org/10.1111/ecog.04290). This function was first translated from R to Julia in August 2024 and is used here for non-commercial scientific research purposes in accordance with Wiley's terms and conditions for use of published content.

The function provides four metrics that are designed to quantify variability at different scales and contexts within a metacommunity:
- Local-scale average species variability (`CV_s_l`)
- Regional-scale average species variability (`CV_s_r`)
- Local-scale average community variability (`CV_c_l`)
- Regional-scale community variability (`CV_c_r`)

## The Function
```@docs
CV_meta
```

## References
- Wang, S., Lamy, T., Hallett, L. M. & Loreau, M. Stability and synchrony across ecological hierarchies in heterogeneous metacommunities: linking theory to data. Ecography 42, 1200-1211 (2019). https://doi.org:https://doi.org/10.1111/ecog.04290

