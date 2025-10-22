# Variability Metrics Functions
```@meta
CurrentModule = MetaCommunityMetrics
```

The variability metrics functions in `MetaCommunityMetrics` are designed to capture changes in dispersal and density-dependent biotic interactions by investigating temporal variability across spatial scales and organizational levels within a metacommunity. These functions are based on Wang et al. (2019), which provides a framework for quantifying variability at different scales and contexts within a community.

## Functionality Overview
In `MetaCommunityMetrics`, the `CV_meta()` function is directly adapted from the R function `var.partition()` in Wang et al. (2019).

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

