# Variability Metrics Functions
```@meta
CurrentModule = MetaCommunityMetrics
```

The variability metrics functions in `MetaCommunityMetrics` are designed to capture changes in dispersal and density-dependent biotic interactions by investigating temporal variability and synchrony across spatial scales and organizational levels within a metacommunity. These functions are based on the work of Wang et al. (2019), which provides a framework for quantifying variability at different scales and contexts within a community.

## An Overview
In `MetaCommunityMetrics`, the `CV_meta_simple` function is directly adapted from the R function `var.partition` in Wang et al. (2019). This function is designed with computational efficiency in mind, particularly for large datasets, by avoiding the calculation of all covariances between species. This approach ensures faster performance while still providing valuable insights into variability across different scales within a metacommunity.

In contrast, the `CV_meta function` extends the analysis by including the calculation of all covariances between species, offering a more detailed and comprehensive examination of variability. This approach captures interactions between species that may be overlooked by more streamlined methods. While the `CV_meta` function provides a richer analysis, the `CV_meta_simple` function remains a valuable tool when computational efficiency is a priority.

These metrics are designed to quantify variability at different scales and contexts within the community:
- Local-scale average species variability (`CV_s_l`)
- Regional-scale average species variability (`CV_s_r`)
- Local-scale average community variability (`CV_c_l`)
- Regional-scale community variability (`CV_c_r`)

The variability metrics are calculated as follows:

The variability metrics `CV_s_l`, `CV_s_r`, `CV_c_l`, and `CV_c_r` are set to zero whenever the mean abundance equals zero at any combination of spatial scales (a patch/all patches) and species number (a species/the whole community). This approach allows us to assess the impact of spatial scale on variability and to understand how different factors influence community dynamics across scales.

## The Function
```@docs
CV_meta
CV_meta_simple
```

## References
- Wang, S., Lamy, T., Hallett, L. M. & Loreau, M. Stability and synchrony across ecological hierarchies in heterogeneous metacommunities: linking theory to data. Ecography 42, 1200-1211 (2019). https://doi.org:https://doi.org/10.1111/ecog.04290

