# Beta Diveristy
This section provides detailed documentation for all the public functions, types, and modules available in `MetaCommunityMetrics`. The documentation is automatically extracted from the source code docstrings.

## Beta Diversity and its components

Calculates beta diversity for given ecological data. This function supports both binary (presence/absence) and quantitative data.

```@docs
beta_diversity
```
## Mean Spatial Beta Diversity
The `mean_spatial_beta_div` function computes the mean of spatial beta diversity and its components using the `beta_diversity` function across different time points. It is used to analyze spatial patterns in ecological data over time.

```@docs
mean_spatial_beta_div
```

## Mean Temporal Beta Diversity
The `mean_temporal_beta_div` function calculates the mean of temporal beta diversity and its components using the `beta_diversity` function across various sites. This function helps in understanding temporal variations in ecological communities.

```@docs
mean_temporal_beta_div
```