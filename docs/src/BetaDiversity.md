# Beta Diversity Functions
```@meta
CurrentModule = MetaCommunityMetrics
```
Beta diversity is a fundamental concept in ecology that quantifies the variation in species composition between different plots, or over time. In the context of community ecology, beta diversity functions help to assess how community composition changes spatially or temporally.

The `beta_diversity()` functions in `MetaCommunityMetrics` are adapted from the `beta.div.comp()` function in the R package `adespatial`. These methods, originally developed by Legendre (2014), are implemented in Julia to provide a more efficient means of computation for large-scale datasets. The functions use indices from the Podani family, Jaccard-based indices, and Ruzicka-based indices to calculate total beta diversity and its components: replacement and richness difference.

## Functionality Overview
- Use `beta_diversity()` for a general, comprehensive measure of beta diversity across your dataset. This function provides an overall assessment of how species composition varies between sites or over time, capturing both replacement (the turnover of species) and richness difference(the difference in species richness). 
- Use `spatial_beta_div()` to compare local communities of a metacommunity across different spatial locations with individual species abundance in each site aggregated across all time points. 
- Use `temporal_beta_div()` to compare metacommunities over time with individual species abundance aggregated across all sites. 

## The Functions
```@docs
beta_diversity
spatial_beta_div
temporal_beta_div
```
## References
- Legendre, P. (2014). Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography, 23(11), 1324-1334. https://doi.org:https://doi.org/10.1111/geb.12207
