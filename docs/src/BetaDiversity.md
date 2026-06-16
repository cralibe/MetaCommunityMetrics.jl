# Beta Diversity Functions
```@meta
CurrentModule = MetaCommunityMetrics
```
Beta diversity is a fundamental concept in ecology that quantifies the variation in species composition between different plots, or over time. In the context of community ecology, beta diversity functions help to assess how community composition changes spatially or temporally.

The `beta_diversity()` functions in `MetaCommunityMetrics` are translated from the `beta.div.comp()` function in the R package `adespatial` 
(https://github.com/adeverse/adespatial), authored by Pierre Legendre (Université de Montréal) and licensed under GPL-3. The underlying methods are described in Legendre (2014). These functions were translated from R to Julia in August 2024 to provide more efficient computation for large-scale datasets, and are redistributed here under GPL-3 in accordance with the terms of the original license. 
The functions use indices from the Podani family, Jaccard-based indices, and Ruzicka-based indices to calculate total beta diversity and its components: replacement and richness difference.

## Functionality Overview
- Use `beta_diversity()` for a general, comprehensive measure of beta diversity across your dataset. This function provides an overall assessment of how species composition varies between sites or over time, capturing both replacement (the turnover of species) and richness difference (or abundance difference when abundance data is used instead of occurrence data). 
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
