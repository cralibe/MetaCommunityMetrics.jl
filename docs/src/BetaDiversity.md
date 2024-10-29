# Beta Diversity Functions
```@meta
CurrentModule = MetaCommunityMetrics
```
Beta diversity is a fundamental concept in ecology that quantifies the variation in species composition between different plots, or over time. In the context of metacommunity analysis, beta diversity functions help to assess how community composition changes spatially across different locations or temporally within a metacommunity.

The `beta_diversity` functions in `MetaCommunityMetrics` are adapted from the `beta.div.comp` function in the R package `adespatial`. These methods, originally developed by Pierre Legendre (2014), are implemented in Julia to provide a more efficient means of computation for large-scale datasets. The functions use indices from the Podani family, Jaccard-based indices, and Ruzicka-based indices to calculate total beta diversity and its components: replacement and richness difference.

## Choosing the Right Function
- Use `beta_diversity` for a general, comprehensive measure of beta diversity across your dataset. This function provides an overall assessment of how species composition varies between sites or over time, capturing both replacement (the turnover of species) and richness difference(the difference in species richness). 
- Use `spatial_beta_div` to comparing diversity between different spatial locations of a metacommunity. 
- Use `temporal_beta_div` to track how diversity changes over time of a metacommunity. 

## The Functions
```@docs
beta_diversity
spatial_beta_div
temporal_beta_div
```
## References
- Guzman, L. M. et al. Accounting for temporal change in multiple biodiversity patterns improves the inference of metacommunity processes. Ecology 103, e3683 (2022). https://doi.org:https://doi.org/10.1002/ecy.3683
- Legendre, P. Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography 23, 1324-1334 (2014). https://doi.org:https://doi.org/10.1111/geb.12207
