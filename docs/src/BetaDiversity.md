# Beta Diversity Functions
The beta diversity functions in `MetaCommunityMetrics` are adapted from the beta.div.comp function in the R package `adespatial`. These methods were originally developed and detailed by Pierre Legendre in his 2014 paper, "Interpreting the replacement and richness difference components of beta diversity" (Global Ecology and Biogeography, 23, 1324-1334). The implementation in Julia provides a more efficient means to compute these metrics, particularly for large-scale datasets, while maintaining the robustness of the original methodology.

Beta diversity is a fundamental concept in ecology that quantifies the variation in species composition between different habitats, plots, or over time. In the context of metacommunity analysis, beta diversity functions help to assess how community composition changes spatially across different locations or temporally within a given location.

## Choosing the Right Function
- Use `beta_diversity` when you want a general overview of diversity across your dataset.
- Opt for `mean_spatial_beta_div` when your focus is on comparing diversity between different spatial locations.
- Select `mean_temporal_beta_div` to track how diversity changes over time within the same location.

```@meta
CurrentModule = MetaCommunityMetrics
```
## The Functions
```@docs
beta_diversity(mat::Matrix; quant::Bool)
mean_spatial_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String}; quant::Bool)
mean_temporal_beta_div(abundance::AbstractVector, time::AbstractVector, patch::Union{AbstractVector, String}, species::Union{AbstractVector, String};quant::Bool)
```
