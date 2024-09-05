# MetaCommunityMetrics.jl 
*A collection of tools and utilities for analyzing meta-communities in Julia.*

Welcome to the documentation for `MetaCommunityMetrics`. Here you can find guides and reference material on how to use the functions.

## An Overview
This package is a comprehensive toolkit designed to characterize the spatiotemporal structure and dynamics of a metacommunity—a network of communities linked by the dispersal of multiple, interacting species, each with unique niche breadths. It includes functions to calculate a range of specific metrics, which have been previously implemented in R and proven valuable for metacommunity analysis. 

However, they come with high computational costs, especially for large species community datasets. To address this issue, MetaCommunityMetrics.jl was developed in Julia, a programming language known for its efficiency in handling computationally intensive tasks. This implementation significantly improves the efficiency of calculating these metrics, making it a powerful tool for metacommunity analysis. 

These metrics include:
- Averaged beta diversity decomposition in space/time: total diversity, species replacement (turnover), and richness differences for both presence-absence and abundance data
- Dispersal-niche continuum index to evaluate the degree to which communities are influenced by dispersal processes and niche breadth
- Niche overlap indices to determine the extent of niche sharing among species within the metacommunity
- The proportion of habitat patches occupied by each species
- The variability of community composition across different spatial and temporal scales


## Getting Started

### Installation

To install MetaCommunityMetrics, use the following command:

```julia
using Pkg

Pkg.add("MetaCommunityMetrics")

using MetaCommunityMetrics
```

## Function Documentation

- [`Beta Diversity`](https://cralibe.github.io/MetaCommunityMetrics.jl/BetaDiversity/)
- [`DNCI`](https://cralibe.github.io/MetaCommunityMetrics.jl/DNCI/)
- [`Niche Overlap Index`](https://cralibe.github.io/MetaCommunityMetrics.jl/NicheOverlapIndex/)
- [`Occupied Patches Proportion`](https://cralibe.github.io/MetaCommunityMetrics.jl/OccupiedPatchesProportion/)
- [`Variability Metrics`](https://cralibe.github.io/MetaCommunityMetrics.jl/VariabilityMetrics/)

## References
1. Legendre, P. Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography 23, 1324-1334 (2014). https://doi.org:https://doi.org/10.1111/geb.12207
2. Clarke, K. R. Non-parametric multivariate analyses of changes in community structure. Australian Journal of Ecology 18, 117-143 (1993). https://doi.org:https://doi.org/10.1111/j.1442-9993.1993.tb00438.x
3. Gibert, C. & Escarguel, G. PER-SIMPER—A new tool for inferring community assembly processes from taxon occurrences. Global Ecology and Biogeography 28, 374-385 (2019). https://doi.org:https://doi.org/10.1111/geb.12859
4. Vilmi, A. et al. Dispersal–niche continuum index: a new quantitative metric for assessing the relative importance of dispersal versus niche processes in community assembly. Ecography 44, 370-379 (2021). https://doi.org:https://doi.org/10.1111/ecog.05356
5. MacArthur, R. & Levins, R. The limiting similarity, convergence, and divergence of coexisting species. The American Naturalist 101, 377-385 (1967). 
6. Pianka, E. R. (1974). "Niche overlap and diffuse competition." Proceedings of the National Academy of Sciences, 71(5), 2141-2145.
7. Pianka, E. R. (1973). "The Structure of Lizard Communities." Annual Review of Ecology and Systematics, 4(1), 53-74.
8. Ehrlén, J., & Eriksson, O. (2000). Dispersal Limitation and Patchy Occupancy in Forest Herbs. Ecology, 81(6), 1667-1674. https://doi.org:https://doi.org/10.1890/0012-9658(2000)081[1667:DLAPOI]2.0.CO;2
9. Wang, S., Lamy, T., Hallett, L. M. & Loreau, M. Stability and synchrony across ecological hierarchies in heterogeneous metacommunities: linking theory to data. Ecography 42, 1200-1211 (2019). https://doi.org:https://doi.org/10.1111/ecog.04290

