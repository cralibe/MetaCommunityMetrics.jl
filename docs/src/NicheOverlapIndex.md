# Niche Overlap Index Function
```@meta
CurrentModule = MetaCommunityMetrics
```
Niche overlap refers to the extent to which different species use the same resources or occupy similar ecological niches. High niche overlap might indicate intense competition, whereas low overlap suggests niche differentiation, allowing species to coexist by minimizing direct competition.

To capture the changes in density-independent abiotic response, also known as niche breadth, this implementation uses Pianka's Niche Overlap Index, as suggested by Pianka (1973). The summary statistics of this index include the mean, maximum, and minimum values across all species pairs, providing a comprehensive understanding of niche sharing within the community.

## Functionality Overview

The Niche Overlap Index functions in MetaCommunityMetrics provide a robust framework for calculating niche overlap based on species abundance or presence-absence data. These functions allow you to evaluate how species share ecological space within a metacommunity, offering valuable insights into community dynamics and species interactions.

## The Function
```@docs
niche_overlap
```
## References
1. MacArthur, R. & Levins, R. The limiting similarity, convergence, and divergence of coexisting species. The American Naturalist 101, 377-385 (1967). 
2. Pianka, E. R. (1974). "Niche overlap and diffuse competition." Proceedings of the National Academy of Sciences, 71(5), 2141-2145.
3. Pianka, E. R. (1973). "The Structure of Lizard Communities." Annual Review of Ecology and Systematics, 4(1), 53-74.

