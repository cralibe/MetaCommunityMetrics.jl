# Niche Overlap Index Function
Niche overlap refers to the extent to which different species use the same resources or occupy similar ecological niches. High niche overlap might indicate intense competition, whereas low overlap suggests niche differentiation, allowing species to coexist by minimizing direct competition.

To capture the changes in density-independent abiotic response, also known as niche breadth, this implementation uses Pianka's Niche Overlap Index, as suggested by Pianka (1973). The summary statistics of this index include the mean, maximum, and minimum values across all species pairs, providing a comprehensive understanding of niche sharing within the community.
Key References

While the specific method implemented here is adapted to the Julia language for efficiency and scalability, the concept of niche overlap has a long history in ecological studies. Researchers interested in the theoretical foundations of niche overlap and its applications in community ecology may refer to the following works:

- MacArthur, R., & Levins, R. (1967). "The Limiting Similarity, Convergence, and Divergence of Coexisting Species." The American Naturalist, 101(921), 377-385.
- Pianka, E. R. (1974). "Niche overlap and diffuse competition." Proceedings of the National Academy of Sciences, 71(5), 2141-2145.
- Pianka, E. R. (1973). "The Structure of Lizard Communities." Annual Review of Ecology and Systematics, 4(1), 53-74.

## Functionality Overview

The Niche Overlap Index functions in MetaCommunityMetrics provide a robust framework for calculating niche overlap based on species abundance or presence-absence data. These functions allow you to evaluate how species share ecological space within a metacommunity, offering valuable insights into community dynamics and species interactions.

## The Function
- `niche_overlap`: This function calculates the Niche Overlap Index for a given set of species across different patches. It provides metrics for mean, minimum, and maximum niche overlap, allowing for a comprehensive assessment of niche sharing within the community.

```@meta
CurrentModule = MetaCommunityMetrics
```
```@docs
niche_overlap(abundance::AbstractVector, species::Union{AbstractVector, String}, patch::Union{AbstractVector, String}, time::AbstractVector)
```

