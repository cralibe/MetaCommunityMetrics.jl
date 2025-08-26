# Niche Overlap Index Function
```@meta
CurrentModule = MetaCommunityMetrics
```
Niche overlap refers to the extent to which different species use the same resources or occupy similar ecological niches. High niche overlap might indicate intense competition, whereas low overlap suggests niche differentiation, allowing species to coexist by minimizing direct competition.

## Functionality Overview
To capture the changes in density-independent abiotic response, also known as niche breadth, this implementation developed based on Pianka's Niche Overlap Index, as suggested by Pianka (1973). The summary statistics of this index include the mean, maximum, and minimum values across all species pairs, providing a comprehensive understanding of niche sharing within the community.

Since most ecological datasets lack direct measurements of species-specific resource utilization, we developed an adaptation of the original niche overlap index that can be applied to abundance data. Specifically, our implementation assumes all species share the same consumption rate and type I functional response, such that species with higher abundance have proportionally higher resource utilization at any given site-time combination. This assumption allows abundance data to serve as a proxy for resource utilization patterns, enabling calculation of niche overlap index from observed abundance data. Keep in mind that this approach will not be able to reflect the fact that some species may be more efficient in consumption than others. Therefore, this approach works best when comparing species with similar body sizes, metabolic rates, or feeding strategies, as these factors strongly influence the abundance-consumption relationship.


## The Function
```@docs
niche_overlap
```
## References
1. MacArthur, R., & Levins, R. (1967). The limiting similarity, convergence, and divergence of coexisting species. The american naturalist, 101(921), 377-385. https://doi.org/10.1086/282505
2. Pianka, E. R. (1974). Niche overlap and diffuse competition. Proceedings of the National Academy of Sciences, 71(5), 2141-2145. https://doi.org/10.1073/pnas.71.5.2141
3. Pianka, E. R. (1973). The structure of lizard communities. Annual review of ecology and systematics, 53-74. https://www.jstor.org/stable/2096804

