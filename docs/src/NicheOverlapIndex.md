# Niche Overlap Index Function
```@meta
CurrentModule = MetaCommunityMetrics
```
Niche overlap refers to the extent to which different species use the same resources or occupy similar ecological niches. High niche overlap might indicate intense competition, whereas low overlap suggests niche differentiation, allowing species to coexist by minimizing direct competition.

## Functionality Overview
To capture the changes in density-independent abiotic response, also known as niche breadth, this implementation developed based on Pianka's Niche Overlap Index, as suggested by Pianka (1973). The summary statistics of this index include the mean, maximum, and minimum values across all species pairs, providing a comprehensive understanding of niche sharing within the community.

Since most ecological datasets lack direct measurements of species-specific resource utilization, we adapted the original niche overlap index using abundance data. Specifically, our implementation assumes all species share the same consumption rate and type I functional response, such that species with higher abundance have proportionally higher resource utilization at any given site-time combination. This assumption allows abundance data to serve as a proxy for resource utilization patterns, enabling calculation of niche overlap index from observed abundance data. Keep in mind that this approach will not be able to reflect the fact that some species may be more efficient in consumption than others. Therefore, this approach works best when comparing species with similar body sizes, metabolic rates, or feeding strategies, as these factors strongly influence the abundance-consumption relationship (Brown et al., 2004;  White et a., 2007; Hudson et al., 2013; Kalinoski et al., 2016; Abrams, 2022).


## The Function
```@docs
niche_overlap
```
## References
1. MacArthur, R., & Levins, R. (1967). The limiting similarity, convergence, and divergence of coexisting species. The american naturalist, 101(921), 377-385. https://doi.org/10.1086/282505
2. Pianka, E. R. (1974). Niche overlap and diffuse competition. Proceedings of the National Academy of Sciences, 71(5), 2141-2145. https://doi.org/10.1073/pnas.71.5.2141
3. Pianka, E. R. (1973). The structure of lizard communities. Annual review of ecology and systematics, 53-74. https://www.jstor.org/stable/2096804
4. Brown, J. H., Gillooly, J. F., Allen, A. P., Savage, V. M., & West, G. B. (2004). Toward a metabolic theory of ecology. Ecology, 85(7), 1771-1789. https://doi.org/https://doi.org/10.1890/03-9000 
5. White, E. P., Ernest, S. M., Kerkhoff, A. J., & Enquist, B. J. (2007). Relationships between body size and abundance in ecology. Trends in ecology & evolution, 22(6), 323-330.
6. Hudson, L. N., Isaac, N. J., & Reuman, D. C. (2013). The relationship between body mass and field metabolic rate among individual birds and mammals. Journal of Animal Ecology, 82(5), 1009-1020.
7. Kalinoski, R. M., & DeLong, J. P. (2016). Beyond body mass: how prey traits improve predictions of functional response parameters. Oecologia, 180(2), 543-550.
8. Abrams, P. A. (2022). Food web functional responses. Frontiers in Ecology and Evolution, 10, 984384.

