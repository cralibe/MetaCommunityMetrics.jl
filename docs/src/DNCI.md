# Dispersal–niche continuum index

Dispersal–niche continuum index (DNCI) estimates whether dispersal or niche processes dominate community assembly in a metacommunity. It is implemented based on the methodology described in the following paper:

> A. Vilmi, C. Gibert, G. Escarguel, K. Happonen, J. Heino, A. Jamoneau, et al.
. (2021). Dispersal–niche continuum index: a new quantitative metric for assessing the relative importance of dispersal versus niche processes in community assembly. Ecography, 44, 370-379. DOI: [doi.org/10.1111/ecog.05356](https://doi.org:https://doi.org/10.1111/ecog.05356)

The R packages, named `[DNCImper]` provides an implementation of the DNCI which has influenced the development of this function in the Julia environment. Since `[DNCImper]` also depends on the R packages `[vegan]`. Therefore, this Julia implementation contains adaptations from these two R packages to optimize performance in Julia or to utilize Julia-specific features. 
. These packages can be found at:
- `[DNCImper]`: [https://github.com/Corentin-Gibert-Paleontology/DNCImper]
- `[vegan]`: [https://cran.r-project.org/web/packages/vegan/index.html]


## Grouping
Grouping of patches are required to calculate DNCI. At least two group are needed, and each group should have at least 5 patches. These groups should have less than 40% and 30% differences in their taxa/species and site numbers, respectively. Please refer to Vilmi(2021) for detailed information. This Julia implementation provide a function to create grouping based on the mention requirement and physical locations of the patches (e.g. Latitdue and Longitude) using agglomerative hierarchical clustering with complete linkage.

```@doc
    create_clusters
```


