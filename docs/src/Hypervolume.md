# Hypervolume Functions
```@meta
CurrentModule = MetaCommunityMetrics
```
Hypervolume is a method originally proposed by Hutchinson (1957). It provide methods to calculate the volume of a niche for a given specice and thus the ncihe overlap between two species. It helps to infer how niche breath of species has contribute to the co-occurance of different species in the same location at the same time. The hypervolume functions provide by this package are adapted from the R package `MVNH` (https://github.com/lvmuyang/MVNH)

## The Functions
```@docs
MVNH_det
MVNH_dissimilarity
average_MVNH_det
average_MVNH_dissimilarity
```