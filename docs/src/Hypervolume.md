# Hypervolume Functions
```@meta
CurrentModule = MetaCommunityMetrics
```
Hypervolume is a method originally proposed by Hutchinson (1957). It provide methods to calculate the volume of a niche for a given specice and thus the ncihe overlap between two species. It helps to infer how niche breath of species has contribute to the co-occurance of different species in the same location at the same time. The hypervolume functions provide by this package are adapted from the R package `MVNH` (https://github.com/lvmuyang/MVNH).

## An Overview
The MVNH framework provides parametric measures for analyzing ecological niches using the multivariate normal distribution model. This framework offers powerful tools for quantifying and comparing the size and dissimilarity of species' niches, with each measure being partitionable into biologically meaningful components.

The framework models a species' niche as a multivariate normal distribution in environmental space, where:
- Each environmental variable represents one dimension of the niche.
- The mean vector represents the niche optimum.
- The covariance matrix represents the niche breadth and shape.

There are four hypervolume functions in this package:
- `MVNH_det` calculates the total hypervolume of a species' niche based on the determinant of the covariance matrix (generalized variance). This measure can be partitioned into:
- `MVNH_dissimilarity` calculates the Bhattacharyya distance between two species' niches, providing a comprehensive measure of niche differentiation. 
- `average_MVNH_det` calculates the mean hypervolume across multiple species in a community, providing an overall measure of niche size at the community level.
- `average_MVNH_dissimilarity` calculates the mean Bhattacharyya distance between all unique pairs of species in a community, providing a measure of overall niche differentiation.

### Practical Considerations
- Statistical assumptions: This framework relies on multivariate normal distribution of environmental data.
    - When encountering skewed variables, apply appropriate transformations to achieve approximately normal distributions
    - Be aware that as you increase the number of variables, you face greater challenges with:
        - Variable interdependence (collinearity) which can drive determinant values toward zero
        - Potential violations of the multivariate normality assumption
    - Address variable interdependence through either:
        - Thoughtful pre-selection of ecologically meaningful variables with direct influence on species distributions
        - Application of dimension reduction methods such as PCA (principal component analysis)
            - Important note: PCA creates orthogonal axes, which forces the correlation component to 1.0, eliminating correlation structure information.
            - For datasets containing multiple distinct groups of related environmental variables (such as climate factors, soil properties, or topographic features), consider using generalized canonical variables to identify the most representative variables within each natural category while preserving the ecological relationships between different variable groups.

- Measurement standardization: Before analysis, standardize all environmental variables to comparable scales to prevent variables with larger numerical ranges from disproportionately influencing results.

## The Functions
```@docs
MVNH_det
MVNH_dissimilarity
average_MVNH_det
average_MVNH_dissimilarity
```
## References
- Lu, Muyang, Kevin Winner, and Walter Jetz. A unifying framework for quantifying and comparing n‐dimensional hypervolumes. Methods in Ecology and Evolution 12.10, 1953-1968 (2021). https://doi.org/10.1111/2041-210X.13665
