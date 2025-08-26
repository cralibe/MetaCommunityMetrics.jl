# Hypervolume Functions
```@meta
CurrentModule = MetaCommunityMetrics
```
Hypervolume is a concept originally proposed by Hutchinson (1957). It provides a framework to calculate the volume of a niche for a given species and thus the niche overlap between two species. It helps to infer how niche breath of species has contributed to the co-occurence of different species in the same location at the same time. Unlike the other metrics in this package which require only abundance or occurrence data, this framework also requires environmental data. The hypervolume functions provide by this package are adapted from the R package `MVNH` (https://github.com/lvmuyang/MVNH).

## An Overview
The MVNH framework provides parametric measures for analyzing ecological niches using the multivariate normal distribution model. This framework offers powerful tools for quantifying and comparing the size and dissimilarity of species' niches. There are four hypervolume functions in this package:
- `MVNH_det` calculates the total hypervolume of a species' niche based on the determinant of the covariance matrix.
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

- Measurement standardization: Before analysis, standardize all environmental variables to comparable scales to prevent variables with larger numerical ranges from disproportionately influencing results.

- Species that only occupied one site should be removed before using these functions to avoid undefined values.

## The Functions
```@docs
MVNH_det
MVNH_dissimilarity
average_MVNH_det
average_MVNH_dissimilarity
```
## References
- Lu, M., Winner, K., & Jetz, W. (2021). A unifying framework for quantifying and comparing n‐dimensional hypervolumes. Methods in Ecology and Evolution, 12(10), 1953-1968. https://doi.org/10.1111/2041-210X.13665
