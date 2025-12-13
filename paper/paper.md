---
title: 'MetaCommunityMetrics.jl: A Julia Package for Spatiotemporal Metacommunity Analysis'
tags:
  - Julia
  - community ecology
  - metacommunity
  - summary statistics
  - spatial and temporal dynamics
authors:
  - name: Yan Yin Jenny Cheung
    orcid: 0000-0002-4080-2951
    affiliation: 1
  - name: Laura Melissa Guzman
    orcid: 0000-0003-2198-9738
    corresponding: true
    affiliation: 1
affiliations:
  - index: 1
    name: Department of Entomology, Cornell University, United States
date: 25 September 2025
bibliography: paper.bib
---

# Summary
Community ecologists are often interested in infering ecological processes from observational data. This involves calculating different metrics such as beta diversity.  Calculating these metrics is easy and fast due to the availability of `R` packages. However, community ecology today is growing beyond analysing single datasets. Large biodiversity databases have expanded the sets of spatio-temporal community data, and there is a growing interest in using simulations to study ecological processes. For these applications, `R` packages become slow to use. 

`MetaCommunityMetrics.jl` is a `Julia` package that helps ecologists to analyze large datasets of species observations to understand how community composition varies across space and time. By implementing
computationally optimized versions of established software and new complementary functions that are designed for handling large spatiotemporal datasets, this package enables ecologists to efficiently analyze large-scale spatiotemporal community datasets that would be computationally prohibitive with existing tools. We designed intuitive interfaces that would feel familiar to ecologists transitioning from `R`. Our function naming conventions closely match their `R` equivalents, and we maintained parameter naming and ordering consistency where possible. Some of the metrics we have implemented include beta diversity
decompositions in space and time, dispersal-niche continuum index, and variability metric, but a full list of metrics is listed in \autoref{functions}. Full documentation and examples, including installation instructions and help files, are available at [https://cralibe.github.io/MetaCommunityMetrics.jl/](https://cralibe.github.io/MetaCommunityMetrics.jl/).

\small

: A summary of the metrics and functions available in `MetaCommunityMetrics.jl`. \label{functions}

| Metrics | Functions | Details |
|:--------|:----------|:--------|
| Beta diversity decompositions in space and time | `beta_diversity()` | Re-implementation of R package `adespatial`[@adespatial]  |
| | `spatial_beta_div()` | Re-implementation of `R` script from @guzman2022accounting |
| | `temporal_beta_div()` | Re-implementation of `R` script from @guzman2022accounting |
| Dispersal-niche continuum index (DNCI) | `create_groups()` | A novel grouping function for DNCI analysis|
| | `plot_groups()` | A novel visualization function for DNCI analysis |
| | `DNCI_multigroup()` | Re-implementation of `R` script from R package `DNCImper` [@DNCImper]|
| Niche overlap index | `niche_overlap()` | Novel implementation |
| Occupied patches proportion | `prop_patche()` | Re-implementation  of `R` script from @guzman2022accounting |
| Variability Metric | `CV_meta()` | Re-implementation of `R` script from @wang2019stability |
| Niche hypervolume measurements | `MVNH_det()` | Re-implementation from R package `MVNH`[@lu2021unifying] |
| | `MVNH_dissimilarity()` | Re-implementation from R package `MVNH`[@lu2021unifying] |
| | `average_MVNH_det()` | Novel extension of `MVNH_det()` |
| | `average_MVNH_dissimilarity()` | Novel extension of `MVNH_dissimilarity()` |

\normalsize

# Statement of need
Many valuable `R` [@R] packages have been developed to aid in community analyses by providing ecological metrics that can summarize processes from biodiversity patterns and have been widely adopted in the field, such as `adespatial` [@adespatial], `codyn` [@codyn], and `vegan` [@vegan]. However, conducting community analyses in R is often computationally expensive, especially when working with large datasets, a limitation that becomes critical in workflows requiring repeated metric calculations. These include simulation studies for method validation, bootstrap resampling for uncertainty estimation [@efron1992bootstrap], null model approaches for statistical testing [@gotelli1996null], and the analysis of large biodiversity databses, such as [BioTime](https://doi.org/10.1111/geb.70003), all of which can create significant computational bottlenecks when processing multiple large datasets iteratively. To improve the efficiency of these computationally intensive analyses, `Julia` becomes a logical choice for re-implementations. `Julia` is known for its just-in-time compilation that optimizes code execution to levels comparable with lower-level languages like C or Fortran [@Julia-2017], yet Julia remains relatively uncommon in ecological workflows. During development, we conducted a comprehensive review of existing `Julia` packages under the [EcoJulia organization](https://ecojulia.org/), a community-driven effort providing tools for ecological and environmental analysis, and the [PoisotLab](https://github.com/PoisotLab), which develop many `Julia` packages focusing on quantitative and computational ecology. We ensure that our package builds beyond the tools currently existing in `Julia` even if they share similar ecological concepts. For example, our beta diversity metrics share a similar concept, beta diversity, with `Diversity.jl` [@DiversityJL]. However, our implementation decomposes beta diversity into species replacement and richness differences (or abundance difference when abundance data is used instead of occurrence data) [@baselga2010partitioning; @baselga2012relationship; @baselga2013separating; @legendre2014interpreting], which differentiates our approach from `Diversity.jl`, which emphasizes alpha, beta, and gamma diversity measures at the metacommunity and community levels [@reeve2014partition]. To our knowledge, no existing Julia packages implement the specific metrics listed in \autoref{functions}.

`MetaCommunityMetrics.jl` delivers substantial computational advantages over existing tools. For example, we re-implemented the `beta.div.comp` function from the `R` package `adespatial` [@adespatial], achieving up to $\sim$ 58× faster execution times while maintaining methodological consistency (see \autoref{speedup}). As shown in \autoref{speedup}, performance gains range from $\sim$ 3× to 58× across all re-implemented functions using [three testing datasets](https://cralibe.github.io/MetaCommunityMetrics.jl/Benchmarking/#Datesets-used-for-this-benchmark) (large: 53,352 observations; medium: 26,676 observations; small: 5,325 observations), positioning `MetaCommunityMetrics.jl` as an efficient solution for diverse ecological research needs, from standard biodiversity analyses to computationally intensive simulation studies and large-scale data processing.

Beyond performance improvements, we extended pre-existing `R` implementations to handle spatiotemporal datasets directly. We added temporal, spatial, and community-level aggregation methods to metrics originally restricted to a single time point, site, and species pair, eliminating extensive data preprocessing. For instance, beta diversity decompositions in `R` operate on a single temporal snapshot of a metacommunity, while niche hypervolume measurements calculate volumes for an individual species and dissimilarity between a single species pair. We extended these with `spatial_beta_div()` and `temporal_beta_div()`, which calculate beta diversity decompositions across time and space, respectively. Similarly, `average_MVNH_det()` calculates average niche volume across all species, and `average_MVNH_dissimilarity()` computes average niche dissimilarity across all species pairs. For DNCI calculations, we developed `create_groups()` and `plot_groups()` to meet the site grouping requirements suggested by @vilmi2021dispersal and to visualize the resulting groupings.

## Additional handling of edge cases when calculating DNCI
We extended our DNCI implementation to handle edge cases common in simulated community datasets (e.g., single-species communities, insufficient permutation variation) by returning status flags that identify when the standard DNCI calculation is not possible. See [documentation](https://cralibe.github.io/MetaCommunityMetrics.jl/DNCI/) for details on all five edge cases.

![A plot showing the speedup of all benchmarked functions across three [datasets](https://cralibe.github.io/MetaCommunityMetrics.jl/Benchmarking/#Benchmarked-using-Large-Dataset): small (5,335 observations), medium (26,676 observations), and large (53,352 observations). Speedup is calculated as the `R` median execution time divided by the `Julia` median execution time. Median speedup and its confidence interval of each comparison is shown. The grey dashed line shows where speedup equals to 1, that is both `R` and `Julia`'s implementations require the same amount of time. The calcualation of beta diversity has two variants: `quant = true` (abundance data is used) and `quant = false`(occurrence data is used). Details about this benchmarking can be found in the documentation under the [benchmark results section](https://cralibe.github.io/MetaCommunityMetrics.jl/Benchmarking/).\label{speedup}](speedup.pdf){width=100%}

## Validation 
All re-implementations were validated against their R equivalents (see [documentation](https://cralibe.github.io/MetaCommunityMetrics.jl/Validation/)).

# Acknowledgement
We are grateful to the authors of the pre-existing `R` packages/implementations that we have implemented in `MetaCommunityMetrics.jl`. We express our gratitude to everyone who has tried out this package and provided feedback on how to improve it. LMG acknowledges funding from USC and Cornell University, YYJC acknowledges funding from the Wrigley Institute.  

# References
