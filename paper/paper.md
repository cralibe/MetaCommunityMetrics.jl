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
Community ecology studies how species coexist and change over time and space [@leibold2004metacommunity]. `MetaCommunityMetrics.jl` is a `Julia` [@Julia-2017] package that helps ecologists to analyze large datasets of species observations to understand how community composition varies across different sites and time points. This package addresses key ecological questions in community ecology from observational data: What drives changes in community composition? How do species coexist locally and regionally over time? How do species dispersal abilities and habitat preferences shape community patterns? By implementing computationally optimized versions of established software and new complementary functions that are designed for handling large spatiotemporal datasets, this package enables ecologists to efficiently analyze large-scale spatiotemporal community datasets that would be computationally prohibitive with existing tools. We designed intuitive interfaces that would feel familiar to ecologists transitioning from `R`. Our function naming conventions closely match their `R` equivalents, and we maintained parameter naming and ordering consistency where possible. For example, if an `R` function accepts parameters named `species` and `site` in that order, our `Julia` implementation follows the same pattern, making it easier for ecologists to use our package. \autoref{functions} shows a summary of that metrics and functions are avaliable in `MetaCommunityMetrics.jl`. Full documentation and examples, including installation instructions and help files, are available at https://cralibe.github.io/MetaCommunityMetrics.jl/.

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
Many valuable `R` [@R] packages have been developed to aid in community analyses by providing ecological metrics that can summarize processes from biodiversity patterns and have been widely adopted in the field, such as `adespatial` [@adespatial], `codyn` [@codyn], and `vegan` [@vegan]. However, conducting community analyses in R is often computationally expensive, especially when working with large datasets, a limitation that becomes critical in workflows requiring repeated metric calculations. These include simulation studies for method validation, bootstrap resampling for uncertainty estimation [@efron1992bootstrap], and null model approaches for statistical testing [@gotelli1996null], all of which can create significant computational bottlenecks when processing multiple large datasets iteratively. To improve the efficiency of these computationally intensive analyses, `Julia` becomes a logical choice for re-implementations. `Julia` is known for its just-in-time compilation that optimizes code execution to levels comparable with lower-level languages like C or Fortran [@Julia-2017], yet Julia remains relatively uncommon in ecological workflows. During development, we conducted a comprehensive review of existing `Julia` packages under the [EcoJulia organization](https://ecojulia.org/), a community-driven effort providing tools for ecological and environmental analysis, and the [PoisotLab](https://github.com/PoisotLab), which develop many `Julia` packages focusing on quantitative and computational ecology. We ensure that our package builds beyond the tools currently existing in `Julia` even if they share similar ecological concepts. For example, our beta diversity metrics share a similar concept, beta diversity, with `Diversity.jl` [@DiversityJL]. However, our implementation decomposes beta diversity into species replacement and richness differences, which differentiates our approach from `Diversity.jl`, which emphasizes alpha, beta, and gamma diversity measures at the metacommunity and community levels. 

`MetaCommunityMetrics.jl` delivers substantial computational advantages over existing tools. For example, we re-implemented the `beta.div.comp` function from the `R` package `adespatial` [@adespatial], achieving up to $\sim$ 45× faster execution times while maintaining methodological consistency (see \autoref{speedup}). The overall performance gains across all re-implementations position `MetaCommunityMetrics.jl` as a versatile solution for diverse ecological research needs, from standard biodiversity analyses to simulation studies and large-scale data processing.

Beyond performance improvements, we extended pre-existing `R` implementations to handle spatiotemporal datasets directly. We added temporal, spatial, and community-level aggregation methods to metrics originally restricted to a single time point, site, and species pair, eliminating extensive data preprocessing. For instance, beta diversity decompositions in `R` operate on a single metacommunity (consisting of all local communities across all sites at a single time point), while niche hypervolume measurements calculate volumes for an individual species and dissimilarity between a single species pair. We extended these with `spatial_beta_div()` and `temporal_beta_div()`, which calculate beta diversity decompositions across time and space, respectively. Similarly, `average_MVNH_det()` calculates average niche volume across all species, and `average_MVNH_dissimilarity()` computes average niche dissimilarity across all species pairs. For DNCI calculations, we developed `create_groups()` and `plot_groups()` to meet the site grouping requirements suggested by Vilmi et al. (2021) and to visualize the resulting groupings.

## Additional handling of edge cases when calculating DNCI
We extended our DNCI implementation to handle edge cases common in simulated community datasets (e.g., single-species communities, insufficient permutation variation) by returning status flags that identify when the standard DNCI calculation is not possible. See [documentation](https://cralibe.github.io/MetaCommunityMetrics.jl/DNCI/) for details on all five edge cases.

![A plot showing the speedup of all benchmarked functions across three [datasets](https://cralibe.github.io/MetaCommunityMetrics.jl/Benchmarking/#Benchmarked-using-Large-Dataset) (small, medium and large). Speedup is calculated as the `R` median execution time divided by the `Julia` median execution time. Median speedup and its confidence interval of each comparison is shown. The grey dashed line shows where speedup equals to 1, that is both `R` and `Julia`'s implementations require the same amount of time. Details about this benchmarking can be found in the documentation under the [benchmark results section](https://cralibe.github.io/MetaCommunityMetrics.jl/Benchmarking/).\label{speedup}](speedup.pdf){width=100%}

## Validation 
All re-implementations were validated against their R equivalents (see [documentation](https://cralibe.github.io/MetaCommunityMetrics.jl/Validation/)).

# Examples
## Using rodent metacommunity data as the sample data
We demonstrate the package using rodent data from the Portal Project [@ernest2018portal], a long-term Chihuahuan desert ecosystem study. The sample dataset includes 21 species across 24 sites from 2010-2023 (117 sampling events), with simulated spatial coordinates and environmental variables (temperature, precipitation). The data can be loaded with:
```julia
julia> load_sample_data()
```
The complete output of `load_sample_data()` is available in the [documentation](https://cralibe.github.io/MetaCommunityMetrics.jl/#Example).

## Illustration using the sample data
Here, we demonstrate a typical community analysis workflow: (1) assessing community stability using variability metrics (`CV_meta()`), (2) decomposing beta diversity to understand spatial and temporal turnover (`spatial_beta_div()`, `temporal_beta_div()`), (3) inferring assembly processes via DNCI analysis (`create_groups()`, `plot_groups()`, `DNCI_multigroup()`), and (4) quantifying niche overlap patterns (`niche_overlap()`, `average_MVNH_det()`, `average_MVNH_dissimilarity()`).

### Assessing temporal stability
We first examine how abundance variability differs across organizational levels (species vs. community) and spatial scales (local vs. regional) using `CV_meta()`:
```julia
julia> using MetaCommunityMetrics, Pipe, DataFrames, Random
julia> df = load_sample_data()
julia> CV_meta(df.Abundance, df.Sampling_date_order, df.plot, 
         df.Species)
1×4 DataFrame
 Row  CV_s_l   CV_s_r    CV_c_l    CV_c_r   
      Float64  Float64   Float64   Float64  

   1  1.48859  0.944937  0.718266  0.580183
```
Variability in abundance decreases from local to regional scales (comparing `CV_s_l` to `CV_s_r` and `CV_c_l` to `CV_c_r`) and from species to community levels (comparing `CV_s_l` to `CV_c_l` and `CV_s_r` to `CV_c_r`), indicating greater stability at broader spatial scales and higher organizational levels.

### Decomposing beta diversity
Next, we quantify compositional turnover and determine whether changes are driven by species replacement or abundance differences using `spatial_beta_div()` and 
`temporal_beta_div()`. The beta diversity decompositions in space:
```julia
julia> spatial_beta_div(df.Abundance, df.Sampling_date_order, 
         df.plot, df.Species; quant = true)
1×3 DataFrame
 Row  spatial_BDtotal  spatial_Repl  spatial_RichDif 
      Float64          Float64       Float64         

   1         0.264822      0.121882         0.142939
```
The beta diversity decompositions in time:
```julia
julia> temporal_beta_div(df.Abundance, df.Sampling_date_order, 
         df.plot, df.Species; quant = true)
1×3 DataFrame
 Row  temporal_BDtotal  temporal_Repl  temporal_RichDif 
      Float64           Float64        Float64          

   1          0.311222      0.0995483          0.211674
```
Spatial composition differences are driven equally by replacement and abundance differences (values are similar), while temporal differences are primarily driven by abundance changes (higher value), suggesting seasonal patterns or disturbance over time [@podani2011new; @legendre2014interpreting].

### Inferring assembly processes
To determine whether communities are predominantly structured by dispersal or environmental filtering, we calculate DNCI using `DNCI_multigroup()`, which quantifies these relative contributions by comparing species-specific contributions to between-group dissimilarity in observed versus null model data. This function automatically filters absent or ubiquitous species and handles empty sites. Before calculating DNCI, we first group sites using `create_groups()` to ensure at least two groups with a minimum of five sites per group, and that variation in the number of species and sites per group does not exceed 40% and 30%, respectively:
```julia
julia> grouping_result = create_groups(df.Sampling_date_order,
         df.Latitude, df.Longitude, df.plot, df.Species, df.Presence)
```
`grouping_result` returns a dictionary indexed by `Sampling_date_order`. Time points not meeting grouping requirements are assigned `missing`. We can access the result at `Sampling_date_order` = 60 with `grouping_result[60]`.

Now, we will use `plot_groups()` to visualize the grouping at `Sampling_date_order()` = 60:
```julia
julia> plot_groups(grouping_result[60].Latitude, 
          grouping_result[60].Longitude, 
          grouping_result[60].Group; output_file = "groups.svg")
```

![A plot showing the location of sites and their groups. Different colors represent different groups at `Sampling_date_order` $=$ 60.\label{groups}](groups.pdf){width=100%}

\autoref{groups} shows the grouping result at `Sampling_date_order` $=$ 60.

We can then join `df` with `grouping_result[60]` to obtain the matrix and group assignment column for calculating DNCI at `Sampling_data_order` $=$ 60 using `DNCI_multigroup()`.
```julia
julia> group_df = @pipe df |>
                  filter(row -> row[:Sampling_date_order] == 60, _) |>
                  select(_, [:plot, :Species, :Presence]) |>
                  innerjoin(_, grouping_result[60], on = [:plot => :Site, 
                  :Species], makeunique = true)|>
                  select(_, [:plot, :Species, :Presence, :Group]) |>
                  unstack(_, :Species, :Presence, fill = 0)
julia> comm= @pipe group_df |>
             select(_, Not([:plot,:Group])) |>
             Matrix(_)       
julia> Random.seed!(1234) 
julia> DNCI_result = DNCI_multigroup(comm, group_df.Group, 1000; 
          Nperm_count = false)
6×6 DataFrame
 Row  group1  group2  DNCI      CI_DNCI  S_DNCI    status 
      Int64   Int64   Float64   Float64  Float64   String 

   1       1       2  -3.41127  2.17348  1.08674   normal
   2       1       3  -2.44866  2.05951  1.02976   normal
   3       1       4  -2.3671   2.45697  1.22848   normal
   4       2       3  -2.65022  2.28931  1.14466   normal
   5       2       4  -3.0168   2.43496  1.21748   normal
   6       3       4  -1.83521  1.9589   0.979449  normal
```
The DNCI values at `Sampling_date_order` $=$ 60 do not differ significantly from zero for group pairs 1-4 and 3-4, while the rest are significantly smaller than zero. This indicates that most local communities in this metacommunity are dominated by dispersal processes, while both dispersal and niche processes contribute evenly to some local communities at `Sampling_date_order` $=$ 60. To draw conclusions about the entire metacommunity over time, we suggest running this function for all available time points and averaging across all group pairs at all time points.

### Quantifying niche overlap
Finally, we examine how species partition niche space by using  `niche_overlap()`, `prop_patches()`, `average_MVNH_det()`, and `average_MVNH_dissimilarity()`. We will start with `niche_overlap()`:
```julia
julia> niche_overlap(df.Abundance, df.Species, df.plot, 
          df.Sampling_date_order)
1×3 DataFrame
 Row  mean_niche_overlap_index  min_niche_overlap_index  max_niche_overlap_index 
      Float64                   Float64                  Float64                 

   1                 0.0923816                      0.0                 0.406837
```
Then, we will use `prop_patches()`:
```julia
julia> prop_patches(df.Presence, df.Species, df.plot)
1×3 DataFrame
 Row  mean_prop_patches  min_prop_patches  max_prop_patches 
      Float64            Float64           Float64          

   1           0.734649         0.0833333               1.0
```
We then calculate the average niche volume and average niche volume dissimilarity across all species using `average_MVNH_det()` and `average_MVNH_dissimilarity()`, respectively. Note that these functions, unlike the others, require environmental data, which is available in our sample dataset (i.e., temperature and precipitation). These environmental variables were standardized and checked for normality before being stored as our sample data. Note that users need to (1) standardized their environmental variables and transform them to normal distribution (if necessary); and (2) remove singletons before using these functions. The following demonstrates how we use these functions:
```julia
julia> df = @pipe load_sample_data() |>
                      groupby(_, :Species) |>
                      filter(row -> sum(row.Presence) > 1, _)|>
                      DataFrame(_)                      
julia> env_data = @pipe df |> 
          select(_, [:standardized_temperature, 
          :standardized_precipitation])       
julia> average_MVNH_det(env_data, df.Presence, df.Species; var_names = 
         ["Temperature", "Precipitation"])
1.2103765096417536
julia> average_MVNH_dissimilarity(env_data, df.Presence, df.Species;
         var_names = ["Temperature", "Precipitation"])     
0.03059942936454443
```
High spatial overlap (`mean_prop_patches` $\sim$ 0.73) combined with low temporal overlap (`mean_niche_overlap_index` $\sim$ 0.09) and high niche similarity (`average_MVNH_dissimilarity` $\sim$ 0.03) suggests temporal niche partitioning, where species occupy similar sites at different times [@albrecht2001spatial; @lear2021temporal].

# Future Directions 
We will apply an approximate Bayesian computation framework to infer community assembly processes from long-term phytoplankton monitoring data. This project will depend on `MetaCommunityMetrics.jl` to calculate summary statistics in its pipeline.

# Acknowledgement
We are grateful to the authors of the pre-existing `R` packages/implementations that we have implemented in `MetaCommunityMetrics.jl`. We express our gratitude to everyone who has tried out this package and provided feedback on how to improve it. LMG acknowledges funding from USC and Cornell University, YYJC acknowledges funding from the Wrigley Institute.  

# References
