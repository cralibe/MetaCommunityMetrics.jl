# Benchmarking Results
```@meta
CurrentModule = MetaCommunityMetrics
```
## Computational Resources
*All benchmarks were performed on the same machine to ensure consistent comparisons.*
- **CPU**: Apple M4
- **Number of Cores**: 10
- **Memory**: 16GB RAM
- **Operating System**: macOS Sequoia 15.6
- **Julia Version**: 1.10.10
- **R Version**: 4.4.2

## Benchmarking Methods
To assess the efficiency of `MetaCommunityMetrics` compared to equivalent `R` implementations, we benchmark our functions against their `R` counterparts, focusing on execution time and memory usage. The following tables summarize the benchmark results based on 100 samples each. 

We tested using datasets of three sizes:
- Large (53,352 observations)
- Medium (26,676 observations, 50% of the size of the large dataset) 
- Small (5,325 observations, 10% of the size of the large dataset)

The large dataset is the sample data included with `MetaCommunityMetrics`, accessible via `load_sample_data()`. The small and medium datasets can be accessed [`here`](https://github.com/cralibe/MetaCommunityMetrics.jl/tree/main/data/data_for_testing).

Each function was benchmarked using 100 samples in both `BenchmarkTools.jl` in `Julia` and `bench::mark()` in `R` to ensure robust statistical sampling. For memory usage comparisons:

- In `Julia`, we report the `memory estimate` from `BenchmarkTools.jl`, which measures bytes allocated for a given expression per sample.

- In `R`, we report the `mem_alloc` metric from `bench::mark()`, which tracks R heap allocation per sample.

According to documentation, the `Julia` metric measures total memory allocation during execution, while the `R` metric specifically tracks heap allocations within the R runtime, excluding "memory allocated outside the R heap, e.g., by `malloc()` or `new` directly." Due to differences in language implementation and measurement methodology, direct numerical comparisons between languages should be interpreted with caution.

## Speedup Summary
*Below is a plot showing the speedup of all benchmarked functions across the three datasets (small, medium and large). Speedup is calculated as the `R` median execution time divided by the `Julia` median execution time.Median speedup and its confidence interval of each comparison is shown. The grey dashed line shows where speedup equals to 1, that is both `R` and `Julia`'s implementations require the same amount of time. The calcualation of beta diversity has two variants: `quant = true` (abundance data is used) and `quant = false`(occurence data is used).*
![Speedup Plot](assets/speedup.png)


## Benchmarking Results in Details
All times are in millisecond (ms), and memory is in mebibytes (MiB). All values are rounded up to 4 decimal places.

### Median Execution Time and Speedup Values
*Bold text indicates the test cases with maximum and minimum speedup values. 95% confidence interval of the speedup is reported.*

| Test Case | Data Size | `Julia` | `R` | Speedup | Lower CI | Upper CI |
| --- | --- | --- | --- | --- | --- | --- |
| Beta Diversity (Abundance, quant=true) | Large | 0.1353 | 2.4174 | 17.8627 | 16.1499 | 19.5014 |
| Beta Diversity (Abundance, quant=true) | Medium | 0.0781 | 1.3655 | 17.4879 | 16.4085 | 18.3717 |
| Beta Diversity (Abundance, quant=true) | Small | 0.0407 | 1.1525 | 28.3121 | 27.2751 | 30.5852 |
| Beta Diversity (Abundance, quant=false) | Large | 0.0169 | 0.2096 | 12.4227 | 11.7001 | 13.475 |
| Beta Diversity (Abundance, quant=false) | Medium | 0.0091 | 0.2756 | 30.2681 | 27.6289 | 32.9273 |
| Beta Diversity (Abundance, quant=false) | Small | 0.0063 | 0.3219 | 50.9991 | 47.7958 | 54.1693 |
| Beta Diversity (Presence, quant=false) | Large | 0.0132 | 0.2101 | 15.929 | 13.0486 | 17.7571 |
| Beta Diversity (Presence, quant=false) | Medium | 0.0112 | 0.309 | 27.47 | 25.373 | 29.6856 |
| **Beta Diversity (Presence, quant=false)** | **Small** | **0.0057** | **0.33** | **58.2335** | **55.2126** | **65.2326** |
| Spatial Beta Diversity (Abundance, quant=true) | Large | 2.1147 | 9.267 | 4.3821 | 4.306 | 4.4893 |
| Spatial Beta Diversity (Abundance, quant=true) | Medium | 1.9895 | 9.2043 | 4.6263 | 4.471 | 4.7898 |
| Spatial Beta Diversity (Abundance, quant=true) | Small | 1.1195 | 8.8448 | 7.9007 | 7.7695 | 8.0108 |
| Spatial Beta Diversity (Abundance, quant=false) | Large | 1.8822 | 6.9715 | 3.7039 | 3.6467 | 3.7513 |
| Spatial Beta Diversity (Abundance, quant=false) | Medium | 1.6694 | 7.1675 | 4.2934 | 4.2002 | 4.399 |
| Spatial Beta Diversity (Abundance, quant=false) | Small | 0.9695 | 6.7469 | 6.9593 | 6.8101 | 7.073 |
| Spatial Beta Diversity (Presence, quant=false) | Large | 1.9022 | 6.9638 | 3.6608 | 3.6186 | 3.73 |
| Spatial Beta Diversity (Presence, quant=false) | Medium | 1.6716 | 7.2182 | 4.3183 | 4.1722 | 4.4475 |
| Spatial Beta Diversity (Presence, quant=false) | Small | 0.9669 | 6.6972 | 6.9262 | 6.8451 | 7.0289 |
| Temporal Beta Diversity (Abundance, quant=true) | Large | 5.8652 | 53.3016 | 9.0878 | 8.7875 | 9.2049 |
| Temporal Beta Diversity (Abundance, quant=true) | Medium | 5.215 | 54.3208 | 10.4162 | 10.0311 | 10.6263 |
| Temporal Beta Diversity (Abundance, quant=true) | Small | 4.5515 | 52.5342 | 11.5422 | 11.294 | 11.7076 |
| Temporal Beta Diversity (Abundance, quant=false) | Large | 2.6689 | 9.769 | 3.6603 | 3.6154 | 3.7156 |
| Temporal Beta Diversity (Abundance, quant=false) | Medium | 2.2572 | 9.672 | 4.285 | 4.1534 | 4.3805 |
| Temporal Beta Diversity (Abundance, quant=false) | Small | 1.5722 | 9.0186 | 5.7364 | 5.5521 | 5.9203 |
| **Temporal Beta Diversity (Presence, quant=false)** | **Large** | **2.7422** | **9.3899** | **3.4242** | **3.3593** | **3.4753** |
| Temporal Beta Diversity (Presence, quant=false) | Medium | 2.3177 | 9.6286 | 4.1543 | 4.0458 | 4.3004 |
| Temporal Beta Diversity (Presence, quant=false) | Small | 1.4952 | 9.0206 | 6.0328 | 5.8424 | 6.135 |
| Dispersal-niche continuum index | Large | 483.221 | 12894.3232 | 26.6841 | 26.4701 | 26.923 |
| Dispersal-niche continuum index | Medium | 391.8913 | 12608.5208 | 32.1735 | 32.0287 | 32.3104 |
| Dispersal-niche continuum index | Small | 102.0653 | 3376.0879 | 33.0777 | 32.8155 | 33.4182 |
| Occupied Patches Proportion | Large | 0.6824 | 8.5893 | 12.5866 | 12.087 | 13.1502 |
| Occupied Patches Proportion | Medium | 0.5913 | 8.3506 | 14.1222 | 13.5761 | 14.9307 |
| Occupied Patches Proportion | Small | 0.4209 | 7.9086 | 18.79 | 16.3675 | 20.2831 |
| Variability Metrics | Large | 13.4022 | 103.596 | 7.7298 | 7.5992 | 9.6017 |
| Variability Metrics | Medium | 7.8724 | 52.2772 | 6.6405 | 6.5977 | 6.7956 |
| Variability Metrics | Small | 2.8329 | 14.0789 | 4.9698 | 4.9394 | 4.9936 |
| Hypervolume Estimation | Large | 0.0042 | 0.0307 | 7.3646 | 7.1158 | 7.5291 |
| Hypervolume Estimation | Medium | 0.0035 | 0.0264 | 7.4546 | 7.3328 | 7.544 |
| Hypervolume Estimation | Small | 0.0034 | 0.0266 | 7.8841 | 7.7032 | 7.9327 |
| Hypervolume Dissimilarity | Large | 0.0087 | 0.1395 | 16.0968 | 15.7783 | 16.5356 |
| Hypervolume Dissimilarity | Medium | 0.0082 | 0.1109 | 13.5755 | 13.4818 | 13.665 |
| Hypervolume Dissimilarity | Small | 0.008 | 0.1125 | 13.9913 | 13.9097 | 14.1604 |

### Memory Usage
#### Benchmarked using Large Dataset
*Bold text indicates the test case with the biggest memory usage difference between `Julia` and `R`.*

| Test Case | `Julia` | `R` |
|----------|------------------|---------------|
| Beta Diversity (Abundance, quant=true) | 0.4341 | 0.3757 |
| Beta Diversity (Abundance, quant=false) | 0.1347 | 0.1252 |
| Beta Diversity (Presence, quant=false) | 0.1347 | 0.1252 |
| Spatial Beta Diversity (Abundance, quant=true) | 3.9178 | 4.1093 |
| Spatial Beta Diversity (Abundance, quant=false) | 3.5170 | 2.6637 |
| Spatial Beta Diversity (Presence, quant=false) | 3.5170 | 2.6637 |
| Temporal Beta Diversity (Abundance, quant=true) | 16.8830 | 16.8787 |
| Temporal Beta Diversity (Abundance, quant=false) | 5.6584 | 5.1676 |
| Temporal Beta Diversity (Presence, quant=false) | 5.6584 | 5.1676 |
| **Dispersal-niche continuum index** | **782.0591** | **77.7877** |
| Occupied Patches Proportion | 1.9095 | 1.9928 |
| Variability Metrics | 12.4573 | 60.2264 |
| Hypervolume Estimation | 0.0122 | 0.0022 |
| Hypervolume Dissimilarity | 0.0168 | 0.0145 |

#### Benchmarked using Medium Dataset
*The test case with the biggest memory usage difference between `Julia` and `R` is highligthed.*

| Test Case | `Julia` | `R` |
|----------|------------------|---------------|
| Beta Diversity (Abundance, quant=true) | 0.2486 | 0.0357 |
| Beta Diversity (Abundance, quant=false) | 0.1086 | 0.0798 |
| Beta Diversity (Presence, quant=false) | 0.1086 | 0.0798 |
| Spatial Beta Diversity (Abundance, quant=true) | 2.3889 | 2.2737 |
| Spatial Beta Diversity (Abundance, quant=false) | 1.9908 | 1.8307 |
| Spatial Beta Diversity (Presence, quant=false) | 1.9908 | 1.8307 |
| Temporal Beta Diversity (Abundance, quant=true) | 15.2761 | 16.2510 |
| Temporal Beta Diversity (Abundance, quant=false) | 4.1321 | 4.5399 |
| Temporal Beta Diversity (Presence, quant=false) | 4.1321 | 4.5399 |
| **Dispersal-niche continuum index** | **543.0698** | **59.2854** |
| Occupied Patches Proportion | 0.9938 | 1.3995 |
| Variability Metrics | 7.6746 | 32.5536 |
| Hypervolume Estimation | 0.0085 | 0.0011 |
| Hypervolume Dissimilarity | 0.0120 | 0.0077 |

#### Benchmarked using Small Dataset
*The test case with the biggest memory usage difference between `Julia` and `R` is highligthed.*

| Test Case | `Julia` | `R` |
|----------|------------------|---------------|
| Beta Diversity (Abundance, quant=true) | 0.1220 | 0.0195 |
| Beta Diversity (Abundance, quant=false) | 0.0883 | 0.0444 |
| Beta Diversity (Presence, quant=false) | 0.0883 | 0.0444 |
| Spatial Beta Diversity (Abundance, quant=true) | 1.1326 | 1.2127 |
| Spatial Beta Diversity (Abundance, quant=false) | 0.7660 | 0.7697 |
| Spatial Beta Diversity (Presence, quant=false) | 0.7660 | 0.7697 |
| Temporal Beta Diversity (Abundance, quant=true) | 12.8803 | 15.4342 |
| Temporal Beta Diversity (Abundance, quant=false) | 2.8934 | 3.7231 |
| Temporal Beta Diversity (Presence, quant=false) | 2.8934 | 3.7231 |
| **Dispersal-niche continuum index** | **192.6010** | **11.1107** |
| Occupied Patches Proportion | 0.2610 | 0.4932 |
| Variability Metrics | 3.8483 | 10.5805 |
| Hypervolume Estimation** | 0.0059 | 0.0003 |
| Hypervolume Dissimilarity | 0.0082 | 0.0014 |

## Datesets used for this benchmark
### Large Dataset
```@jildoctest
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  standardized_temperature  standardized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0                0.829467              -1.4024
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5               -1.12294               -0.0519895
     3 │  2010      1     16                    1      4  BA               0         0      35.0     -108.5               -0.409808              -0.803663
     4 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               -1.35913               -0.646369
     5 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0                0.0822                 1.09485
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 53348 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               -0.571565              -0.836345
 53349 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               -2.33729               -0.398522
 53350 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5                0.547169               1.03257
 53351 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               -0.815015               0.95971
 53352 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0                0.48949               -1.59416
                                                                                                                                            53342 rows omitted
```
### Medium Dataset
```@jildoctest
26676×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  standardized_temperature  standardized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2015      7     18                   56     12  PM               0         0      35.5     -107.5              -1.30965                  0.967859
     2 │  2016      8      6                   66     13  SF               0         0      36.0     -110.0               1.45692                  1.77253
     3 │  2017      2     25                   71     21  SF               0         0      36.5     -109.0              -1.50086                  0.993311
     4 │  2018      5     19                   82     16  PB               0         0      36.0     -108.5              -1.2202                   0.684295
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 26673 │  2010     11      6                    8     13  BA               0         0      36.0     -110.0              -1.05336                 -0.250324
 26674 │  2013      9     14                   36     12  NA               0         0      35.5     -107.5               0.213222                 0.12
 26675 │  2023      2     18                  116     12  DS               0         0      35.5     -107.5              -0.217475                 0.042571
 26676 │  2014     11     22                   49     13  PF               0         0      36.0     -110.0               0.613491                -1.17076
                                                                                                                                            26668 rows omitted
```
### Small Dataset
```@jildoctest
5335×12 DataFrame
  Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  standardized_temperature  standardized_precipitation 
      │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
──────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    1 │  2015      7     18                   56     12  PM               0         0      35.5     -107.5               -1.30965                 0.967859
    2 │  2016      8      6                   66     13  SF               0         0      36.0     -110.0                1.45692                 1.77253
    3 │  2017      2     25                   71     21  SF               0         0      36.5     -109.0               -1.50086                 0.993311
    4 │  2018      5     19                   82     16  PB               0         0      36.0     -108.5               -1.2202                  0.684295
  ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 5332 │  2013     11      2                   37     13  RF               0         0      36.0     -110.0                2.07639                 2.72841
 5333 │  2018     11     10                   88      6  PH               0         0      35.0     -107.5               -0.197006                0.842547
 5334 │  2018      3     24                   80      2  DS               0         0      35.0     -109.5               -0.920093               -0.276074
 5335 │  2014      4     26                   42      7  PH               0         0      35.5     -110.0                0.848755               -0.247144
                                                                                                                                            5327 rows omitted
```

## Remarks
- For `DNCI_multigroup_result`, 100 permutations per sample are used in both the `Julia` and `R` implementation, and `parallelComputing` was set to be `TRUE` when benchmarking `DNCImper:::DNCI_multigroup()` in `R`.

## The Scripts Used for Benchmarking
- [`Julia`](https://github.com/cralibe/MetaCommunityMetrics.jl/blob/main/benchmarks/benchmark_julia.jl)
- [`R`](https://github.com/cralibe/MetaCommunityMetrics.jl/blob/main/benchmarks/benchmark_r/benchmark_r.R)

## Packages used for benchmarking
- [`bench`](https://github.com/r-lib/bench)
- [`BenchmarkTools.jl`](https://github.com/JuliaCI/BenchmarkTools.jl)