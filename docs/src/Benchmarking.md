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
*The test cases with the maximum and minimum speedup values are highligthed. 95% confidence interval of the speedup is reported.*

| TestCase                                         | DataSize | `Julia`      | `R`        | Speedup        | Lower CI | Upper CI |
|--------------------------------------------------|----------|--------------|------------|----------------|----------|----------|
| Beta Diversity (Abundance, quant=true)           | Large    | 0.1244       | 2.1471     | 17.2602        | 16.6071  | 17.7147  |
| Beta Diversity (Abundance, quant=true)           | Medium   | 0.0778       | 1.3732     | 17.6518        | 17.3768  | 18.1861  |
| Beta Diversity (Abundance, quant=true)           | Small    | 0.0394       | 0.9179     | 23.3000        | 22.0561  | 25.3414  |
| Beta Diversity (Abundance, quant=false)          | Large    | 0.0146       | 0.1997     | 13.6934        | 13.2358  | 16.8835  |
| Beta Diversity (Abundance, quant=false)          | Medium   | 0.0104       | 0.2138     | 20.5226        | 19.6683  | 22.2222  |
| ==Beta Diversity (Abundance, quant=false)==          | Small    | 0.0054       | 0.2414     | 44.5725        | 40.8472  | 47.8959  |
| Beta Diversity (Presence, quant=false)           | Large    | 0.0139       | 0.1962     | 14.1653        | 13.5828  | 14.7816  |
| Beta Diversity (Presence, quant=false)           | Medium   | 0.0113       | 0.1994     | 17.5968        | 17.1127  | 18.1192  |
| Beta Diversity (Presence, quant=false)           | Small    | 0.0057       | 0.2513     | 44.1863        | 39.5111  | 46.2503  |
| Spatial Beta Diversity (Abundance, quant=true)   | Large    | 2.1688       | 9.1234     | 4.2068         | 4.1387   | 4.3877   |
| Spatial Beta Diversity (Abundance, quant=true)   | Medium   | 1.4725       | 9.5022     | 6.4529         | 6.0358   | 6.5497   |
| Spatial Beta Diversity (Abundance, quant=true)   | Small    | 1.5839       | 9.4987     | 5.9971         | 5.7505   | 6.6586   |
| Spatial Beta Diversity (Abundance, quant=false)  | Large    | 1.9444       | 7.3179     | 3.7636         | 3.6124   | 3.8332   |
| Spatial Beta Diversity (Abundance, quant=false)  | Medium   | 1.6618       | 7.1957     | 4.3301         | 4.2245   | 4.4562   |
| Spatial Beta Diversity (Abundance, quant=false)  | Small    | 1.5207       | 6.9109     | 4.5447         | 4.3966   | 4.8593   |
| Spatial Beta Diversity (Presence, quant=false)   | Large    | 2.0419       | 8.0388     | 3.9370         | 3.8449   | 4.0554   |
| Spatial Beta Diversity (Presence, quant=false)   | Medium   | 1.7067       | 7.3284     | 4.2939         | 4.1892   | 4.4127   |
| Spatial Beta Diversity (Presence, quant=false)   | Small    | 1.4723       | 7.0079     | 4.7599         | 4.5737   | 5.1147   |
| Temporal Beta Diversity (Abundance, quant=true)  | Large    | 5.3096       | 53.9455    | 10.1600        | 10.0780  | 10.2358  |
| Temporal Beta Diversity (Abundance, quant=true)  | Medium   | 5.3015       | 53.3531    | 10.0638        | 9.8600   | 10.2946  |
| Temporal Beta Diversity (Abundance, quant=true)  | Small    | 4.7431       | 53.3732    | 11.2528        | 10.7984  | 11.5056  |
| ==Temporal Beta Diversity (Abundance, quant=false)== | Large    | 2.6548       | 9.5602     | 3.6012         | 3.5442   | 3.6613   |
| Temporal Beta Diversity (Abundance, quant=false) | Medium   | 2.2395       | 9.4459     | 4.2179         | 4.0805   | 4.3567   |
| Temporal Beta Diversity (Abundance, quant=false) | Small    | 2.2169       | 9.2761     | 4.1843         | 4.0034   | 4.3031   |
| Temporal Beta Diversity (Presence, quant=false)  | Large    | 2.6877       | 10.7795    | 4.0106         | 3.8444   | 4.1896   |
| Temporal Beta Diversity (Presence, quant=false)  | Medium   | 2.2363       | 9.7801     | 4.3734         | 4.2446   | 4.5393   |
| Temporal Beta Diversity (Presence, quant=false)  | Small    | 2.2682       | 9.2270     | 4.0680         | 3.9668   | 4.3038   |
| Dispersal-niche continuum index                  | Large    | 588.7658     | 12574.8461 | 21.3580        | 21.2157  | 21.4509  |
| Dispersal-niche continuum index                  | Medium   | 505.2023     | 12268.2182 | 24.2838        | 24.1870  | 24.3857  |
| Dispersal-niche continuum index                  | Small    | 122.3677     | 3140.2943  | 25.6628        | 25.4323  | 25.8697  |
| Occupied Patches Proportion                      | Large    | 0.8469       | 7.2455     | 8.5554         | 8.2737   | 8.8168   |
| Occupied Patches Proportion                      | Medium   | 0.6504       | 7.3115     | 11.2416        | 10.8683  | 11.8014  |
| Occupied Patches Proportion                      | Small    | 0.3657       | 7.0934     | 19.3974        | 18.6676  | 20.1553  |
| Variability Metrics                              | Large    | 13.8979      | 98.9990    | 7.1233         | 7.0156   | 7.1537   |
| Variability Metrics                              | Medium   | 7.9131       | 52.7903    | 6.6713         | 6.6224   | 6.7430   |
| Variability Metrics                              | Small    | 2.7518       | 14.2607    | 5.1824         | 5.1289   | 5.2272   |
| Hypervolume Estimation                           | Large    | 0.0035       | 0.0328     | 9.3831         | 8.4226   | 18.6374  |
| Hypervolume Estimation                           | Medium   | 0.0032       | 0.0276     | 8.7127         | 8.5896   | 8.8098   |
| Hypervolume Estimation                           | Small    | 0.0029       | 0.0264     | 9.2643         | 9.0658   | 9.4181   |
| Hypervolume Dissimilarity                        | Large    | 0.0064       | 0.1147     | 17.9885        | 17.6258  | 18.2162  |
| Hypervolume Dissimilarity                        | Medium   | 0.0056       | 0.1188     | 21.3605        | 21.1646  | 21.6136  |
| Hypervolume Dissimilarity                        | Small    | 0.0050       | 0.1852     | 37.3449        | 34.6406  | 38.0277  |


### Memory Usage
#### Benchmarked using Large Dataset
*The test case with the biggest memory usage difference between `Julia` and `R` is highligthed.*

| TestCase                                         | `Julia`      | `R`      |
|--------------------------------------------------|--------------|----------|
| Beta Diversity (Abundance, quant=true)           | 0.4346       | 0.0566   |
| Beta Diversity (Abundance, quant=false)          | 0.1347       | 0.1252   |
| Beta Diversity (Presence, quant=false)           | 0.1347       | 0.1252   |
| Spatial Beta Diversity (Abundance, quant=true)   | 3.9185       | 3.1142   |
| Spatial Beta Diversity (Abundance, quant=false)  | 3.5173       | 2.6712   |
| Spatial Beta Diversity (Presence, quant=false)   | 3.5173       | 2.6712   |
| Temporal Beta Diversity (Abundance, quant=true)  | 16.8838      | 16.8863  |
| Temporal Beta Diversity (Abundance, quant=false) | 5.6586       | 5.1752   |
| Temporal Beta Diversity (Presence, quant=false)  | 5.6586       | 5.1752   |
| ==Dispersal-niche continuum index==              | 781.9220     | 71.9833  |
| Occupied Patches Proportion                      | 1.9098       | 1.8857   |
| Variability Metrics                              | 12.4573      | 60.2217  |
| Hypervolume Estimation                           | 0.0122       | 0.0022   |
| Hypervolume Dissimilarity                        | 0.0168       | 0.0145   |

#### Benchmarked using Medium Dataset
*The test case with the biggest memory usage difference between `Julia` and `R` is highligthed.*

| TestCase                                         | `Julia`      | `R`      |
|--------------------------------------------------|--------------|----------|
| Beta Diversity (Abundance, quant=true)           | 0.2491       | 0.0357   |
| Beta Diversity (Abundance, quant=false)          | 0.1086       | 0.0798   |
| Beta Diversity (Presence, quant=false)           | 0.1086       | 0.0798   |
| Spatial Beta Diversity (Abundance, quant=true)   | 2.3901       | 2.2736   |
| Spatial Beta Diversity (Abundance, quant=false)  | 1.9915       | 1.8306   |
| Spatial Beta Diversity (Presence, quant=false)   | 1.9915       | 1.8306   |
| Temporal Beta Diversity (Abundance, quant=true)  | 15.2773      | 16.2509  |
| Temporal Beta Diversity (Abundance, quant=false) | 4.1328       | 4.5397   |
| Temporal Beta Diversity (Presence, quant=false)  | 4.1328       | 4.5397   |
| ==Dispersal-niche continuum index==              | 542.9327     | 59.2854  |
| Occupied Patches Proportion                      | 0.9940       | 1.3994   |
| Variability Metrics                              | 7.6746       | 32.5536  |
| Hypervolume Estimation                           | 0.0085       | 0.0011   |
| Hypervolume Dissimilarity                        | 0.0120       | 0.0077   |

#### Benchmarked using Small Dataset
| TestCase                                         | `Julia`      | `R`      |
|--------------------------------------------------|--------------|----------|
| Beta Diversity (Abundance, quant=true)           | 0.1225       | 0.0195   |
| Beta Diversity (Abundance, quant=false)          | 0.0883       | 0.0444   |
| Beta Diversity (Presence, quant=false)           | 0.0883       | 0.0444   |
| Spatial Beta Diversity (Abundance, quant=true)   | 1.1334       | 1.2126   |
| Spatial Beta Diversity (Abundance, quant=false)  | 0.7663       | 0.7696   |
| Spatial Beta Diversity (Presence, quant=false)   | 0.7663       | 0.7696   |
| Temporal Beta Diversity (Abundance, quant=true)  | 12.8811      | 15.4341  |
| Temporal Beta Diversity (Abundance, quant=false) | 2.8937       | 3.7230   |
| Temporal Beta Diversity (Presence, quant=false)  | 2.8937       | 3.7230   |
| ==Dispersal-niche continuum index==              | 192.5782     | 11.1101  |
| Occupied Patches Proportion                      | 0.2612       | 0.4931   |
| Variability Metrics                              | 3.8483       | 10.5805  |
| Hypervolume Estimation                           | 0.0059       | 0.0003   |
| Hypervolume Dissimilarity                        | 0.0082       | 0.0014   |

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