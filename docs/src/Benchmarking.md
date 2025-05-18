# Benchmarking Results
```@meta
CurrentModule = MetaCommunityMetrics
```
## Computational Resources
*All benchmarks were performed on the same machine to ensure consistent comparisons.*
- **CPU**: Apple M4
- **Number of Cores**: 10
- **Memory**: 16GB RAM
- **Operating System**: macOS Sequoia 15.5
- **Julia Version**: 1.10.9
- **R Version**: 4.4.2

## Bechmarking Methods
To assess the efficiency of `MetaCommunityMetrics`compared to equivalent `R` implementations, we benchmark our functions against equivalent implementations in `R`, focusing on comparing the time and memory usage. The following tables summarized the benchmarking results. We tested using datasets of three sizes: small (5325 observations), medium (26676 observations), and large (53352 observations). The larger dataset is the sample data of `MetaCommunityMetrics`, which can be accessed by using `load_sample_data()`. The small and medium dataset can be accessed [`here`](https://github.com/cralibe/MetaCommunityMetrics.jl/tree/main/data/data_for_testing). Each function was benchmarked using 100 samples in both `BenchmarkTools.jl` in `Julia` and `bench::mark()` in `R` to ensure robust statistical sampling. For memory usage comparisons, we report Julia's `memory estimate` from `BenchmarkTools.jl`. According to the `BenchmarkTools.j` documentation, this measurement returns the number of bytes allocated when executing a given expression, corresponding to the trial with the minimum elapsed time measured during benchmarking. For `R`, we report the `mem_alloc` metric from `bench::mark()`, which specifically tracks R heap allocations, with the documentation noting that it excludes 'memory allocated outside the R heap, e.g., by `malloc()` or `new` directly. Both metrics measure heap memory allocations within their respective language runtimes. Due to differences in language implementation and measurement methodology, direct numerical comparisons between languages should be interpreted with caution. 

### Benchmarking using the large dataset
*All times are in millisecond (ms), and memory is in mebibytes (MiB). All values are rounded up to 4 decimal places.*

#### Julia
| TestCase                                         | Time_minimum | Time_median | Time_mean | Time_maximum | Time_std | memory   |
|--------------------------------------------------|--------------|-------------|-----------|--------------|----------|----------|
| Beta Diversity (Abundance, quant=true)           | 0.1047       | 0.1359      | 0.2862    | 2.3685       | 0.3502   | 0.4346   |
| Beta Diversity (Abundance, quant=false)          | 0.0108       | 0.0170      | 0.1825    | 8.1686       | 1.0067   | 0.1347   |
| Beta Diversity (Presence, quant=false)           | 0.0102       | 0.0115      | 0.1738    | 10.0957      | 1.0664   | 0.1347   |
| Spatial Beta Diversity (Abundance, quant=true)   | 1.0407       | 1.1159      | 1.2525    | 4.1448       | 0.5467   | 3.9185   |
| Spatial Beta Diversity (Abundance, quant=false)  | 0.9028       | 1.0091      | 1.0973    | 3.7046       | 0.4355   | 3.5173   |
| Spatial Beta Diversity (Presence, quant=false)   | 0.8751       | 1.0448      | 1.1481    | 3.7535       | 0.4430   | 3.5173   |
| Temporal Beta Diversity (Abundance, quant=true)  | 4.0484       | 4.4040      | 5.1821    | 15.0618      | 1.6165   | 16.8838  |
| Temporal Beta Diversity (Abundance, quant=false) | 1.6742       | 1.8266      | 2.1790    | 9.2889       | 1.1900   | 5.6586   |
| Temporal Beta Diversity (Presence, quant=false)  | 1.6390       | 1.8623      | 2.0460    | 4.4411       | 0.5831   | 5.6586   |
| Dispersal-niche continuum index                  | 143.1201     | 147.7740    | 149.8734  | 206.0497     | 7.7826   | 408.4675 |
| Occupied Patches Proportion                      | 0.7545       | 0.9106      | 1.0019    | 4.1220       | 0.4708   | 2.5230   |
| Variability Metrics                              | 11.9627      | 12.5781     | 12.8409   | 15.1580      | 0.7824   | 12.4573  |
| Hypervolume Estimation                           | 0.0033       | 0.0038      | 0.0043    | 0.0537       | 0.0050   | 0.0122   |
| Hypervolume Dissimilarity                        | 0.0057       | 0.0064      | 0.0075    | 0.0872       | 0.0082   | 0.0168   |

#### R
| TestCase                                         | Time_minimum | Time_median | Time_mean  | Time_maximum | Time_std  | memory    |
|--------------------------------------------------|--------------|-------------|------------|--------------|-----------|-----------|
| Beta Diversity (Abundance, quant=true)           | 2.7006       | 4.4541      | 10.4123    | 541.1395     | 53.6438   | 0.0566    |
| Beta Diversity (Abundance, quant=false)          | 0.2130       | 0.2601      | 0.4028     | 2.3464       | 0.3122    | 0.1252    |
| Beta Diversity (Presence, quant=false)           | 0.2256       | 0.2670      | 0.3392     | 1.4041       | 0.1947    | 0.1252    |
| Spatial Beta Diversity (Abundance, quant=true)   | 15.1884      | 21.9025     | 22.7474    | 82.2930      | 7.8687    | 3.1113    |
| Spatial Beta Diversity (Abundance, quant=false)  | 10.8790      | 16.5718     | 20.5865    | 285.4170     | 27.7484   | 2.6683    |
| Spatial Beta Diversity (Presence, quant=false)   | 10.3769      | 14.6347     | 15.7660    | 54.5412      | 6.1221    | 2.6683    |
| Temporal Beta Diversity (Abundance, quant=true)  | 95.3134      | 114.8494    | 127.8791   | 513.7155     | 60.3508   | 16.8833   |
| Temporal Beta Diversity (Abundance, quant=false) | 14.2315      | 18.4020     | 22.2465    | 332.8670     | 31.5405   | 5.1722    |
| Temporal Beta Diversity (Presence, quant=false)  | 14.2623      | 18.3448     | 19.2711    | 43.8397      | 4.0817    | 5.1722    |
| Dispersal-niche continuum index                  | 16150.8938   | 23881.0763  | 23694.5085 | 33782.7078   | 2608.4650 | 8347.5709 |
| Occupied Patches Proportion                      | 11.4863      | 16.1346     | 16.6754    | 26.2392      | 3.2114    | 1.8827    |
| Variability Metrics                              | 223.6846     | 314.4044    | 413.6490   | 1656.9433    | 289.0475  | 60.0667   |
| Hypervolume Estimation                           | 0.0558       | 0.1042      | 0.1054     | 0.2251       | 0.0256    | 0.0022    |
| Hypervolume Dissimilarity                        | 0.1825       | 0.2540      | 0.3457     | 0.8766       | 0.1756    | 0.0145    |

#### Speedup
| TestCase                                         | `Time_median_julia` | `Time_median_r` | `Speedup` |
|--------------------------------------------------|---------------------|-----------------|-----------|
| Beta Diversity (Abundance, quant=true)           | 0.1359              | 4.4541          | 32.7757   |
| Beta Diversity (Abundance, quant=false)          | 0.0170              | 0.2601          | 15.3026   |
| Beta Diversity (Presence, quant=false)           | 0.0115              | 0.2670          | 23.1789   |
| Spatial Beta Diversity (Abundance, quant=true)   | 1.1159              | 21.9025         | 19.6281   |
| Spatial Beta Diversity (Abundance, quant=false)  | 1.0091              | 16.5718         | 16.4223   |
| Spatial Beta Diversity (Presence, quant=false)   | 1.0448              | 14.6347         | 14.0068   |
| Temporal Beta Diversity (Abundance, quant=true)  | 4.4040              | 114.8494        | 26.0784   |
| Temporal Beta Diversity (Abundance, quant=false) | 1.8266              | 18.4020         | 10.0744   |
| Temporal Beta Diversity (Presence, quant=false)  | 1.8623              | 18.3448         | 9.8507    |
| Dispersal-niche continuum index                  | 147.7740            | 23881.0763      | 161.6054  |
| Occupied Patches Proportion                      | 0.9106              | 16.1346         | 17.7182   |
| Variability Metrics                              | 12.5781             | 314.4044        | 24.9962   |
| Hypervolume Estimation                           | 0.0038              | 0.1042          | 27.7980   |
| Hypervolume Dissimilarity                        | 0.0064              | 0.2540          | 39.7179   |

### Benchmarking using the medium dataset
*All times are in millisecond (ms), and memory is in mebibytes (MiB). All values are rounded up to 4 decimal places.*

#### Julia
| TestCase                                         | Time_minimum | Time_median | Time_mean | Time_maximum | Time_std | memory   |
|--------------------------------------------------|--------------|-------------|-----------|--------------|----------|----------|
| Beta Diversity (Abundance, quant=true)           | 0.0647       | 0.0853      | 0.1210    | 0.5125       | 0.0796   | 0.2491   |
| Beta Diversity (Abundance, quant=false)          | 0.0078       | 0.0160      | 0.0154    | 0.0355       | 0.0044   | 0.1086   |
| Beta Diversity (Presence, quant=false)           | 0.0067       | 0.0083      | 0.0090    | 0.0310       | 0.0030   | 0.1086   |
| Spatial Beta Diversity (Abundance, quant=true)   | 0.7606       | 0.9022      | 1.0497    | 4.5912       | 0.5313   | 2.3897   |
| Spatial Beta Diversity (Abundance, quant=false)  | 0.6500       | 0.7594      | 0.8366    | 2.8606       | 0.3347   | 1.9910   |
| Spatial Beta Diversity (Presence, quant=false)   | 0.6063       | 0.7095      | 0.8539    | 7.7896       | 0.8100   | 1.9910   |
| Temporal Beta Diversity (Abundance, quant=true)  | 3.7435       | 4.0110      | 4.4374    | 9.5260       | 0.9203   | 15.2769  |
| Temporal Beta Diversity (Abundance, quant=false) | 1.1183       | 1.3314      | 1.5423    | 6.5418       | 0.7577   | 4.1324   |
| Temporal Beta Diversity (Presence, quant=false)  | 1.0437       | 1.2106      | 1.3106    | 4.0138       | 0.4640   | 4.1324   |
| Dispersal-niche continuum index                  | 51.6221      | 53.4701     | 53.6385   | 57.3854      | 0.9387   | 143.7355 |
| Occupied Patches Proportion                      | 0.5801       | 0.6862      | 0.7330    | 1.7060       | 0.1706   | 1.3019   |
| Variability Metrics                              | 7.2410       | 7.4703      | 7.7076    | 10.3640      | 0.6414   | 7.6746   |
| Hypervolume Estimation                           | 0.0030       | 0.0034      | 0.0039    | 0.0504       | 0.0047   | 0.0085   |
| Hypervolume Dissimilarity                        | 0.0048       | 0.0054      | 0.0062    | 0.0773       | 0.0072   | 0.0120   |

#### R
| TestCase                                         | Time_minimum | Time_median | Time_mean | Time_maximum | Time_std  | memory    |
|--------------------------------------------------|--------------|-------------|-----------|--------------|-----------|-----------|
| Beta Diversity (Abundance, quant=true)           | 1.8467       | 2.7219      | 3.3254    | 11.5697      | 1.8266    | 0.0357    |
| Beta Diversity (Abundance, quant=false)          | 0.2459       | 0.3299      | 0.4395    | 2.0706       | 0.3433    | 0.0798    |
| Beta Diversity (Presence, quant=false)           | 0.2383       | 0.3195      | 0.4647    | 3.1369       | 0.4921    | 0.0798    |
| Spatial Beta Diversity (Abundance, quant=true)   | 14.7974      | 20.6098     | 29.2914   | 385.6328     | 50.7938   | 2.2707    |
| Spatial Beta Diversity (Abundance, quant=false)  | 9.6693       | 13.5177     | 15.0492   | 49.3984      | 5.8203    | 1.8277    |
| Spatial Beta Diversity (Presence, quant=false)   | 10.2190      | 14.7190     | 15.0275   | 28.9880      | 3.1351    | 1.8277    |
| Temporal Beta Diversity (Abundance, quant=true)  | 81.9918      | 114.5111    | 149.0064  | 1049.1010    | 169.9126  | 16.2479   |
| Temporal Beta Diversity (Abundance, quant=false) | 12.0927      | 17.6840     | 25.9879   | 672.2291     | 65.5459   | 4.5368    |
| Temporal Beta Diversity (Presence, quant=false)  | 14.0173      | 18.3145     | 18.5485   | 24.8850      | 2.0702    | 4.5368    |
| Dispersal-niche continuum index                  | 7890.5657    | 9162.7688   | 9536.9359 | 22655.8338   | 1805.3635 | 3573.2095 |
| Occupied Patches Proportion                      | 10.2833      | 14.9776     | 17.9419   | 89.8005      | 10.7927   | 1.3965    |
| Variability Metrics                              | 97.7768      | 141.9484    | 172.3243  | 744.6694     | 132.1003  | 32.5536   |
| Hypervolume Estimation                           | 0.0292       | 0.0402      | 0.0568    | 0.4308       | 0.0518    | 0.0011    |
| Hypervolume Dissimilarity                        | 0.1382       | 0.2190      | 0.4250    | 3.6017       | 0.5595    | 0.0077    |

#### Speedup
| TestCase                                         | `Time_median_julia` | `Time_median_r` | `Speedup` |
|--------------------------------------------------|---------------------|-----------------|-----------|
| Beta Diversity (Abundance, quant=true)           | 0.0853              | 2.7219          | 31.8970   |
| Beta Diversity (Abundance, quant=false)          | 0.0160              | 0.3299          | 20.5928   |
| Beta Diversity (Presence, quant=false)           | 0.0083              | 0.3195          | 38.3334   |
| Spatial Beta Diversity (Abundance, quant=true)   | 0.9022              | 20.6098         | 22.8437   |
| Spatial Beta Diversity (Abundance, quant=false)  | 0.7594              | 13.5177         | 17.8001   |
| Spatial Beta Diversity (Presence, quant=false)   | 0.7095              | 14.7190         | 20.7449   |
| Temporal Beta Diversity (Abundance, quant=true)  | 4.0110              | 114.5111        | 28.5494   |
| Temporal Beta Diversity (Abundance, quant=false) | 1.3314              | 17.6840         | 13.2821   |
| Temporal Beta Diversity (Presence, quant=false)  | 1.2106              | 18.3145         | 15.1279   |
| Dispersal-niche continuum index                  | 53.4701             | 9162.7688       | 171.3625  |
| Occupied Patches Proportion                      | 0.6862              | 14.9776         | 21.8266   |
| Variability Metrics                              | 7.4703              | 141.9484        | 19.0017   |
| Hypervolume Estimation                           | 0.0034              | 0.0402          | 11.7786   |
| Hypervolume Dissimilarity                        | 0.0054              | 0.2190          | 40.7483   |

### Benchmarking using the small dataset
*All times are in millisecond (ms), and memory is in mebibytes (MiB). All values are rounded up to 4 decimal places.*

#### Julia
| TestCase                                         | Time_minimum | Time_median | Time_mean | Time_maximum | Time_std | memory  |
|--------------------------------------------------|--------------|-------------|-----------|--------------|----------|---------|
| Beta Diversity (Abundance, quant=true)           | 0.0350       | 0.0423      | 0.0428    | 0.0793       | 0.0067   | 0.1225  |
| Beta Diversity (Abundance, quant=false)          | 0.0050       | 0.0115      | 0.0112    | 0.0360       | 0.0047   | 0.0883  |
| Beta Diversity (Presence, quant=false)           | 0.0048       | 0.0115      | 0.0113    | 0.0346       | 0.0045   | 0.0883  |
| Spatial Beta Diversity (Abundance, quant=true)   | 0.6158       | 0.6942      | 0.7670    | 3.3913       | 0.3912   | 1.1334  |
| Spatial Beta Diversity (Abundance, quant=false)  | 0.5008       | 0.5395      | 0.5717    | 3.3988       | 0.2870   | 0.7663  |
| Spatial Beta Diversity (Presence, quant=false)   | 0.4988       | 0.5296      | 0.5668    | 3.6272       | 0.3108   | 0.7663  |
| Temporal Beta Diversity (Abundance, quant=true)  | 3.4157       | 3.5487      | 3.8326    | 5.5600       | 0.6124   | 12.8811 |
| Temporal Beta Diversity (Abundance, quant=false) | 0.8895       | 1.0069      | 1.1532    | 4.7410       | 0.5860   | 2.8937  |
| Temporal Beta Diversity (Presence, quant=false)  | 0.8746       | 0.9875      | 1.1156    | 5.7435       | 0.6782   | 2.8937  |
| Dispersal-niche continuum index                  | 12.3592      | 14.7842     | 13.9412   | 16.2559      | 1.2164   | 37.7604 |
| Occupied Patches Proportion                      | 0.1983       | 0.2222      | 0.2643    | 1.6241       | 0.1684   | 0.3249  |
| Variability Metrics                              | 2.4110       | 2.4428      | 2.5369    | 4.3091       | 0.3546   | 3.8483  |
| Hypervolume Estimation                           | 0.0028       | 0.0031      | 0.0035    | 0.0426       | 0.0040   | 0.0059  |
| Hypervolume Dissimilarity                        | 0.0048       | 0.0056      | 0.0063    | 0.0720       | 0.0067   | 0.0082  |

#### R
| TestCase                                         | Time_minimum | Time_median | Time_mean | Time_maximum | Time_std  | memory   |
|--------------------------------------------------|--------------|-------------|-----------|--------------|-----------|----------|
| Beta Diversity (Abundance, quant=true)           | 1.1577       | 1.5902      | 2.1087    | 10.1530      | 1.4055    | 0.0195   |
| Beta Diversity (Abundance, quant=false)          | 0.2662       | 0.3371      | 0.5438    | 1.8151       | 0.3854    | 0.0444   |
| Beta Diversity (Presence, quant=false)           | 0.2606       | 0.3569      | 0.5201    | 1.1688       | 0.2659    | 0.0444   |
| Spatial Beta Diversity (Abundance, quant=true)   | 11.9829      | 16.8506     | 23.9044   | 622.2835     | 60.6749   | 1.2097   |
| Spatial Beta Diversity (Abundance, quant=false)  | 8.4718       | 12.0729     | 12.5712   | 23.3589      | 2.9096    | 0.7667   |
| Spatial Beta Diversity (Presence, quant=false)   | 8.4567       | 11.8902     | 11.9560   | 20.2556      | 2.5048    | 0.7667   |
| Temporal Beta Diversity (Abundance, quant=true)  | 74.1172      | 98.1867     | 124.9184  | 848.6431     | 123.4590  | 15.4311  |
| Temporal Beta Diversity (Abundance, quant=false) | 12.3936      | 17.6707     | 22.9253   | 411.9853     | 40.0099   | 3.7200   |
| Temporal Beta Diversity (Presence, quant=false)  | 11.4340      | 16.1048     | 16.5250   | 22.9118      | 2.8710    | 3.7200   |
| Dispersal-niche continuum index                  | 3940.5316    | 4616.0343   | 4861.0437 | 16822.0133   | 1369.1769 | 912.5856 |
| Occupied Patches Proportion                      | 10.3292      | 15.5675     | 24.8453   | 821.6664     | 80.6197   | 0.4902   |
| Variability Metrics                              | 27.9684      | 43.4103     | 58.6905   | 820.6562     | 92.9202   | 10.5805  |
| Hypervolume Estimation                           | 0.0279       | 0.0465      | 0.0652    | 0.2666       | 0.0523    | 0.0003   |
| Hypervolume Dissimilarity                        | 0.1456       | 0.1991      | 0.3145    | 1.4536       | 0.2514    | 0.0014   |

#### Speedup
| TestCase                                         | `Time_median_julia` | `Time_median_r` | `Speedup` |
|--------------------------------------------------|---------------------|-----------------|-----------|
| Beta Diversity (Abundance, quant=true)           | 0.0423              | 1.5902          | 37.6006   |
| Beta Diversity (Abundance, quant=false)          | 0.0115              | 0.3371          | 29.3155   |
| Beta Diversity (Presence, quant=false)           | 0.0115              | 0.3569          | 31.0334   |
| Spatial Beta Diversity (Abundance, quant=true)   | 0.6942              | 16.8506         | 24.2731   |
| Spatial Beta Diversity (Abundance, quant=false)  | 0.5395              | 12.0729         | 22.3779   |
| Spatial Beta Diversity (Presence, quant=false)   | 0.5296              | 11.8902         | 22.4529   |
| Temporal Beta Diversity (Abundance, quant=true)  | 3.5487              | 98.1867         | 27.6683   |
| Temporal Beta Diversity (Abundance, quant=false) | 1.0069              | 17.6707         | 17.5490   |
| Temporal Beta Diversity (Presence, quant=false)  | 0.9875              | 16.1048         | 16.3087   |
| Dispersal-niche continuum index                  | 14.7842             | 4616.0343       | 312.2269  |
| Occupied Patches Proportion                      | 0.2222              | 15.5675         | 70.0713   |
| Variability Metrics                              | 2.4428              | 43.4103         | 17.7706   |
| Hypervolume Estimation                           | 0.0031              | 0.0465          | 14.8715   |
| Hypervolume Dissimilarity                        | 0.0056              | 0.1991          | 35.6574   |

## Remarks
- Speedup is calculated as the `R` median execution time divided by the `Julia` median execution time. 
- For `DNCI_multigroup_result`, 100 permutations are used in both the julia and R function.

## The Scripts Used for Benchmarking
- [`Julia`](https://github.com/cralibe/MetaCommunityMetrics.jl/blob/main/benchmarks/benchmark_julia.jl)
- [`R`](https://github.com/cralibe/MetaCommunityMetrics.jl/blob/main/benchmarks/benchmark_r/benchmark_r.R)

## Packages used for benchmarking
- [`bench`](https://github.com/r-lib/bench)
- [`BenchmarkTools.jl`](https://github.com/JuliaCI/BenchmarkTools.jl)