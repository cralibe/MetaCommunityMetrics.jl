# Benchmarking Results
```@meta
CurrentModule = MetaCommunityMetrics
```
## Computational Resources
*All benchmarks were performed on the same machine to ensure consistent comparisons.*
- **CPU**: Apple M1 Pro
- **Number of Cores**: 10
- **Memory**: 16GB RAM
- **Operating System**: macOS Sonoma 14.7
- **Julia Version**: 1.9.3
- **R Version**: 4.3.1


Below is a comparison of the benchmarking results between my Julia package and functions/equivalent implmentation in R.

## Direct Comparison (When an equivalent fuction in R is avaliable)

| Test Case               | Julia Execution Time | R Execution Time  | Speedup (Julia/R) | Memory Usage (Julia) | Memory Usage (R) | R Function and Package used              | 
|-------------------------|----------------------|-------------------|-------------------|----------------------|------------------|------------------------------------------|
| beta_diversity_1        | 0.118                | 2.248             | 19x               | 0.132                | 0.057            | `beta.div.comp` from `adespatial`        |
| beta_diversity_2        | 0.039                | 0.270             | 7x                | 0.133                | 0.125            | `beta.div.comp` from `adespatial`        |
| beta_diversity_3        | 0.038                | 0.280             | 7x                | 0.133                | 0.125            | `beta.div.comp` from `adespatial`        |
| spatial_beta_div_1      | 0.954                | 19.530            | 20x               | 2.300                | 3.451            | implementation from Guzman et al. (2022) |
| spatial_beta_div_2      | 0.631                | 17.106            | 27x               | 1.900                | 3.009            | implementation from Guzman et al. (2022) | 
| spatial_beta_div_3      | 0.586                | 16.742            | 29x               | 1.900                | 3.009            | implementation from Guzman et al. (2022) | 
| temporal_beta_div_1     | 9.012                | 86.396            | 10x               | 14.760               | 17.342           | implementation from Guzman et al. (2022) | 
| temporal_beta_div_2     | 1.439                | 20.516            | 14x               | 3.540                | 5.630            | implementation from Guzman et al. (2022) | 
| temporal_beta_div_3     | 1.438                | 21.006            | 15x               | 3.540                | 5.630            | implementation from Guzman et al. (2022) | 
| DNCI_multigroup_result  | 261.033              | 45849.450         | 176x              | 407.07               | 10630.390        | `DNCI_multigroup` from `DNCImper`        |       
| prop_patches_result     | 1.491                | 18.726            | 13x               | 2.320                | 2.412            | implementation from Guzman et al. (2022) |
| CV_meta_simple_result   | 105.155              | 152.730           | 1x                | 35.41                | 55.435           | implementation from Wang et al. (2019)   |



*Note: All times are in millisecond, and memory is in mebibytes (MiB). Only mean execution times are presented here. Speedup is calculated as the ratio of R execution time to Julia execution time and rounded to the nearest integer. The test case names are assigned according to the object names used for the benchmarking results. The same object name is used for the same test case between Julia and R for consistency. The same data input are used for each test case.*

## Self-Benchmarking (For Julia-Only Functions)
*These functions are unique to the Julia package and do not have equivalents in R. We provide benchmarks in terms of execution time and memory usage as a point of comparison for future work or potential replication in other languages.*

| Test Case                | Julia Execution Time | Memory Usage (Julia) | Data Size (row, column)|
|--------------------------|----------------------|----------------------|------------------------|
| cluster_result           | 1.601                | 916.880              | 2545, 5                |
| plot_clusters_result     | 3.656                | 817.930              | 23, 3                  |
| niche_overlap_result     | 6823.000             | 118.280              | 48735, 4               |
| CV_meta_result           | 161600.000           | 96317.44             | 48735, 4               |

*Note: All times are in millisecond, and memory is in mebibytes (MiB).*

## Remarks
For `CV_meta_simple_result` and `CV_meta_result`, they are both benchmarked using the same dataset (48735, 4).
For `DNCI_multigroup_result`, 100 permutations are used in both the julia and R function.
