# Benchmarking Results
```@meta
CurrentModule = MetaCommunityMetrics
```
## Computational Resources
*All benchmarks were performed on the same machine to ensure consistent comparisons.*
- **CPU**: Apple M1 Pro
- **Number of Cores**: 10
- **Memory**: 16GB RAM
- **Operating System**: macOS Sonoma Version 14.2.1
- **Julia Version**: 1.9.3
- **R Version**: 4.3.1


Below is a comparison of the benchmarking results between my Julia package and functions/equivalent implmentation in R.

## Direct Comparison (When an equivalent fuction in R is avaliable)

| Test Case              | Julia Execution Time | R Execution Time  | Speedup (Julia/R) | Memory Usage (Julia) | Memory Usage (R) |
|------------------------|----------------------|-------------------|-------------------|----------------------|------------------|
| beta_diversity_1       | 0.118                | 2.248             | 19x               | 0.132                | 0.057            |
| beta_diversity_2       | 0.039                | 0.270             | 7x                | 0.133                | 0.125            |
| beta_diversity_3       | 0.038                | 0.280             | 7x                | 0.133                | 0.125            |
| DNCI_multigroup_result | 188.369              |                   |                   | 816.60        

*Note: All times are in millisecond, and memory is in mebibytes (MiB). Speedup is calculated as the ratio of R execution time to Julia execution time and rounded to the nearest integer. The test case names are assigned according to the object names used for the benchmarking results. The same object name is used for the same test case between Julia and R for consistency. The same data input are used for each test case.*

## Self-Benchmarking (For Julia-Only Functions)
*These functions are unique to the Julia package and do not have equivalents in R. We provide benchmarks in terms of execution time and memory usage as a point of comparison for future work or potential replication in other languages.*

| Test Case                | Julia Execution Time | Memory Usage (Julia) | Data Size (row, column)|
|--------------------------|----------------------|----------------------|------------------------|
| mean_spatial_beta_div_1  | 632.693              | 239.660              | 48735, 10              |
| mean_spatial_beta_div_2  | 625.321              | 213.450              | 48735, 10              |
| mean_spatial_beta_div_3  | 623.162              | 213.450              | 48735, 10              |
| mean_temporal_beta_div_1 | 623.162              | 213.450              | 48735, 10              |
| mean_temporal_beta_div_2 | 141.001              | 76.520               | 48735, 10              |
| mean_temporal_beta_div_3 | 141.477              | 76.520               | 48735, 10              |
| cluster_result           | 1.491                | 860.280              | 2565, 5                |
| plot_clusters_result     | 2.852                | 816.600              | 15, 3                  |
| niche_overlap_result     | 6823.000             | 118.280              | 48735, 4               |







*Note: All times are in millisecond, and memory is in mebibytes (MiB).*


