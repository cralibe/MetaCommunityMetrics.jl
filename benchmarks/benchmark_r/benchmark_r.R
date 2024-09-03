# benchmarks/benchmark_r.R

# Load necessary libraries
library(microbenchmark)
library(adespatial)  # Assuming this is the package with beta.div.comp
library(readr)  # For reading CSV files

# Load the sample data
sample_matrix_abundance <- read_csv("path/to/data.csv")

# Define the R function you want to benchmark
r_results <- microbenchmark(
  beta_div_comp = beta.div.comp(sample_matrix_abundance),
  times = 100L  # Run the benchmark 100 times
)

# Print the benchmarking results
print(r_results)

# Save the results to a file
write.csv(as.data.frame(r_results), "r_benchmark_results.csv")
