# benchmarks/data_for_testing/data_generation.jl
using Pkg
Pkg.activate()

using CSV
using MetaCommunityMetrics
using DataFrames
using Pipe
using Random

#Read in the sample data
df = load_sample_data()

# Define the percentages of data to subset for small and medium datasets
small_percent=0.1
medium_percent=0.5

small_size = round(Int, nrow(df) * small_percent)
medium_size = round(Int, nrow(df) * medium_percent)

# Shuffle the data to ensure random sampling
Random.seed!(123)  # Set seed for reproducibility

# Shuffle the DataFrame
shuffled_df = df[randperm(nrow(df)), :]

# Create datasets of different sizes
small_df = shuffled_df[1:small_size, :]
medium_df = shuffled_df[1:medium_size, :]

# Save the datasets to CSV files
small_path = joinpath("data", "data_for_testing", "small_dataset.csv")
medium_path = joinpath("data", "data_for_testing", "medium_dataset.csv")

CSV.write(small_path, small_df)
CSV.write(medium_path, medium_df)