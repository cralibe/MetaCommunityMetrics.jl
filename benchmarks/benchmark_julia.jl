# benchmarks/benchmark_julia_vs_r.jl
using Pkg
Pkg.activate()
using BenchmarkTools
using CSV
using MetaCommunityMetrics
using DataFrames
using Pipe

#Read in the sample data
df = load_sample_data()
    
# Benchmarking Julia functions

## Benchmark the beta diversity functions
### The abundance matrix
matrix_with_abundance =  @pipe df |>
           filter(row -> row[:Sampling_date_order] == 50, _) |> 
           select(_, Not(:Presence)) |>
           unstack(_, :Species, :Abundance, fill=0) |>
           select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> 
           Matrix(_) |> 
           x -> x[:, sum(x, dims=1)[1, :] .!= 0] |> 
           x -> x[vec(sum(x, dims=2)) .!= 0, :]
### The binary matrix
matrix_with_presence = @pipe df |> 
filter(row -> row[:Sampling_date_order] == 50, _) |> 
select(_, Not(:Abundance)) |>
unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> 
Matrix(_) |>     
x -> x[:, sum(x, dims=1)[1, :] .!= 0]

#### Benchmark the `beta_diversity` function
beta_diversity_1 = @benchmark beta_diversity(matrix_with_abundance; quant=true)
beta_diversity_2 = @benchmark beta_diversity(matrix_with_abundance; quant=false)
beta_diversity_3 = @benchmark beta_diversity(matrix_with_presence; quant=false)

# Benchmark the mean_spatial_beta_div function
mean_spatial_beta_div_1 = @benchmark mean_spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
mean_spatial_beta_div_2 = @benchmark mean_spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
mean_spatial_beta_div_3 = @benchmark mean_spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)
# Benchmark the mean_temporal_beta_div function
@benchmark mean_temporal_beta_div(metacomm_df.Abundance, 
        metacomm_df.Sampling_date_order, 
        metacomm_df.plot, 
        metacomm_df.Species; quant=true)
    
# Test the create_clusters function
cluster_result = @benchmark create_clusters(total_richness_df.Sampling_date_order, 
        total_richness_df.Latitude, 
        total_richness_df.Longitude, 
        total_richness_df.plot,
        total_richness_df.Total_Richness) 


# Test the plot_clusters function
cluster_result = create_clusters(total_richness_df.Sampling_date_order, 
        total_richness_df.Latitude, 
        total_richness_df.Longitude, 
        total_richness_df.plot,
        total_richness_df.Total_Richness) 

@benchmark plot_clusters(cluster_result[1].Latitude, cluster_result[1].Longitude, cluster_result[1].Group)

# Test the DNCI_multigroup function

    
@benchmark DNCI_multigroup(comm, presence_df.Group; count = false) 

# Test the niche_overlap function
@benchmark niche_overlap(metacomm_df.Abundance, 
                        metacomm_df.Species, 
                        metacomm_df.plot, 
                        metacomm_df.Sampling_date_order)
                                                                        
# Test the prop_patches function
@benchmark prop_patches(metacomm_df.Presence, metacomm_df.Species, metacomm_df.plot)

# Test the CV_meta function
CV_test_df = @pipe metacomm_df |>
filter(row -> row[:Sampling_date_order] < 20, _)

@benchmark CV_meta(CV_test_df.Abundance, 
                    CV_test_df.Sampling_date_order,
                    CV_test_df.plot, 
                    CV_test_df.Species)

# Test the CV_meta_simple function
@benchmark CV_meta_simple(CV_test_df.Abundance, 
                            CV_test_df.Sampling_date_order,
                            CV_test_df.plot, 
                            CV_test_df.Species)    


# Save results to a file
open("julia_benchmark_results.txt", "w") do f
    write(f, string(julia_results))
end

