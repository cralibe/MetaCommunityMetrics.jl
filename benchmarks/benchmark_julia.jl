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
mean_temporal_beta_div_1 = @benchmark mean_temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
mean_temporal_beta_div_2 = @benchmark mean_temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
mean_temporal_beta_div_3 = @benchmark mean_temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)

# Preparing the data to benchmark the DNCI functions
total_presence_df=@pipe df|>
                        groupby(_,[:Species,:Sampling_date_order])|>
                        combine(_,:Presence=>sum=>:Total_Presence) |>
                        filter(row -> row[:Total_Presence] > 1, _)

total_richness_df= @pipe df|>
innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
combine(_,:Presence=>sum=>:Total_Richness)


comm= @pipe df|>
                  innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                  filter(row -> row[:Sampling_date_order] == 1, _) |>
                  select(_, [:plot, :Species, :Presence]) |>
                  unstack(_, :Species, :Presence, fill=0) |>
                  select(_, Not(:plot)) |>
                  Matrix(_)

# Benchmark  the create_clusters function
cluster_result = @benchmark create_clusters(total_richness_df.Sampling_date_order, 
                        total_richness_df.Latitude, 
                        total_richness_df.Longitude, 
                        total_richness_df.plot, 
                        total_richness_df.Total_Richness)

# Benchmark the plot_clusters function
cluster_list = create_clusters(total_richness_df.Sampling_date_order, 
        total_richness_df.Latitude, 
        total_richness_df.Longitude, 
        total_richness_df.plot,
        total_richness_df.Total_Richness) 

plot_clusters_result = @benchmark plot_clusters(cluster_list[1].Latitude, cluster_list[1].Longitude, cluster_list[1].Group)

# Benchmark the DNCI_multigroup function
DNCI_multigroup_result = @benchmark DNCI_multigroup(comm, cluster_list[1].Group; count = false) 

# Benchmark the niche_overlap function
niche_overlap_result = @benchmark niche_overlap(df.Abundance, 
                                df.Species, 
                                df.plot, 
                                df.Sampling_date_order)
                                                                        
# Benchmark the prop_patches function
prop_patches_result = @benchmark prop_patches(metacomm_df.Presence, metacomm_df.Species, metacomm_df.plot)

# Benchmark the CV_meta function
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


proportion_use_df=
@pipe df[:,[:Abundance, :Species, :plot, :Sampling_date_order]] |>#select column N, Species, Patch, Time and env
groupby(_, [:Species]) |>
transform(_, :Abundance => (x -> x ./ sum(x)) => :relativ_N) |> #relative abundance (proportion use) of species i in patch k at time t across all sites and times
groupby(_, [:Species]) |>
transform(_, :Abundance => sum => :total_N)|> #total abundance of species i in all sites and times for cross checking
select(_, [:Species,:relativ_N, :plot ,:Sampling_date_order]) |>#select columns Species, total_N, relativ_N, env
unstack(_, :Species,:relativ_N) |> #pivot wider
_[!, Not(:plot, :Sampling_date_order)] |> # only retain the Proportional use values for each species
permutedims(_)