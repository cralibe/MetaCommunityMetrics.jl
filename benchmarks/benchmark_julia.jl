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
spatial_beta_div_1 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
spatial_beta_div_2 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
spatial_beta_div_3 = @benchmark spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)

# Benchmark the mean_temporal_beta_div function
temporal_beta_div_1 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true)
temporal_beta_div_2 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false)
temporal_beta_div_3 = @benchmark temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false)

## Preparing the data to benchmark the DNCI functions
total_presence_df=@pipe df|>
                        groupby(_,[:Species,:Sampling_date_order])|>
                        combine(_,:Presence=>sum=>:Total_Presence) |>
                        filter(row -> row[:Total_Presence] > 1, _)

total_richness_df= @pipe df|>
innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
combine(_,:Presence=>sum=>:Total_Richness)|>
filter(row -> row[:Total_Richness] > 0, _) 


comm= @pipe df|>
            innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
            innerjoin(_,  total_richness_df, on = [:plot, :Sampling_date_order], makeunique = true) |>
            filter(row -> row[:Sampling_date_order] == 50, _) |>
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

plot_clusters_result = @benchmark plot_clusters(cluster_list[50].Latitude, cluster_list[50].Longitude, cluster_list[50].Group)

# Save the wrangled data to CSV files for the R benchmarks
comm_for_R_t50= @pipe df|>
            innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
            innerjoin(_,  total_richness_df, on = [:plot, :Sampling_date_order], makeunique = true) |>
            filter(row -> row[:Sampling_date_order] == 50, _) |>
            select(_, [:plot, :Species, :Presence]) |>
            unstack(_, :Species, :Presence, fill=0) |>
            select(_, Not(:plot)) 
CSV.write("benchmarks/benchmark_r/data/DNCI_comm_t50.csv", comm_for_R_t50)
CSV.write("benchmarks/benchmark_r/data/cluster_list_t50.csv", cluster_list[50])

# Benchmark the DNCI_multigroup function
DNCI_multigroup_result = @benchmark DNCI_multigroup(comm, cluster_list[50].Group, 100; count = false) 

## Benchmark the niche_overlap function
niche_overlap_result = @benchmark niche_overlap(df.Abundance, 
                                df.Species, 
                                df.plot, 
                                df.Sampling_date_order)
                                                                        
## Benchmark the prop_patches function
prop_patches_result = @benchmark prop_patches(df.Presence, df.Species, df.plot)

## Benchmark the variability metrics function
# Benchmark the CV_meta function
CV_meta_result = @benchmark CV_meta(df.Abundance, 
                    df.Sampling_date_order,
                    df.plot, 
                    df.Species)

# Benchmark the CV_meta_simple function
CV_meta_simple_result = @benchmark CV_meta_simple(df.Abundance, 
                            df.Sampling_date_order,
                            df.plot, 
                            df.Species)    
                            