# benchmarks/benchmark_julia_vs_r.jl
using Pkg
Pkg.activate(".")
using BenchmarkTools
using CSV
using MetaCommunityMetrics
using DataFrames
using Pipe

#Read in the sample data (small, medium and full sizes)
full_df = load_sample_data()
medium_df = CSV.read("data/data_for_testing/medium_dataset.csv",  DataFrame)
small_df = CSV.read("data/data_for_testing/small_dataset.csv",  DataFrame)
    
# Data wrangling 
##Full dataset
df=full_df
### The abundance matrix
matrix_with_abundance = @pipe df |> 
           select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Presence)) |> 
           filter(row -> row[:Sampling_date_order] == 50, _) |> 
           unstack(_, :Species, :Abundance, fill=0) |>  
           select(_, Not(:Sampling_date_order, :plot)) |> 
           Matrix(_) |> #convert the dataframe to a matrix
           _[:, sum(_, dims=1)[1, :] .!= 0] |> 
           _[sum(_, dims=2)[:, 1] .!= 0,:] 
### The binary matrix
matrix_with_presence =  @pipe df |> 
           select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Abundance)) |> 
           filter(row -> row[:Sampling_date_order] == 50, _) |>
           unstack(_, :Species, :Presence, fill=0) |> 
           select(_, Not(:Sampling_date_order, :plot)) |> 
           Matrix(_) |> #convert the dataframe to a matrix      
           _[:, sum(_, dims=1)[1, :] .!= 0] |> 
           _[sum(_, dims=2)[:, 1] .!= 0,:]

#### Benchmark the `beta_diversity` function
beta_diversity_1 = @benchmark beta_diversity(matrix_with_abundance; quant=true) samples=100 evals=1
beta_diversity_2 = @benchmark beta_diversity(matrix_with_abundance; quant=false) samples=100 evals=1
beta_diversity_3 = @benchmark beta_diversity(matrix_with_presence; quant=false) samples=100 evals=1

@ballocated beta_diversity(matrix_with_abundance; quant=false) samples=100 evals=1

# Benchmark the mean_spatial_beta_div function
spatial_beta_div_1 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true) samples=100 evals=1
spatial_beta_div_2 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1
spatial_beta_div_3 = @benchmark spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1

# Benchmark the mean_temporal_beta_div function
temporal_beta_div_1 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true) samples=100 evals=1 seconds=1800
temporal_beta_div_2 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1
temporal_beta_div_3 = @benchmark temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1

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
                        total_richness_df.Total_Richness) samples=100 evals=1 seconds=1800

# Benchmark the plot_clusters function
cluster_list = create_clusters(total_richness_df.Sampling_date_order, 
        total_richness_df.Latitude, 
        total_richness_df.Longitude, 
        total_richness_df.plot,
        total_richness_df.Total_Richness)

plot_clusters_result = @benchmark plot_clusters(cluster_list[50].Latitude, cluster_list[50].Longitude, cluster_list[50].Group) samples=100 evals=1

# Save the wrangled data to CSV files for the R benchmarks
comm_full_df= @pipe full_df|>
            innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
            innerjoin(_,  total_richness_df, on = [:plot, :Sampling_date_order], makeunique = true) |>
            filter(row -> row[:Sampling_date_order] == 50, _) |>
            select(_, [:plot, :Species, :Presence]) |>
            unstack(_, :Species, :Presence, fill=0) |>
            select(_, Not(:plot)) 
#Save the data to CSV files for the R benchmarks
CSV.write("data/data_for_testing/comm_full_df.csv", comm_full_df)
CSV.write("data/data_for_testing/groups_full_df.csv", cluster_list[50])

# Benchmark the DNCI_multigroup function
DNCI_multigroup_result = @benchmark DNCI_multigroup(comm, cluster_list[50].Group, 100; count = false) samples=100 evals=1 seconds=1800

## Benchmark the niche_overlap function
niche_overlap_result = @benchmark niche_overlap(df.Abundance, 
                                df.Species, 
                                df.plot, 
                                df.Sampling_date_order) samples=100 evals=1 seconds=1800

save_object("benchmarks/result/niche_overlap_full_df_result.jld2", niche_overlap_result)

                                                                        
## Benchmark the prop_patches function
prop_patches_result = @benchmark prop_patches(df.Presence, df.Species, df.plot) samples=100 evals=1

## Benchmark the variability metrics function
# Benchmark the CV_meta function
CV_meta_result = @benchmark CV_meta(df.Abundance, 
                    df.Sampling_date_order,
                    df.plot, 
                    df.Species) samples=100 evals=1 seconds=1800


# Benchmark the hypervolume functions
data_1 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])

data_2 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "SH", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])

hypervolume_det_result = @benchmark MVNH_det(data_1; var_names=["Temperature", "Precipitation"]) samples=100 evals=1 seconds=1800

hypervolume_dis_result = @benchmark MVNH_dissimilarity(data_1, data_2; var_names=["Temperature", "Precipitation"]) samples=100 evals=1 seconds=1800

# Save the results to a DataFrame
benchmark_result_full_df = DataFrame(TestCase = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                        "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                        "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                        "cluster_result", "plot_clusters_result",
                                        "DNCI_multigroup_result", "niche_overlap_result",
                                        "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result"],
                                        Time_minimum = [minimum(beta_diversity_1.times), minimum(beta_diversity_2.times), minimum(beta_diversity_3.times),
                                        minimum(spatial_beta_div_1.times), minimum(spatial_beta_div_2.times), minimum(spatial_beta_div_3.times),
                                        minimum(temporal_beta_div_1.times), minimum(temporal_beta_div_2.times), minimum(temporal_beta_div_3.times),
                                        minimum(cluster_result.times), minimum(plot_clusters_result.times),
                                        minimum(DNCI_multigroup_result.times), minimum(niche_overlap_result.times),
                                        minimum(prop_patches_result.times), minimum(CV_meta_result.times), minimum(hypervolume_det_result.times), minimum(hypervolume_dis_result.times)],
                                        Time_median = [median(beta_diversity_1.times), median(beta_diversity_2.times), median(beta_diversity_3.times),
                                        median(spatial_beta_div_1.times), median(spatial_beta_div_2.times), median(spatial_beta_div_3.times),
                                        median(temporal_beta_div_1.times), median(temporal_beta_div_2.times), median(temporal_beta_div_3.times),
                                        median(cluster_result.times), median(plot_clusters_result.times),
                                        median(DNCI_multigroup_result.times), median(niche_overlap_result.times),
                                        median(prop_patches_result.times), median(CV_meta_result.times), median(hypervolume_det_result.times), median(hypervolume_dis_result.times)],
                                        Time_mean = [mean(beta_diversity_1.times), mean(beta_diversity_2.times), mean(beta_diversity_3.times),
                                        mean(spatial_beta_div_1.times), mean(spatial_beta_div_2.times), mean(spatial_beta_div_3.times),
                                        mean(temporal_beta_div_1.times), mean(temporal_beta_div_2.times), mean(temporal_beta_div_3.times),
                                        mean(cluster_result.times), mean(plot_clusters_result.times),
                                        mean(DNCI_multigroup_result.times), mean(niche_overlap_result.times),
                                        mean(prop_patches_result.times), mean(CV_meta_result.times), mean(hypervolume_det_result.times), mean(hypervolume_dis_result.times)],
                                        Time_maximum = [maximum(beta_diversity_1.times), maximum(beta_diversity_2.times), maximum(beta_diversity_3.times),
                                        maximum(spatial_beta_div_1.times), maximum(spatial_beta_div_2.times), maximum(spatial_beta_div_3.times),
                                        maximum(temporal_beta_div_1.times), maximum(temporal_beta_div_2.times), maximum(temporal_beta_div_3.times),
                                        maximum(cluster_result.times), maximum(plot_clusters_result.times),
                                        maximum(DNCI_multigroup_result.times), maximum(niche_overlap_result.times),
                                        maximum(prop_patches_result.times), maximum(CV_meta_result.times), maximum(hypervolume_det_result.times), maximum(hypervolume_dis_result.times)],
                                        Time_std = [std(beta_diversity_1.times), std(beta_diversity_2.times), std(beta_diversity_3.times),
                                        std(spatial_beta_div_1.times), std(spatial_beta_div_2.times), std(spatial_beta_div_3.times),
                                        std(temporal_beta_div_1.times), std(temporal_beta_div_2.times), std(temporal_beta_div_3.times),
                                        std(cluster_result.times), std(plot_clusters_result.times),
                                        std(DNCI_multigroup_result.times), std(niche_overlap_result.times),
                                        std(prop_patches_result.times), std(CV_meta_result.times), std(hypervolume_det_result.times), std(hypervolume_dis_result.times)],
                                        memory = [beta_diversity_1.memory, beta_diversity_2.memory, beta_diversity_3.memory,
                                        spatial_beta_div_1.memory, spatial_beta_div_2.memory, spatial_beta_div_3.memory,
                                        temporal_beta_div_1.memory, temporal_beta_div_2.memory, temporal_beta_div_3.memory,
                                        cluster_result.memory, plot_clusters_result.memory,
                                        DNCI_multigroup_result.memory, niche_overlap_result.memory,
                                        prop_patches_result.memory, CV_meta_result.memory, hypervolume_det_result.memory, hypervolume_dis_result.memory])   

# Divide the time by 1e6 to convert from nanoseconds to milliseconds
benchmark_result_full_df[:, 2:6] .= benchmark_result_full_df[:, 2:6] ./ 1e6
# Divide the memory by 1024^2 to convert from bytes to mebibytes (Mib)
# convert column type to Float64 before division
benchmark_result_full_df[!, 7] = convert.(Float64, benchmark_result_full_df[:, 7])
# perform the division
benchmark_result_full_df[:, 7] .= benchmark_result_full_df[:, 7] ./ (1024^2)

CSV.write("benchmarks/result/benchmark_result_full_df_julia.csv", benchmark_result_full_df)

##Medium dataset
df=medium_df
### The abundance matrix
matrix_with_abundance = @pipe df |> 
                select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Presence)) |> 
                filter(row -> row[:Sampling_date_order] == 50, _) |> 
                unstack(_, :Species, :Abundance, fill=0) |>  
           select(_, Not(:Sampling_date_order, :plot)) |> 
           Matrix(_) |> #convert the dataframe to a matrix
           _[:, sum(_, dims=1)[1, :] .!= 0] |> 
           _[sum(_, dims=2)[:, 1] .!= 0,:] 
### The binary matrix
matrix_with_presence =  @pipe df |> 
           select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Abundance)) |> 
           filter(row -> row[:Sampling_date_order] == 50, _) |>
           unstack(_, :Species, :Presence, fill=0) |> 
           select(_, Not(:Sampling_date_order, :plot)) |> 
           Matrix(_) |> #convert the dataframe to a matrix      
           _[:, sum(_, dims=1)[1, :] .!= 0] |> 
           _[sum(_, dims=2)[:, 1] .!= 0,:]

#### Benchmark the `beta_diversity` function
beta_diversity_1 = @benchmark beta_diversity(matrix_with_abundance; quant=true) samples=100 evals=1
beta_diversity_2 = @benchmark beta_diversity(matrix_with_abundance; quant=false) samples=100 evals=1
beta_diversity_3 = @benchmark beta_diversity(matrix_with_presence; quant=false) samples=100 evals=1

# Benchmark the mean_spatial_beta_div function
spatial_beta_div_1 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true) samples=100 evals=1
spatial_beta_div_2 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1
spatial_beta_div_3 = @benchmark spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1

# Benchmark the mean_temporal_beta_div function
temporal_beta_div_1 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true) samples=100 evals=1 seconds=1800
temporal_beta_div_2 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1
temporal_beta_div_3 = @benchmark temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1

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
                        total_richness_df.Total_Richness) samples=100 evals=1 seconds=1800

# Benchmark the plot_clusters function
cluster_list = create_clusters(total_richness_df.Sampling_date_order, 
        total_richness_df.Latitude, 
        total_richness_df.Longitude, 
        total_richness_df.plot,
        total_richness_df.Total_Richness)

plot_clusters_result = @benchmark plot_clusters(cluster_list[50].Latitude, cluster_list[50].Longitude, cluster_list[50].Group) samples=100 evals=1

# Save the wrangled data to CSV files for the R benchmarks
comm_medium_df= @pipe df|>
            innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
            innerjoin(_,  total_richness_df, on = [:plot, :Sampling_date_order], makeunique = true) |>
            filter(row -> row[:Sampling_date_order] == 50, _) |>
            select(_, [:plot, :Species, :Presence]) |>
            unstack(_, :Species, :Presence, fill=0) |>
            select(_, Not(:plot)) 
CSV.write("data/data_for_testing/comm_medium_df.csv", comm_medium_df)
CSV.write("data/data_for_testing/groups_medium_df.csv", cluster_list[50])

# Benchmark the DNCI_multigroup function
DNCI_multigroup_result = @benchmark DNCI_multigroup(comm, cluster_list[50].Group, 100; count = false) samples=100 evals=1 seconds=1800 

## Benchmark the niche_overlap function
niche_overlap_result = @benchmark niche_overlap(df.Abundance, 
                                df.Species, 
                                df.plot, 
                                df.Sampling_date_order) samples=100 evals=1 seconds=1800

save_object("benchmarks/result/niche_overlap_medium_df_result.jld2", niche_overlap_result)
                                
## Benchmark the prop_patches function
prop_patches_result = @benchmark prop_patches(df.Presence, df.Species, df.plot) samples=100 evals=1

## Benchmark the variability metrics function
# Benchmark the CV_meta function
CV_meta_result = @benchmark CV_meta(df.Abundance, 
                    df.Sampling_date_order,
                    df.plot, 
                    df.Species) samples=100 evals=1 seconds=1800

# Benchmark the hypervolume functions
data_1 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])

data_2 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "SH", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])

hypervolume_det_result = @benchmark MVNH_det(data_1; var_names=["Temperature", "Precipitation"]) samples=100 evals=1 seconds=1800

hypervolume_dis_result = @benchmark MVNH_dissimilarity(data_1, data_2; var_names=["Temperature", "Precipitation"]) samples=100 evals=1 seconds=1800

                    # Save the results to a DataFrame
benchmark_result_medium_df = DataFrame(TestCase = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                                "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                                "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                                "cluster_result", "plot_clusters_result",
                                                "DNCI_multigroup_result", "niche_overlap_result",
                                                "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result"],
                                                Time_minimum = [minimum(beta_diversity_1.times), minimum(beta_diversity_2.times), minimum(beta_diversity_3.times),
                                                minimum(spatial_beta_div_1.times), minimum(spatial_beta_div_2.times), minimum(spatial_beta_div_3.times),
                                                minimum(temporal_beta_div_1.times), minimum(temporal_beta_div_2.times), minimum(temporal_beta_div_3.times),
                                                minimum(cluster_result.times), minimum(plot_clusters_result.times),
                                                minimum(DNCI_multigroup_result.times), minimum(niche_overlap_result.times),
                                                minimum(prop_patches_result.times), minimum(CV_meta_result.times), minimum(hypervolume_det_result.times), minimum(hypervolume_dis_result.times)],
                                                Time_median = [median(beta_diversity_1.times), median(beta_diversity_2.times), median(beta_diversity_3.times),
                                                median(spatial_beta_div_1.times), median(spatial_beta_div_2.times), median(spatial_beta_div_3.times),
                                                median(temporal_beta_div_1.times), median(temporal_beta_div_2.times), median(temporal_beta_div_3.times),
                                                median(cluster_result.times), median(plot_clusters_result.times),
                                                median(DNCI_multigroup_result.times), median(niche_overlap_result.times),
                                                median(prop_patches_result.times), median(CV_meta_result.times), median(hypervolume_det_result.times), median(hypervolume_dis_result.times)],
                                                Time_mean = [mean(beta_diversity_1.times), mean(beta_diversity_2.times), mean(beta_diversity_3.times),
                                                mean(spatial_beta_div_1.times), mean(spatial_beta_div_2.times), mean(spatial_beta_div_3.times),
                                                mean(temporal_beta_div_1.times), mean(temporal_beta_div_2.times), mean(temporal_beta_div_3.times),
                                                mean(cluster_result.times), mean(plot_clusters_result.times),
                                                mean(DNCI_multigroup_result.times), mean(niche_overlap_result.times),
                                                mean(prop_patches_result.times), mean(CV_meta_result.times), mean(hypervolume_det_result.times), mean(hypervolume_dis_result.times)],
                                                Time_maximum = [maximum(beta_diversity_1.times), maximum(beta_diversity_2.times), maximum(beta_diversity_3.times),
                                                maximum(spatial_beta_div_1.times), maximum(spatial_beta_div_2.times), maximum(spatial_beta_div_3.times),
                                                maximum(temporal_beta_div_1.times), maximum(temporal_beta_div_2.times), maximum(temporal_beta_div_3.times),
                                                maximum(cluster_result.times), maximum(plot_clusters_result.times),
                                                maximum(DNCI_multigroup_result.times), maximum(niche_overlap_result.times),
                                                maximum(prop_patches_result.times), maximum(CV_meta_result.times), maximum(hypervolume_det_result.times), maximum(hypervolume_dis_result.times)],
                                                Time_std = [std(beta_diversity_1.times), std(beta_diversity_2.times), std(beta_diversity_3.times),
                                                std(spatial_beta_div_1.times), std(spatial_beta_div_2.times), std(spatial_beta_div_3.times),
                                                std(temporal_beta_div_1.times), std(temporal_beta_div_2.times), std(temporal_beta_div_3.times),
                                                std(cluster_result.times), std(plot_clusters_result.times),
                                                std(DNCI_multigroup_result.times), std(niche_overlap_result.times),
                                                std(prop_patches_result.times), std(CV_meta_result.times), std(hypervolume_det_result.times), std(hypervolume_dis_result.times)],
                                                memory = [beta_diversity_1.memory, beta_diversity_2.memory, beta_diversity_3.memory,
                                                spatial_beta_div_1.memory, spatial_beta_div_2.memory, spatial_beta_div_3.memory,
                                                temporal_beta_div_1.memory, temporal_beta_div_2.memory, temporal_beta_div_3.memory,
                                                cluster_result.memory, plot_clusters_result.memory,
                                                DNCI_multigroup_result.memory, niche_overlap_result.memory,
                                                prop_patches_result.memory, CV_meta_result.memory, hypervolume_det_result.memory, hypervolume_dis_result.memory])   
# Save the results to a CSV file        
# Divide the time by 1e6 to convert from nanoseconds to milliseconds
benchmark_result_medium_df[:, 2:6] .= benchmark_result_medium_df[:, 2:6] ./ 1e6
# Divide the memory by 1024^2 to convert from bytes to mebibytes (Mib)
# convert column type to Float64 before division
benchmark_result_medium_df[!, 7] = convert.(Float64, benchmark_result_medium_df[:, 7])
# perform the division
benchmark_result_medium_df[:, 7] .= benchmark_result_medium_df[:, 7] ./ (1024^2)

CSV.write("benchmarks/result/benchmark_result_medium_df_julia.csv", benchmark_result_medium_df)


##Full dataset
df=small_df
### The abundance matrix
matrix_with_abundance = @pipe df |> 
           select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Presence)) |> 
           filter(row -> row[:Sampling_date_order] == 55, _) |> 
           unstack(_, :Species, :Abundance, fill=0) |>  
           select(_, Not(:Sampling_date_order, :plot)) |> 
           Matrix(_) |> #convert the dataframe to a matrix
           _[:, sum(_, dims=1)[1, :] .!= 0] |> 
           _[sum(_, dims=2)[:, 1] .!= 0,:] 
### The binary matrix
matrix_with_presence =  @pipe df |> 
           select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Abundance)) |> 
           filter(row -> row[:Sampling_date_order] == 55, _) |>
           unstack(_, :Species, :Presence, fill=0) |> 
           select(_, Not(:Sampling_date_order, :plot)) |> 
           Matrix(_) |> #convert the dataframe to a matrix      
           _[:, sum(_, dims=1)[1, :] .!= 0] |> 
           _[sum(_, dims=2)[:, 1] .!= 0,:]

#### Benchmark the `beta_diversity` function
beta_diversity_1 = @benchmark beta_diversity(matrix_with_abundance; quant=true) samples=100 evals=1
beta_diversity_2 = @benchmark beta_diversity(matrix_with_abundance; quant=false) samples=100 evals=1
beta_diversity_3 = @benchmark beta_diversity(matrix_with_presence; quant=false) samples=100 evals=1

# Benchmark the mean_spatial_beta_div function
spatial_beta_div_1 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true) samples=100 evals=1
spatial_beta_div_2 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1
spatial_beta_div_3 = @benchmark spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1

# Benchmark the mean_temporal_beta_div function
temporal_beta_div_1 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true) samples=100 evals=1 seconds=1800
temporal_beta_div_2 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1
temporal_beta_div_3 = @benchmark temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1

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
            filter(row -> row[:Sampling_date_order] == 55, _) |>
            select(_, [:plot, :Species, :Presence]) |>
            unstack(_, :Species, :Presence, fill=0) |>
            select(_, Not(:plot)) |>
            Matrix(_)

# Benchmark  the create_clusters function
cluster_result = @benchmark create_clusters(total_richness_df.Sampling_date_order, 
                        total_richness_df.Latitude, 
                        total_richness_df.Longitude, 
                        total_richness_df.plot, 
                        total_richness_df.Total_Richness) samples=100 evals=1 seconds=1800

# Benchmark the plot_clusters function
cluster_list = create_clusters(total_richness_df.Sampling_date_order, 
        total_richness_df.Latitude, 
        total_richness_df.Longitude, 
        total_richness_df.plot,
        total_richness_df.Total_Richness)

plot_clusters_result = @benchmark plot_clusters(cluster_list[55].Latitude, cluster_list[55].Longitude, cluster_list[55].Group) samples=100 evals=1

# Save the wrangled data to CSV files for the R benchmarks
comm_small_df= @pipe full_df|>
            innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
            innerjoin(_,  total_richness_df, on = [:plot, :Sampling_date_order], makeunique = true) |>
            filter(row -> row[:Sampling_date_order] == 55, _) |>
            select(_, [:plot, :Species, :Presence]) |>
            unstack(_, :Species, :Presence, fill=0) |>
            select(_, Not(:plot)) 
#Save the data to CSV files for the R benchmarks
CSV.write("data/data_for_testing/comm_small_df.csv", comm_small_df)
CSV.write("data/data_for_testing/groups_small_df.csv", cluster_list[55])

# Benchmark the DNCI_multigroup function
DNCI_multigroup_result = @benchmark DNCI_multigroup(comm, cluster_list[55].Group, 100; count = false) samples=100 evals=1 seconds=1800 

## Benchmark the niche_overlap function
niche_overlap_result = @benchmark niche_overlap(df.Abundance, 
                                df.Species, 
                                df.plot, 
                                df.Sampling_date_order) samples=100 evals=1 seconds=1800

save_object("benchmarks/result/niche_overlap_small_df_result.jld2", niche_overlap_result)
                                                                       
## Benchmark the prop_patches function
prop_patches_result = @benchmark prop_patches(df.Presence, df.Species, df.plot) samples=100 evals=1

## Benchmark the variability metrics function
# Benchmark the CV_meta function
CV_meta_result = @benchmark CV_meta(df.Abundance, 
                    df.Sampling_date_order,
                    df.plot, 
                    df.Species) samples=100 evals=1 seconds=1800


# Benchmark the hypervolume functions
data_1 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "BA", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])

data_2 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "SH", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])

hypervolume_det_result = @benchmark MVNH_det(data_1; var_names=["Temperature", "Precipitation"]) samples=100 evals=1 seconds=1800

hypervolume_dis_result = @benchmark MVNH_dissimilarity(data_1, data_2; var_names=["Temperature", "Precipitation"]) samples=100 evals=1 seconds=1800

# Save the results to a DataFrame
benchmark_result_small_df = DataFrame(TestCase = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                                "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                                "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                                "cluster_result", "plot_clusters_result",
                                                "DNCI_multigroup_result", "niche_overlap_result",
                                                "prop_patches_result", "CV_meta_result", "hypervolume_det_result", "hypervolume_dis_result"],
                                                Time_minimum = [minimum(beta_diversity_1.times), minimum(beta_diversity_2.times), minimum(beta_diversity_3.times),
                                                minimum(spatial_beta_div_1.times), minimum(spatial_beta_div_2.times), minimum(spatial_beta_div_3.times),
                                                minimum(temporal_beta_div_1.times), minimum(temporal_beta_div_2.times), minimum(temporal_beta_div_3.times),
                                                minimum(cluster_result.times), minimum(plot_clusters_result.times),
                                                minimum(DNCI_multigroup_result.times), minimum(niche_overlap_result.times),
                                                minimum(prop_patches_result.times), minimum(CV_meta_result.times), minimum(hypervolume_det_result.times), minimum(hypervolume_dis_result.times)],
                                                Time_median = [median(beta_diversity_1.times), median(beta_diversity_2.times), median(beta_diversity_3.times),
                                                median(spatial_beta_div_1.times), median(spatial_beta_div_2.times), median(spatial_beta_div_3.times),
                                                median(temporal_beta_div_1.times), median(temporal_beta_div_2.times), median(temporal_beta_div_3.times),
                                                median(cluster_result.times), median(plot_clusters_result.times),
                                                median(DNCI_multigroup_result.times), median(niche_overlap_result.times),
                                                median(prop_patches_result.times), median(CV_meta_result.times), median(hypervolume_det_result.times), median(hypervolume_dis_result.times)],
                                                Time_mean = [mean(beta_diversity_1.times), mean(beta_diversity_2.times), mean(beta_diversity_3.times),
                                                mean(spatial_beta_div_1.times), mean(spatial_beta_div_2.times), mean(spatial_beta_div_3.times),
                                                mean(temporal_beta_div_1.times), mean(temporal_beta_div_2.times), mean(temporal_beta_div_3.times),
                                                mean(cluster_result.times), mean(plot_clusters_result.times),
                                                mean(DNCI_multigroup_result.times), mean(niche_overlap_result.times),
                                                mean(prop_patches_result.times), mean(CV_meta_result.times), mean(hypervolume_det_result.times), mean(hypervolume_dis_result.times)],
                                                Time_maximum = [maximum(beta_diversity_1.times), maximum(beta_diversity_2.times), maximum(beta_diversity_3.times),
                                                maximum(spatial_beta_div_1.times), maximum(spatial_beta_div_2.times), maximum(spatial_beta_div_3.times),
                                                maximum(temporal_beta_div_1.times), maximum(temporal_beta_div_2.times), maximum(temporal_beta_div_3.times),
                                                maximum(cluster_result.times), maximum(plot_clusters_result.times),
                                                maximum(DNCI_multigroup_result.times), maximum(niche_overlap_result.times),
                                                maximum(prop_patches_result.times), maximum(CV_meta_result.times), maximum(hypervolume_det_result.times), maximum(hypervolume_dis_result.times)],
                                                Time_std = [std(beta_diversity_1.times), std(beta_diversity_2.times), std(beta_diversity_3.times),
                                                std(spatial_beta_div_1.times), std(spatial_beta_div_2.times), std(spatial_beta_div_3.times),
                                                std(temporal_beta_div_1.times), std(temporal_beta_div_2.times), std(temporal_beta_div_3.times),
                                                std(cluster_result.times), std(plot_clusters_result.times),
                                                std(DNCI_multigroup_result.times), std(niche_overlap_result.times),
                                                std(prop_patches_result.times), std(CV_meta_result.times), std(hypervolume_det_result.times), std(hypervolume_dis_result.times)],
                                                memory = [beta_diversity_1.memory, beta_diversity_2.memory, beta_diversity_3.memory,
                                                spatial_beta_div_1.memory, spatial_beta_div_2.memory, spatial_beta_div_3.memory,
                                                temporal_beta_div_1.memory, temporal_beta_div_2.memory, temporal_beta_div_3.memory,
                                                cluster_result.memory, plot_clusters_result.memory,
                                                DNCI_multigroup_result.memory, niche_overlap_result.memory,
                                                prop_patches_result.memory, CV_meta_result.memory, hypervolume_det_result.memory, hypervolume_dis_result.memory])   
# Divide the time by 1e6 to convert from nanoseconds to milliseconds
benchmark_result_small_df[:, 2:6] .= benchmark_result_small_df[:, 2:6] ./ 1e6
# Divide the memory by 1024^2 to convert from bytes to mebibytes (Mib)
# convert column type to Float64 before division
benchmark_result_small_df[!, 7] = convert.(Float64, benchmark_result_small_df[:, 7])
# perform the division
benchmark_result_small_df[:, 7] .= benchmark_result_small_df[:, 7] ./ (1024^2)

CSV.write("benchmarks/result/benchmark_result_small_df_julia.csv", benchmark_result_small_df)
