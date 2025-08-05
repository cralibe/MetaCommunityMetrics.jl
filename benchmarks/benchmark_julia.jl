# benchmarks/benchmark_julia_vs_r.jl
using Pkg
Pkg.activate("benchmarks")
using BenchmarkTools
using CSV
using MetaCommunityMetrics
using DataFrames
using Pipe

#Read in the sample data (small, medium and full sizes)
full_df = load_sample_data()

medium_df = CSV.read("data/data_for_testing/medium_dataset.csv",  DataFrame)
medium_df = select(medium_df, All())   #convert all columns back to regular vectors

small_df = CSV.read("data/data_for_testing/small_dataset.csv",  DataFrame)
small_df = select(small_df, All())   #convert all columns back to regular vectors


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

# Benchmark the mean_spatial_beta_div function
spatial_beta_div_1 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true) samples=100 evals=1
spatial_beta_div_2 = @benchmark spatial_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1
spatial_beta_div_3 = @benchmark spatial_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1

# Benchmark the mean_temporal_beta_div function
temporal_beta_div_1 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=true) samples=100 evals=1 seconds=1800
temporal_beta_div_2 = @benchmark temporal_beta_div(df.Abundance, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1
temporal_beta_div_3 = @benchmark temporal_beta_div(df.Presence, df.Sampling_date_order, df.plot, df.Species; quant=false) samples=100 evals=1

## Benchmark the DNCI functions
clustering_result = create_clusters(df.Sampling_date_order, 
                                            df.Latitude, 
                                            df.Longitude,                                      
                                            df.plot, 
                                            df.Species, 
                                            df.Presence)

# Save the wrangled data to CSV files for the R benchmarks
group_df = @pipe df |>
                filter(row -> row[:Sampling_date_order] == 60, _) |>
                select(_, [:plot, :Species, :Presence]) |>
                innerjoin(_, clustering_result[60], on = [:plot => :Site, :Species], makeunique = true)|>
                select(_, [:plot, :Species, :Presence, :Group]) |>
                unstack(_, :Species, :Presence, fill=0)

#Save the data to CSV files for the R benchmarks
CSV.write("data/data_for_testing/comm_full_df.csv", group_df)

# Prepare the community matrix for the DNCI_multigroup function
comm= @pipe group_df |>
            select(_, Not([:plot,:Group]))|>
               Matrix(_)  #convert the dataframe to a matrix

# Benchmark the DNCI_multigroup function
DNCI_multigroup_result = @benchmark DNCI_multigroup(comm, group_df.Group, 100; Nperm_count = false) samples=100 evals=1 seconds=1800
save_object("benchmarks/result/DNCI_multigroup_full_df_result.jld2", DNCI_multigroup_result)
                                                                        
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

#Save all the time samples to a DataFrame
# Initialize empty DataFrame
all_time_full = DataFrame()

test_case_list = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                "DNCI_multigroup_result", "prop_patches_result", "CV_meta_result", 
                "hypervolume_det_result", "hypervolume_dis_result"]

for testcase in test_case_list
    # Convert string to symbol to access the variable
    test_case_sym = Symbol(testcase)
    
    # Get the actual times (assuming they're in a field called 'time')
    # Multiply by 1000 to convert to milliseconds
    times = eval(test_case_sym).times ./ 1e6
    
    # Create DataFrame for this test case
    data = DataFrame(
        TestCase = fill(testcase, length(times)),
        Time = times
    )
    
    # Append to the full DataFrame
    all_time_full = vcat(all_time_full, data)
end
# Save the time samples to a CSV file
CSV.write("benchmarks/result/all_time_full_julia.csv", all_time_full)

# Save the results to a DataFrame
benchmark_result_full_df = DataFrame(TestCase = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                        "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                        "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                        "DNCI_multigroup_result", 
                                        "prop_patches_result", 
                                        "CV_meta_result", 
                                        "hypervolume_det_result", 
                                        "hypervolume_dis_result"],
                                        Time_minimum = [minimum(beta_diversity_1.times), minimum(beta_diversity_2.times), minimum(beta_diversity_3.times),
                                        minimum(spatial_beta_div_1.times), minimum(spatial_beta_div_2.times), minimum(spatial_beta_div_3.times),
                                        minimum(temporal_beta_div_1.times), minimum(temporal_beta_div_2.times), minimum(temporal_beta_div_3.times),
                                        minimum(DNCI_multigroup_result.times),
                                        minimum(prop_patches_result.times), 
                                        minimum(CV_meta_result.times), 
                                        minimum(hypervolume_det_result.times), minimum(hypervolume_dis_result.times)],
                                        Time_median = [median(beta_diversity_1.times), median(beta_diversity_2.times), median(beta_diversity_3.times),
                                        median(spatial_beta_div_1.times), median(spatial_beta_div_2.times), median(spatial_beta_div_3.times),
                                        median(temporal_beta_div_1.times), median(temporal_beta_div_2.times), median(temporal_beta_div_3.times),
                                        median(DNCI_multigroup_result.times), 
                                        median(prop_patches_result.times), 
                                        median(CV_meta_result.times), 
                                        median(hypervolume_det_result.times), median(hypervolume_dis_result.times)],
                                        Time_mean = [mean(beta_diversity_1.times), mean(beta_diversity_2.times), mean(beta_diversity_3.times),
                                        mean(spatial_beta_div_1.times), mean(spatial_beta_div_2.times), mean(spatial_beta_div_3.times),
                                        mean(temporal_beta_div_1.times), mean(temporal_beta_div_2.times), mean(temporal_beta_div_3.times),
                                        mean(DNCI_multigroup_result.times), 
                                        mean(prop_patches_result.times), 
                                        mean(CV_meta_result.times), 
                                        mean(hypervolume_det_result.times), mean(hypervolume_dis_result.times)],
                                        Time_maximum = [maximum(beta_diversity_1.times), maximum(beta_diversity_2.times), maximum(beta_diversity_3.times),
                                        maximum(spatial_beta_div_1.times), maximum(spatial_beta_div_2.times), maximum(spatial_beta_div_3.times),
                                        maximum(temporal_beta_div_1.times), maximum(temporal_beta_div_2.times), maximum(temporal_beta_div_3.times),
                                        maximum(DNCI_multigroup_result.times), 
                                        maximum(prop_patches_result.times), 
                                        maximum(CV_meta_result.times), 
                                        maximum(hypervolume_det_result.times), maximum(hypervolume_dis_result.times)],
                                        Time_std = [std(beta_diversity_1.times), std(beta_diversity_2.times), std(beta_diversity_3.times),
                                        std(spatial_beta_div_1.times), std(spatial_beta_div_2.times), std(spatial_beta_div_3.times),
                                        std(temporal_beta_div_1.times), std(temporal_beta_div_2.times), std(temporal_beta_div_3.times),
                                        std(DNCI_multigroup_result.times),
                                        std(prop_patches_result.times), 
                                        std(CV_meta_result.times), 
                                        std(hypervolume_det_result.times), std(hypervolume_dis_result.times)],
                                        memory = [beta_diversity_1.memory, beta_diversity_2.memory, beta_diversity_3.memory,
                                        spatial_beta_div_1.memory, spatial_beta_div_2.memory, spatial_beta_div_3.memory,
                                        temporal_beta_div_1.memory, temporal_beta_div_2.memory, temporal_beta_div_3.memory,
                                        DNCI_multigroup_result.memory, 
                                        prop_patches_result.memory, 
                                        CV_meta_result.memory, 
                                        hypervolume_det_result.memory, hypervolume_dis_result.memory])   

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

## Benchmark the DNCI functions
clustering_result = create_clusters(df.Sampling_date_order, 
                                            df.Latitude, 
                                            df.Longitude,                                      
                                            df.plot, 
                                            df.Species, 
                                            df.Presence)

# Save the wrangled data to CSV files for the R benchmarks
group_df = @pipe df |>
                filter(row -> row[:Sampling_date_order] == 60, _) |>
                select(_, [:plot, :Species, :Presence]) |>
                innerjoin(_, clustering_result[60], on = [:plot => :Site, :Species], makeunique = true)|>
                select(_, [:plot, :Species, :Presence, :Group]) |>
                unstack(_, :Species, :Presence, fill=0)

#Save the data to CSV files for the R benchmarks
CSV.write("data/data_for_testing/comm_medium_df.csv", group_df)

# Prepare the community matrix for the DNCI_multigroup function
comm= @pipe group_df |>
            select(_, Not([:plot,:Group]))|>
               Matrix(_)  #convert the dataframe to a matrix

# Benchmark the DNCI_multigroup function
DNCI_multigroup_result = @benchmark DNCI_multigroup(comm, group_df.Group, 100; Nperm_count = false) samples=100 evals=1 seconds=1800
save_object("benchmarks/result/DNCI_multigroup_medium_df_result.jld2", DNCI_multigroup_result)

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

#Save all the time samples to a DataFrame
# Initialize empty DataFrame
all_time_medium = DataFrame()

test_case_list = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                "DNCI_multigroup_result", "prop_patches_result", "CV_meta_result", 
                "hypervolume_det_result", "hypervolume_dis_result"]

for testcase in test_case_list
    # Convert string to symbol to access the variable
    test_case_sym = Symbol(testcase)
    
    # Get the actual times (assuming they're in a field called 'time')
    # Multiply by 1000 to convert to milliseconds
    times = eval(test_case_sym).times ./ 1e6
    
    # Create DataFrame for this test case
    data = DataFrame(
        TestCase = fill(testcase, length(times)),
        Time = times
    )
    
    # Append to the full DataFrame
    all_time_medium = vcat(all_time_medium, data)
end
# Save the time samples to a CSV file
CSV.write("benchmarks/result/all_time_medium_julia.csv", all_time_medium)

                    # Save the results to a DataFrame
benchmark_result_medium_df = DataFrame(TestCase = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                                "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                                "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                                "DNCI_multigroup_result", 
                                                "prop_patches_result", 
                                                "CV_meta_result", 
                                                "hypervolume_det_result", 
                                                "hypervolume_dis_result"],
                                                Time_minimum = [minimum(beta_diversity_1.times), minimum(beta_diversity_2.times), minimum(beta_diversity_3.times),
                                                minimum(spatial_beta_div_1.times), minimum(spatial_beta_div_2.times), minimum(spatial_beta_div_3.times),
                                                minimum(temporal_beta_div_1.times), minimum(temporal_beta_div_2.times), minimum(temporal_beta_div_3.times),
                                                minimum(DNCI_multigroup_result.times),
                                                minimum(prop_patches_result.times), 
                                                minimum(CV_meta_result.times), 
                                                minimum(hypervolume_det_result.times), minimum(hypervolume_dis_result.times)],
                                                Time_median = [median(beta_diversity_1.times), median(beta_diversity_2.times), median(beta_diversity_3.times),
                                                median(spatial_beta_div_1.times), median(spatial_beta_div_2.times), median(spatial_beta_div_3.times),
                                                median(temporal_beta_div_1.times), median(temporal_beta_div_2.times), median(temporal_beta_div_3.times),
                                                median(DNCI_multigroup_result.times), 
                                                median(prop_patches_result.times), 
                                                median(CV_meta_result.times), 
                                                median(hypervolume_det_result.times), median(hypervolume_dis_result.times)],
                                                Time_mean = [mean(beta_diversity_1.times), mean(beta_diversity_2.times), mean(beta_diversity_3.times),
                                                mean(spatial_beta_div_1.times), mean(spatial_beta_div_2.times), mean(spatial_beta_div_3.times),
                                                mean(temporal_beta_div_1.times), mean(temporal_beta_div_2.times), mean(temporal_beta_div_3.times),
                                                mean(DNCI_multigroup_result.times), 
                                                mean(prop_patches_result.times), 
                                                mean(CV_meta_result.times), 
                                                mean(hypervolume_det_result.times), mean(hypervolume_dis_result.times)],
                                                Time_maximum = [maximum(beta_diversity_1.times), maximum(beta_diversity_2.times), maximum(beta_diversity_3.times),
                                                maximum(spatial_beta_div_1.times), maximum(spatial_beta_div_2.times), maximum(spatial_beta_div_3.times),
                                                maximum(temporal_beta_div_1.times), maximum(temporal_beta_div_2.times), maximum(temporal_beta_div_3.times),
                                                maximum(DNCI_multigroup_result.times), 
                                                maximum(prop_patches_result.times), 
                                                maximum(CV_meta_result.times), 
                                                maximum(hypervolume_det_result.times), maximum(hypervolume_dis_result.times)],
                                                Time_std = [std(beta_diversity_1.times), std(beta_diversity_2.times), std(beta_diversity_3.times),
                                                std(spatial_beta_div_1.times), std(spatial_beta_div_2.times), std(spatial_beta_div_3.times),
                                                std(temporal_beta_div_1.times), std(temporal_beta_div_2.times), std(temporal_beta_div_3.times),
                                                std(DNCI_multigroup_result.times),
                                                std(prop_patches_result.times), 
                                                std(CV_meta_result.times), 
                                                std(hypervolume_det_result.times), std(hypervolume_dis_result.times)],
                                                memory = [beta_diversity_1.memory, beta_diversity_2.memory, beta_diversity_3.memory,
                                                spatial_beta_div_1.memory, spatial_beta_div_2.memory, spatial_beta_div_3.memory,
                                                temporal_beta_div_1.memory, temporal_beta_div_2.memory, temporal_beta_div_3.memory,
                                                DNCI_multigroup_result.memory, 
                                                prop_patches_result.memory, 
                                                CV_meta_result.memory, 
                                                hypervolume_det_result.memory, hypervolume_dis_result.memory])   

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

## Benchmark the DNCI functions
clustering_result = create_clusters(df.Sampling_date_order, 
                                            df.Latitude, 
                                            df.Longitude,                                      
                                            df.plot, 
                                            df.Species, 
                                            df.Presence)

# Save the wrangled data to CSV files for the R benchmarks
group_df = @pipe df |>
                filter(row -> row[:Sampling_date_order] == 60, _) |>
                select(_, [:plot, :Species, :Presence]) |>
                innerjoin(_, clustering_result[60], on = [:plot => :Site, :Species], makeunique = true)|>
                select(_, [:plot, :Species, :Presence, :Group]) |>
                unstack(_, :Species, :Presence, fill=0)

#Save the data to CSV files for the R benchmarks
CSV.write("data/data_for_testing/comm_small_df.csv", group_df)

# Prepare the community matrix for the DNCI_multigroup function
comm= @pipe group_df |>
            select(_, Not([:plot,:Group]))|>
               Matrix(_)  #convert the dataframe to a matrix

# Benchmark the DNCI_multigroup function
DNCI_multigroup_result = @benchmark DNCI_multigroup(comm, group_df.Group, 100; Nperm_count = false) samples=100 evals=1 seconds=1800 
save_object("benchmarks/result/DNCI_multigroup_small_df_result.jld2", DNCI_multigroup_result)
                                                                       
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

#Save all the time samples to a DataFrame
# Initialize empty DataFrame
all_time_small = DataFrame()

test_case_list = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                "DNCI_multigroup_result", "prop_patches_result", "CV_meta_result", 
                "hypervolume_det_result", "hypervolume_dis_result"]

for testcase in test_case_list
    # Convert string to symbol to access the variable
    test_case_sym = Symbol(testcase)
    
    # Get the actual times (assuming they're in a field called 'time')
    # Multiply by 1000 to convert to milliseconds
    times = eval(test_case_sym).times ./ 1e6
    
    # Create DataFrame for this test case
    data = DataFrame(
        TestCase = fill(testcase, length(times)),
        Time = times
    )
    
    # Append to the full DataFrame
    all_time_small = vcat(all_time_small, data)
end
# Save the time samples to a CSV file
CSV.write("benchmarks/result/all_time_small_julia.csv", all_time_small)

# Save the results to a DataFrame
benchmark_result_small_df = DataFrame(TestCase = ["beta_diversity_1", "beta_diversity_2", "beta_diversity_3",
                                                "spatial_beta_div_1", "spatial_beta_div_2", "spatial_beta_div_3",
                                                "temporal_beta_div_1", "temporal_beta_div_2", "temporal_beta_div_3",
                                                "DNCI_multigroup_result", 
                                                "prop_patches_result", 
                                                "CV_meta_result", 
                                                "hypervolume_det_result", 
                                                "hypervolume_dis_result"],
                                                Time_minimum = [minimum(beta_diversity_1.times), minimum(beta_diversity_2.times), minimum(beta_diversity_3.times),
                                                minimum(spatial_beta_div_1.times), minimum(spatial_beta_div_2.times), minimum(spatial_beta_div_3.times),
                                                minimum(temporal_beta_div_1.times), minimum(temporal_beta_div_2.times), minimum(temporal_beta_div_3.times),
                                                minimum(DNCI_multigroup_result.times),
                                                minimum(prop_patches_result.times), 
                                                minimum(CV_meta_result.times), 
                                                minimum(hypervolume_det_result.times), minimum(hypervolume_dis_result.times)],
                                                Time_median = [median(beta_diversity_1.times), median(beta_diversity_2.times), median(beta_diversity_3.times),
                                                median(spatial_beta_div_1.times), median(spatial_beta_div_2.times), median(spatial_beta_div_3.times),
                                                median(temporal_beta_div_1.times), median(temporal_beta_div_2.times), median(temporal_beta_div_3.times),
                                                median(DNCI_multigroup_result.times), 
                                                median(prop_patches_result.times), 
                                                median(CV_meta_result.times), 
                                                median(hypervolume_det_result.times), median(hypervolume_dis_result.times)],
                                                Time_mean = [mean(beta_diversity_1.times), mean(beta_diversity_2.times), mean(beta_diversity_3.times),
                                                mean(spatial_beta_div_1.times), mean(spatial_beta_div_2.times), mean(spatial_beta_div_3.times),
                                                mean(temporal_beta_div_1.times), mean(temporal_beta_div_2.times), mean(temporal_beta_div_3.times),
                                                mean(DNCI_multigroup_result.times), 
                                                mean(prop_patches_result.times), 
                                                mean(CV_meta_result.times), 
                                                mean(hypervolume_det_result.times), mean(hypervolume_dis_result.times)],
                                                Time_maximum = [maximum(beta_diversity_1.times), maximum(beta_diversity_2.times), maximum(beta_diversity_3.times),
                                                maximum(spatial_beta_div_1.times), maximum(spatial_beta_div_2.times), maximum(spatial_beta_div_3.times),
                                                maximum(temporal_beta_div_1.times), maximum(temporal_beta_div_2.times), maximum(temporal_beta_div_3.times),
                                                maximum(DNCI_multigroup_result.times), 
                                                maximum(prop_patches_result.times), 
                                                maximum(CV_meta_result.times), 
                                                maximum(hypervolume_det_result.times), maximum(hypervolume_dis_result.times)],
                                                Time_std = [std(beta_diversity_1.times), std(beta_diversity_2.times), std(beta_diversity_3.times),
                                                std(spatial_beta_div_1.times), std(spatial_beta_div_2.times), std(spatial_beta_div_3.times),
                                                std(temporal_beta_div_1.times), std(temporal_beta_div_2.times), std(temporal_beta_div_3.times),
                                                std(DNCI_multigroup_result.times),
                                                std(prop_patches_result.times), 
                                                std(CV_meta_result.times), 
                                                std(hypervolume_det_result.times), std(hypervolume_dis_result.times)],
                                                memory = [beta_diversity_1.memory, beta_diversity_2.memory, beta_diversity_3.memory,
                                                spatial_beta_div_1.memory, spatial_beta_div_2.memory, spatial_beta_div_3.memory,
                                                temporal_beta_div_1.memory, temporal_beta_div_2.memory, temporal_beta_div_3.memory,
                                                DNCI_multigroup_result.memory, 
                                                prop_patches_result.memory, 
                                                CV_meta_result.memory, 
                                                hypervolume_det_result.memory, hypervolume_dis_result.memory])   
# Save the results to a CSV file
# Divide the time by 1e6 to convert from nanoseconds to milliseconds
benchmark_result_small_df[:, 2:6] .= benchmark_result_small_df[:, 2:6] ./ 1e6
# Divide the memory by 1024^2 to convert from bytes to mebibytes (Mib)
# convert column type to Float64 before division
benchmark_result_small_df[!, 7] = convert.(Float64, benchmark_result_small_df[:, 7])
# perform the division
benchmark_result_small_df[:, 7] .= benchmark_result_small_df[:, 7] ./ (1024^2)

CSV.write("benchmarks/result/benchmark_result_small_df_julia.csv", benchmark_result_small_df)