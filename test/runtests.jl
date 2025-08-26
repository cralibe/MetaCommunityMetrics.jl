# test/runtests.jl
#This is a script to test the functions in the MetaCommunityMetrics.jl package
using Test
using MetaCommunityMetrics
using Random
using Pipe: @pipe
using CSV
using DataFrames
using .MetaCommunityMetrics.Internal 

@testset "MetaCommunityMetrics.jl" begin

    #Example Data
    #Load sample data included in the package
    df = load_sample_data()

    #Preparing the matric for calculating beta diveristy
    #matrix with the species abundance
    sample_matrix_abundance = @pipe df |> 
    select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :standardized_temperature, :standardized_precipitation, :Presence)) |> # remove the unused four columns
    filter(row -> row[:Sampling_date_order] == 1, _) |> # filter rows where Sampling_date_order == 1
    unstack(_, :Species, :Abundance) |> #convert it back to the wide format 
    select(_, Not(:Sampling_date_order, :plot)) |> # remove the unused columns
    Matrix(_) |> #convert the dataframe to a matrix
    _[:, sum(_, dims=1)[1, :] .!= 0] |> # Remove columns with sum of zero
    _[sum(_, dims=2)[:, 1] .!= 0,:] # Remove rows with sum of zero


    #matrix with the species presence-absence data
    sample_matrix_presence = @pipe df |> 
    select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :standardized_temperature, :standardized_precipitation, :Abundance)) |> # remove the unused four columns
    filter(row -> row[:Sampling_date_order] == 1, _) |># filter rows where Sampling_date_order == 1
    unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
    select(_, Not(:Sampling_date_order, :plot)) |> # remove the unused columns
    Matrix(_) |> #convert the dataframe to a matrix      
    _[:, sum(_, dims=1)[1, :] .!= 0] |> # Remove columns with sum of zero
    _[sum(_, dims=2)[:, 1] .!= 0,:] # Remove rows with sum of zero

    #Preparing the data for the DNCI analysis
    #Remove singletons and empty sites
    total_presence_df=@pipe df|>
                    groupby(_,[:Species,:Sampling_date_order])|>
                    combine(_,:Presence=>sum=>:Total_Presence) |>
                    filter(row -> row[:Total_Presence] > 0, _) |>
                    select(_, [:Species, :Sampling_date_order])

    non_empty_site_df = @pipe df|>
                    innerjoin(_,  total_presence_df, on = [:Species], makeunique = true)|>
                    groupby(_, [:plot]) |>
                    combine(_, :Presence=>sum=>:Total_N) |>
                    filter(row -> row[:Total_N] > 0, _) |>
                    select(_, [:plot])

    ubiquitous_species_df = @pipe df |>
                    groupby(_, [:Species, :Sampling_date_order]) |>
                    combine(_, :Presence => sum => :Total_Presence) |>
                    filter(row -> row[:Total_Presence] < length(unique(non_empty_site_df.plot)), _) |>
                    select(_, [:Species, :Sampling_date_order])

    filtered_df = @pipe df|>
                    innerjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
                    innerjoin(_,  non_empty_site_df, on = [:plot], makeunique = true) |>
                    innerjoin(_, ubiquitous_species_df, on = [:Species, :Sampling_date_order], makeunique = true) 
    

    # Test the beta_diversity function
    ## For abundance data
    @test isapprox(beta_diversity(sample_matrix_abundance; quant=true), DataFrame(BDtotal = 0.390317, 
                                                                            Repl = 0.2678, 
                                                                            RichDif = 0.122517),
                                                                            atol = 1e-5)
    ## For treating abundance data as presence-absence data using the function
    @test isapprox(beta_diversity(sample_matrix_abundance; quant=false), DataFrame(BDtotal = 0.357143, 
                                                                            Repl = 0.284127, 
                                                                            RichDif = 0.0730159),
                                                                            atol = 1e-5)
    ## For presence-absence data
    @test isapprox(beta_diversity(sample_matrix_presence; quant=false), DataFrame(BDtotal = 0.357143, 
                                                                            Repl = 0.284127, 
                                                                            RichDif = 0.0730159),
                                                                            atol = 1e-5)
    # Test the mean_spatial_beta_div function
    @test isapprox(spatial_beta_div(df.Abundance, 
        df.Sampling_date_order, 
        df.plot, 
        df.Species; quant=true), DataFrame( spatial_BDtotal = 0.26482181105238306, 
                                                        spatial_Repl = 0.12188231356454636, 
                                                        spatial_RichDif = 0.14293949748783713),
                                                        atol = 1e-8)

    # Test the mean_temporal_beta_div function
    @test isapprox(temporal_beta_div(df.Abundance, 
        df.Sampling_date_order, 
        df.plot, 
        df.Species; quant=true), DataFrame( temporal_BDtotal = 0.31122221160244523, 
                                                        temporal_Repl =  0.09954831011523949, 
                                                        temporal_RichDif = 0.21167390148720572),
                                                        atol = 1e-8)
    
    # Test the create_clusters function
    clustering_result = create_groups(filtered_df.Sampling_date_order, 
                                            filtered_df.Latitude, 
                                            filtered_df.Longitude,                                      
                                            filtered_df.plot, 
                                            filtered_df.Species, 
                                            filtered_df.Presence)

    @test size(clustering_result[3]) == (144, 7)
    @test names(clustering_result[3]) == ["Time", "Latitude", "Longitude", "Site", "Species", "Presence", "Group"]
    
    # Test specific rows
    @test clustering_result[3][1, :] == (Time=3, Latitude=35.0, Longitude=-110.0, Site=1, Species="DM", Presence=1, Group=1)
    @test clustering_result[3][25, :] == (Time=3, Latitude=35.0, Longitude=-110.0, Site=1, Species="DO", Presence=0, Group=1)
    @test clustering_result[3][54, :] == (Time=3, Latitude=35.5, Longitude=-108.0, Site=11, Species="OL", Presence=0, Group=2)
    @test clustering_result[3][72, :] == (Time=3, Latitude=36.5, Longitude=-107.5, Site=24, Species="OL", Presence=0, Group=2)
    
    # Test column properties
    @test all(clustering_result[3].Time .== 3)
    @test length(unique(clustering_result[3].Species)) == 6
    @test Set(unique(clustering_result[3].Group)) == Set([1, 2, 3])
    @test length(unique(clustering_result[3].Site)) == 24
    
    # Test species counts
    @test count(==("BA"), clustering_result[3].Species) == 0
    @test count(==("DM"), clustering_result[3].Species) == 24
    @test count(==("SH"), clustering_result[3].Species) == 0
    
    # Test presence patterns
    @test sum(clustering_result[3][clustering_result[3].Species .== "BA", :Presence]) == 0
    @test sum(clustering_result[3][clustering_result[3].Species .== "DM", :Presence]) == 7
    @test sum(clustering_result[3][clustering_result[3].Species .== "SH", :Presence]) == 0


    # Test the niche_overlap function
    @test isapprox(niche_overlap(df.Abundance, 
                        df.Species, 
                        df.plot, 
                        df.Sampling_date_order), DataFrame(mean_niche_overlap_index = 0.0923816,
                                                                        min_niche_overlap_index = 0.0,
                                                                        max_niche_overlap_index =  0.406837),
                                                                        atol = 1e-5)      
                                                                        
    # Test the prop_patches function
    @test isapprox(prop_patches(df.Presence, df.Species, df.plot), DataFrame(mean_prop_patches = 0.7346491228070174,
                                                                                                    min_prop_patches = 0.08333333333333333,
                                                                                                    max_prop_patches = 1.0),
                                                                                                    atol = 1e-8)      


    # Test the CV_meta function

    @test isapprox(CV_meta(df.Abundance, 
                    df.Sampling_date_order,
                    df.plot, 
                    df.Species), DataFrame(CV_s_l = 1.48859,
                                            CV_s_r = 0.944937,
                                            CV_c_l = 0.718266,
                                            CV_c_r = 0.580183),
                                            atol = 1e-4)      


    # Test the MVNH_det function
    data = @pipe df |> 
    filter(row -> row[:Presence] > 0, _) |>
    filter(row -> row[:Species] == "BA", _) |>
    select(_, [:standardized_temperature, :standardized_precipitation])

    @test isapprox(MVNH_det(data; var_names=["Temperature", "Precipitation"]), DataFrame(total = 1.15268,
                                            correlation = 0.999732,
                                            Temperature = 0.962495,
                                            Precipitation =  1.19792),
                                            atol = 1e-5) 
                                            
    # Test the MVNH_dissimilarity function
    data_2 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "SH", _) |>
            select(_, [:standardized_temperature, :standardized_precipitation])

    result_df = MVNH_dissimilarity(data, data_2; var_names = ["Temperature", "Precipitation"])
    expected_df = DataFrame(
                metric = ["Bhattacharyya_distance", "Mahalanobis_distance", "Determinant_ratio"],
                total = [0.00980771, 0.00664862, 0.00315908],
                correlation = [0.00015205 ,5.06232e-6, 0.000146988],
                Temperature = [0.00388058, 0.00234902, 0.00153156],
                Precipitation = [0.0057750, 0.00429454, 0.00148054])
    
    @test result_df.metric == expected_df.metric
    @test isapprox(result_df.total, expected_df.total, atol=1e-5)
    @test isapprox(result_df.correlation, expected_df.correlation, atol=1e-5)
    @test isapprox(result_df.Temperature, expected_df.Temperature, atol=1e-5)
    @test isapprox(result_df.Precipitation, expected_df.Precipitation, atol=1e-5)                            
    
    # Test the average_MVNH_det function    
    data = @pipe df |> 
    select(_, [:standardized_temperature, :standardized_precipitation])

    @test isapprox(average_MVNH_det(data, Vector{Int64}(df.Presence), df.Species; var_names = ["Temperature", "Precipitation"]), 
    1.2103765096417536, atol = 1e-5)

    # Test the average_MVNH_dissimilarity function
    @test isapprox(average_MVNH_dissimilarity(data, Vector{Int64}(df.Presence), df.Species; var_names = ["Temperature", "Precipitation"]), 
    0.03059942936454443, atol = 1e-5)

end

