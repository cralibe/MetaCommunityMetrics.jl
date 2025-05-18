# test/runtests.jl
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
    select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Presence)) |> # remove the unused four columns
    filter(row -> row[:Sampling_date_order] == 1, _) |> # filter rows where Sampling_date_order == 1
    unstack(_, :Species, :Abundance) |> #convert it back to the wide format 
    select(_, Not(:Sampling_date_order, :plot)) |> # remove the unused columns
    Matrix(_) |> #convert the dataframe to a matrix
    _[:, sum(_, dims=1)[1, :] .!= 0] |> # Remove columns with sum of zero
    _[sum(_, dims=2)[:, 1] .!= 0,:] # Remove rows with sum of zero


    #matrix with the species presence-absence data
    sample_matrix_presence = @pipe df |> 
    select(_, Not(:Year, :Month, :Day, :Longitude, :Latitude, :normalized_temperature, :normalized_precipitation, :Abundance)) |> # remove the unused four columns
    filter(row -> row[:Sampling_date_order] == 1, _) |># filter rows where Sampling_date_order == 1
    unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
    select(_, Not(:Sampling_date_order, :plot)) |> # remove the unused columns
    Matrix(_) |> #convert the dataframe to a matrix      
    _[:, sum(_, dims=1)[1, :] .!= 0] |> # Remove columns with sum of zero
    _[sum(_, dims=2)[:, 1] .!= 0,:] # Remove rows with sum of zero

    #Preparing the data for the DNCI analysis
    total_presence_df=@pipe df|>
    groupby(_,[:Species,:Sampling_date_order])|>
    combine(_,:Presence=>sum=>:Total_Presence)

    #Remove singletons 
    total_richness_df= @pipe df|>
    leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
    filter(row -> row[:Total_Presence] > 1, _)|> 
    groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
    combine(_,:Presence=>sum=>:Total_Richness)
    

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
    cluster_result=create_clusters(total_richness_df.Sampling_date_order, 
        total_richness_df.Latitude, 
        total_richness_df.Longitude, 
        total_richness_df.plot,
        total_richness_df.Total_Richness) 

    @test cluster_result[1] == DataFrame(
        Time = repeat([1], 15),
        Latitude = [35.0, 35.0, 35.5, 35.5, 35.5, 36.0, 36.0, 36.5, 35.0, 36.0, 36.0, 36.5, 36.5, 35.0, 36.5],
        Longitude = [-110.0, -109.5, -109.5, -109.0, -108.0, -109.5, -108.0, -108.5, -107.5, -110.0, -109.0, -109.5, -109.0, -108.0, -108.0],
        Patch = [1, 2, 8, 9, 11, 14, 17, 22, 6, 13, 15, 20, 21, 5, 23],
        Total_Richness = [1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 1],
        Group = [1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2])
    # Test the plot_clusters function
    #p = plot_clusters(cluster_result[1].Latitude, cluster_result[1].Longitude, cluster_result[1].Group; output_file="clusters.svg")

    # Test the DNCI_multigroup function
    comm= @pipe df|>
    filter(row -> row[:Sampling_date_order] == 50, _) |>
    select(_, [:plot, :Species, :Presence]) |>
    unstack(_, :Species, :Presence, fill=0) |>
    select(_, Not(:plot)) |>
    Matrix(_)
    
    presence_df = @pipe df|>
    leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
    filter(row -> row[:Total_Presence] > 1, _) |>
    filter(row -> row[:Sampling_date_order] == 50, _) |>
    select(_, [:plot, :Species, :Presence]) |>
    unstack(_, :Species, :Presence, fill=0) |>
    leftjoin(_, cluster_result[50], on = [:plot => :Patch], makeunique = true) 
    
    #Random.seed!(1234)
    #DNCI_result = DNCI_multigroup(comm, presence_df.Group, 1000; count = false)
    #@test isapprox(DNCI_result, DataFrame(
    #    group1 = [1, 1, 2],
    #    group2 = [2, 3, 3],
    #    DNCI = [ -2.96216, -2.93839, -1.45788],
    #    CI_DNCI = [1.81732, 2.24484, 2.27338],
    #    S_DNCI = [0.908661, 1.12242, 1.13669]),
    #atol = 1e-4)

    # Test the niche_overlap function
    @test isapprox(niche_overlap(df.Abundance, 
                        df.Species, 
                        df.plot, 
                        df.Sampling_date_order), DataFrame(mean_niche_overlap_index = 0.8277387279283689,
                                                                        min_niche_overlap_index = 0.5918363543,
                                                                        max_niche_overlap_index = 1.0),
                                                                        atol = 1e-8)      
                                                                        
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
    select(_, [:normalized_temperature, :normalized_precipitation])

    @test isapprox(MVNH_det(data; var_names=["Temperature", "Precipitation"]), DataFrame(total = 1.15268,
                                            correlation = 0.999732,
                                            Temperature = 0.962495,
                                            Precipitation =  1.19792),
                                            atol = 1e-5) 
                                            
    # Test the MVNH_dissimilarity function
    data_2 = @pipe df |> 
            filter(row -> row[:Presence] > 0, _) |>
            filter(row -> row[:Species] == "SH", _) |>
            select(_, [:normalized_temperature, :normalized_precipitation])

    result_df = MVNH_dissimilarity(data, data_2; var_names = ["Temperature", "Precipitation"])
    expected_df = DataFrame(
                metric = ["Bhattacharyya_distance", "Mahalanobis_distance", "Determinant_ratio"],
                total = [0.0739592, 0.0322218, 0.0417374],
                correlation = [0.000155556,-0.000360636, 0.000516192],
                Temperature = [0.0478506, 0.0269479, 0.0209027],
                Precipitation = [0.0259531, 0.00563454, 0.0203185])
    
    @test result_df.metric == expected_df.metric
    @test isapprox(result_df.total, expected_df.total, atol=1e-5)
    @test isapprox(result_df.correlation, expected_df.correlation, atol=1e-5)
    @test isapprox(result_df.Temperature, expected_df.Temperature, atol=1e-5)
    @test isapprox(result_df.Precipitation, expected_df.Precipitation, atol=1e-5)                            
    
    # Test the average_MVNH_det function    
    data = @pipe df |> 
    select(_, [:normalized_temperature, :normalized_precipitation])

    @test isapprox(average_MVNH_det(data, df.Presence, String.(df.Species); var_names = ["Temperature", "Precipitation"]), 
    1.2103765096417536, atol = 1e-5)

    # Test the average_MVNH_dissimilarity function
    @test isapprox(average_MVNH_dissimilarity(data, df.Presence, String.(df.Species); var_names = ["Temperature", "Precipitation"]), 
    0.03059942936454443, atol = 1e-5)

end

