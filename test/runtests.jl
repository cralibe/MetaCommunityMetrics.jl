# test/runtests.jl
using Test
using MetaCommunityMetrics
using Random
using Pipe: @pipe
using CSV
using DataFrames
using Plots
using .MetaCommunityMetrics.Internal 

@testset "MetaCommunityMetrics.jl" begin

    #Example Data
    #Load sample data included in the package
    df = load_sample_data()

    #Preparing the matric for calculating beta diveristy
    #matrix with the species abundance
    sample_matrix_abundance = @pipe df |> 
    select(_, Not(:Presence)) |>
    filter(row -> row[:Sampling_date_order] == 50, _) |> # filter rows where Sampling_date_order == 50
    unstack(_, :Species, :Abundance, fill=0) |> #convert it back to the wide format 
    select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> # remove the first four columns
    Matrix(_) |> #convert the dataframe to a matrix
    _[:, sum(_, dims=1)[1, :] .!= 0] # Remove columns with sum of zero
    #matrix with the species presence-absence data
    sample_matrix_presence = @pipe df |> 
    select(_, Not(:Abundance)) |>
    filter(row -> row[:Sampling_date_order] == 50, _) |># filter rows where Sampling_date_order == 50
    unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
    select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> # remove the first four columns
    Matrix(_) |> #convert the dataframe to a matrix      
    _[:, sum(_, dims=1)[1, :] .!= 0] # Remove columns with sum of zero

    #Preparing the data for the DNCI analysis
    total_presence_df=@pipe df|>
    groupby(_,[:Species,:Sampling_date_order])|>
    combine(_,:Presence=>sum=>:Total_Presence)
    #filter(row -> row[:Total_Presence] <= 1, _)|> 

    #Remove singletons 
    total_richness_df= @pipe df|>
    leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
    filter(row -> row[:Total_Presence] > 1, _)|> 
    groupby(_,[:plot,:Sampling_date_order,:Longitude, :Latitude])|>
    combine(_,:Presence=>sum=>:Total_Richness)
    

    # Test the beta_diversity function
    # For abundance data
    @test isapprox(beta_diversity(sample_matrix_abundance; quant=true), DataFrame(BDtotal = 0.413476082224629, 
                                                                            Repl = 0.23722016292087425, 
                                                                            RichDif = 0.176255919303755),
                                                                            atol = 1e-8)
    # For presence-absence data
    @test isapprox(beta_diversity(sample_matrix_abundance; quant=false), DataFrame(BDtotal = 0.38877781199915573, 
                                                                            Repl = 0.2419032219427474, 
                                                                            RichDif = 0.14687459005640824),
                                                                            atol = 1e-8)

    @test isapprox(beta_diversity(sample_matrix_presence; quant=false), DataFrame(BDtotal = 0.38877781199915573, 
                                                                            Repl = 0.2419032219427474, 
                                                                            RichDif = 0.14687459005640824),
                                                                            atol = 1e-8)
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
    p = plot_clusters(cluster_result[1].Latitude, cluster_result[1].Longitude, cluster_result[1].Group)
    @test p isa Plots.Plot

    # Test the DNCI_multigroup function
    comm= @pipe df|>
    leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
    filter(row -> row[:Total_Presence] > 1, _) |>
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
    
    DNCI_result = DNCI_multigroup(comm, 
    presence_df.Group; count = false) 
    
    #=@test isapprox(DNCI_result, DataFrame(
        group1 = [1, 1, 2],
        group2 = [2, 3, 3],
        DNCI = [-2.6398129641535233, -1.2754356031283491, -1.5357372549747799],
        CI_DNCI = [1.792594117783981, 2.646493538408896, 2.321281196527955],
        S_DNCI = [0.8962970588919905, 1.323246769204448, 1.1606405982639776]),
    atol = 1e-8)=#

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
    CV_test_df = @pipe df |>
    filter(row -> row[:Sampling_date_order] < 20, _)

    @test isapprox(CV_meta(CV_test_df.Abundance, 
                    CV_test_df.Sampling_date_order,
                    CV_test_df.plot, 
                    CV_test_df.Species), DataFrame(CV_s_l=1.05607236064782,
                                            CV_s_r=0.8133239050934045,
                                            CV_c_l=0.6561837274590009,
                                            CV_c_r=0.5372991135672851),
                                            atol = 1e-8)      
    
    # Test the CV_meta_simple function
    @test isapprox(CV_meta_simple(CV_test_df.Abundance, 
                            CV_test_df.Sampling_date_order,
                            CV_test_df.plot, 
                            CV_test_df.Species), DataFrame(CV_s_l= 1.2080844086744866,
                                            CV_s_r= 0.9188619724139119,
                                            CV_c_l= 0.8012160438000538,
                                            CV_c_r= 0.6650823635388619),
                                            atol = 1e-8)      

end
