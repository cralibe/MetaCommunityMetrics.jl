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
    sample_df = CSV.read(joinpath(@__DIR__, "..", "data", "rodent_abundance_data.csv"), DataFrame, skipto=2)#read in the sample data
    sample_df = stack(sample_df, Not(:Sampling_date_order, :Year, :Month, :Day, :plot), variable_name = :Species, value_name = :Abundance) #stack the dataframe

    #A dataframe that idetifies the persisting species (The species that existed at least once in all time steps)
    persist_sp_df=
    @pipe sample_df |> 
    groupby(_, :Species) |> #group the dataframe by every species
    combine(_, :Abundance => sum => :Total_Abundance) |> #calculate the total abundance of every spesice within the whole sampling period
    filter(row -> row[:Total_Abundance] > 0, _) #select rows that has a total abundance >0.

    #A dataframe that idetifies non empty patch at every time step
    non_empty_site = 
    @pipe sample_df[in(persist_sp_df.Species).(sample_df.Species), :] |>
    groupby(_, [:Sampling_date_order, :plot]) |>
    combine(_, :Abundance => sum => :Total_Abundance) |>
    filter(row -> row[:Total_Abundance] > 0, _)

    #The data that contains the abundace and presence-absence of persisting species and non empty site at every time step
    metacomm_df=
    @pipe sample_df|>
    _[in(persist_sp_df.Species).(_.Species), :]|> #keeping species with a total abundance >0 across all patches all time steps. 
    innerjoin(_, non_empty_site, on = [:Sampling_date_order, :plot]) |> #keeping non empty site at every time step. 
    select(_, Not(:Total_Abundance)) |> #remove the Total_Abundance column to avoid confusion
    transform(_, :Abundance => ByRow(x->ifelse.(x> 0, 1, 0)) => :Presence) #create a presence-absence column from the abundance column

    #Simulating coordinates for the sample data. This step is not nessary if: (1) you have the coordinates data, (2) you don't need to caluculating DNCI, or (3) you are doing the groupings for the DNCI by yourself.
    #   Set grid dimensions
    num_rows = 4
    num_columns = 6
    num_plots = 24  # Total number of plots
    # Define base coordinates (for example purposes, these are arbitrary)
    base_latitude = 35.0
    base_longitude = -110.0
    # Define distance between plots in degrees (arbitrary small values for demonstration)
    lat_spacing = 0.5
    lon_spacing = 0.5
    # Initialize DataFrame
    patch_coord_df = DataFrame(plot = Int[], Latitude = Float64[], Longitude = Float64[])

    # Generate grid coordinates
    plot_count = 1
    for i in 1:num_rows
        for j in 1:num_columns
            if plot_count > num_plots
                break
            end
            lat = base_latitude + (i - 1) * lat_spacing
            lon = base_longitude + (j - 1) * lon_spacing
            push!(patch_coord_df, (plot_count, lat, lon))
            plot_count += 1
        end
    end

    #Adding the coordinates to the metacomm_df
    metacomm_df = @pipe metacomm_df |>
    innerjoin(_, patch_coord_df, on = :plot) #joining the metacomm_df with the patch_coord_df

    #Preparing the matric for calculating beta diveristy
    #matrix with the species abundance
    sample_matrix_abundance = @pipe metacomm_df |> 
    select(_, Not(:Presence)) |>
    filter(row -> row[:Sampling_date_order] == 50, _) |> # filter rows where Sampling_date_order == 50
    unstack(_, :Species, :Abundance, fill=0) |> #convert it back to the wide format 
    select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> # remove the first four columns
    Matrix(_) |> #convert the dataframe to a matrix
    _[:, sum(_, dims=1)[1, :] .!= 0] # Remove columns with sum of zero
    #matrix with the species presence-absence data
    sample_matrix_presence = @pipe metacomm_df |> 
    select(_, Not(:Abundance)) |>
    filter(row -> row[:Sampling_date_order] == 50, _) |># filter rows where Sampling_date_order == 50
    unstack(_, :Species, :Presence, fill=0) |> #convert it back to the wide format 
    select(_, Not(:Year, :Month, :Day, :Sampling_date_order, :plot, :Longitude, :Latitude)) |> # remove the first four columns
    Matrix(_) |> #convert the dataframe to a matrix      
    _[:, sum(_, dims=1)[1, :] .!= 0] # Remove columns with sum of zero

    #Preparing the data for the DNCI analysis
    total_presence_df=@pipe metacomm_df|>
    groupby(_,[:Species,:Sampling_date_order])|>
    combine(_,:Presence=>sum=>:Total_Presence)
    #filter(row -> row[:Total_Presence] <= 1, _)|> 

    #Remove singletons 
    total_richness_df= @pipe metacomm_df|>
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
    @test isapprox(mean_spatial_beta_div(metacomm_df.Abundance, 
        metacomm_df.Sampling_date_order, 
        metacomm_df.plot, 
        metacomm_df.Species; quant=true), DataFrame( mean_spatial_BDtotal = 0.3538123660657936, 
                                                        mean_spatial_Repl = 0.1685842743545843, 
                                                        mean_spatial_RichDif = 0.18522809171120944),
                                                        atol = 1e-8)

    # Test the mean_temporal_beta_div function
    @test isapprox(mean_temporal_beta_div(metacomm_df.Abundance, 
        metacomm_df.Sampling_date_order, 
        metacomm_df.plot, 
        metacomm_df.Species; quant=true), DataFrame( mean_temporal_BDtotal = 0.3718602494286472, 
                                                        mean_temporal_Repl = 0.15241610198410346, 
                                                        mean_temporal_RichDif = 0.21944414744454374),
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
    comm= @pipe metacomm_df|>
    leftjoin(_,  total_presence_df, on = [:Species, :Sampling_date_order], makeunique = true) |>
    filter(row -> row[:Total_Presence] > 1, _) |>
    filter(row -> row[:Sampling_date_order] == 50, _) |>
    select(_, [:plot, :Species, :Presence]) |>
    unstack(_, :Species, :Presence, fill=0) |>
    select(_, Not(:plot)) |>
    Matrix(_)
    
    presence_df = @pipe metacomm_df|>
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
    @test isapprox(niche_overlap(metacomm_df.Abundance, 
                        metacomm_df.Species, 
                        metacomm_df.plot, 
                        metacomm_df.Sampling_date_order), DataFrame(mean_niche_overlap_index = 0.8277387279283689,
                                                                        min_niche_overlap_index = 0.5918363543,
                                                                        max_niche_overlap_index = 1.0),
                                                                        atol = 1e-8)      
                                                                        
    # Test the prop_patches function
    @test isapprox(prop_patches(metacomm_df.Presence, metacomm_df.Species, metacomm_df.plot), DataFrame(mean_prop_patches = 0.7346491228070174,
                                                                                                    min_prop_patches = 0.08333333333333333,
                                                                                                    max_prop_patches = 1.0),
                                                                                                    atol = 1e-8)      


    # Test the CV_meta function
    CV_test_df = @pipe metacomm_df |>
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
