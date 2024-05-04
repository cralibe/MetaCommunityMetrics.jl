using Pkg
Pkg.activate(".")

using Distributions
using Plots
using LinearAlgebra
using GaussianRandomFields
using Distances
using DataFrames
using CSV
using StatsBase
using DataStructures
using ProgressMeter
using Pipe: @pipe
using Combinatorics
###Beta Diveristy
#a function to calculate beta diversity, only for binary data
function beta_diversity(beta_div_matrix)
    if sum((beta_div_matrix .- mean(beta_div_matrix, dims=1)).^2) == 0
        println("The data matrix has no variation in beta space")
        component=DataFrames.DataFrame(
        BDtotal=0, 
        Repl=0, 
        RichDif=0)    
    else
            
        a=beta_div_matrix * transpose(beta_div_matrix)
        b=beta_div_matrix * (1 .- transpose(beta_div_matrix))
        c=(1 .- beta_div_matrix) * transpose(beta_div_matrix)
        min_bc=min.(b, c)
        repl=2 .* min_bc
        rich=abs.(b - c)
        

        #Jaccard-based components
        repl = repl ./ (a + b + c)
        rich = rich ./ (a + b + c)
        D = (b + c) ./ (a + b + c)

        n=size(beta_div_matrix,1) #Count the number of rows
        total_div = (sum(D)/2) / (n * (n - 1))
        repl_div = (sum(repl)/2) / (n * (n - 1))
        rich_div = (sum(rich)/2) / (n * (n - 1))

        component=DataFrames.DataFrame(
            BDtotal=total_div, 
            Repl=repl_div, 
            RichDif=rich_div)
    end 
end

#temporal mean of spatial beta-diversity,  beta-diversity among sites averaged across time
function mean_spatial_beta_div(dynamic_df)
    #Prepare a presence-absence matrix for beta diversity calculation
    beta_matrix = dynamic_df[:,[:Presence, :Species, :Time, :Patch]]#select column N, Species, Time, Patch    
    mean_spatial_beta_div_dict = Dict{Int, Matrix}()  # a dictionary to store data frames for each time points with all the sites
    mean_spatial_beta_div_df = DataFrames.DataFrame()  #a data frame to store beta diversity components from evey time point

    for t in unique(beta_matrix[:, :Time])
        #printIn("Time$t")
        #subset_df=beta_matrix[:,[:Presence, :Species, :Time, :Patch]]
        subset_df=filter(row -> row[:Time]==t, beta_matrix)
        df_wide =unstack(subset_df, :Patch, :Species, :Presence) #pivot wider
        #df_wide .= coalesce.(df_wide, 0) #Fill the missing values will zeros
        df_wide=df_wide[:,2:end] #Remove the patch column
        df_wide=Matrix(df_wide)
        df_wide .= coalesce.(df_wide, 0) #replace missing values with zeros
        #Add the data frame to the dictionary with the time point as the key
        mean_spatial_beta_div_dict[t] = df_wide
    end

    for i in eachindex(mean_spatial_beta_div_dict)
        component=beta_diversity(mean_spatial_beta_div_dict[i])
        mean_spatial_beta_div_df= [mean_spatial_beta_div_df; component]
    end

    mean_spatial_beta_div_summary = DataFrames.DataFrame(
        mean_spatial_BDtotal = mean(mean_spatial_beta_div_df.BDtotal),
        mean_spatial_Repl = mean(mean_spatial_beta_div_df.Repl),
        mean_spatial_RichDif = mean(mean_spatial_beta_div_df.RichDif))
end
#spatial mean of temporal beta-diversity, beta-diversity along all time points averaged across site
function mean_temporal_beta_div(dynamic_df)
    #Prepare a presence-absence matrix for beta diversity calculation
    beta_matrix = dynamic_df[:,[:Presence, :Species, :Time, :Patch]]#select column N, Species, Time, Patch

    mean_temporal_beta_div_dict = Dict{Int, Matrix}()  # a dictionary to store data frames for every site with all time points
    mean_temporal_beta_div_df = DataFrames.DataFrame()  #a data frame to store beta diversity components at all patches

    for p in unique(beta_matrix[:,:Patch])
        subset_df = filter(row -> row[:Patch]==p, beta_matrix)
        df_wide = unstack(subset_df, :Time, :Species, :Presence) #pivot wider
        df_wide = df_wide[:,2:end] #Remove the patch column
        df_wide = Matrix(df_wide)
        df_wide .= coalesce.(df_wide, 0)
        #Add the data frame to the dictionary with the site as the key
        mean_temporal_beta_div_dict[p] = df_wide
    end

    for i in eachindex(mean_temporal_beta_div_dict)
        component = beta_diversity(mean_temporal_beta_div_dict[i])
        mean_temporal_beta_div_df = [mean_temporal_beta_div_df; component]
    end
    mean_temporal_beta_div_summary = DataFrames.DataFrame(
        mean_temporal_BDtotal = mean(mean_temporal_beta_div_df.BDtotal),
        mean_temporal_Repl = mean(mean_temporal_beta_div_df.Repl),
        mean_temporal_RichDif = mean(mean_temporal_beta_div_df.RichDif))
end

