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


###Niche overlap index
    #we do not have to remove absences from the data frame because of how the niche overlap index is calculated
    #the sum in the nominator and denominator will skip the zeros when calculating the index
    function niche_overlap(dynamic_df)

        proportion_use_df=
        @pipe dynamic_df[:,[:N, :Species, :Patch, :Time]] |>#select column N, Species, Patch, Time and env
        groupby(_, [:Species]) |>
        transform(_, :N => (x -> x ./ sum(x)) => :relativ_N) |> #relative abundance (proportion use) of species i in patch k at time t across all sites and times
        groupby(_, [:Species]) |>
        transform(_, :N => sum => :total_N)|> #total abundance of species i in all sites and times for cross checking
        select(_, [:Species,:relativ_N, :Patch ,:Time]) |>#select columns Species, total_N, relativ_N, env
        unstack(_, :Species,:relativ_N) |> #pivot wider
        _[!, Not(:Patch, :Time)] |> # only retain the Proportional use values for each species
        permutedims(_) #transpose the data frame 
    
    
        combs=collect(combinations(1:size(proportion_use_df,1), 2)) # Generate all combinations of 2 elements (rows) from the indices 1 to n (the number of rows in the data frame)
    
        pairwise_df=hcat(combs, zeros(Float64, size(combs, 1))) #Initialize an array with zeros to contain the niche overlpa index of each pair of species
    
        for i in 1:size(pairwise_df,1)
    
            first_sp=values(proportion_use_df[pairwise_df[i][1],:])
            second_sp=values(proportion_use_df[pairwise_df[i][2],:])
    
            final_numerator = 0.0
            denominator_species_j = 0.0
            denominator_species_k = 0.0
    
            for j in 1:length(first_sp)
                if first_sp[j] !=0.0 && second_sp[j] !=0.0
                    final_numerator += first_sp[j]*second_sp[j]
                    denominator_species_j += first_sp[j]^2
                    denominator_species_k += second_sp[j]^2
                else
                    final_numerator += 0
                    denominator_species_j += 0
                    denominator_species_k += 0
                end
            end
    
            final_denominator = sqrt(denominator_species_j * denominator_species_k)
            pairwise_df[i, 2] = final_numerator / final_denominator
        end
    
        #@save ("/home/jenny/phyto/niche_overlap_index_output/niche_overlap_matrix.jld2") pairwise_df
        
        if isempty(pairwise_df) #check if pairwise_df is empty
        
            niche_overlap_index=DataFrames.DataFrame(mean_niche_overlap_index=0,
                min_niche_overlap_index=0,
                max_niche_overlap_index=0)
            
        else
            niche_overlap_index=DataFrames.DataFrame(mean_niche_overlap_index=mean(filter(!isnan, round.(pairwise_df[:, 2], digits=10))),
                min_niche_overlap_index=minimum(filter(!isnan, round.(pairwise_df[:, 2], digits=10))),
                max_niche_overlap_index=maximum(filter(!isnan, round.(pairwise_df[:, 2], digits=10)))) #rounding to avoid floating-point precision errors
            
        end
        #histogram(pairwise_df[:,2], bins=20, xlabel="Value", ylabel="Frequency", title="Histogram")
        #scatter(1:nrow(test_df), test_df[:,5], xlabel="Row ID", ylabel="Column Value", title="Scatter Plot")
    end