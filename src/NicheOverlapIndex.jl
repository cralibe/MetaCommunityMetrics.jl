# src/NicheOverlapIndex.jl


"""
    niche_overlap(abundance::AbstractVector, species::Union{AbstractVector, String}, patch::Union{AbstractVector, String}, time::AbstractVector) -> DataFrame

Calculates the overall mean, maximum, and minimum values of the niche overlap index from all species pairs in the provided data.
# Arguments
- `abundance::AbstractVector`: Vector representing the abundance of species.
- `species::Union{AbstractVector, String}`: Vector or string representing species names or IDs.
- `patch::Union{AbstractVector, String}`: Vector or string representing patch names or IDs.
- `time::AbstractVector`: Vector representing the time points.

# Description
The niche overlap index is calculated based on the method suggested by Pianka (1973), with the assumption that the proportional use of a species at a specific site and time equals its relative abundance at that site and time. To determine relative abundance, the abundance of each species in a particular patch is divided by the total abundance of that species across all patches and times.

# Returns
- `DataFrame`: A DataFrame containing the overall mean, maximum, and minimum values of the niche overlap index from all species pairs.

# Example
```@jildoctest
julia> using MetaCommunityMetrics

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  normalized_temperature  normalized_precipitation 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64    Float64                 Float64                  
───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0                0.829467              -1.4024
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5               -1.12294               -0.0519895
     3 │  2010      1     16                    1      4  BA               0         0      35.0     -108.5               -0.409808              -0.803663
     4 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5               -1.35913               -0.646369
     5 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0                0.0822                 1.09485
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮                ⋮                        ⋮
 53348 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0               -0.571565              -0.836345
 53349 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5               -2.33729               -0.398522
 53350 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5                0.547169               1.03257
 53351 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5               -0.815015               0.95971
 53352 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0                0.48949               -1.59416
                                                                                                                                            53342 rows omitted
                                                                                          
julia> result = niche_overlap(df.Abundance, df.Species, df.plot, df.Sampling_date_order)
1×3 DataFrame
 Row │ mean_niche_overlap_index  min_niche_overlap_index  max_niche_overlap_index 
     │ Float64                   Float64                  Float64                 
─────┼────────────────────────────────────────────────────────────────────────────
   1 │                 0.827739                 0.591836                      1.0
```
"""
function niche_overlap(abundance::AbstractVector, species::Union{AbstractVector, String}, patch::Union{AbstractVector, String}, time::AbstractVector)
    
    df = DataFrames.DataFrame(
        N=abundance,
        Species=species,
        Time=time,
        Patch=patch
        )
    
    proportion_use_df=
    @pipe df[:,[:N, :Species, :Patch, :Time]] |>#select column N, Species, Patch, Time and env
    groupby(_, [:Species]) |>
    transform(_, :N => (x -> x ./ sum(x)) => :relativ_N) |> #relative abundance (proportion use) of species i in patch k at time t across all sites and times
    groupby(_, [:Species]) |>
    transform(_, :N => sum => :total_N)|> #total abundance of species i in all sites and times for cross checking
    select(_, [:Species,:relativ_N, :Patch ,:Time]) |>#select columns Species, total_N, relativ_N, env
    unstack(_, :Species,:relativ_N, fill = 0) |> #pivot wider
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