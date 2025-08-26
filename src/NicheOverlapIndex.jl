# src/NicheOverlapIndex.jl


"""
    niche_overlap(abundance::AbstractVector, species::AbstractVector, site::AbstractVector, time::AbstractVector) -> DataFrame

Calculates the overall mean, maximum, and minimum values of the niche overlap index from all species pairs in the provided data.
# Arguments
- `abundance::AbstractVector`: Vector representing the abundance of species.
- `species::AbstractVector`: Vector representing species names or IDs. 
- `site::AbstractVector`: Vector representing site names or IDs. 
- `time::AbstractVector`: Vector representing sampling dates. 

# Description
The niche overlap index is calculated based on the method suggested by Pianka (1973), with the assumption that the proportional use of resources by a species at a specific site and time is equivalent to its relative abundance at that location and time period across all sampled sites and times.

# Returns
- `DataFrame`: A DataFrame containing the overall mean, maximum, and minimum values of the niche overlap index from all species pairs.

# Example
```@jildoctest
julia> using MetaCommunityMetrics

julia> df = load_sample_data()
53352×12 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude  standardized_temperature  standardized_precipitation 
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
   1 │                0.0923816                      0.0                 0.406837
```
"""
function niche_overlap(abundance::AbstractVector, species::AbstractVector, site::AbstractVector, time::AbstractVector)
        
    # Input validation
    if !(length(abundance) == length(species) == length(site) == length(time))
        throw(ArgumentError("All input vectors must have the same length"))
    end

    if isempty(abundance)  # Check before DataFrame creation
        throw(ArgumentError("Input vectors cannot be empty"))
    end

    # Construct working dataframe
    df = DataFrame(N = abundance, Species = species, Time = time, Patch = site)
    
    # Remove species with zero total abundance
    species_totals = combine(groupby(df, :Species), :N => sum => :total_N)
    valid_species = species_totals[species_totals.total_N .> 0, :Species]
    df = filter(row -> row.Species in valid_species, df)
    
    # Compute relative abundance per species across all patches/times
    proportion_use_df = @pipe df[:, [:N, :Species, :Patch, :Time]] |>
        groupby(_, :Species) |>
        transform(_, :N => (x -> x ./ sum(x)) => :relative_N) |>
        select(_, [:Species, :relative_N, :Patch, :Time]) |>
        unstack(_, :Species, :relative_N, fill = 0) |>
        select(_, Not([:Patch, :Time])) |>
        permutedims(_)
    
    # If fewer than 2 species, return zeros (no pairwise comparisons possible)
    if size(proportion_use_df, 1) < 2
        return DataFrame(
            mean_niche_overlap_index = 0.0,
            min_niche_overlap_index  = 0.0,
            max_niche_overlap_index  = 0.0
        )
    end
    
    # Generate all unique pairwise combinations
    combs = collect(combinations(1:size(proportion_use_df, 1), 2))
    pairwise_scores = Float64[]
    
    for (i, j) in combs
        sp1 = values(proportion_use_df[i, :])
        sp2 = values(proportion_use_df[j, :])
        
        # Calculate Pianka's overlap index
        numerator = sum(sp1[k] * sp2[k] for k in 1:length(sp1))
        denom1 = sum(sp1[k]^2 for k in 1:length(sp1))
        denom2 = sum(sp2[k]^2 for k in 1:length(sp2))
        denom = sqrt(denom1 * denom2)
        
        push!(pairwise_scores, numerator / denom)

    end
    
    # Round to avoid floating point precision errors
    rounded_scores = round.(pairwise_scores, digits = 10)
    
    return DataFrame(
        mean_niche_overlap_index = mean(rounded_scores),
        min_niche_overlap_index  = minimum(rounded_scores),
        max_niche_overlap_index  = maximum(rounded_scores)
    )
end