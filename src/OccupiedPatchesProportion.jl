# src/OccupiedPatchesProportion.jl

"""
    prop_patches(presence::AbstractVector, species::Union{AbstractVector, String}, patch::Union{AbstractVector, String}) -> DataFrame

Calculate the proportion of patches occupied by each species and summarize the results.

This function takes three vectors: `presence`, `species`, and `patch`, and performs the following steps:

Arguments
- `presence::AbstractVector`: A vector indicating the presence (1) or absence (0) of a species in a patch.
- `species::Union{AbstractVector, String}`: A vector of species names.
- `patch::Union{AbstractVector, String}`: A vector of patch identifiers.

Returns
- `DataFrame`: A DataFrame containing the mean, minimum, and maximum proportion of patches 
               occupied across all species.

Example
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

                                                                                          
julia> prop_patches(df.Presence, df.Species, df.plot)
1×3 DataFrame
 Row │ mean_prop_patches  min_prop_patches  max_prop_patches 
     │ Float64            Float64           Float64          
─────┼───────────────────────────────────────────────────────
   1 │          0.734649         0.0833333               1.0
```
"""
function prop_patches(presence::AbstractVector, species::Union{AbstractVector, String}, patch::Union{AbstractVector, String})

    df = DataFrames.DataFrame(
        Presence=presence,
        Species=species,
        Patch=patch,
       )

    species_patches_df= 
        @pipe df[:,[:Presence, :Species, :Patch]]|>#select column N, Species, Patch
        groupby(_, [:Patch, :Species])|>
    combine(_, :Presence => sum => :Total_N)|>
    transform(_, :Total_N=> ByRow(x->ifelse(x> 0, 1.0, 0.0))=> :Presence)|>
    groupby(_, [:Species])


    patches_occupied_matrix=hcat(unique(df.Species), zeros(Float64, size(unique(df.Species), 1))) 
    for group in 1:size(patches_occupied_matrix,1)
        patches_occupied_matrix[group,2]= sum(species_patches_df[group].Presence)/nrow(species_patches_df[group])
    end

    prop_patches_df=DataFrames.DataFrame(mean_prop_patches = mean(patches_occupied_matrix[:,2]),
        min_prop_patches = minimum(patches_occupied_matrix[:,2]),
        max_prop_patches = maximum(patches_occupied_matrix[:,2]))   
end
