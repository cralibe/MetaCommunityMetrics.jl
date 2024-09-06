# src/OccupiedPatchesProportion.jl

"""
    prop_patches(presence::AbstractVector, species::Union{AbstractVector, String}, patch::Union{AbstractVector, String}) -> DataFrame

Calculate the proportion of patches occupied by each species and summarize the results.

This function takes three vectors: `presence`, `species`, and `patch`, and performs the following steps:

Arguments:
    presence::AbstractVector: A vector indicating the presence (1) or absence (0) of a species in a patch.
    species::Union{AbstractVector, String}: A vector of species names.
    patch::Union{AbstractVector, String}: A vector of patch identifiers.

Returns:
    DataFrame: A DataFrame containing the mean, minimum, and maximum proportion of patches 
               occupied across all species.

Example:
```@repl
    presence = [1, 0, 1, 0, 1]
    species = ["A", "A", "B", "B", "C"]
    patch = [1, 2, 1, 2, 1]
    prop_patches(presence, species, patch) 
````
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
