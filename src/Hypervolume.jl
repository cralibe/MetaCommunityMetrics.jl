# src/Hypervolume.jl


"""
    MVNH_det(data::DataFrame; var_names::Vector{String}=String[]) -> DataFrame

Calculate the niche volume for a species in a given biodiversity data. 

Arguments
- `data::DataFrame`: A dataframe where each row represents an obsevation of species and each column represent the environmental variables for that observation.
- `var_names::Vector{String}`: A vector of strings with the names of the environmental variables. Default is an empty vector.
Returns
- `DataFrame`: A DataFrame with the following columns:
    - `total`: The total hypervolume of the niche.
    - `cor`: The correlation component of the hypervolume, also known as the determinant of the correlation matrix.
    - `variable_i`: The variance of the environmental variable i.

Details
- Evirnomental variables are assumed to be normally distributed.
- Evirnomental variables are encouraged to be normalized to avoid bias due to different scaling before using this function.
- This function is a translation/adaptation of the MVNH_det function from the R package `MVNH`,licensed under GPL-3.
- Original package and documentation available at: https://github.com/lvmuyang/MVNH

Example
```jildoctest
julia> using MetaCommunityMetrics, Pipe, DataFrames, Statistics

julia> df = load_sample_data()
48735×10 DataFrame
   Row │ Year   Month  Day    Sampling_date_order  plot   Species  Abundance  Presence  Latitude  Longitude 
       │ Int64  Int64  Int64  Int64                Int64  String3  Int64      Int64     Float64   Float64   
───────┼────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │  2010      1     16                    1      1  BA               0         0      35.0     -110.0
     2 │  2010      1     16                    1      2  BA               0         0      35.0     -109.5
     3 │  2010      1     16                    1      8  BA               0         0      35.5     -109.5
     4 │  2010      1     16                    1      9  BA               0         0      35.5     -109.0
     5 │  2010      1     16                    1     11  BA               0         0      35.5     -108.0
   ⋮   │   ⋮      ⋮      ⋮             ⋮             ⋮       ⋮         ⋮         ⋮         ⋮          ⋮
 48731 │  2023      3     21                  117      9  SH               0         0      35.5     -109.0
 48732 │  2023      3     21                  117     10  SH               0         0      35.5     -108.5
 48733 │  2023      3     21                  117     12  SH               1         1      35.5     -107.5
 48734 │  2023      3     21                  117     16  SH               0         0      36.0     -108.5
 48735 │  2023      3     21                  117     23  SH               0         0      36.5     -108.0
                                                                                          48725 rows omitted


julia> result = MVNH_det(df; var_names=["Latitude", "Longitude"])

```

"""

function MVNH_det(data::DataFrame; var_names::Vector{String}=String[]) 
    # If variable names are not provided, generate default names
    if isempty(var_names)
        var_names = ["variable" * string(i) for i in 1:size(data, 2)]
    end

    # Convert DataFrame to Matrix (drop non-numeric columns if any)
    data_matrix = Matrix(data[:, :])

    # Compute the covariance matrix from the data
    COV = cov(data_matrix)

    # Calculate the variance (diagonal of the covariance matrix)
    s = diag(COV)

    # Calculate the correlation component (determinant ratio)
    rho = det(COV) / prod(s)

    # Prepare the data for creating the DataFrame
    # Start with the determinant and correlation component
    df_data = Dict("total" => det(COV), "Correlation" => rho)

    # Add each variance with corresponding variable name
    for i in 1:length(var_names)
        df_data[var_names[i]] = s[i]
    end

    # Convert to a DataFrame
    result = DataFrame(df_data)

    return result
end

"""
    MVNH_dissimilarity(data_1::DataFrame, data_2::DataFrame; var_names::Vector{String}=String[]) -> Dict

Calculate the niche overlap between two species in a given biodiversity data. 

"""
function MVNH_dissimilarity(data_1::DataFrame, data_2::DataFrame; var_names::Vector{String}=String[]) 
    # If variable names are not provided, generate default names
    if isempty(var_names)
        var_names = ["variable" * string(i) for i in 1:size(db1, 2)]
    end

    # Compute covariance matrices for both datasets
    db1 = Matrix(data_1[:, :])
    db2 = Matrix(data_2[:, :])
    S1 = cov(db1)
    S2 = cov(db2)
    S_avg = (S1 + S2) / 2  # Average covariance matrix

    # Compute the means for both datasets and convert them to vectors
    u1 = vec(mean(db1, dims=1))  # Convert 1-row matrix to vector
    u2 = vec(mean(db2, dims=1))  # Convert 1-row matrix to vector

    # Calculate the Mahalanobis component
    diff = u1 .- u2  # Element-wise subtraction of vectors
    mahalanobis_dist = (1 / 8) * (diff' * inv(S_avg) * diff)
    mahalanobis_1d = (1 / 8) * (diff.^2) ./ diag(S_avg)
    mahalanobis_cor = mahalanobis_dist[1] - sum(mahalanobis_1d)

    # Calculate the determinant ratio component
    determinant_ratio = (1 / 2) * log(det(S_avg) / sqrt(det(S1) * det(S2)))
    determinant_1d = (1 / 2) * log.(diag(S_avg) ./ sqrt.(diag(S1) .* diag(S2)))
    determinant_cor = determinant_ratio - sum(determinant_1d)

    # Total Bhattacharyya distance
    bhattacharyya_dist = mahalanobis_dist[1] + determinant_ratio
    bhattacharyya_1d = mahalanobis_1d + determinant_1d
    bhattacharyya_cor = mahalanobis_cor + determinant_cor

    # Create Dictionary for each component
    mahalanobis_results = Dict("total" => mahalanobis_dist[1], "correlation" => mahalanobis_cor)
    determinant_results = Dict("total" => determinant_ratio, "correlation" => determinant_cor)
    bhattacharyya_results = Dict("total" => bhattacharyya_dist, "correlation" => bhattacharyya_cor)

    for i in 1:length(var_names)
        mahalanobis_results[var_names[i]] = mahalanobis_1d[i]
        determinant_results[var_names[i]] = determinant_1d[i]
        bhattacharyya_results[var_names[i]] = bhattacharyya_1d[i]
    end

    # Convert to DataFrames
    mahalanobis_df = DataFrame(mahalanobis_results)
    determinant_df = DataFrame(determinant_results)
    bhattacharyya_df = DataFrame(bhattacharyya_results)

    # Set the desired column order: ["total", "correlation", var_names...]
    column_order = vcat(["total", "correlation"], var_names)

    # Reorder columns in the DataFrames
    mahalanobis_df = select(mahalanobis_df, column_order)
    determinant_df = select(determinant_df, column_order)
    bhattacharyya_df = select(bhattacharyya_df, column_order)

    # Return all three DataFrames in a dictionary
    return Dict(
        "Bhattacharyya_distance" => bhattacharyya_df,
        "Mahalanobis_distance" => mahalanobis_df,
        "Determinant_ratio" => determinant_df
    )

end

"""
    MVNH_dissimilarity(data_1::DataFrame, data_2::DataFrame; var_names::Vector{String}=String[]) -> Dict

Calculate the avaraged niche volume across all species in a given biodiversity data.

"""
function average_MVNH_det(data::DataFrame, presence_absence::Vector{Int}, species::Vector{String}; var_names::Vector{String}=String[])
    # Ensure the presence_absence and species vectors match the DataFrame size
    @assert length(presence_absence) == nrow(data) "Mismatch in data length"
    @assert length(species) == nrow(data) "Mismatch in species length"

    # Add presence_absence and species as columns to the DataFrame
    data = hcat(data, DataFrame(presence_absence=presence_absence, species=species))

    # Initialize final result DataFrame
    final_result = DataFrame(species=String[], hypervolume=Float64[])

    # Loop over each unique species
    for sp in unique(species)
        # Filter rows for the current species with presence > 0
        filtered_data = @pipe data |>
            filter(row -> row.species == sp && row.presence_absence == 1, _) |>
            select(_, Not(:species,:presence_absence))

        # Skip if no valid rows for the species
        if nrow(filtered_data) == 0
            continue
        end     
            
        # Compute the MVNH_det result
        result = MVNH_det(filtered_data; var_names)

        # Add the species and its total hypervolume to the final result DataFrame
        push!(final_result, (species=sp, hypervolume=result.total[1]))

    end

    # Compute the average hypervolume across species
    averaged_hypervolume = mean(final_result.hypervolume)

    return averaged_hypervolume

end

"""
    MVNH_dissimilarity(data_1::DataFrame, data_2::DataFrame; var_names::Vector{String}=String[]) -> Dict

Calculate the avearaged niche overlap across all species pairs in a given biodiversity data. 

"""
function average_MVNH_dissimilarity(data::DataFrame, presence_absence::Vector{Int}, species::Vector{String}; var_names::Vector{String}=String[])
    @assert length(presence_absence) == nrow(data) "Mismatch in data length"
    @assert length(species) == nrow(data) "Mismatch in species length"

    # Add presence_absence and species as columns to the DataFrame
    data = hcat(data, DataFrame(presence_absence=presence_absence, species=species))

    # Initialize final result DataFrame
    final_result = DataFrame(species_1=String[], species_2=String[], hypervolume_dissimilarity=Float64[])

     # Loop over each unique species pair
    for sp1 in unique(species)
        for sp2 in unique(species)
            # Skip if the same species or if the pair has already been processed or if the pair is the reverse
            if sp1 == sp2 || nrow(filter(row -> row.species_1 == sp1 && row.species_2 == sp2, final_result)) > 0 || nrow(filter(row -> row.species_1 == sp2 && row.species_2 == sp1, final_result)) > 0
                continue
            end

            # Filter rows for the current species pair with presence > 0
            filtered_data_1 = @pipe data |>
                filter(row -> row.species == sp1 && row.presence_absence == 1, _) |>
                select(_, Not(:species,:presence_absence))

            filtered_data_2 = @pipe data |>
                filter(row -> row.species == sp2 && row.presence_absence == 1, _) |>
                select(_, Not(:species,:presence_absence))

            # Skip if no valid rows for the species pair
            if nrow(filtered_data_1) == 0 || nrow(filtered_data_2) == 0
                continue
            end

            # Compute the MVNH_dissimilarity result
            result = MVNH_dissimilarity(filtered_data_1, filtered_data_2; var_names)

            # Add the species pair and its hypervolume dissimilarity to the final result DataFrame
            push!(final_result, (species_1=sp1, species_2=sp2, hypervolume_dissimilarity=result["Bhattacharyya_distance"].total[1]))
        
        end
    end

    mean_hypervolume_dissimilarity = mean(final_result.hypervolume_dissimilarity)
            
    return mean_hypervolume_dissimilarity
    
end